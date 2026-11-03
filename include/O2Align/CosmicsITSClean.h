// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Selection of the ITS road hits of a B = 0 cosmic (TrackCosmicsExtended) which lie on it.
//
// The ITS road of the collector keeps, per half of the cosmic and layer, the cluster closest to the trajectory predicted by the TPC legs
// (in the inner barrel the 5 closest), from every ROF in the cosmic's time window. Without a TOF time this window is the TPC one (tens of
// mus), hence the large ROF spread and the z tolerance of tens of cm: most of the IB road hits and 5-20 % of the OB ones are not on the
// cosmic. The `leg` of a road hit is the half of the cosmic (global y below / above the closest approach to the beam line), not a
// half-barrel and not a validated association.
//
// clean(): a straight line (principal axis) through the OB hits (r > 10 cm); while the worst OB hit is farther than `cut` from the line in
// its sensor plane and more than 3 OB hits are left, it is dropped and the line refitted; then every hit, IB included, is kept only
// within `cut` of that line. Straight lines only: B = 0.
// gold(): clean() + a TOF time (the road is then a single ROF), no duplicate, >= minOBPerHalf kept OB hits on each half, kept hits
// within rofWindow ROFs of the majority ROF.
#pragma once

#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITSBase/GeometryTGeo.h"

namespace o2::alignrs::cosmics_its_clean
{

/// global position of a stored ITS cluster: the stored pattern, else the dictionary, else the anchor pixel's centre
inline o2::math_utils::Point3D<float> globalPosition(const o2::dataformats::TrackCosmicsExtended& cosm, const o2::dataformats::CosmicITSCluster& c,
                                                     const o2::itsmft::TopologyDictionary* dict, const o2::its::GeometryTGeo& geom)
{
  const o2::itsmft::CompClusterExt compCluster(c.row, c.col, c.pattID, c.chipID);
  o2::math_utils::Point3D<float> local;
  if (c.pattEntry >= 0) {
    o2::itsmft::ClusterPattern pattern;
    auto pattIt = cosm.itsPatterns.begin() + c.pattEntry;
    pattern.acquirePattern(pattIt);
    local = o2::itsmft::TopologyDictionary::getClusterCoordinates<float>(compCluster, pattern, c.pattID != o2::itsmft::CompCluster::InvalidPatternID);
  } else if (dict && c.pattID != o2::itsmft::CompCluster::InvalidPatternID) {
    local = dict->getClusterCoordinates<float>(compCluster);
  } else {
    o2::itsmft::SegmentationAlpide::detectorToLocalUnchecked(c.row, c.col, local);
  }
  return geom.getMatrixL2G(c.chipID) * local;
}

/// keep[i] = cosm.clITS[i] is on the cosmic's straight line; false (keep all false) without 3 OB hits on one line within cut [cm].
/// geom: GeometryTGeo with the L2G matrices filled (fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G))).
inline bool clean(const o2::dataformats::TrackCosmicsExtended& cosm, const o2::itsmft::TopologyDictionary* dict, const o2::its::GeometryTGeo& geom,
                  std::vector<bool>& keep, double cut = 0.03)
{
  const int n = (int)cosm.clITS.size();
  keep.assign(n, false);
  std::vector<Eigen::Vector3d> global(n);
  std::vector<double> alpha(n);
  std::vector<bool> inLine(n);
  int nLine = 0;
  for (int i = 0; i < n; ++i) {
    const auto p = globalPosition(cosm, cosm.clITS[i], dict, geom);
    global[i] = Eigen::Vector3d(p.X(), p.Y(), p.Z());
    alpha[i] = geom.getSensorRefAlpha(cosm.clITS[i].chipID);
    inLine[i] = std::hypot(p.X(), p.Y()) > 10.;
    nLine += inLine[i];
  }
  if (nLine < 3) {
    return false;
  }
  std::vector<double> residual(n, 0.);
  while (true) {
    Eigen::Vector3d centre = Eigen::Vector3d::Zero();
    for (int i = 0; i < n; ++i) {
      if (inLine[i]) {
        centre += global[i];
      }
    }
    centre /= nLine;
    Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n; ++i) {
      if (inLine[i]) {
        const Eigen::Vector3d d = global[i] - centre;
        scatter += d * d.transpose();
      }
    }
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(scatter);
    const Eigen::Vector3d direction = solver.eigenvectors().col(2); // eigenvalues ascending: the line's direction
    int worst = -1;
    double worstResidual = 0.;
    for (int i = 0; i < n; ++i) {
      // residual in the sensor plane: the line intersected with the plane x_tracking = const of the hit's sensor
      const double cosAlpha = std::cos(alpha[i]);
      const double sinAlpha = std::sin(alpha[i]);
      const double normalDirection = cosAlpha * direction.x() + sinAlpha * direction.y();
      if (std::abs(normalDirection) < 1e-3) {
        residual[i] = std::numeric_limits<double>::infinity();
      } else {
        const double xHit = cosAlpha * global[i].x() + sinAlpha * global[i].y();
        const double t = (xHit - (cosAlpha * centre.x() + sinAlpha * centre.y())) / normalDirection;
        const Eigen::Vector3d r = global[i] - (centre + t * direction);
        residual[i] = std::hypot(-sinAlpha * r.x() + cosAlpha * r.y(), r.z());
      }
      if (inLine[i] && residual[i] > worstResidual) {
        worstResidual = residual[i];
        worst = i;
      }
    }
    if (worstResidual <= cut) {
      break;
    }
    if (nLine <= 3) {
      return false;
    }
    inLine[worst] = false;
    --nLine;
  }
  for (int i = 0; i < n; ++i) {
    keep[i] = residual[i] <= cut;
  }
  return true;
}

struct GoldSettings {
  bool requireTOF = true; ///< the cosmic has a TOF time (TrackCosmicsExtended::hasTOFTime())
  double cut = 0.03;      ///< clean(): max. distance of a hit from the OB line in its sensor plane [cm]
  int minOBPerHalf = 2;   ///< kept OB hits (r > 10 cm) on each half of the cosmic (CosmicITSCluster::leg)
  int rofWindow = 1;      ///< kept hits within this many ROFs of the majority ROF of the kept hits (< 0: no ROF check)
  int rofLengthBC = 198;  ///< ITS ROF length [BC] of the run (o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::roFrameLengthInBC)
};

/// keep[i] = cosm.clITS[i] is a gold hit; false (keep all false) if the cosmic is not gold
inline bool gold(const o2::dataformats::TrackCosmicsExtended& cosm, const o2::itsmft::TopologyDictionary* dict, const o2::its::GeometryTGeo& geom,
                 std::vector<bool>& keep, const GoldSettings& settings = {})
{
  const int n = (int)cosm.clITS.size();
  keep.assign(n, false);
  if ((settings.requireTOF && !cosm.hasTOFTime()) || cosm.duplicateOf >= 0 || !clean(cosm, dict, geom, keep, settings.cut)) {
    keep.assign(n, false);
    return false;
  }
  if (settings.rofWindow >= 0) {
    std::map<int, int> rofs;
    for (int i = 0; i < n; ++i) {
      if (keep[i]) {
        rofs[cosm.clITS[i].rofBC]++;
      }
    }
    int majority = 0;
    int majorityRof = 0;
    for (const auto& [rof, count] : rofs) {
      if (count > majority) {
        majority = count;
        majorityRof = rof;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (keep[i] && std::abs(cosm.clITS[i].rofBC - majorityRof) > settings.rofWindow * settings.rofLengthBC) {
        keep[i] = false;
      }
    }
  }
  int nOB[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    if (keep[i] && geom.getLayer(cosm.clITS[i].chipID) >= 3) {
      nOB[cosm.clITS[i].leg & 1]++;
    }
  }
  if (nOB[0] < settings.minOBPerHalf || nOB[1] < settings.minOBPerHalf) {
    keep.assign(n, false);
    return false;
  }
  return true;
}

} // namespace o2::alignrs::cosmics_its_clean
