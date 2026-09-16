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

#include <cmath>
#include <cstdlib>
#include <format>
#include <memory>

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsTPC/Constants.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/WorkflowHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/DetectorTPC.h"
#include "O2Align/Params.h"
#include "O2Align/SensorTPC.h"
#include "O2Align/TrackFit.h"
#include "TPCFastTransformPOD.h"

#include "GPUO2ExternalUser.h"
#include "GPUParam.h"
#include "GPUParam.inc"

namespace o2::alignrs
{

void DetectorTPC::prepareData(o2::globaltracking::RecoContainer* /*recoData*/)
{
  // Nothing can be precomputed per TF: the cluster transformation needs the track time and both
  // the supercluster building and the cluster errors need the track angles, hence everything is
  // done in prepareTrack.
}

Volume::Ptr DetectorTPC::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  // The TPC has no alignable entries in the geometry: both the envelope and the sector volumes
  // are fictitious and define their frames themselves.
  uint32_t gLbl{0};
  const uint32_t det = mDetIdx;
  auto root = std::make_unique<Volume>("TPC_envelope", gLbl++, det, false, true);
  mSensors.assign(o2::tpc::constants::MAXSECTOR, nullptr);
  for (int isec = 0; isec < o2::tpc::constants::MAXSECTOR; ++isec) {
    const Label lbl(det, isec, true);
    const auto symName = std::format("TPC/sec{:02d}", isec);
    auto* sector = root->addChild<SensorTPC>(symName.c_str(), lbl, true);
    sector->setSensorId(isec);
    sensorMap[lbl] = sector;
    mSensors[isec] = sector;
  }
  return root;
}

bool DetectorTPC::prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack)
{
  const auto& params = Params::Instance();
  const auto gid = ids[GTrackID::TPC];
  if (!gid.isIndexSet()) {
    return false;
  }
  if (!mCorrMaps || !mTPCParam) {
    LOGP(fatal, "TPC correction maps ({}) and GPUParam ({}) must be set before processing TPC tracks",
         (void*)mCorrMaps, (void*)mTPCParam);
  }
  const auto& trk = recoData->getTrack<o2::tpc::TrackTPC>(gid);
  const int nClus = trk.getNClusters();
  if (nClus < params.minTPCClusters) {
    return false;
  }
  const auto clusterIdxStruct = recoData->getTPCTracksClusterRefs();
  const auto clusterNativeAccess = recoData->inputsTPCclusters->clusterIndex;
  const auto shMap = recoData->clusterShMapTPC;

  float trackTime{0.f}, trackTimeErr{0.f};
  recoData->getTrackTime(resTrack.gid, trackTime, trackTimeErr);
  const float tOffset = trackTime / (o2::constants::lhc::LHCBunchSpacingMUS * 8);

  auto prop = o2::base::PropagatorD::Instance();
  // we continue the fit of the seed prepared by the inner detectors outward, w/o resetting its cov.matrix
  auto trkParam = resTrack.track;
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trkParam);
  }
  float chi2 = 0.f;

  constexpr float TAN10 = 0.17632698f; // tan of the half sector opening angle
  constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
  const size_t nPointsIni = resTrack.info.size();
  // the min number of points is relaxed by each merging of clusters into a supercluster
  int npntCut = params.minTPCClusters;
  int npoints = 0;
  bool stopLoop = false;
  const o2::tpc::ClusterNative* cl = nullptr;
  uint8_t sector = 0, row = 0, currentSector = 0, currentRow = 0;
  int16_t clusterState = 0, nextState = 0;

  // the clusters are ordered from the outermost to the innermost, we traverse them in the outward direction
  for (int i = nClus - 1; i >= 0; i -= cl ? 0 : 1) {
    float x{0.f}, y{0.f}, z{0.f}, xTmp{0.f}, yTmp{0.f}, zTmp{0.f}, charge{0.f};
    int clusters = 0;
    double combRow = 0;

    while (true) {
      if (!cl) {
        const auto* clTmp = &trk.getCluster(clusterIdxStruct, i, clusterNativeAccess, sector, row);
        if (row > params.maxTPCPadRow) { // outward refit: all following clusters will have a larger padrow
          stopLoop = true;
          break;
        }
        if (row < params.minTPCPadRow) { // the following clusters still have a chance to be accepted
          break;
        }
        if (params.discardEdgePadrows > 0 && getDistanceToStackEdge(row) < params.discardEdgePadrows) {
          if (i != 0) {
            --i;
            continue;
          }
          stopLoop = true;
          break;
        }
        mCorrMaps->Transform(sector, row, clTmp->getPad(), clTmp->getTime(), xTmp, yTmp, zTmp, tOffset);
        if (params.discardSectorEdgeDepth > 0 && std::abs(yTmp) + params.discardSectorEdgeDepth > xTmp * TAN10) {
          if (i != 0) {
            --i;
            continue;
          }
          stopLoop = true;
          break;
        }
        cl = clTmp;
        nextState = shMap[cl - clusterNativeAccess.clustersLinear];
      }
      if (clusters == 0 || (sector == currentSector && std::abs(row - currentRow) < params.maxTPCRowsCombined)) {
        if (clusters == 1) { // charge-weight the 1st cluster before merging the 2nd one
          x *= charge;
          y *= charge;
          z *= charge;
          combRow *= charge;
        }
        if (clusters == 0) { // start a new supercluster
          x = xTmp;
          y = yTmp;
          z = zTmp;
          currentRow = row;
          currentSector = sector;
          charge = cl->getQtot();
          clusterState = nextState;
          combRow = row;
        } else { // merge to the supercluster started at currentRow
          x += xTmp * cl->getQtot();
          y += yTmp * cl->getQtot();
          z += zTmp * cl->getQtot();
          combRow += row * cl->getQtot();
          charge += cl->getQtot();
          clusterState |= nextState;
          --npntCut;
        }
        cl = nullptr;
        ++clusters;
        if (i != 0) {
          --i;
          continue;
        }
      }
      break;
    }
    if (stopLoop) {
      break;
    }
    if (clusters == 0) {
      continue;
    }
    if (clusters > 1) {
      x /= charge;
      y /= charge;
      z /= charge;
      currentRow = static_cast<uint8_t>(combRow / charge);
    }

    const double alpha = o2::math_utils::detail::sector2Angle<double>(currentSector % NSectorsPerSide);
    if (!prop->propagateToAlphaX(trkParam, refLin, alpha, x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      break;
    }
    std::array<float, 3> cov{0.f, 0.f, 0.f};
    // TODO: this disables the occupancy / charge components of the error estimation
    mTPCParam->GetClusterErrors2(currentSector, currentRow, z, trkParam.getSnp(), trkParam.getTgl(), -1.f, 0.f, 0.f, cov[0], cov[2]);
    mTPCParam->UpdateClusterError2ByState(clusterState, cov[0], cov[2]);
    const int nrComb = std::abs(row - currentRow) + 1;
    if (nrComb > 1) {
      const float fact = 1.f / std::sqrt(static_cast<float>(nrComb));
      cov[0] *= fact;
      cov[2] *= fact;
    }
    cov[0] += params.extraClsErrYTPC * params.extraClsErrYTPC;
    cov[2] += params.extraClsErrZTPC * params.extraClsErrZTPC;

    const std::array<double, 2> pos{y, z};
    const std::array<double, 3> covD{cov[0], cov[1], cov[2]};
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pos, covD));
    if (!trkParam.update(pos, covD)) {
      break;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(pos[0]);
      refLin->setZ(pos[1]);
    }

    auto& pnt = resTrack.info.emplace_back();
    pnt.lr = static_cast<int8_t>(getStack(currentRow));
    pnt.label = Label(mDetIdx, currentSector, true);
    pnt.x = x;
    pnt.alpha = static_cast<float>(alpha);
    pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(currentSector), x, y, z, cov[0], cov[2], cov[1]);
    ++npoints;
  }

  if (npoints < npntCut) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  resTrack.track = trkParam; // the seed is replaced only by a successful fit
  resTrack.kfFit.chi2 += chi2;
  return true;
}

} // namespace o2::alignrs
