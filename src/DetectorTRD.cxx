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

#include <array>
#include <cmath>
#include <format>
#include <memory>

#include <TGeoManager.h>
#include <TMath.h>

#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TrackTRD.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/DetectorTRD.h"
#include "O2Align/Params.h"
#include "O2Align/SensorTRD.h"
#include "O2Align/TrackFit.h"
#include "TRDBase/Geometry.h"
#include "TRDBase/TrackletTransformer.h"

namespace o2::alignrs
{

Volume::Ptr DetectorTRD::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  auto geo = o2::trd::Geometry::instance();
  geo->createPadPlaneArray();
  geo->createClusterMatrixArray(); // ideal T2L matrices

  uint32_t gLbl{0};
  const uint32_t det = mDetIdx;
  // the TRD envelope has no alignable entry in the geometry
  auto root = std::make_unique<Volume>("TRD_envelope", gLbl++, det, false, true);
  mSensors.assign(o2::trd::constants::MAXCHAMBER, nullptr);

  for (int isector = 0; isector < o2::trd::constants::NSECTOR; ++isector) {
    const auto symSector = std::format("TRD/sm{:02d}", isector);
    if (!gGeoManager->GetAlignableEntry(symSector.c_str())) { // disabled supermodule
      continue;
    }
    Volume* volSector = nullptr;
    for (int ilr = 0; ilr < o2::trd::constants::NLAYER; ++ilr) {
      for (int istack = 0; istack < o2::trd::constants::NSTACK; ++istack) {
        const auto symChamber = std::format("TRD/sm{:02d}/st{}/pl{}", isector, istack, ilr);
        if (!gGeoManager->GetAlignableEntry(symChamber.c_str())) { // missing chamber (e.g. PHOS hole)
          continue;
        }
        if (!volSector) { // create the supermodule volume only once, and only if it has chambers
          volSector = root->addChild(symSector.c_str(), gLbl++, det, false);
        }
        const int sid = o2::trd::Geometry::getDetector(ilr, istack, isector);
        const Label lbl(det, sid, true);
        auto* chamber = volSector->addChild<SensorTRD>(symChamber.c_str(), lbl);
        chamber->setSensorId(sid);
        sensorMap[lbl] = chamber;
        mSensors[sid] = chamber;
      }
    }
  }
  return root;
}

void DetectorTRD::prepareData(o2::globaltracking::RecoContainer* recoData)
{
  // The transformation of the raw tracklet to the calibrated LOCAL frame position does not depend
  // on the track, hence it is done once per TF. The tilt correction and the covariance depend on
  // the track and are applied in prepareTrack.
  if (!mTransformer) {
    LOGP(fatal, "TRD tracklet transformer must be set before processing TRD data");
  }
  const auto trackletsRaw = recoData->getTRDTracklets();
  mTrackletsLoc.clear();
  mTrackletsLoc.reserve(trackletsRaw.size());
  for (const auto& tracklet : trackletsRaw) {
    mTrackletsLoc.push_back(mTransformer->transformTracklet(tracklet, false));
  }
}

bool DetectorTRD::prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack)
{
  const auto& params = Params::Instance();
  const auto gid = ids[GTrackID::TRD];
  if (!gid.isIndexSet()) {
    return false;
  }
  const auto& trk = recoData->getTrack<o2::trd::TrackTRD>(gid);
  if (trk.getNtracklets() < params.minTRDTracklets) {
    return false;
  }
  auto prop = o2::base::PropagatorD::Instance();
  if (!mRecoParamInit) {
    mRecoParam.init(prop->getNominalBz());
    mRecoParamInit = true;
  }
  const auto trackletsRaw = recoData->getTRDTracklets();
  const auto* geo = o2::trd::Geometry::instance();
  const size_t nPointsIni = resTrack.info.size();

  // we refit the outer param inward to get the tracklet coordinates accounting for the tilt
  o2::track::TrackParD trkParam = convertTrack<double>(trk.getOuterParam());

  for (int il = o2::trd::constants::NLAYER; il--;) {
    const int trkltId = trk.getTrackletIndex(il);
    if (trkltId < 0) {
      continue;
    }
    const auto& trackletRaw = trackletsRaw[trkltId];
    const int trkltDet = trackletRaw.getDetector();
    auto* sensor = mSensors[trkltDet];
    if (!sensor) {
      LOGP(error, "TRD chamber {} referred by a track has no alignable volume", trkltDet);
      resTrack.info.resize(nPointsIni);
      return false;
    }
    // the calibrated tracklet is in the LOCAL frame, T2L of the (possibly pre-aligned) volume brings it to TRK
    const auto& trkltLoc = mTrackletsLoc[trkltId];
    const double locXYZ[3] = {trkltLoc.getX(), trkltLoc.getY(), trkltLoc.getZ()};
    double traXYZ[3]{};
    sensor->getT2L().MasterToLocal(locXYZ, traXYZ);

    const int trkltSec = o2::trd::Geometry::getSector(trkltDet);
    const double alpSens = o2::math_utils::detail::sector2Angle<double>(trkltSec);
    if (trkltSec != o2::math_utils::angle2Sectord(trkParam.getAlpha()) ||
        !trkParam.rotateParam(alpSens) ||
        // we don't need high precision here
        !prop->propagateTo(trkParam, traXYZ[0], false, o2::base::PropagatorD::MAX_SIN_PHI, 10., o2::base::PropagatorD::MatCorrType::USEMatCorrNONE)) {
      resTrack.info.resize(nPointsIni);
      return false;
    }

    const o2::trd::PadPlane* pad = geo->getPadPlane(trkltDet);
    const double tilt = std::tan(TMath::DegToRad() * pad->getTiltingAngle()); // tilt is signed and returned in degrees
    const double padLength = pad->getRowSize(trackletRaw.getPadRow());
    double tiltCorr = tilt * (traXYZ[2] - trkParam.getZ());
    if (std::abs(traXYZ[2] - trkParam.getZ()) < padLength) { // RS do we need this?
      tiltCorr = 0.;
    }
    double posY = traXYZ[1] - tiltCorr;
    const double posZ = trkParam.getZ() + params.TRDNonRCCorrDzDtgl * trkParam.getTgl();
    // correction for DVT, equivalent to a shift in X at which Y is evaluated: dY = tg_phi * dvt
    if (std::abs(params.TRDCorrDVT) > 1e-9f) {
      const double snp = trkParam.getSnp();
      posY += snp / std::sqrt((1. - snp) * (1. + snp)) * params.TRDCorrDVT;
    }
    std::array<float, 3> cov{};
    mRecoParam.recalcTrkltCov(static_cast<float>(tilt), static_cast<float>(trkParam.getSnp()), static_cast<float>(padLength), cov);

    auto& pnt = resTrack.info.emplace_back();
    pnt.lr = static_cast<int8_t>(o2::trd::Geometry::getLayer(trkltDet));
    pnt.label = Label(mDetIdx, trkltDet, true);
    pnt.x = static_cast<float>(traXYZ[0]);
    pnt.alpha = static_cast<float>(alpSens);
    pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(trkltDet), static_cast<float>(traXYZ[0]),
                                         static_cast<float>(posY), static_cast<float>(posZ),
                                         cov[0] + params.extraClsErrYTRD * params.extraClsErrYTRD,
                                         cov[2] + params.extraClsErrZTRD * params.extraClsErrZTRD, cov[1]);
  }

  if (static_cast<int>(resTrack.info.size() - nPointsIni) < params.minTRDTracklets) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  return true;
}

} // namespace o2::alignrs
