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
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"
#include <cstring>
#include "O2Align/SensorTRD.h"
#include "O2Align/TrackFit.h"
#include "TRDBase/Geometry.h"

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
  // w/o a geometry counterpart the envelope cannot be shifted, but it may own calibration DOFs
  root->setRigidBodyAllowed(false);
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

void DetectorTRD::initCalib()
{
  if (!mTransformer) {
    mTransformer = std::make_unique<o2::trd::TrackletTransformer>();
    if (Params::Instance().applyXORTRD) {
      mTransformer->setApplyXOR();
    }
    const char* prodType = std::getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
    if (prodType && std::strcmp(prodType, "MC") == 0) {
      // apply artificial pad shift in case non-ideal alignment is used to compensate for shift in current alignment from real data
      mTransformer->setApplyShift(false);
    }
    mTransformer->init();
    mTransformer->setCalVdriftExB(mCalVdriftExB);
  }
  mRecoParam.init(o2::base::PropagatorD::Instance()->getNominalBz());
}

void DetectorTRD::prepareData(o2::globaltracking::RecoContainer* recoData)
{
  const auto trackletsRaw = recoData->getTRDTracklets();
  mTrackletsLoc.clear();
  mTrackletsLoc.reserve(trackletsRaw.size());
  for (const auto& tracklet : trackletsRaw) {
    mTrackletsLoc.push_back(mTransformer->transformTracklet(tracklet, false));
  }
}

// The calibrated tracklet is in the LOCAL frame of its chamber, T2L of the (possibly pre-aligned) volume brings it to TRK.
// False if the chamber has no volume in the hierarchy.
bool DetectorTRD::trackletToTrackingFrame(const o2::trd::Tracklet64& trackletRaw, const o2::trd::CalibratedTracklet& trkltLoc, std::array<double, 3>& traXYZ, double& alpha) const
{
  const int trkltDet = trackletRaw.getDetector();
  auto* sensor = mSensors[trkltDet];
  if (!sensor) {
    return false;
  }
  const double locXYZ[3] = {trkltLoc.getX(), trkltLoc.getY(), trkltLoc.getZ()};
  sensor->getT2L().MasterToLocal(locXYZ, traXYZ.data());
  alpha = o2::math_utils::detail::sector2Angle<double>(o2::trd::Geometry::getSector(trkltDet));
  return true;
}

// Point of a tracklet given the track state at its X: the track provides the Z and the angles needed for the
// tilt correction and the covariance
FrameInfoExt DetectorTRD::makeTrackletPoint(const o2::trd::Tracklet64& trackletRaw, const std::array<double, 3>& traXYZ, double alpSens, const o2::track::TrackParD& trkParam) const
{
  const auto& params = Params::Instance();
  const auto* geo = o2::trd::Geometry::instance();
  const int trkltDet = trackletRaw.getDetector();
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

  FrameInfoExt pnt;
  pnt.detIdx = mDetIdx;
  pnt.lr = static_cast<int8_t>(o2::trd::Geometry::getLayer(trkltDet));
  pnt.label = Label(mDetIdx, trkltDet, true);
  pnt.x = static_cast<float>(traXYZ[0]);
  pnt.alpha = static_cast<float>(alpSens);
  pnt.zFromTrack = true; // posZ is the track Z at the tracklet: it does not depend on the chamber position
  pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(trkltDet), static_cast<float>(traXYZ[0]), static_cast<float>(posY), static_cast<float>(posZ),
                                       cov[0] + params.extraClsErrYTRD * params.extraClsErrYTRD,
                                       cov[2] + params.extraClsErrZTRD * params.extraClsErrZTRD, cov[1]);
  return pnt;
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
  const auto trackletsRaw = recoData->getTRDTracklets();
  const size_t nPointsIni = resTrack.info.size();

  // we continue the fit of the seed prepared by the inner detectors outward, w/o resetting its cov.matrix;
  // the fitted track also provides the angles and the Z needed for the tilt correction and the covariance
  auto trkParam = resTrack.track;
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trkParam);
  }
  float chi2 = 0.f;

  for (int il = 0; il < o2::trd::constants::NLAYER; ++il) {
    const int trkltId = trk.getTrackletIndex(il);
    if (trkltId < 0) {
      continue;
    }
    const auto& trackletRaw = trackletsRaw[trkltId];
    std::array<double, 3> traXYZ{};
    double alpSens = 0.;
    if (!trackletToTrackingFrame(trackletRaw, mTrackletsLoc[trkltId], traXYZ, alpSens)) {
      LOGP(error, "TRD chamber {} referred by a track has no alignable volume", trackletRaw.getDetector());
      resTrack.info.resize(nPointsIni);
      return false;
    }
    if (!prop->propagateToAlphaX(trkParam, refLin, alpSens, traXYZ[0], false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      resTrack.info.resize(nPointsIni);
      return false;
    }
    auto& pnt = resTrack.info.emplace_back(makeTrackletPoint(trackletRaw, traXYZ, alpSens, trkParam));
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pnt.cluster));
    if (!trkParam.update(pnt.cluster)) {
      resTrack.info.resize(nPointsIni);
      return false;
    }
    if (refLin) { // displace the reference to the last updated tracklet
      refLin->setY(pnt.cluster.getY());
      refLin->setZ(pnt.cluster.getZ());
    }
  }
  int nAdded = static_cast<int>(resTrack.info.size() - nPointsIni);
  if (params.verbose > 0) {
    LOGP(info, "Prepare {}: NTrackletsIni={} -> Npoints={}, chi2 = {}", gid.asString(), trk.getNtracklets(), nAdded, chi2);
  }
  if (nAdded < params.minTRDTracklets) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  resTrack.track = trkParam; // the seed is replaced only by a successful fit
  resTrack.kfFit.chi2 += chi2;
  return true;
}

// Per layer the tracklet of the leg (the one matched to its track and those found on the road) closest to the reference track is taken
int DetectorTRD::prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames)
{
  const auto& params = Params::Instance();
  std::array<FrameInfoExt, o2::trd::constants::NLAYER> best{};
  std::array<double, o2::trd::constants::NLAYER> bestChi2{};
  bestChi2.fill(params.cosmMaxChi2Cl);
  for (const auto& trkl : cosm.trdTracklets) {
    if (trkl.leg != leg) {
      continue;
    }
    const o2::trd::Tracklet64 trackletRaw(trkl.word);
    const auto trkltLoc = mTransformer->transformTracklet(trackletRaw, false);
    std::array<double, 3> traXYZ{};
    double alpSens = 0.;
    if (!trackletToTrackingFrame(trackletRaw, trkltLoc, traXYZ, alpSens)) {
      continue;
    }
    o2::track::TrackParCovD trk;
    if (!propagateCosmicRef(ref, alpSens, traXYZ[0], trk)) {
      continue;
    }
    const auto pnt = makeTrackletPoint(trackletRaw, traXYZ, alpSens, trk);
    const double chi2 = trk.getPredictedChi2Quiet(pnt.cluster);
    if (chi2 < bestChi2[pnt.lr]) {
      bestChi2[pnt.lr] = chi2;
      best[pnt.lr] = pnt;
    }
  }
  int nAdded = 0;
  for (const auto& pnt : best) {
    if (pnt.isValid()) {
      frames.push_back(pnt);
      ++nAdded;
    }
  }
  return nAdded;
}

int DetectorTRD::getMinCosmicMeasurements() const
{
  return Params::Instance().minTRDTracklets;
}

} // namespace o2::alignrs
