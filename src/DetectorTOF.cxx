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

#include <format>
#include <memory>

#include <TGeoManager.h>

#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsTOF/Cluster.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "O2Align/DetectorTOF.h"
#include "O2Align/Params.h"
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"
#include "DetectorsBase/Propagator.h"
#include "O2Align/SensorTOF.h"
#include "TOFBase/Geo.h"

namespace o2::alignrs
{

Volume::Ptr DetectorTOF::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  uint32_t gLbl{0};
  const uint32_t det = mDetIdx;
  // the TOF envelope has no alignable entry in the geometry
  auto root = std::make_unique<Volume>("TOF_envelope", gLbl++, det, false, true);
  // w/o a geometry counterpart the envelope cannot be shifted, but it may own calibration DOFs
  root->setRigidBodyAllowed(false);
  mSensors.assign(o2::tof::Geo::NSECTORS * o2::tof::Geo::NSTRIPXSECTOR, nullptr);

  int cnt = 0;
  for (int isect = 0; isect < o2::tof::Geo::NSECTORS; ++isect) {
    const auto symSector = std::format("TOF/sm{:02d}", isect);
    Volume* volSector = nullptr;
    for (int istr = 1; istr <= o2::tof::Geo::NSTRIPXSECTOR; ++istr) {
      const int sid = cnt++; // continuous strip number, must follow the numbering used at the geometry creation
      const auto symStrip = std::format("TOF/sm{:02d}/strip{:02d}", isect, istr);
      if (!gGeoManager->GetAlignableEntry(symStrip.c_str())) { // missing strip (disabled SM or TOF hole)
        continue;
      }
      if (!volSector) { // create the supermodule volume only once, and only if it has strips
        if (!gGeoManager->GetAlignableEntry(symSector.c_str())) {
          LOGP(fatal, "TOF supermodule {} has strips but no alignable entry", symSector);
        }
        volSector = root->addChild(symSector.c_str(), gLbl++, det, false);
      }
      const Label lbl(det, sid, true);
      auto* strip = volSector->addChild<SensorTOF>(symStrip.c_str(), lbl);
      strip->setSensorId(sid);
      sensorMap[lbl] = strip;
      mSensors[sid] = strip;
    }
  }
  return root;
}

// Point of a TOF cluster with the given main contributing channel and errors, in the tracking frame of its strip;
// invalid if the strip is not in the hierarchy. The position is fully defined by the channel.
FrameInfoExt DetectorTOF::makeClusterPoint(int channel, float sigY2, float sigZ2, float sigYZ) const
{
  const auto& params = Params::Instance();
  FrameInfoExt pnt;
  int detInd[5] = {};
  o2::tof::Geo::getVolumeIndices(channel, detInd);
  const int sector = channel / o2::tof::Geo::NPADSXSECTOR;
  const int sid = o2::tof::Geo::getStripNumberPerSM(detInd[1], detInd[2]) + sector * o2::tof::Geo::NSTRIPXSECTOR;
  if (sid < 0 || sid >= static_cast<int>(mSensors.size()) || !mSensors[sid]) {
    LOGP(debug, "TOF channel {} refers to strip {} which has no alignable volume", channel, sid);
    return pnt;
  }
  const auto* sensor = mSensors[sid];
  const double loc[3] = {(detInd[4] + 0.5) * o2::tof::Geo::XPAD - o2::tof::Geo::XHALFSTRIP, 0., (detInd[3] - 0.5) * o2::tof::Geo::ZPAD};
  double tra[3]{};
  sensor->getT2L().MasterToLocal(loc, tra); // LOC -> TRK

  pnt.lr = 0;
  pnt.detIdx = mDetIdx;
  pnt.label = Label(mDetIdx, sid, true);
  pnt.x = static_cast<float>(tra[0]);
  pnt.alpha = static_cast<float>(o2::math_utils::detail::sector2Angle<double>(sector));
  pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(sid), static_cast<float>(tra[0]), static_cast<float>(tra[1]), static_cast<float>(tra[2]),
                                       sigY2 + params.extraClsErrYTOF * params.extraClsErrYTOF, sigZ2 + params.extraClsErrZTOF * params.extraClsErrZTOF, sigYZ);
  return pnt;
}

void DetectorTOF::prepareData(o2::globaltracking::RecoContainer* recoData)
{
  // The TOF cluster position is fully defined by its main contributing channel, hence the
  // conversion to the tracking frame does not depend on the track and is done once per TF.
  const auto clusters = recoData->getTOFClusters();
  mTOFPointsInfo.clear();
  mTOFPointsInfo.resize(clusters.size()); // unusable clusters are left with the default lr = Invalid
  for (size_t ic = 0; ic < clusters.size(); ++ic) {
    const auto& clus = clusters[ic];
    mTOFPointsInfo[ic] = makeClusterPoint(clus.getMainContributingChannel(), clus.getSigmaY2(), clus.getSigmaZ2(), clus.getSigmaYZ());
  }
}

bool DetectorTOF::prepareTrack(o2::globaltracking::RecoContainer* /*recoData*/, const GlobalIDSet& ids, Track& resTrack)
{
  const auto gid = ids[GTrackID::TOF];
  if (!gid.isIndexSet()) {
    return false;
  }
  const int idx = gid.getIndex();
  if (idx < 0 || idx >= static_cast<int>(mTOFPointsInfo.size())) {
    LOGP(error, "TOF cluster index {} is out of the range of {} points prepared for this TF", idx, mTOFPointsInfo.size());
    return false;
  }
  const auto& pnt = mTOFPointsInfo[idx];
  if (pnt.lr == FrameInfoExt::Invalid) { // the strip of this cluster is not in the hierarchy
    return false;
  }
  const size_t nPointsIni = resTrack.info.size();
  resTrack.info.push_back(pnt);
  // continue the fit of the seed prepared by the inner detectors outward, w/o resetting its cov.matrix
  return resTrack.continueFitOutward(static_cast<int>(nPointsIni));
}

// The raw TOF clusters of the cosmic do not carry the errors: those of a pad are used. Among the clusters of the leg (the one matched to
// its track and those found on the road) the one closest to the reference track is taken.
int DetectorTOF::prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames)
{
  const auto& params = Params::Instance();
  constexpr float SigY2 = o2::tof::Geo::XPAD * o2::tof::Geo::XPAD / 12.f;
  constexpr float SigZ2 = o2::tof::Geo::ZPAD * o2::tof::Geo::ZPAD / 12.f;
  auto prop = o2::base::PropagatorD::Instance();
  FrameInfoExt bestPnt;
  bool found = false;
  double bestChi2 = params.cosmMaxChi2Cl;
  for (const auto& cl : cosm.clTOF) {
    if (cl.leg != leg) {
      continue;
    }
    auto pnt = makeClusterPoint(cl.channel, SigY2, SigZ2, 0.f);
    if (!pnt.isValid()) {
      continue;
    }
    auto trk = ref;
    if (!prop->propagateToAlphaX(trk, nullptr, pnt.alpha, pnt.x, false, params.maxSnp, params.maxStep, 1, o2::base::PropagatorD::MatCorrType::USEMatCorrNONE)) {
      continue;
    }
    const double chi2 = trk.getPredictedChi2Quiet(pnt.cluster);
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestPnt = pnt;
      found = true;
    }
  }
  if (!found) {
    return 0;
  }
  frames.push_back(bestPnt);
  return 1;
}

} // namespace o2::alignrs
