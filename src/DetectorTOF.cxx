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

void DetectorTOF::prepareData(o2::globaltracking::RecoContainer* recoData)
{
  // The TOF cluster position is fully defined by its main contributing channel, hence the
  // conversion to the tracking frame does not depend on the track and is done once per TF.
  const auto& params = Params::Instance();
  const auto clusters = recoData->getTOFClusters();
  mTOFPointsInfo.clear();
  mTOFPointsInfo.resize(clusters.size()); // unusable clusters are left with the default lr = -1

  const float extraErrY2 = params.extraClsErrYTOF * params.extraClsErrYTOF;
  const float extraErrZ2 = params.extraClsErrZTOF * params.extraClsErrZTOF;

  for (size_t ic = 0; ic < clusters.size(); ++ic) {
    const auto& clus = clusters[ic];
    int detInd[5] = {};
    o2::tof::Geo::getVolumeIndices(clus.getMainContributingChannel(), detInd);
    const int sid = o2::tof::Geo::getStripNumberPerSM(detInd[1], detInd[2]) + clus.getSector() * o2::tof::Geo::NSTRIPXSECTOR;
    if (sid < 0 || sid >= static_cast<int>(mSensors.size()) || !mSensors[sid]) {
      LOGP(debug, "TOF cluster {} refers to strip {} which has no alignable volume", ic, sid);
      continue;
    }
    const auto* sensor = mSensors[sid];
    const double loc[3] = {(detInd[4] + 0.5) * o2::tof::Geo::XPAD - o2::tof::Geo::XHALFSTRIP, 0., (detInd[3] - 0.5) * o2::tof::Geo::ZPAD};
    double tra[3]{};
    sensor->getT2L().MasterToLocal(loc, tra); // LOC -> TRK

    auto& pnt = mTOFPointsInfo[ic];
    pnt.lr = 0;
    pnt.label = Label(mDetIdx, sid, true);
    pnt.x = static_cast<float>(tra[0]);
    pnt.alpha = static_cast<float>(o2::math_utils::detail::sector2Angle<double>(clus.getSector()));
    pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(sid), static_cast<float>(tra[0]), static_cast<float>(tra[1]), static_cast<float>(tra[2]),
                                         clus.getSigmaY2() + extraErrY2, clus.getSigmaZ2() + extraErrZ2, clus.getSigmaYZ());
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
  if (pnt.lr < 0) { // the strip of this cluster is not in the hierarchy
    return false;
  }
  resTrack.info.push_back(pnt);
  return true;
}

} // namespace o2::alignrs
