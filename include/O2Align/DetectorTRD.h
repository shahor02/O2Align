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

#ifndef O2_ALIGN_DETECTORTRD_H
#define O2_ALIGN_DETECTORTRD_H

#include <vector>

#include "DataFormatsTRD/CalibratedTracklet.h"
#include "GPUTRDRecoParam.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::trd
{
class TrackletTransformer;
}

namespace o2::alignrs
{

class DetectorTRD final : public Detector
{
 public:
  DetectorTRD() : Detector(DetTRD) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;

  /// must be provided by the caller before the 1st prepareData call
  void setTransformer(const o2::trd::TrackletTransformer* tr) { mTransformer = tr; }

 private:
  const o2::trd::TrackletTransformer* mTransformer{nullptr};
  o2::gpu::GPUTRDRecoParam mRecoParam; // parameters required for the tracklet covariance
  bool mRecoParamInit{false};
  std::vector<o2::trd::CalibratedTracklet> mTrackletsLoc; // calibrated tracklets of the TF in the LOCAL frame
  std::vector<Volume*> mSensors;                          // chamber volumes, indexed by TRD chamber ID
};

} // namespace o2::alignrs

#endif
