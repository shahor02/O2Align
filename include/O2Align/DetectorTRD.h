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
#include "DataFormatsTRD/CalVdriftExB.h"
#include "TRDBase/TrackletTransformer.h"
#include "GPUTRDRecoParam.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::alignrs
{

class DetectorTRD final : public Detector
{
 public:
  DetectorTRD() : Detector(DetTRD) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;
  int prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames) final;
  int getMinCosmicMeasurements() const final;
  void initCalib();
  void setCalVdriftExB(const o2::trd::CalVdriftExB* cal)
  {
    mCalVdriftExB = cal; // may arrive before or after initCalib()
    if (mTransformer) {
      mTransformer->setCalVdriftExB(cal);
    }
  }

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;

 private:
  bool trackletToTrackingFrame(const o2::trd::Tracklet64& trackletRaw, const o2::trd::CalibratedTracklet& trkltLoc, std::array<double, 3>& traXYZ, double& alpha) const;
  FrameInfoExt makeTrackletPoint(const o2::trd::Tracklet64& trackletRaw, const std::array<double, 3>& traXYZ, double alpSens, const o2::track::TrackParD& trkParam) const;

  std::unique_ptr<o2::trd::TrackletTransformer> mTransformer; // created at the 1st prepareData, once the geometry is loaded
  const o2::trd::CalVdriftExB* mCalVdriftExB{nullptr};
  o2::gpu::GPUTRDRecoParam mRecoParam; // parameters required for the tracklet covariance, set in prepareData
  float mRecoParamField{-999.f};       // field mRecoParam was initialised for
  std::vector<o2::trd::CalibratedTracklet> mTrackletsLoc; // calibrated tracklets of the TF in the LOCAL frame
  std::vector<Volume*> mSensors;                          // chamber volumes, indexed by TRD chamber ID
};

} // namespace o2::alignrs

#endif
