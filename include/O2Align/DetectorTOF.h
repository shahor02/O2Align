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

#ifndef O2_ALIGN_DETECTORTOF_H
#define O2_ALIGN_DETECTORTOF_H

#include <vector>

#include "O2Align/AlignmentTypes.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::alignrs
{

class DetectorTOF final : public Detector
{
 public:
  DetectorTOF() : Detector(DetTOF) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;

  const std::vector<FrameInfoExt>& getPointsInfo() const { return mTOFPointsInfo; }

 private:
  std::vector<FrameInfoExt> mTOFPointsInfo; // point per TOF cluster of the TF, lr<0 if unusable
  std::vector<Volume*> mSensors;            // strip volumes, indexed by continuous strip number
};

} // namespace o2::alignrs

#endif
