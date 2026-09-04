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

#ifndef O2_ALIGN_DETECTORITS_H
#define O2_ALIGN_DETECTORITS_H

#include <vector>

#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::itsmft
{
class TopologyDictionary;
}

namespace o2::its3
{
class TopologyDictionary;
}

namespace o2::alignrs
{

class DetectorITS : public Detector
{
 public:
  using OVL = o2::itsmft::ChipMappingITS::Overlaps;
  enum EdgeFlags : int8_t {
    NONE = -1,
    LowRow = OVL::LowRow,
    HighRow = OVL::HighRow,
    Biased = 2
  };

  explicit DetectorITS(bool isITS3) : mIsITS3(isITS3) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final 
  {
    return mIsITS3 ? buildHierarchyIT3(sensorMap) : buildHierarchyITS(sensorMap);
  }

  void setTopologyDictionaries(const o2::itsmft::TopologyDictionary* itsDict, const o2::its3::TopologyDictionary* its3Dict)
  {
    mITSDict = itsDict;
    mIT3Dict = its3Dict;
  }
  const std::vector<FrameInfoExt>& getPointsInfo() const { return mITSPointsInfo; }

 private:
  Volume::Ptr buildHierarchyITS(Volume::SensorMapping& sensorMap);
  Volume::Ptr buildHierarchyIT3(Volume::SensorMapping& sensorMap);
 
  bool mIsITS3{false};
  const o2::itsmft::TopologyDictionary* mITSDict{nullptr};
  const o2::its3::TopologyDictionary* mIT3Dict{nullptr};
  std::vector<o2::itsmft::ChipMappingITS::Overlaps> mOverlaps;
  std::vector<int> mITSOvlCandidateID;
  std::vector<int> mITSOvlClusRef;
  std::vector<FrameInfoExt> mITSPointsInfo;
};

} // namespace o2::alignrs

#endif
