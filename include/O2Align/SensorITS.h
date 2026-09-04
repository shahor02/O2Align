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

#ifndef O2_ALIGN_SENSORITS_H
#define O2_ALIGN_SENSORITS_H

#include <vector>

#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Volume.h"

namespace o2::globaltracking
{
class RecoContainer;
}

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

Volume::Ptr buildHierarchyITS(Volume::SensorMapping& sensorMap);
Volume::Ptr buildHierarchyIT3(Volume::SensorMapping& sensorMap);

class SensorITS final : public Volume
{
 public:
  using Volume::Volume;
  using OVL = o2::itsmft::ChipMappingITS::Overlaps;
  enum EdgeFlags : int8_t {
    NONE = -1,
    LowRow = OVL::LowRow,
    HighRow = OVL::HighRow,
    Biased = 2
  };

  explicit SensorITS(bool isITS3) : Volume("alignment-data", 0, 0, false), mIsITS3(isITS3) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData);
  void setTopologyDictionaries(const o2::itsmft::TopologyDictionary* itsDict, const o2::its3::TopologyDictionary* its3Dict)
  {
    mITSDict = itsDict;
    mIT3Dict = its3Dict;
  }
  const std::vector<FrameInfoExt>& getPointsInfo() const { return mITSPointsInfo; }

  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
  void computeJacobianL2T(const double* pos, Matrix66& jac) const final;

 private:
  bool mIsITS3{false};
  const o2::itsmft::TopologyDictionary* mITSDict{nullptr};
  const o2::its3::TopologyDictionary* mIT3Dict{nullptr};
  std::vector<o2::itsmft::ChipMappingITS::Overlaps> mOverlaps;
  std::vector<int> mITSOvlCandidateID;
  std::vector<int> mITSOvlClusRef;
  std::vector<FrameInfoExt> mITSPointsInfo;
};

class SensorIT3 final : public Volume
{
  using Volume::Volume;
  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
  void computeJacobianL2T(const double* pos, Matrix66& jac) const final;
};

} // namespace o2::alignrs

#endif
