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

#include "DataFormatsITSMFT/CompCluster.h"
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

  explicit DetectorITS(bool isITS3) : Detector(DetITS), mIsITS3(isITS3) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;
  int prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames) final;
  int getMinCosmicMeasurements() const final;

  void setTopologyDictionaries(const o2::itsmft::TopologyDictionary* itsDict, const o2::its3::TopologyDictionary* its3Dict)
  {
    if (itsDict) {
      mITSDict = itsDict;
    }
    if (its3Dict) {
      mIT3Dict = its3Dict;
    }
  }

  /// Independent KF refit of the ITS track gidITS with its clusters prepared by prepareData, starting
  /// from the outer parameters with the covariance reset. `track` keeps its PID. Used for the PV refit.
  bool refitITSTrack(o2::globaltracking::RecoContainer* recoData, GTrackID gidITS, o2::track::TrackParCov& track) const;

  const std::vector<FrameInfoExt>& getPointsInfo() const { return mITSPointsInfo; }

  std::string reportDOFSummary() const final;

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;

 private:
  template <typename PattIt>
  FrameInfoExt makeClusterPoint(const o2::itsmft::CompClusterExt& cls, PattIt& pattIt, bool findOverlaps) const;
  bool isOverlapCandidate(const o2::BaseCluster<float>& cl, const o2::BaseCluster<float>& oClus, int side) const;
  bool isOverlapPair(const FrameInfoExt& a, const FrameInfoExt& b) const;

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
