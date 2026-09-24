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

#ifndef O2_ALIGN_DETECTOR_H
#define O2_ALIGN_DETECTOR_H

#include <array>
#include <cstdint>

#include "O2Align/Volume.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"


namespace o2::alignrs
{
struct Track;

using GTrackID = o2::dataformats::GlobalTrackID;
using GlobalIDSet = std::array<GTrackID, GTrackID::NSources>;

class Detector
{
 public:
  /// detector index, stored in the DET bits of the Millepede Label
  enum DetIdx : uint32_t {
    DetPVT = 0, // PVT is a virtual detector, used for the primary vertex
    DetITS,
    DetTPC,
    DetTRD,
    DetTOF,
    NDetectors
  };
  static constexpr const char* DetName[NDetectors] = {"PVT", "ITS", "TPC", "TRD", "TOF"};
  static_assert(NDetectors <= Label::DET_GLOBAL, "Detector index clashes with the code reserved for the root of the hierarchy");

  explicit Detector(DetIdx det) : mDetIdx(det) {}
  virtual ~Detector() = default;
  virtual void prepareData(o2::globaltracking::RecoContainer* recoData) = 0;
  virtual bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) = 0;

  /// build the hierarchy of this detector and graft it under the common root of all detectors.
  /// Returns the top volume of the detector, nullptr if the detector has no volumes at all.
  Volume* attachTo(Volume* root, Volume::SensorMapping& sensorMap);
  /// top volume of this detector within the common hierarchy, nullptr before attachTo
  Volume* getTopVolume() const noexcept { return mTopVolume; }

  DetIdx getDetIdx() const noexcept { return mDetIdx; }
  const char* getDetName() const noexcept { return DetName[mDetIdx]; }

 protected:
  /// create the stand-alone hierarchy of this detector, its top volume being its own root.
  /// Called by attachTo, which makes it a branch of the common hierarchy.
  virtual Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) = 0;

  DetIdx mDetIdx{DetPVT};
  Volume* mTopVolume{nullptr}; // top volume of the detector, owned by the common hierarchy
};

} // namespace o2::alignrs

#endif
