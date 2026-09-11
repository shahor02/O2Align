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
class Track;

using GTrackID = o2::dataformats::GlobalTrackID;
using GlobalIDSet = std::array<GTrackID, GTrackID::NSources>;

class Detector
{
 public:
  /// detector index, stored in the DET bits of the Millepede Label (only 4 values fit there!)
  enum DetIdx : uint32_t {
    DetITS = 0,
    DetTPC,
    DetTRD,
    DetTOF,
    NDetectors
  };
  static constexpr const char* DetName[NDetectors] = {"ITS", "TPC", "TRD", "TOF"};

  explicit Detector(DetIdx det) : mDetIdx(det) {}
  virtual ~Detector() = default;
  virtual void prepareData(o2::globaltracking::RecoContainer* recoData) = 0;
  virtual bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) = 0;
  virtual Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) = 0;

  DetIdx getDetIdx() const noexcept { return mDetIdx; }
  const char* getDetName() const noexcept { return DetName[mDetIdx]; }

 protected:
  DetIdx mDetIdx{DetITS};
};

} // namespace o2::alignrs

#endif
