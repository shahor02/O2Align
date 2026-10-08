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
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "O2Align/TimeSlotsSet.h"
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
  enum DetIdx : int8_t {
    DetPVT = 0, // PVT is a virtual detector, used for the primary vertex
    DetITS,
    DetTPC,
    DetTRD,
    DetTOF,
    NDetectors
  };
  static constexpr const char* DetName[NDetectors] = {"PVT", "ITS", "TPC", "TRD", "TOF"};
  static constexpr int O2DetID[NDetectors] = {-1, o2::detectors::DetID::ITS, o2::detectors::DetID::TPC, o2::detectors::DetID::TRD, o2::detectors::DetID::TOF};
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
  /// o2::detectors::DetID of the detector, -1 for the virtual PVT
  int getO2DetID() const noexcept { return O2DetID[mDetIdx]; }
  const char* getDetName() const noexcept { return DetName[mDetIdx]; }

  /// \name Time-sliced calibration
  /// A detector may calibrate some of its DOFs independently in consecutive time intervals: the
  /// DOF set is configured once, on the volume built for the slot 0, while the emitted Millepede
  /// labels carry the ID of the slot the processed TF belongs to. A detector which does not use
  /// this leaves getTimeSlotsJson() empty and stays in the single slot 0.
  ///@{
  /// json file with the calibration intervals of this detector, empty if it has none. Detector
  /// specific, since it names its own Params field.
  virtual const std::string& getTimeSlotsJson() const;
  /// read getTimeSlotsJson() into mTimeSlots. Must be called before buildHierarchy, which may need
  /// the list of the slots to prepare the per-slot labels of its volumes.
  void loadTimeSlots();
  /// calibration intervals of this detector, nullptr if it has none
  const TimeSlotsSet* getTimeSlots() const noexcept { return mTimeSlots.get(); }
  /// ID of the calibration slot of the processed TF, 0 if the detector has no slots or no TF was
  /// processed yet
  int getSlotID() const noexcept { return mSlotID < 0 ? 0 : mSlotID; }
  /// Select the calibration slot covering the timestamp (in ms) and report whether anything changed
  /// for this detector, i.e. whether the caller must refresh what depends on it. Fatal if the slots
  /// are defined and the timestamp is covered by none of them unles ignoreMismatch, e.g. for MilleRes (results extraction) mode
  virtual bool setTimeStamp(long tsMS, bool ignoreMismatch = false);
  ///@}

  /// Convert the rigid-body corrections fitted by Millepede for the branch of this detector into
  /// the AlignParam objects to be applied to the IDEAL geometry (GeometryManager::applyAlignment),
  /// i.e. the fitted corrections combined with the initial alignment the fit started from.
  /// \param labelToValue fitted values of the free parameters, see Volume::readMillepedeResults
  /// \param initial      alignment of this detector only, as provided by the CCDB, whose
  ///                     application to the ideal geometry gave the geometry the fit was done on.
  ///                     Its objects may be in the global or in the local delta convention.
  /// \param writeLocal   write all output objects as local deltas (true) or as global ones (false),
  ///                     whatever the convention of the initial objects
  /// \return             updated alignment of this detector: every initial object rewritten (a new
  ///                     correction of a parent changes the global delta of all its daughters),
  ///                     plus the volumes which were not in the initial vector but are moved by the
  ///                     fit, ordered by the geometry level.
  /// Must be called after attachTo, with the geometry the fit was done on loaded; the ideal geometry
  /// itself is not needed. The derivation (thesis A.7, A.8) is documented in Detector.cxx.
  virtual std::vector<o2::detectors::AlignParam> MP2AlignParams(const std::map<uint32_t, double>& labelToValue, bool writeLocal) const;

  /// Compact, human-readable summary of how many volumes of this detector's branch have at least
  /// one free rigid-body and/or calibration DOF, one line per hierarchy level: "<Nrb>/<Ncal>/<N>"
  /// counts, with N the total number of volumes of that kind. The base implementation groups by
  /// the generic tree depth below getTopVolume(); a detector whose levels have a meaningful name
  /// (ITS: half-barrel/stave/half-stave/module/sensor) should override it. Empty before attachTo.
  virtual std::string reportDOFSummary() const;

 protected:
  /// create the stand-alone hierarchy of this detector, its top volume being its own root.
  /// Called by attachTo, which makes it a branch of the common hierarchy.
  virtual Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) = 0;

  /// hook invoked by setTimeStamp when the calibration slot changes, for the detector to re-point
  /// its slot-dependent labels and priors. It is invoked also for the 1st TF, whatever its slot, so
  /// that the slot-dependent state (e.g. the mean vertex prior) is always initialised by it.
  virtual void onSlotChange(int /*slotID*/) {}

  DetIdx mDetIdx{DetPVT};
  Volume* mTopVolume{nullptr};              // top volume of the detector, owned by the common hierarchy
  std::unique_ptr<TimeSlotsSet> mTimeSlots; // calibration intervals of this detector, if any
  int mSlotID{-1};                          // calibration slot of the processed TF, -1 before the 1st one  
};

} // namespace o2::alignrs

#endif
