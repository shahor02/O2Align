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

#ifndef O2_ALIGN_DETECTORPVT_H
#define O2_ALIGN_DETECTORPVT_H

#include <string>
#include <vector>

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::alignrs
{

/// Virtual detector of the mean interaction point. It provides no measurements and no geometry of
/// its own: its single dummy volume only declares the position of the mean vertex as a global
/// alignment parameter. The prior value of this position, as well as the sigmas of the luminous
/// region constraining the vertex of every collision to it, are those of the MeanVertexObject and
/// are used directly by the track fit.
/// The detector owns the whole mean-vertex calibration state: the prior of the current calibration
/// slot, the last object delivered by the CCDB and the Millepede labels of the slot being aligned.
class DetectorPVT final : public Detector
{
 public:
  DetectorPVT();

  /// the mean vertex is not extracted from the TF data
  void prepareData(o2::globaltracking::RecoContainer* /*recoData*/) final {}
  /// the mean vertex contributes no measured point to an individual track
  bool prepareTrack(o2::globaltracking::RecoContainer* /*recoData*/, const GlobalIDSet& /*ids*/, Track& /*resTrack*/) final { return false; }

  /// label of the mean vertex volume. Its position is a calibration DOF set (MeanVertexDOFSet),
  /// whose labels carry the ID of the calibration slot, see Volume::getActiveCalibLabel
  static Label getVertexLabel() { return Label(DetPVT, 0, true); }
  Volume* getVertexVolume() const { return mVertexVolume; }

  const std::string& getTimeSlotsJson() const final;
  /// A CCDB update within an ongoing calibration slot is ignored, the slot being the unit of the
  /// calibration; w/o slots the prior follows the CCDB object. Reports whether the prior changed.
  bool setTimeStamp(long tsMS) final;

  /// hand over a MeanVertexObject delivered by the CCDB. It becomes the prior of the next
  /// calibration slot, or immediately the prior if no slots are defined (see setTimeStamp)
  void setMeanVertexCCDB(const o2::dataformats::MeanVertexObject& mv);
  const o2::dataformats::MeanVertexObject& getMeanVertexCCDB() const { return mMeanVtxCCDB; }
  /// prior of the primary vertex of every collision of the current calibration slot: the CCDB
  /// object valid at the start of the slot, frozen for its whole duration
  const o2::dataformats::MeanVertexObject& getMeanVertexPrior() const { return mMeanVtxSlot; }

  /// Millepede labels of the mean vertex position in the current calibration slot, ordered as
  /// (X, Y, Z) in the global frame. Empty if the position is not a free parameter.
  const std::vector<int>& getPositionLabels() const { return mPositionLabels; }
  /// (re)build the labels of the current slot; call once the DOF configuration is applied
  void updatePositionLabels();

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;
  void onSlotChange(int slotID) final;

 private:
  Volume* mVertexVolume{nullptr};                   // the dummy volume of the mean vertex
  o2::dataformats::MeanVertexObject mMeanVtxCCDB{}; // last mean vertex object received from the CCDB
  o2::dataformats::MeanVertexObject mMeanVtxSlot{}; // prior of the current calibration slot
  bool mMeanVtxCCDBUpdated{false};                  // a CCDB object arrived and was not consumed yet
  std::vector<int> mPositionLabels;                 // labels of the vertex position in the current slot
};

} // namespace o2::alignrs

#endif
