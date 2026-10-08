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

  /// this is not a real detector, it has no geometry and no measurements, hence it produces no AlignParam objects
  std::vector<o2::detectors::AlignParam> MP2AlignParams(const std::map<uint32_t, double>& labelToValue, bool writeLocal) const final { return {}; }

  /// label of the mean vertex volume. Its position is a calibration DOF set (MeanVertexDOFSet),
  /// whose labels carry the ID of the calibration slot, see Volume::getActiveCalibLabel
  static Label getVertexLabel() { return Label(DetPVT, 0, true); }
  Volume* getVertexVolume() const { return mVertexVolume; }

  const std::string& getTimeSlotsJson() const final;
  /// A CCDB update within an ongoing calibration slot is ignored, the slot being the unit of the
  /// calibration; w/o slots the prior follows the CCDB object. Reports whether the prior changed.
  bool setTimeStamp(long tsMS, bool ignoreMismath = true) final;

  /// hand over a MeanVertexObject delivered by the CCDB. It becomes the prior of the next
  /// calibration slot, or immediately the prior if no slots are defined (see setTimeStamp)
  void setMeanVertexCCDB(const o2::dataformats::MeanVertexObject& mv);
  const o2::dataformats::MeanVertexObject& getMeanVertexCCDB() const { return mMeanVtxCCDB; }
  /// prior of the primary vertex of every collision of the current calibration slot: the CCDB
  /// object valid at the start of the slot, frozen for its whole duration
  const o2::dataformats::MeanVertexObject& getMeanVertexPrior() const { return mMeanVtxSlot; }

  /// Millepede labels of the free DOFs (any subset of X, Y, Z, SlopeX, SlopeY) of the mean vertex in
  /// the current calibration slot, each independently free or fixed. Empty if none is free.
  const std::vector<int>& getFreeDOFLabels() const { return mFreeDOFLabels; }
  /// MeanVertexDOFSet::MeanVertexDOF index of every label in getFreeDOFLabels(), same order, used to
  /// pick the matching columns of the (Y,Z) vs (X,Y,Z,SlopeX,SlopeY) derivative matrix
  const std::vector<int>& getFreeDOFIndices() const { return mFreeDOFIndices; }
  /// (re)build the free-DOF labels/indices of the current slot; call once the DOF configuration is applied
  void updateFreeDOFs();

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;
  void onSlotChange(int slotID) final;

 private:
  Volume* mVertexVolume{nullptr};                   // the dummy volume of the mean vertex
  o2::dataformats::MeanVertexObject mMeanVtxCCDB{}; // last mean vertex object received from the CCDB
  o2::dataformats::MeanVertexObject mMeanVtxSlot{}; // prior of the current calibration slot
  bool mMeanVtxCCDBUpdated{false};                  // a CCDB object arrived and was not consumed yet
  std::vector<int> mFreeDOFLabels;                  // labels of the free DOFs in the current slot
  std::vector<int> mFreeDOFIndices;                 // MeanVertexDOF index of each label above, same order
};

} // namespace o2::alignrs

#endif
