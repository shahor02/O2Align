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

#include <vector>

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
class DetectorPVT final : public Detector
{
 public:
  DetectorPVT() : Detector(DetPVT) {}

  /// the mean vertex is not extracted from the TF data
  void prepareData(o2::globaltracking::RecoContainer* /*recoData*/) final {}
  /// the mean vertex contributes no measured point to an individual track
  bool prepareTrack(o2::globaltracking::RecoContainer* /*recoData*/, const GlobalIDSet& /*ids*/, Track& /*resTrack*/) final { return false; }
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;

  /// label of the mean vertex volume, its DOFs TX/TY/TZ are the global parameters of its position
  static Label getVertexLabel() { return Label(DetPVT, 0, true); }
  Volume* getVertexVolume() const { return mVertexVolume; }

  /// Millepede labels of the mean vertex position, ordered as (X, Y, Z) in the global frame.
  /// Empty if the position is not a free parameter.
  std::vector<int> getPositionLabels() const;

  int getMVSlotID() const { return mMVSlotID; }
  void setMVSlotID(int slotID) { mMVSlotID = slotID; }

 private:
  Volume* mVertexVolume{nullptr}; // the dummy volume of the mean vertex
  int mMVSlotID = 0; // the intervalID of the mean vertex calibration slot, 0 if no slots are defined
};

} // namespace o2::alignrs

#endif
