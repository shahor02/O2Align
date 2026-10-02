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

#include <memory>

#include "Framework/Logger.h"
#include "O2Align/DetectorPVT.h"
#include "O2Align/Params.h"

namespace o2::alignrs
{

Volume::Ptr DetectorPVT::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  const auto lbl = getVertexLabel();
  // The mean vertex is the only volume of this detector, hence it is directly the top volume of its
  // branch of the hierarchy, with no envelope above it. It has no alignable entry in the geometry
  // and carries no geometry of its own: the frame in which the vertex is measured is defined by each
  // track separately and the prior position is used directly by the fit, not via a matrix.
  auto vtx = std::make_unique<Volume>("PVT/meanVertex", lbl, true);
  // being a point w/o orientation, only its position is alignable. It changes with time and is
  // therefore a calibration (there are no time-sliced rigid-body DOFs): the position of every slot
  // is an independent set of global parameters. A calib rule of the DOF configuration may replace
  // this default, all-free set, e.g. to fix the position.
  vtx->setRigidBodyAllowed(false);
  vtx->setCalib(std::make_unique<MeanVertexDOFSet>());
  if (mTimeSlots) {
    vtx->setCalibSlots(*mTimeSlots);
  }
  vtx->setSensorId(0);
  sensorMap[lbl] = vtx.get();
  mVertexVolume = vtx.get();
  return vtx;
}

void DetectorPVT::updatePositionLabels()
{
  mPositionLabels.clear();
  const auto* dofs = mVertexVolume ? mVertexVolume->getCalib() : nullptr;
  if (!dofs) {
    return;
  }
  if (dofs->type() != DOFSet::Type::MeanVertex) {
    LOGP(fatal, "The calibration DOFs of {} must be of the meanvertex type, check the DOF configuration", mVertexVolume->getSymName());
  }
  // the labels of the calibration slot of the processed TF
  const auto& lbl = mVertexVolume->getActiveCalibLabel();
  for (auto dof : {MeanVertexDOFSet::X, MeanVertexDOFSet::Y, MeanVertexDOFSet::Z}) {
    if (!dofs->isFree(dof)) { // a fixed position imposes the prior w/o being fitted
      mPositionLabels.clear();
      return;
    }
    mPositionLabels.push_back(lbl.rawGBL(dof));
  }
}

const std::string& DetectorPVT::getTimeSlotsJson() const
{
  return Params::Instance().MVTimeSlotsJson;
}

void DetectorPVT::setMeanVertexCCDB(const o2::dataformats::MeanVertexObject& mv)
{
  mMeanVtxCCDB = mv;
  mMeanVtxCCDBUpdated = true;
  LOGP(info, "New CCDB MeanVertex: {}", mMeanVtxCCDB.asString());
}

void DetectorPVT::onSlotChange(int slotID)
{
  // the slot is the unit of the mean vertex calibration: its prior is the CCDB object valid at the
  // start of the slot, an eventual later CCDB update within the same slot is ignored
  mMeanVtxSlot = mMeanVtxCCDB;
  mMeanVtxCCDBUpdated = false;
  if (mVertexVolume) { // the vertex of every slot is aligned via its own global parameters
    mVertexVolume->setActiveCalibSlot(slotID);
    updatePositionLabels();
  }
  LOGP(info, "Mean vertex prior for time slot {}: {}", slotID, mMeanVtxSlot.asString());
}

bool DetectorPVT::setTimeStamp(long tsMS)
{
  if (Detector::setTimeStamp(tsMS)) { // onSlotChange has refreshed the prior and the labels
    return true;
  }
  if (!mMeanVtxCCDBUpdated) {
    return false;
  }
  mMeanVtxCCDBUpdated = false;
  if (mTimeSlots) { // the prior of the ongoing slot is kept
    LOGP(info, "Ignoring the new CCDB MeanVertex within the time slot {}, its prior stays {}", getSlotID(), mMeanVtxSlot.asString());
    return false;
  }
  mMeanVtxSlot = mMeanVtxCCDB; // w/o calibration slots the prior follows the CCDB object
  LOGP(info, "Mean vertex prior at timestamp {}: {}", tsMS, mMeanVtxSlot.asString());
  return true;
}

} // namespace o2::alignrs
