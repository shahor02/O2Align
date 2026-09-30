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
  // only the position of the mean vertex is alignable: being a point, it has no orientation
  auto dofs = std::make_unique<RigidBodyDOFSet>();
  dofs->setAllFree(false);
  dofs->setFree(RigidBodyDOFSet::TX, true);
  dofs->setFree(RigidBodyDOFSet::TY, true);
  dofs->setFree(RigidBodyDOFSet::TZ, true);
  vtx->setRigidBody(std::move(dofs));
  vtx->setSensorId(0);
  sensorMap[lbl] = vtx.get();
  mVertexVolume = vtx.get();
  return vtx;
}

void DetectorPVT::updatePositionLabels()
{
  mPositionLabels.clear();
  const auto* dofs = mVertexVolume ? mVertexVolume->getRigidBody() : nullptr;
  if (!dofs) {
    return;
  }
  // the free/fixed DOFs are those configured for the volume (built for the slot 0), but the labels
  // are those of the current calibration slot, which has its own set of global parameters
  const auto lbl = getVertexLabel(mSlotID);
  for (auto dof : {RigidBodyDOFSet::TX, RigidBodyDOFSet::TY, RigidBodyDOFSet::TZ}) {
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
  updatePositionLabels(); // the vertex of every slot is aligned via its own global parameters
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
    LOGP(info, "Ignoring the new CCDB MeanVertex within the time slot {}, its prior stays {}", mSlotID, mMeanVtxSlot.asString());
    return false;
  }
  mMeanVtxSlot = mMeanVtxCCDB; // w/o calibration slots the prior follows the CCDB object
  LOGP(info, "Mean vertex prior at timestamp {}: {}", tsMS, mMeanVtxSlot.asString());
  return true;
}

} // namespace o2::alignrs
