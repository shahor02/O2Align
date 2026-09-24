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

std::vector<int> DetectorPVT::getPositionLabels() const
{
  std::vector<int> labels;
  if (!mVertexVolume) {
    return labels;
  }
  const auto* dofs = mVertexVolume->getRigidBody();
  if (!dofs) {
    return labels;
  }
  // the free/fixed DOFs are those configured for the volume (built for the slot 0), but the labels
  // are those of the current calibration slot, which has its own set of global parameters
  const auto lbl = getVertexLabel(mMVSlotID);
  for (auto dof : {RigidBodyDOFSet::TX, RigidBodyDOFSet::TY, RigidBodyDOFSet::TZ}) {
    if (!dofs->isFree(dof)) { // a fixed position imposes the prior w/o being fitted
      return {};
    }
    labels.push_back(lbl.rawGBL(dof));
  }
  return labels;
}

} // namespace o2::alignrs
