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

#include "Framework/Logger.h"
#include "O2Align/Detector.h"

namespace o2::alignrs
{

Volume* Detector::attachTo(Volume* root, Volume::SensorMapping& sensorMap)
{
  if (root == nullptr) {
    LOGP(fatal, "Cannot attach the hierarchy of {} to a null root", getDetName());
  }
  auto top = buildHierarchy(sensorMap);
  if (!top) {
    LOGP(fatal, "Detector {} provided no hierarchy", getDetName());
  }
  // a branch which provides no measurement at all is not attached: this happens when the geometry
  // has no volume of this detector (all its supermodules disabled), leaving an envelope alone.
  // A detector consisting of a single sensor (the mean vertex) is perfectly legitimate: it needs no
  // envelope above it, its only volume being directly the top of its branch.
  bool hasSensor = false;
  top->traverse([&hasSensor](Volume* vol) { hasSensor = hasSensor || vol->getLabel().sens(); });
  if (!hasSensor) {
    LOGP(warn, "Detector {} has no sensor volume, its branch {} is discarded", getDetName(), top->getSymName());
    return nullptr;
  }
  mTopVolume = root->adoptChild(std::move(top));
  LOGP(info, "Attached the branch {} of {} to {}", mTopVolume->getSymName(), getDetName(), root->getSymName());
  return mTopVolume;
}

// The KF (re)fit of the prepared track moved to Track::fitTrack / Track::continueFitOutward

} // namespace o2::alignrs
