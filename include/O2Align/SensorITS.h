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

#ifndef O2_ALIGN_SENSORITS_H
#define O2_ALIGN_SENSORITS_H

#include <map>
#include <cstdint>
#include <nlohmann/json_fwd.hpp>

#include "O2Align/Volume.h"

namespace o2::alignrs
{

/// ITS2 chip: the measurement plane is the effective sensitive layer, not the chip volume
class SensorITS final : public Volume
{
 public:
  using Volume::Volume;
  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
  /// writes the Legendre or Inextensional calibration block, moved here from
  /// Volume::writeMillepedeResults since the calibration DOFSet layout is detector-specific
  bool MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const final;
};

/// ITS3 tile: a segment of a cylinder whose local origin is the centre of the circle. The tracking
/// frame of every cluster is the one normal to the surface at the cluster, hence it varies over the
/// tile: getT2L() refers to the frame at the centre of the tile only.
class SensorIT3 final : public Volume
{
 public:
  SensorIT3(const char* symName, Label label, bool virt = false) : Volume(symName, label, virt) { mFixedTrackingFrame = false; }
  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
};

} // namespace o2::alignrs

#endif
