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

#ifndef O2_ALIGN_SENSORTPC_H
#define O2_ALIGN_SENSORTPC_H

#include "O2Align/Volume.h"

namespace o2::alignrs
{

/// A TPC "sensor" is a readout sector. It has no counterpart in the geometry, hence the volume
/// is always created as virtual and both its frames are derived from the sector angle alone:
/// the local frame is the global one rotated by the sector angle and the tracking frame coincides
/// with the local one. The sensor ID is the sector number [0 : MAXSECTOR).
class SensorTPC final : public Volume
{
 public:
  using Volume::Volume;
  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
};

} // namespace o2::alignrs

#endif
