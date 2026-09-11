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

#ifndef O2_ALIGN_SENSORTRD_H
#define O2_ALIGN_SENSORTRD_H

#include "O2Align/Volume.h"

namespace o2::alignrs
{

/// A TRD "sensor" is a readout chamber. Its TGeo local frame has extra X,Y rotations wrt the
/// standard convention, so L2G is redefined to include them: the frame in which the alignment
/// parameters are defined is therefore the pseudo-local (rotated) chamber frame.
/// The sensor ID is the TRD chamber number [0 : MAXCHAMBER).
class SensorTRD final : public Volume
{
 public:
  using Volume::Volume;
  void defineMatrixL2G() final;
  void defineMatrixT2L() final;
};

} // namespace o2::alignrs

#endif
