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

#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TMath.h>

#include "MathUtils/Utils.h"
#include "O2Align/SensorTOF.h"
#include "TOFBase/Geo.h"

namespace o2::alignrs
{

void SensorTOF::defineMatrixL2G()
{
  mL2G = *mPN->GetMatrix();
}

void SensorTOF::defineMatrixT2L()
{
  setT2LFromAlpha(o2::math_utils::detail::sector2Angle<double>(getSensorId() / o2::tof::Geo::NSTRIPXSECTOR));
}

} // namespace o2::alignrs
