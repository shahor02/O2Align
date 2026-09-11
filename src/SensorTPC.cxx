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
#include <TMath.h>

#include "DataFormatsTPC/Constants.h"
#include "MathUtils/Utils.h"
#include "O2Align/SensorTPC.h"

namespace o2::alignrs
{

void SensorTPC::defineMatrixL2G()
{
  // there is no TGeo volume for the sector: the local frame is the global one rotated to the sector
  constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
  const double alp = o2::math_utils::detail::sector2Angle<double>(getSensorId() % NSectorsPerSide);
  TGeoHMatrix rotZ; // identity
  rotZ.RotateZ(alp * TMath::RadToDeg());
  mL2G = rotZ;
}

void SensorTPC::defineMatrixT2L()
{
  mT2L = TGeoHMatrix(); // the tracking and local frames of a TPC sector coincide
}

} // namespace o2::alignrs
