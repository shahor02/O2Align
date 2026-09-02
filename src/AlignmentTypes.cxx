// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <string>
#include <format>

#include "O2Align/AlignmentTypes.h"

namespace o2::alignrs
{

std::string FrameInfoExt::asString() const
{
  return std::format("Sensor={} Layer={} X={} Alpha={}\n\tMEAS: y={} z={}", sens, lr, x, alpha, positionTrackingFrame[0], positionTrackingFrame[1]);
}

// RSTODO temporary here
std::vector<double> legendrePols(int order, double x)
{
  std::vector<double> p(order + 1);
  p[0] = 1.;
  if (order > 0) {
    p[1] = x;
  }
  for (int n = 1; n < order; ++n) {
    p[n + 1] = ((2 * n + 1) * x * p[n] - n * p[n - 1]) / (n + 1);
  }
  return p;
}

} // namespace o2::alignrs
