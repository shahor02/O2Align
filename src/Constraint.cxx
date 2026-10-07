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


#include <format>
#include <ostream>

#include "O2Align/Constraint.h"

namespace o2::alignrs
{

void Constraint::write(std::ostream& os) const
{
  os << "!!! " << mName << '\n';
  if (isMeasurement()) {
    os << std::format("Measurement {:.12g} {:.12g}\n", mValue, mSigma);
  } else {
    os << std::format("Constraint {:.12g}\n", mValue);
  }
  for (size_t i{0}; i < mLabels.size(); ++i) {
    // full precision: the coefficients mix rotations with lever arms of up to a few metres
    os << std::format("{} {:.12g}\n", mLabels[i], mCoeffs[i]);
  }
  os << '\n';
}

} // namespace o2::alignrs
