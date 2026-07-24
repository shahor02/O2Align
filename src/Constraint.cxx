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
#include <fstream>
#include <sstream>
#include <cmath>

#include "Align/Constraint.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"

namespace o2::alignrs
{

void Constraint::write(std::ostream& os) const
{
  os << "!!! " << mName << '\n';
  os << "Constraint " << mValue << '\n';
  for (size_t i{0}; i < mLabels.size(); ++i) {
    os << mLabels[i] << " " << mCoeffs[i] << '\n';
  }
  os << '\n';
}

}

