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

#ifndef O2_ALIGN_COSTRAINT_H
#define O2_ALIGN_COSTRAINT_H

#include <Rtypes.h>
#include <vector>

namespace o2::alignrs
{

class Constraint
{
 public:
  Constraint(const std::string& name, float value) : mName(name), mValue(value) {}
  void add(uint32_t lab, float coeff)
  {
    mLabels.push_back(lab);
    mCoeffs.push_back(coeff);
  }
  auto& getName() const { return mName; }
  auto& getLabels() const { return mLabels; }
  auto& getCoeffs() const { return mCoeffs; }
  
  void write(std::ostream& os) const;
  auto getSize() const noexcept { return mLabels.size(); }

 private:
  std::string mName;             // name of the constraint
  float mValue{0.0f};            // constraint value
  std::vector<uint32_t> mLabels; // parameter labels
  std::vector<float> mCoeffs;    // their coefficients

  ClassDefNV(Constraint, 1);
};


};

#endif  // O2_ALIGN_COSTRAINT_H
