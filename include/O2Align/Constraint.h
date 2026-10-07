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


#ifndef O2_ALIGN_CONSTRAINT_H
#define O2_ALIGN_CONSTRAINT_H

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#include <Rtypes.h>

namespace o2::alignrs
{

/// linear equality constraint sum_i c_i * p_i = value of the Millepede steering file or, if a
/// sigma > 0 is given, the "Measurement" sum_i c_i * p_i = value +- sigma, i.e. a soft constraint
class Constraint
{
 public:
  Constraint(const std::string& name, double value, double sigma = -1.) : mName(name), mValue(value), mSigma(sigma) {}
  void add(uint32_t lab, double coeff)
  {
    mLabels.push_back(lab);
    mCoeffs.push_back(coeff);
  }
  auto& getName() const { return mName; }
  auto& getLabels() const { return mLabels; }
  auto& getCoeffs() const { return mCoeffs; }

  void write(std::ostream& os) const;
  auto getSize() const noexcept { return mLabels.size(); }
  bool isMeasurement() const noexcept { return mSigma > 0.; }

 private:
  std::string mName;             // name of the constraint
  double mValue{0.};             // constraint value
  double mSigma{-1.};            // its uncertainty if it is a measurement, otherwise <= 0 (exact constraint)
  std::vector<uint32_t> mLabels; // parameter labels
  std::vector<double> mCoeffs;   // their coefficients

  ClassDefNV(Constraint, 3);
};

} // namespace o2::alignrs

#endif // O2_ALIGN_CONSTRAINT_H
