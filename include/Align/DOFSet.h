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

#ifndef O2_ALIGN_DOFSET_H
#define O2_ALIGN_DOFSET_H

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace o2::alignrs
{

struct DerivativeContext {
  int sensorID{-1};
  int layerID{-1};
  float measX{0.};
  float measAlpha{0.};
  float measZ{0.};
  float trkY{0.};
  float trkZ{0.};
  float snp{0.};
  float tgl{0.};
  float dydx{0.};
  float dzdx{0.};
  
  ClassDefNV(DerivativeContext);
};

// Generic set of DOF
class DOFSet
{
 public:
  enum class Type : uint8_t {
    RigidBody,
    Legendre,
    Inextensional
  };
  
  virtual ~DOFSet() = default;
  virtual Type getType() const = 0;
  virtual std::string getDOFName(int idx) const = 0;
  virtual void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const = 0;
  
  int getNDOFs() const { return static_cast<int>(mFree.size()); }
  bool isFree(int idx) const { return mFree[idx]; }
  void setFreeStatus(int idx, bool f) { mFree[idx] = f; }
  void setAllFreeStatus(bool f) { std::fill(mFree.begin(), mFree.end(), f); }
  int getNFreeDOFs() const {
    int n = 0;
    for (bool f : mFree) {
      n += f;
    }
    return n;
  }

  void setGlobalsPointers(float* pars, float* errs) {
    gParVals = pars;
    gParErrs = errs;
  }

  void setFirstEntry(int i) { mFirstEntry = i; }
  int getFirstEntry() const { return mFirstEntry; }

  float getParVal(int par) const { return getParVals()[par]; }
  float getParErr(int par) const { return getParErrs()[par]; }

  void setParVal(int par, double v = 0) { getParVals()[par] = v; }
  void setParErr(int par, double e = 0) { getParErrs()[par] = e; }

protected:  
  DOFSet(int n) : mFree(n, true) {}
  float* getParVals() { return gParVals + mFirstEntry; }
  float* getParErrs() { return gParErrs + mFirstEntry; }
  void validateDerivativeOutput(const DOFSet& dofSet, Eigen::Ref<Eigen::MatrixXd> out);

  std::vector<bool> mFree;         // status of each DOF
  int mFirstEntry = -1;            // ID of the 1st parameter in the global results array

  static float* gParVals = nullptr; // start of global parameters array
  static float* gParErrs = nullptr; // start of global parameters errors array
  
  ClassDef(DOFSet,1);
};

// Rigid body set (rotations and offset)
class RigidBodyDOFSet final : public DOFSet
{
 public:
  enum RigidBodyDOF : uint8_t { TX = 0, TY, TZ, RX, RY, RZ, NDO };   // indices for rigid body parameters in LOC frame
  static constexpr const char* RigidBodyDOFNames[RigidBodyDOF::NDOF] = {"TX", "TY", "TZ", "RX", "RY", "RZ"};

  RigidBodyDOFSet() : DOFSet(NDOF) {}
  
  // mask: bitmask of free DOFs (bit i = DOF i is free)
  explicit RigidBodyDOFSet(uint8_t mask) : DOFSet(NDOF)
  {
    for (int i = 0; i < NDOF; ++i) {
      mFree[i] = (mask >> i) & 1;
    }
  }
  Type getType() const override { return Type::RigidBody; }
  std::string dofName(int idx) const override { return RigidBodyDOFNames[idx]; }
  void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const override;
  uint8_t mask() const
  {
    uint8_t m = 0;
    for (int i = 0; i < NDOF; ++i) {
      m |= (uint8_t(mFree[i]) << i);
    }
    return m;
  }
};

  
} // namespace o2::alignrs

#endif
