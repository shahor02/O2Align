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
#include <array>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Rtypes.h>

namespace o2::alignrs
{

struct DerivativeContext {
  int detID{-1};    // Detector::DetIdx of the measurement
  int volID{-1};    // ID of the measured volume within its detector (TPC: sector)
  int sensorID{-1}; // ITS3 sensor ID, -1 for any other detector
  int layerID{-1};  // ITS3 layer ID, -1 for any other detector
  double measX{0.};
  double measAlpha{0.};
  double measZ{0.};
  double trkY{0.};
  double trkZ{0.};
  double snp{0.};
  double tgl{0.};
  double dydx{0.};
  double dzdx{0.};
  
  ClassDefNV(DerivativeContext,1);
};

// Generic set of DOF
class DOFSet
{
 public:
  enum class Type : uint8_t {
    RigidBody,
    Legendre,
    Inextensional,
    TPCVDrift,
    MeanVertex
  };
  
  virtual ~DOFSet() = default;
  virtual Type type() const = 0;
  virtual std::string dofName(int idx) const = 0;
  virtual void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const = 0;
  
  int nDOFs() const { return static_cast<int>(mFree.size()); }
  bool isFree(int idx) const { return mFree[idx]; }
  void setFree(int idx, bool f) { mFree[idx] = f; }
  void setAllFree(bool f) { std::fill(mFree.begin(), mFree.end(), f); }
  int nFreeDOFs() const {
    int n = 0;
    for (bool f : mFree) {
      n += f;
    }
    return n;
  }

  void setGlobalsPointers(double* pars, double* errs) {
    gParVals = pars;
    gParErrs = errs;
  }

  void setFirstEntry(int i) { mFirstEntry = i; }
  int getFirstEntry() const { return mFirstEntry; }

  double getParVal(int par) const { return getParVals()[par]; }
  double getParErr(int par) const { return getParErrs()[par]; }

  void setParVal(int par, double v = 0) { getParVals()[par] = v; }
  void setParErr(int par, double e = 0) { getParErrs()[par] = e; }

protected:  
  DOFSet(int n) : mFree(n, true) {}
  double* getParVals() { return gParVals + mFirstEntry; }
  double* getParErrs() { return gParErrs + mFirstEntry; }
  double* getParVals() const { return gParVals + mFirstEntry; }
  double* getParErrs() const { return gParErrs + mFirstEntry; }
  void validateDerivativeOutput(Eigen::Ref<Eigen::MatrixXd> out) const;

  std::vector<bool> mFree;         // status of each DOF
  int mFirstEntry = -1;            // ID of the 1st parameter in the global results array

  static double* gParVals; // start of global parameters array
  static double* gParErrs; // start of global parameters errors array
  
  ClassDef(DOFSet,1);
};

// Rigid body set (rotations and offset)
class RigidBodyDOFSet final : public DOFSet
{
 public:
  enum RigidBodyDOF : uint8_t { TX = 0, TY, TZ, RX, RY, RZ, NDOF };   // indices for rigid body parameters in LOC frame
  static constexpr const char* RigidBodyDOFNames[RigidBodyDOF::NDOF] = {"TX", "TY", "TZ", "RX", "RY", "RZ"};

  RigidBodyDOFSet() : DOFSet(NDOF) {}
  
  // mask: bitmask of free DOFs (bit i = DOF i is free)
  explicit RigidBodyDOFSet(uint8_t mask) : DOFSet(NDOF)
  {
    for (int i = 0; i < NDOF; ++i) {
      mFree[i] = (mask >> i) & 1;
    }
  }
  Type type() const override { return Type::RigidBody; }
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
  ClassDefOverride(RigidBodyDOFSet,1);
};

// Legendre DOFs
// Describing radial misplacement
class LegendreDOFSet final : public DOFSet
{
 public:
  explicit LegendreDOFSet(int order) : DOFSet((order + 1) * (order + 2) / 2), mOrder(order) {}
  Type type() const override { return Type::Legendre; }
  int order() const { return mOrder; }
  std::string dofName(int idx) const override;
  
  void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const override;

  static double getSensorPhiWidth(int sensorID, double radius);
  static std::pair<double, double> computeUV(double gloX, double gloY, double gloZ, int sensorID, double radius);
  /// phi range [phi1, phi2] (in [0, 2pi)) of the ITS3 half-barrel sensor, excluding the equatorial gap
  static std::pair<double, double> getSensorPhiBorders(int sensorID, double radius);

  /// 8-point Gauss-Legendre quadrature on [-1, 1] (Numerical Recipes 4.6), used to integrate the
  /// radial deformation along the arc
  static constexpr std::array<double, 8> GaussNodes = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
  static constexpr std::array<double, 8> GaussWeights = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
  
 private:
  int mOrder;
  ClassDefOverride(LegendreDOFSet, 1);
};

// In-extensional deformation DOFs for cylindrical half-shells
// Fourier modes n=2..N: 4 params each (a_n, b_n, c_n, d_n)
// Plus 2 non-periodic modes (alpha, beta) for the half-cylinder open edges
// Total: 4*(N-1) + 2
class InextensionalDOFSet final : public DOFSet
{
 public:
  explicit InextensionalDOFSet(int maxOrder) : DOFSet((4 * (maxOrder - 1)) + 2), mMaxOrder(maxOrder)
  {
    if (maxOrder < 2) {
      // the rest is eq. to rigid body
      throw std::invalid_argument("InextensionalDOFSet requires maxOrder >= 2");
    }
  }
  Type type() const override { return Type::Inextensional; }
  int maxOrder() const { return mMaxOrder; }

  // number of periodic DOFs (before alpha, beta)
  int nPeriodic() const { return 4 * (mMaxOrder - 1); }

  // flat index layout: [a_2, b_2, c_2, d_2, a_3, b_3, c_3, d_3, ..., alpha, beta]
  // index of first DOF for mode n
  static int modeOffset(int n) { return 4 * (n - 2); }

  // indices of the non-periodic modes
  int alphaIdx() const { return nPeriodic(); }
  int betaIdx() const { return nPeriodic() + 1; }

  std::string dofName(int idx) const override;
  void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const override;

 private:
  int mMaxOrder;
  ClassDefOverride(InextensionalDOFSet,1);
};

/// Calibration of the TPC drift, owned by the TPC envelope volume. The cluster Z is reconstructed
/// from the drift time as z = s * (Lz - l), l = (t - T0 - tOffset) * vDrift being the drift length
/// and s = +1 (-1) on the A (C) side, hence a relative correction d of vDrift and an offset o of
/// the drift length displace the measurement along the local (= tracking frame) Z by
///     dz = -s * l * d - s * o,      l = Lz - s * z
/// The two are separated by the lever arm in l. Being a pure Z shift of the measurement, the
/// derivatives need no transformation to the local frame and are scaled with the same convention
/// as the TZ column of RigidBodyDOFSet.
class TPCVDriftDOFSet final : public DOFSet
{
 public:
  enum VDriftDOF : uint8_t { VDRIFT = 0, DRIFTOFF, NDOF }; // relative vDrift correction, drift length offset [cm]
  static constexpr const char* VDriftDOFNames[VDriftDOF::NDOF] = {"VDRIFT", "DRIFTOFF"};

  TPCVDriftDOFSet() : DOFSet(NDOF) {}
  Type type() const override { return Type::TPCVDrift; }
  std::string dofName(int idx) const override { return VDriftDOFNames[idx]; }
  void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const override;

  /// length of the TPC drift volume, TPCFastTransformGeo::getTPCzLength(). Known only once the
  /// correction maps are delivered, i.e. after the DOF set is created by the DOF configuration.
  void setZLength(double zLength) { mZLength = zLength; }
  double getZLength() const { return mZLength; }

 private:
  double mZLength{-1.}; // negative until set from the correction maps
  ClassDefOverride(TPCVDriftDOFSet, 1);
};

/// Position of the mean interaction vertex, owned by the single volume of DetectorPVT and fitted
/// independently in every time slot of the mean vertex calibration. The DOFs are the corrections
/// to the X, Y, Z of the MeanVertexObject in the global frame. Being the position of the measurement
/// (the mean-vertex prior) rather than of a sensor, its derivatives are not produced by the generic
/// chain but directly by the vertex constraint (AlignmentSpec::addMeanVertexPrior).
class MeanVertexDOFSet final : public DOFSet
{
 public:
  enum MeanVertexDOF : uint8_t { X = 0, Y, Z, NDOF };
  static constexpr const char* MeanVertexDOFNames[MeanVertexDOF::NDOF] = {"X", "Y", "Z"};

  MeanVertexDOFSet() : DOFSet(NDOF) {}
  Type type() const override { return Type::MeanVertex; }
  std::string dofName(int idx) const override { return MeanVertexDOFNames[idx]; }
  /// not applicable: throws, see the class description
  void fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const override;

  ClassDefOverride(MeanVertexDOFSet, 1);
};

} // namespace o2::alignrs

#endif
