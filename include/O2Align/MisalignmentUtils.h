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

#ifndef O2_ALIGN_MISALIGNMENTUTILS_H
#define O2_ALIGN_MISALIGNMENTUTILS_H

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "O2Align/AlignmentTypes.h"
#include "MathUtils/LegendrePols.h"

namespace o2::alignrs
{
class Volume;

struct InextensionalMisalignment {
  std::map<int, std::array<double, 4>> modes; // n -> (a_n, b_n, c_n, d_n)
  double alpha{0.};
  double beta{0.};
};

/// non-linear (exact) model of the deformation of an ITS3 half-barrel; the rigid-body part is always injected linearly
struct SensorMisalignment {
  o2::math_utils::Legendre2DPolynominal legendre;
  bool hasLegendre{false};
  InextensionalMisalignment inextensional;
  bool hasInextensional{false};

  bool empty() const noexcept { return !hasLegendre && !hasInextensional; }
};

struct MisalignmentFrame {
  int sensorID{-1};
  int layerID{-1};
  double x{0.};     // tracking-frame X / nominal radius at the measurement
  double alpha{0.}; // tracking-frame alpha
  double z{0.};     // tracking-frame measurement z
};

struct MisalignmentShift {
  double dy{0.};
  double dz{0.};
  bool accepted{true};

  MisalignmentShift& operator+=(const MisalignmentShift& other)
  {
    dy += other.dy;
    dz += other.dz;
    accepted = accepted && other.accepted;
    return *this;
  }
};

MisalignmentShift evaluateLegendreShift(const SensorMisalignment& sensor, const MisalignmentFrame& frame, const TrackSlopes& slopes);
MisalignmentShift evaluateInextensionalShift(const SensorMisalignment& sensor, const MisalignmentFrame& frame, const TrackSlopes& slopes);

/// Misalignment to inject, defined by the rule file AlignParams.misAlgJson (see doc/MisalignmentTest.md).
/// The rules name volumes (DOF-config style regex), their rigid-body DOFs and calibration DOFs and the
/// values to inject. They are validated against the configured hierarchy: a DOF which is not free in the
/// DOF configuration cannot be recovered by the fit and is fatal.
/// The values are injected linearly through the derivatives of the fit itself (getLinearValue by the raw
/// label of the DOF), except for the calibration clauses with "model": "exact", the non-linear ITS3
/// Legendre / inextensional deformation evaluated by evaluateLegendreShift / evaluateInextensionalShift.
class MisalignmentConfig
{
 public:
  /// parse the file and validate it against the DOF configuration of the finalised hierarchy, fatal on any
  /// mismatch. exactLegendreAllowed: the Monte Carlo truth needed by the exact Legendre model is available
  static MisalignmentConfig load(Volume* root, const std::string& jsonPath, bool exactLegendreAllowed);

  bool empty() const noexcept { return mAllLabels.empty(); }
  /// value injected linearly (through the fit derivatives) for a raw label, 0 if none
  double getLinearValue(uint32_t rawLabel) const
  {
    const auto it = mLinear.find(rawLabel);
    return it != mLinear.end() ? it->second : 0.;
  }
  /// exact-model misalignment of an ITS3 half-barrel sensor volume, nullptr if none
  const SensorMisalignment* getExact(const Volume* sensorVol) const
  {
    const auto it = mExact.find(sensorVol);
    return it != mExact.end() ? &it->second : nullptr;
  }
  /// all injected values (linear and exact), by raw label: the reference subtracted from the fit result
  const std::map<uint32_t, double>& getInjectedLabels() const noexcept { return mAllLabels; }
  /// log the injected volumes, DOFs and values
  void print() const;

 private:
  void addLinear(uint32_t label, const std::string& vol, const std::string& dof, double value);
  void addCalib(uint32_t label, const std::string& vol, const std::string& dof, double value, bool exact);
  struct Entry { // one injected DOF, for print()
    std::string volume;
    std::string dof;
    double value;
    bool exact;
  };
  std::unordered_map<uint32_t, double> mLinear;
  std::unordered_map<const Volume*, SensorMisalignment> mExact;
  std::map<uint32_t, double> mAllLabels;
  std::vector<Entry> mEntries;
};

} // namespace o2::alignrs

#endif
