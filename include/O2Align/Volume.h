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

#ifndef O2_ALIGN_VOLUME_H
#define O2_ALIGN_VOLUME_H

#include <memory>
#include <utility>
#include <vector>
#include <ostream>
#include <string>
#include <map>
#include <array>
#include <algorithm>

#include <Eigen/Dense>
#include <nlohmann/json_fwd.hpp>

#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include "O2Align/Label.h"
#include "O2Align/DOFSet.h"

namespace o2::alignrs
{
using Matrix36 = Eigen::Matrix<double, 3, 6>;
using Matrix66 = Eigen::Matrix<double, 6, 6>;

class Volume
{
 protected:
  /// injected misalignment loaded from a closure-test JSON, indexed by sensor ID; used by
  /// writeMillepedeResults to subtract the known input from the fitted values. Protected (rather
  /// than private) so the type is nameable in the MP2JSON_Calib/MP2ROOT_Calib overrides. Declared
  /// first since it must be visible when used as a parameter type by the public methods below.
  struct InjectedMisalignment {
    std::map<int, std::vector<double>> rigidBody;
    std::map<int, std::vector<std::vector<double>>> matrix;
    struct Inextensional {
      std::map<int, std::array<double, 4>> modes;
      double alpha{0.};
      double beta{0.};
    };
    std::map<int, Inextensional> inextensional;
  };

 public:
  using Ptr = std::unique_ptr<Volume>;
  using SensorMapping = std::map<Label, Volume*>;

  Volume(const Volume&) = delete;
  Volume(Volume&&) = delete;
  Volume& operator=(const Volume&) = delete;
  Volume& operator=(Volume&&) = delete;
  /// \param virt if true the volume has no counterpart in the geometry (fictitious envelope,
  ///             TPC readout sector, ...); no alignable entry is looked up and the volume is
  ///             expected to define its own L2G matrix in defineMatrixL2G()
  Volume(const char* symName, uint32_t label, uint32_t det, bool sens, bool virt = false);
  Volume(const char* symName, Label label, bool virt = false);
  virtual ~Volume() = default;

  static void applyDOFConfig(Volume* root, const std::string& jsonPath);
  static void writeMillepedeResults(Volume* root, const std::string& milleResPath, const std::string& outJsonPath, const std::string& injectedJsonPath = "");

  /// create the common root of the whole hierarchy: a fictitious volume with the identity L2G, whose
  /// children are the top volumes of the individual detectors. It is not rigid-body alignable, hence
  /// the detector top volumes are not subjected to any automatic mutual constraint.
  static Ptr makeRoot(const char* symName = "ALICE");

  void finalise(uint8_t level = 0);

  // steering file output
  void writeRigidBodyConstraints(std::ostream& os) const;
  void writeParameters(std::ostream& os) const;
  void writeTree(std::ostream& os, int indent = 0) const;

  // tree-like
  auto getLevel() const noexcept { return mLevel; }
  bool isRoot() const noexcept { return mParent == nullptr; }
  bool isLeaf() const noexcept { return mChildren.empty(); }
  
  template <class T = Volume>
    requires std::derived_from<T, Volume>
  Volume* addChild(const char* symName, uint32_t label, uint32_t det, bool sens, bool virt = false)
  {
    auto c = std::make_unique<T>(symName, label, det, sens, virt);
    return setParent(std::move(c));
  }

  template <class T = Volume>
    requires std::derived_from<T, Volume>
  Volume* addChild(const char* symName, Label lbl, bool virt = false)
  {
    auto c = std::make_unique<T>(symName, lbl, virt);
    return setParent(std::move(c));
  }

  /// graft an externally built subtree (e.g. the hierarchy of a single detector) under this volume
  Volume* adoptChild(Ptr c) { return setParent(std::move(c)); }

  // bfs traversal
  void traverse(const std::function<void(Volume*)>& visitor)
  {
    visitor(this);
    for (auto& c : mChildren) {
      c->traverse(visitor);
    }
  }

  std::string getSymName() const noexcept { return mSymName; }
  Label getLabel() const noexcept { return mLabel; }
  Volume* getParent() const { return mParent; }
  int getNChildren() const noexcept { return static_cast<int>(mChildren.size()); }

  // DOF management
  void setRigidBody(std::unique_ptr<DOFSet> rb) { mRigidBody = std::move(rb); }
  void setCalib(std::unique_ptr<DOFSet> cal) { mCalib = std::move(cal); }
  DOFSet* getRigidBody() const { return mRigidBody.get(); }
  DOFSet* getCalib() const { return mCalib.get(); }
  void setPseudo(bool p) noexcept { mIsPseudo = p; }
  bool isPseudo() const noexcept { return mIsPseudo; }
  bool isVirtual() const noexcept { return mVirtual; }
  /// A volume which is not rigid-body alignable never receives rigid-body DOFs (a matching
  /// rigidBody rule of the DOF config is ignored) and, having none, generates no rigid-body
  /// constraint over its children: its branches are mutually unconstrained. Calibration DOFs are
  /// not affected: an envelope w/o its own geometry may still own a calibration DOF set.
  void setRigidBodyAllowed(bool v) noexcept { mRBAllowed = v; }
  bool isRigidBodyAllowed() const noexcept { return mRBAllowed; }
  void setSensorId(int id) noexcept { mSensorId = id; }
  int getSensorId() const noexcept { return mSensorId; }
  // true if this volume participates in the hierarchy (has DOFs or is pseudo)
  bool isActive() const noexcept { return mRigidBody != nullptr || mIsPseudo; }

  // transformation matrices
  virtual void defineMatrixL2G() {}
  virtual void defineMatrixT2L() {}
  /// jacobian of the (LOC)->(TRK) rigid-body parameter transformation at the local point posLoc
  virtual void computeJacobianL2T(const double* posLoc, Matrix66& jac) const;
  /// (LOC)->(GLO) matrix: for sensors and fictitious volumes it is defined by the volume itself,
  /// for the rest it is taken from the (possibly pre-aligned) geometry
  const TGeoHMatrix& getMatrixL2G() const;
  const TGeoHMatrix& getL2P() const { return mL2P; }
  const TGeoHMatrix& getT2L() const { return mT2L; }
  const Matrix66& getJL2P() const { return mJL2P; }
  const Matrix66& getJP2L() const { return mJP2L; }

  /// write the rigid-body block of the closure-test JSON output for this volume, subtracting the
  /// injected reference value (if inj is not nullptr) from the fitted one. Generic over any
  /// RigidBodyDOFSet layout, hence implemented once in the base class. The caller
  /// (writeMillepedeResults) only calls this when the volume has free rigid-body DOFs.
  virtual bool MP2JSON_RB(const std::map<uint32_t, double>& labelToValue, const std::vector<double>* inj, nlohmann::json& entry) const;
  /// write the calibration block of the closure-test JSON output for this volume. Detector-specific
  /// (the layout of the calibration DOFSet varies), hence the base implementation only warns that
  /// it is not implemented for this volume. inj may be nullptr if no misalignment was injected. The
  /// caller only calls this when the volume has free calibration DOFs.
  virtual bool MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const;

  /// ROOT-output counterparts of MP2JSON_RB/MP2JSON_Calib, same signature, not yet implemented.
  virtual bool MP2ROOT_RB(const std::map<uint32_t, double>& labelToValue, const std::vector<double>* inj, nlohmann::json& entry) const;
  virtual bool MP2ROOT_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const;

 protected:
  /// matrices
  Volume* mParent{nullptr}; // parent
  TGeoPNEntry* mPNE{nullptr};        // physical entry
  TGeoPhysicalNode* mPN{nullptr};    // physical node
  TGeoHMatrix mL2G;                  // (LOC) -> (GLO)
  TGeoHMatrix mL2P;                  // (LOC) -> (PAR)
  Matrix66 mJL2P;                    // jac (LOC) -> (PAR)
  Matrix66 mJP2L;                    // jac (PAR) -> (LOC)
  TGeoHMatrix mT2L;                  // (TRK) -> (LOC)

 private:
  static InjectedMisalignment loadInjectedMisalignment(const std::string& injectedJsonPath);

  std::string mSymName;
  Label mLabel;
  bool mVirtual{false}; // no counterpart in the geometry
  uint8_t mLevel{0};
  bool mIsPseudo{false};
  bool mRBAllowed{true}; // rigid-body DOFs may be assigned to this volume
  int mSensorId{-1}; // RS check if needed
  std::unique_ptr<DOFSet> mRigidBody;
  std::unique_ptr<DOFSet> mCalib;

  Volume* setParent(Ptr c)
  {
    c->mParent = this;
    mChildren.push_back(std::move(c));
    return mChildren.back().get();
  }
  std::vector<Ptr> mChildren; // children

  void init();
};

} // namespace o2::alignment

#endif
