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
class TimeSlotsSet;
class Constraint;
}

namespace o2::alignrs
{
using Matrix26 = Eigen::Matrix<double, 2, 6>; // d(Y,Z)/d(rigid-body DOFs)
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
  /// parse millepede.res into label -> fitted value, skipping the fixed parameters
  static std::map<uint32_t, double> readMillepedeResults(const std::string& milleResPath);
  static void writeMillepedeResults(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outJsonPath, const std::string& injectedJsonPath = "");

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

  /// \name Time-sliced calibration DOFs
  /// The calibration DOF set is configured once, but its parameters may be fitted independently in
  /// consecutive time slots, the slot ID being encoded in the volume ID of the calibration label
  /// (see Detector::setTimeStamp). By default a volume has the single label of the slot 0.
  /// The rigid-body DOFs are never time-sliced: a time-dependent position (the mean vertex) is a
  /// calibration DOF set.
  ///@{
  /// declare the calibration slots of this volume. The slots are not owned, they are kept to report
  /// the validity of the fitted values.
  void setCalibSlots(const TimeSlotsSet& slots);
  /// select the slot the processed TF belongs to. Fatal if the slot was not declared.
  void setActiveCalibSlot(int slotID);
  /// calibration label of the processed TF, the one the derivatives are attributed to
  const Label& getActiveCalibLabel() const { return mCalibLabels[mActiveCalibSlot]; }
  /// calibration labels of all the declared slots, for the steering file output
  const std::vector<Label>& getCalibLabels() const { return mCalibLabels; }
  ///@}
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
  // true if this volume participates in the hierarchy (has a free RB DOF or is pseudo)
  bool isActive() const noexcept { return (mRigidBody && mRigidBody->nFreeDOFs() > 0) || mIsPseudo; }

  // transformation matrices
  virtual void defineMatrixL2G() {}
  virtual void defineMatrixT2L() {}
  /// Jacobian of the (LOC)->(TRK) rigid-body parameter transformation (thesis A.64), for the TRK
  /// rotations taken about the point posLoc (in LOC) instead of the TRK origin: the translation
  /// block is then -R*[posLoc]x, which equals [t']x*R of A.64 with t' the LOC origin wrt the pivot.
  /// t2l is the (TRK)->(LOC) matrix, R = its rotation transposed.
  static void computeJacobianL2T(const TGeoHMatrix& t2l, const double* posLoc, Matrix66& jac);
  /// (TRK)->(LOC) matrix of the tracking frame rotated by alpha from the global one; equals getT2L()
  /// for the alpha of the volume. Needed when the tracking frame depends on the measured point (ITS3).
  /// Valid for the leaves (sensors) only, once their L2G is defined.
  TGeoHMatrix computeT2L(double alpha) const;
  /// the tracking frame is the same for all the points of this sensor, i.e. getT2L() is valid for all
  bool hasFixedTrackingFrame() const noexcept { return mFixedTrackingFrame; }
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
  /// fitted value of a parameter, 0 if it was not fitted (fixed or absent from millepede.res)
  static double getFittedValue(const std::map<uint32_t, double>& labelToValue, uint32_t rawLabel)
  {
    const auto it = labelToValue.find(rawLabel);
    return it != labelToValue.end() ? it->second : 0.;
  }
  /// write the calibration block of the closure-test JSON output for this volume. The base
  /// implementation writes, under "calib", the absolute fitted values of the free DOFs of every
  /// calibration slot (see calibSlotsToJSON); a detector-specific layout or the subtraction of an
  /// injected calibration (inj, nullptr if none) needs an override. The caller only calls this when
  /// the volume has free calibration DOFs.
  virtual bool MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const;

  /// ROOT-output counterparts of MP2JSON_RB/MP2JSON_Calib, same signature, not yet implemented.
  virtual bool MP2ROOT_RB(const std::map<uint32_t, double>& labelToValue, const std::vector<double>* inj, nlohmann::json& entry) const;
  virtual bool MP2ROOT_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const;

  uint32_t getDataCounter() const noexcept { return mStat; }
  void incDataCounter() noexcept { ++mStat; }
  void setDataCounter(uint32_t c) noexcept { mStat = c; }

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
  TGeoHMatrix mG2L;                  // (GLO) -> (LOC), the inverse of getMatrixL2G(), for leaves only
  bool mFixedTrackingFrame{true};    // see hasFixedTrackingFrame

  /// set mT2L for the standard ALICE tracking frame, i.e. the global one rotated by alpha around Z
  void setT2LFromAlpha(double alpha) { mT2L = computeT2L(alpha); }
  /// one record per calibration slot with at least one fitted DOF: the slot ID and validity (if the
  /// slots were declared) and the absolute fitted values of the free DOFs, by DOF name
  nlohmann::json calibSlotsToJSON(const std::map<uint32_t, double>& labelToValue) const;

 private:
  static InjectedMisalignment loadInjectedMisalignment(const std::string& injectedJsonPath);
  /// apply the "fixed"/"free"/"fix" clauses of a calib rule, common to all calibration DOF types
  static void applyFreeFixConfig(DOFSet& dofSet, const nlohmann::json& cal);
  /// DOF sets described by the rigidBody / calib clauses of a rule of the DOF configuration
  static std::unique_ptr<DOFSet> makeRigidBodyDOFSet(const nlohmann::json& rb, const std::string& pattern);
  static std::unique_ptr<DOFSet> makeCalibDOFSet(const nlohmann::json& cal, const std::string& pattern, const std::string& sym);
  /// apply the "pinChildrenMean" clause of a rule: flags rigid-body DOFs (free or fixed) for which
  /// writeRigidBodyConstraints must still force the weighted mean of the active children's movement
  /// to vanish, mirroring the "replace, don't merge" semantics of the rigidBody/calib clauses
  void applyPinChildrenMeanConfig(const nlohmann::json& pin, const std::string& pattern);
  /// constraints of this volume over its children: the vanishing mean movement of the children in
  /// every free or pinned rigid-body DOF of this volume, see writeRigidBodyConstraints
  void writeChildrenMeanConstraints(std::ostream& os) const;
  /// a free-DOF descendant contributing to the mean constraints of this volume, with the jacobian
  /// transporting its DOFs to the frame of this volume
  struct MeanContributor {
    const Volume* vol;
    Matrix66 jToThis;
  };
  /// collect the nearest descendants having a free RB DOF: a child w/o one (no DOF set or all fixed)
  /// is looked through, its own children being taken with its getJP2L() folded in
  void collectMeanContributors(const Matrix66& jToThis, std::vector<MeanContributor>& out) const;
  /// add to con the free DOFs of the contributors feeding the DOF iDOF of this volume, scaled by weight
  static void addChildrenMeanTerms(Constraint& con, int iDOF, const std::vector<MeanContributor>& contributors, double weight);

  std::string mSymName;
  Label mLabel;
  bool mVirtual{false}; // no counterpart in the geometry
  uint8_t mLevel{0};
  bool mIsPseudo{false};
  bool mRBAllowed{true}; // rigid-body DOFs may be assigned to this volume
  int mSensorId{-1}; // RS check if needed
  uint32_t mStat{0}; // counter for tracks contributing to this volume
  std::unique_ptr<DOFSet> mRigidBody;
  std::unique_ptr<DOFSet> mCalib;
  /// per rigid-body DOF: force the weighted mean of the active children's movement (transported to
  /// this volume's frame) to vanish even though the DOF itself is fixed here, or the volume has no
  /// rigid-body DOF set at all. Set by the
  /// "pinChildrenMean" clause of the DOF config; see writeRigidBodyConstraints.
  std::array<bool, RigidBodyDOFSet::NDOF> mPinChildrenMean{};
  std::vector<Label> mCalibLabels{};            // calibration label per time slot, filled with the slot 0 in the ctor
  size_t mActiveCalibSlot{0};                   // index in mCalibLabels of the slot being processed
  const TimeSlotsSet* mCalibSlots{nullptr};     //! calibration slots, if declared; not owned

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
