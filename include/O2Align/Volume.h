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
#include <limits>
#include <algorithm>
#include <functional>

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

/// How the coefficients of the "globalModes"/"groupModes" equations are built, see
/// Volume::writeGlobalModeConstraints and doc/DOFConfig_rules.md
struct ModeEquationOptions {
  /// coefficient of a DOF = projection of the mode on the displacement this DOF gives to the
  /// sensor-box corners of its subtree; false: the mode evaluated at the volume's origin (legacy)
  bool sensorProjection = true;
  /// "groupModes": constrain the mode relative to the free rigid motion of each group volume, i.e.
  /// remove from the mode its least-squares rigid part over the volume's sensors (needs sensorProjection)
  bool relativeToGroupRoots = true;
};

class Volume
{
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
  /// parse millepede.res into label -> fitted value, fixed parameters included (at their preset value)
  static std::map<uint32_t, double> readMillepedeResults(const std::string& milleResPath);
  static void writeMillepedeResults(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outJsonPath);

  /// create the common root of the whole hierarchy: a fictitious volume with the identity L2G, whose
  /// children are the top volumes of the individual detectors. It is not rigid-body alignable, hence
  /// the detector top volumes are not subjected to any automatic mutual constraint.
  static Ptr makeRoot(const char* symName = "ALICE");

  void finalise(uint8_t level = 0);

  // steering file output
  /// statWeighted: default weighting of the hierarchical mean-of-children constraint (see
  /// writeChildrenMeanConstraints) for every DOF of every volume whose own "weightedMean" clause did
  /// not override it; modeOpt: see ModeEquationOptions
  void writeRigidBodyConstraints(std::ostream& os, bool statWeighted = false, const ModeEquationOptions& modeOpt = {}) const;
  void writeParameters(std::ostream& os) const;
  void writeTree(std::ostream& os, int indent = 0) const;

  /// \name Global deformation-mode constraints
  /// The hierarchical constraint of writeChildrenMeanConstraints / pinChildrenMean is a *rigid-body*
  /// constraint: 6 scalar equations per parent, removing the net translation and rotation of that
  /// parent's set of free descendants and nothing else. A coherent deformation of the whole barrel --
  /// radial scaling, elliptical squashing, twist, ... -- is not a net rigid motion of any single
  /// parent, carries no net rigid component, and therefore passes straight through it. Such modes are
  /// also (very nearly) null directions of the track chi2, so Millepede takes large amplitudes in them
  /// for a negligible chi2 gain. The "globalModes" clause of the DOF config flags the volumes whose
  /// free rigid-body DOFs take part in one global equation per mode, forcing the total amplitude of
  /// that mode over the whole flagged set to vanish. This is a different operator from the mean
  /// constraint, not a different setting of it.
  ///@{
  enum GlobalMode : int {
    RadialScaling,      ///< dr = eps*r       : uniform transverse dilation, absorbed by the 1/pt scale
    ZScaling,           ///< dz = eps*z       : longitudinal scale
    Curl,               ///< r*dphi = eps*r   : rigid rotation about the beam
    Twist,              ///< r*dphi = eps*z*r : torsion along the beam
    Elliptical,         ///< dr = eps*cos(2phi) and dr = eps*sin(2phi): TWO equations
    Bowing,             ///< dr = eps*z^2     : radial sagitta along the beam
    RadialHarmonic,     ///< dr = eps*r*cos(l*phi) / sin(l*phi), l = 1..GlobalHarmonicLMax
    TangentialHarmonic, ///< r*dphi = eps*r*cos(l*phi) / sin(l*phi), l = 1..GlobalHarmonicLMax
    GlobalTilt,         ///< rigid rotation of the whole branch about the global X / Y axis: TWO equations
    RadialPolar,        ///< dr     = eps*r^p*t^k, t = z/r: radial displacement set by the polar angle
    TangentialPolar,    ///< r*dphi = eps*r^p*t^k
    LongitudinalPolar,  ///< dz     = eps*r^p*t^k
    GlobalShift,        ///< rigid translation along the global X / Y / Z axis: THREE equations
    NGlobalModes
  };
  /// highest azimuthal order of the RadialHarmonic / TangentialHarmonic families
  static constexpr int GlobalHarmonicLMax = 4;
  /// highest power k of t = z/r in the polar families
  static constexpr int GlobalPolarKMax = 2;
  static constexpr bool isPolarMode(int mode) { return mode == RadialPolar || mode == TangentialPolar || mode == LongitudinalPolar; }
  /// The polar families are the displacement fields a track from the vertex cannot see: at radius r it
  /// is at z ~ r*t (t = tgLambda), so a field r*f(t) is absorbed by its phi (rphi), tgLambda (z) or,
  /// with B != 0, curvature (r), and a field f(t) by its d0 / z0, i.e. broken only by the vertex.
  /// They are invisible as long as every free DOF sits at z ~ 0 (stave level), and become accessible
  /// once DOFs distributed in z (modules, chips) are freed. The members (p, k), p = 0, 1, k = 0..KMax,
  /// which are already modes of their own are left out: radial (1,0) = RadialScaling, tangential
  /// (1,0) = Curl, longitudinal (1,1) = ZScaling; so is longitudinal (0,0), a rigid z shift.
  static constexpr bool isPolarDuplicate(int mode, int p, int k)
  {
    return (p == 1 && k == 0 && (mode == RadialPolar || mode == TangentialPolar)) ||
           (mode == LongitudinalPolar && ((p == 0 && k == 0) || (p == 1 && k == 1)));
  }
  static constexpr int polarNEq(int mode)
  {
    int n = 0;
    for (int p = 0; p < 2; ++p) {
      for (int k = 0; k <= GlobalPolarKMax; ++k) {
        n += isPolarDuplicate(mode, p, k) ? 0 : 1;
      }
    }
    return n;
  }
  /// powers (p, k) of the polar-family equation 'phase', in the order p-major, k-minor
  static constexpr std::pair<int, int> polarPK(int mode, int phase)
  {
    for (int p = 0; p < 2; ++p) {
      for (int k = 0; k <= GlobalPolarKMax; ++k) {
        if (!isPolarDuplicate(mode, p, k) && phase-- == 0) {
          return {p, k};
        }
      }
    }
    return {-1, -1};
  }
  /// number of scalar equations of a mode: 2 for Elliptical (phases cos2phi, sin2phi) and for
  /// GlobalTilt (about X, about Y), 3 for GlobalShift (along X, Y, Z), 2*GlobalHarmonicLMax for the
  /// harmonic families, polarNEq for the polar families, 1 for the rest
  static constexpr int globalModeNEq(int mode)
  {
    if (isPolarMode(mode)) {
      return polarNEq(mode);
    }
    if (mode == GlobalShift) {
      return 3;
    }
    return (mode == Elliptical || mode == GlobalTilt) ? 2 : ((mode == RadialHarmonic || mode == TangentialHarmonic) ? 2 * GlobalHarmonicLMax : 1);
  }
  /// flat index of the first equation of a mode, and the total number of equations of all the modes.
  /// The flat index is what the DOF-config flags, what writeGlobalModeConstraints loops over and what
  /// indexes the columns of the WeakModes diagnostics.
  static constexpr int globalEqFirst(int mode)
  {
    int n = 0;
    for (int m = 0; m < mode; ++m) {
      n += globalModeNEq(m);
    }
    return n;
  }
  /// Total number of mode equations. Spelled out rather than computed with globalEqFirst() because a
  /// constexpr member function may not be called in a constant expression inside its own class body
  /// (rootcling rejects it); the static_assert right after the class checks the two agree.
  static constexpr int NGlobalEq = 5 /* the one-equation modes */ + 2 /* Elliptical */ + 2 /* GlobalTilt */ + 2 * (2 * GlobalHarmonicLMax) +
                                   3 * 2 * (GlobalPolarKMax + 1) - 4 /* polar families minus their duplicates */ + 3 /* GlobalShift */;
  static constexpr int globalEqIndex(int mode, int phase) { return globalEqFirst(mode) + phase; }
  static int globalEqMode(int eq);
  static int globalEqPhase(int eq);
  /// azimuthal order and cos/sin selection encoded in the phase of a harmonic-family equation
  static constexpr int globalHarmonicL(int phase) { return phase / 2 + 1; }
  static constexpr bool globalHarmonicIsSin(int phase) { return phase % 2 != 0; }
  /// name of a mode as used in the "globalModes" clause
  static const char* globalModeName(int mode);
  /// name of one equation: the mode name for a single-equation mode, otherwise qualified by the phase,
  /// e.g. "elliptical_cos2phi", "radialHarmonic3_sin3phi". Accepted by the "globalModes" clause too.
  static std::string globalEqName(int eq);
  /// displacement field u (translation) and local rotation omega = curl(u)/2 of the unit-amplitude mode
  /// at the global position pos; phase selects the equation within the mode; false if undefined there
  static bool globalModeField(int mode, int phase, const double* pos, double* u, double* om);
  /// one Millepede Constraint per requested mode and phase, spanning the free rigid-body DOFs of every
  /// volume flagged by a "globalModes" rule. A mode equation is global, so this is written once for the
  /// whole tree, not per volume: call it on the root (writeRigidBodyConstraints already does).
  /// It also writes the "groupModes" equations: one set per group, each over the free rigid-body DOFs
  /// of the whole subtrees of the group's volumes, see writeGroupModeConstraints.
  void writeGlobalModeConstraints(std::ostream& os, const ModeEquationOptions& opt = {}) const;
  ///@}

  /// \name Ray-mode constraints
  /// A track from the vertex crosses every layer at the same azimuth phi and t = z/r (its tgLambda),
  /// so a sensor displacement r*F(phi, t)*e, e = phiHat, zHat or rHat, is absorbed by the track's phi,
  /// tgLambda or (B != 0) curvature for *any* function F: a function space of null directions, not a
  /// short list of global modes. It becomes accessible as soon as free DOFs are distributed in z
  /// (modules, chips). The "rayModes" clause bins the sensors of the flagged branch in (phi, t) cells
  /// and writes, per cell and component, sum_s (e . d_s) / r_s = 0 over the sensors s of the cell,
  /// d_s being the linearised displacement of the sensor centre by every free rigid-body DOF on its
  /// path to the root. This removes the common-ray motion of each cell and nothing finer than the cell.
  ///@{
  struct RayModesConfig {
    int nPhi{24};                     ///< azimuthal cells over [phi0, phi0 + 2pi)
    double phi0{0.};                  ///< lower edge of the first azimuthal cell
    int nT{16};                       ///< polar cells over [-tMax, tMax], sensors beyond go to the edge cells
    double tMax{4.};
    double rSpanMin{1.};              ///< cm: a cell whose sensors span a smaller radial range (one layer) is skipped
    std::vector<double> rEdges{};     ///< cm: optional radial boundaries, each radial range gets its own cells
    double rMin{0.};                  ///< cm: sensors at smaller radius are not binned (see writeRayModeConstraints)
    double rMax{std::numeric_limits<double>::max()}; ///< cm: nor those at larger radius
    std::array<bool, 3> comps{true, true, true}; ///< components r*dphi, dz, dr
  };
  static constexpr int NRayComps = 3;
  static const char* rayCompName(int comp);
  /// one Constraint per (radial range, phi cell, t cell, component) for every branch flagged by a
  /// "rayModes" rule; written once for the whole tree from the root (writeRigidBodyConstraints does)
  void writeRayModeConstraints(std::ostream& os) const;
  ///@}

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

  /// visit the direct children
  void forEachChild(const std::function<void(const Volume*)>& visitor) const
  {
    for (const auto& c : mChildren) {
      visitor(c.get());
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

  /// write the rigid-body block of the JSON output (result.json) for this volume, the fitted values
  /// as given in labelToValue. Generic over any
  /// RigidBodyDOFSet layout, hence implemented once in the base class. The caller
  /// (writeMillepedeResults) only calls this when the volume has free rigid-body DOFs.
  virtual bool MP2JSON_RB(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const;
  /// local delta (LOC->LOC, l' = delta*l) built from the fitted rigid-body parameters of this volume.
  /// Returns false if it has no free rigid-body DOF or all its fitted values are 0 (fixed parameters
  /// are absent from labelToValue, hence count as 0).
  bool getFittedLocalDelta(const std::map<uint32_t, double>& labelToValue, TGeoHMatrix& delta) const;
  /// fitted value of a parameter, 0 if it was not fitted (fixed or absent from millepede.res)
  static double getFittedValue(const std::map<uint32_t, double>& labelToValue, uint32_t rawLabel)
  {
    const auto it = labelToValue.find(rawLabel);
    return it != labelToValue.end() ? it->second : 0.;
  }
  /// write the calibration block of the closure-test JSON output for this volume. The base
  /// implementation writes, under "calib", the absolute fitted values of the free DOFs of every
  /// calibration slot (see calibSlotsToJSON); a detector-specific layout needs an override. The caller only calls this when
  /// the volume has free calibration DOFs.
  virtual bool MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const;

  /// ROOT-output counterparts of MP2JSON_RB/MP2JSON_Calib, same signature, not yet implemented.
  virtual bool MP2ROOT_RB(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const;
  virtual bool MP2ROOT_Calib(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const;

  uint32_t getDataCounter() const noexcept { return mStat; }
  void incDataCounter() noexcept { ++mStat; }
  void setDataCounter(uint32_t c) noexcept { mStat = c; }

  /// true if this volume or any descendant has an explicit "weightedMean" rule requesting statistics
  /// weighting for at least one DOF, regardless of the statWeighted default passed to
  /// writeRigidBodyConstraints. Used to decide whether per-volume statistics must be available (loaded
  /// from Params::statInFile in the standalone MilleSteer case) even if Params::useStatWeightedMeanConstraints
  /// itself is off.
  bool hasAnyWeightedMean() const;

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
  /// apply the "fixed"/"free"/"fix" clauses of a calib rule, common to all calibration DOF types
  static void applyFreeFixConfig(DOFSet& dofSet, const nlohmann::json& cal);
  /// DOF sets described by the rigidBody / calib clauses of a rule of the DOF configuration
  static std::unique_ptr<DOFSet> makeRigidBodyDOFSet(const nlohmann::json& rb, const std::string& pattern);
  static std::unique_ptr<DOFSet> makeCalibDOFSet(const nlohmann::json& cal, const std::string& pattern, const std::string& sym);
  /// apply the "pinChildrenMean" (noCross = false) or "pinChildrenMeanNoCross" (noCross = true) clause
  /// of a rule: flags rigid-body DOFs (free or fixed) for which writeRigidBodyConstraints must still
  /// force the weighted mean of the active children's movement to vanish, mirroring the
  /// "replace, don't merge" semantics of the rigidBody/calib clauses. The two clauses have independent
  /// flags and may both be set on the same volume (they are different equations).
  void applyPinChildrenMeanConfig(const nlohmann::json& pin, const std::string& pattern, bool noCross);
  /// apply the "weightedMean" clause of a rule: a bool (all DOFs), "all", or an array of DOF names.
  /// Per DOF of this volume, overrides whether its mean-of-children constraint (the automatic one for
  /// a free DOF, or the one written for pinChildrenMean/pinChildrenMeanNoCross) weights each
  /// contributor by its accumulated getDataCounter() instead of splitting the mean equally (1/N) --
  /// see writeChildrenMeanConstraints. A DOF never touched by any matching "weightedMean" rule keeps
  /// falling back to the statWeighted default passed to writeRigidBodyConstraints; once a rule does
  /// match, that DOF is set explicitly (true or false) and no longer falls back, exactly as the other
  /// clauses "replace, don't merge" on a later match.
  void applyWeightedMeanConfig(const nlohmann::json& cfg, const std::string& pattern);
  /// apply the "measurement" clause of a rule: {"TX": [value, sigma], "SlopeX": {"value": v, "sigma": s}, ...}
  /// external measurements of rigid-body AND/OR calibration DOFs of this volume, resolved by name
  /// against whichever of mRigidBody/mCalib has a matching dofName() (the two never share a name),
  /// written as Millepede "Measurement" records; a later matching rule replaces all the measurements
  /// of an earlier one (of either DOF set), like the other clauses
  void applyMeasurementConfig(const nlohmann::json& meas, const std::string& pattern);
  /// apply the "globalModes" clause of a rule: flags this volume as a participant of the global
  /// deformation-mode equations, with the same "replace, don't merge" semantics as the other clauses
  void applyGlobalModesConfig(const nlohmann::json& modes, const std::string& pattern);
  /// apply the "groupModes" clause of a rule: a mode list as for "globalModes" (this volume is a group
  /// of its own), or {"modes": <list>, "group": "<name>"} (all volumes given the same name form one
  /// group); see writeGroupModeConstraints
  void applyGroupModesConfig(const nlohmann::json& cfg, const std::string& pattern);
  /// the mode equations selected by a "globalModes"/"groupModes" list: a bool, "all", or an array of
  /// mode, harmonic-order or equation names
  static std::array<bool, NGlobalEq> parseGlobalModeList(const nlohmann::json& modes, const std::string& clause, const std::string& pattern);
  /// accumulate the terms of one mode equation over this volume (if flagged for that mode, or always
  /// with allInSubtree) and its whole subtree, as (label, coefficient) pairs still carrying their
  /// physical lever arms, counting in nFlagged the volumes which took part with a free rigid-body DOF
  void collectGlobalModeTerms(int eq, std::vector<std::pair<uint32_t, double>>& terms, int& nFlagged, bool allInSubtree = false) const;
  /// displacement field of a mode at a global position (0 where the mode is undefined)
  using ModeField = std::function<void(const double* pos, double* u)>;
  static ModeField modeField(int eq);
  /// sensor-projected version of collectGlobalModeTerms (ModeEquationOptions::sensorProjection): the
  /// coefficient of a DOF is sum_p u(p).d(p), p running over the sample points of the sensors of the
  /// volume's subtree and d(p) the displacement of p for a unit value of the DOF. Returns the sums
  /// {U = sum u(p), M = sum p x u(p)} over this subtree, from which every coefficient follows.
  std::array<double, 6> collectProjectedModeTerms(const ModeField& field, std::vector<std::pair<uint32_t, double>>& terms, int& nFlagged, bool allInSubtree,
                                                  int eq) const;
  /// sample points (global frame) of the sensors of this subtree: the 8 corners of the box of each
  /// leaf, placed in its L2G frame, or its origin alone if it has no box
  void collectSensorSamplePoints(std::vector<std::array<double, 3>>& pts) const;
  /// field minus its least-squares projection, over the sensor sample points of this subtree, on the
  /// rigid motions of this volume's free DOFs; sum |u|^2 over the points before and after in norm2
  ModeField removeFreeRigidMotion(const ModeField& field, std::array<double, 2>& norm2) const;
  /// the volumes of this subtree carrying a "groupModes" clause, by group name
  void collectGroupModeRoots(std::map<std::string, std::vector<const Volume*>>& groups) const;
  /// one Millepede Constraint per group and requested mode equation, see the "groupModes" clause
  void writeGroupModeConstraints(std::ostream& os, uint32_t anchorLabel, const ModeEquationOptions& opt) const;
  /// normalise one collected mode equation, drop its rounding residues and write it with the pede
  /// grouping anchor; false (with a warning) if fewer than 2 terms survive
  static bool writeModeEquation(std::ostream& os, const std::string& name, const std::vector<std::pair<uint32_t, double>>& terms, int nFlagged, uint32_t anchorLabel);
  /// smallest Millepede label of a free rigid-body DOF in this subtree, UINT32_MAX if there is none
  uint32_t minFreeRigidBodyLabel() const;
  /// workaround for the constraint grouping of pede, see writeGlobalModeConstraints
  static void addPedeGroupingAnchor(Constraint& con, uint32_t anchorLabel);
  /// apply the "rayModes" clause of a rule: true (defaults), false, or an object overriding the
  /// fields of RayModesConfig ("nPhi", "phi0", "nT", "tMax", "rSpanMin", "rEdges", "components")
  void applyRayModesConfig(const nlohmann::json& cfg, const std::string& pattern);
  /// the sensors (leaves which are neither virtual nor pseudo) of this subtree
  void collectSensors(std::vector<const Volume*>& out) const;
  /// one Millepede "Measurement" per free rigid-body DOF of this volume having one configured
  void writeRigidBodyDOFMeasurements(std::ostream& os) const;
  /// one Millepede "Measurement" per free calibration DOF having one configured, for every declared
  /// calibration slot (the same value/sigma applies to every slot)
  void writeCalibDOFMeasurements(std::ostream& os) const;
  /// constraints of this volume over its children: the vanishing mean movement of the children in
  /// every free or pinned rigid-body DOF of this volume, see writeRigidBodyConstraints
  void writeChildrenMeanConstraints(std::ostream& os, bool statWeighted) const;
  /// a free-DOF descendant contributing to the mean constraints of this volume, with the jacobian
  /// transporting its DOFs to the frame of this volume
  struct MeanContributor {
    const Volume* vol;
    Matrix66 jToThis;
  };
  /// collect the nearest descendants having a free RB DOF: a child w/o one (no DOF set or all fixed)
  /// is looked through, its own children being taken with its getJP2L() folded in
  void collectMeanContributors(const Matrix66& jToThis, std::vector<MeanContributor>& out) const;
  /// Add to con the contributors' free DOFs feeding the DOF iDOF of this volume, scaled by weight.
  /// noCross = false: every free DOF jDOF of a contributor enters with the jacobian element
  /// jToThis(iDOF, jDOF), i.e. the contributor's *pose* transported into this volume's frame -- the
  /// correct operator when this volume's DOF is a fit parameter, since that is exactly the degeneracy.
  /// noCross = true: only the contributor's own iDOF enters, with the weight alone. Use this to pin a
  /// collective drift: the transported form mixes translations and rotations with lever arms of up to
  /// tens of cm, so it is cheaply satisfied by a weakly-observable rotation (an ITS half-stave RZ of
  /// 7.5 mrad cancelling a 97 um mean TY) instead of by the displacement one meant to constrain.
  /// weights has one entry per contributor, in the same order, already normalised (e.g. equal 1/N, or
  /// each contributor's share of the accumulated statistics) -- see writeChildrenMeanConstraints.
  static void addChildrenMeanTerms(Constraint& con, int iDOF, const std::vector<MeanContributor>& contributors, const std::vector<double>& weights, bool noCross);
  /// Largest |jacobian element| coupling the DOF iDOF of this volume to a free DOF of any contributor.
  static double maxChildrenCoupling(int iDOF, const std::vector<MeanContributor>& contributors);
  /// Below this coupling the contributors cannot realise the parent DOF (only frame misrotations
  /// leak into it), so its transported mean equation is a rescaled copy of another one.
  static constexpr double MinChildrenCoupling = 1e-2;

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
  /// as mPinChildrenMean, but the equation sums only the children's *same-named* DOF, without the
  /// jacobian cross-terms. Set by the "pinChildrenMeanNoCross" clause; see addChildrenMeanTerms.
  std::array<bool, RigidBodyDOFSet::NDOF> mPinChildrenMeanNoCross{};
  /// per rigid-body DOF: whether this volume's mean-of-children constraint is weighted by each
  /// contributor's accumulated statistics instead of an equal split; Unset (the default) falls back to
  /// the statWeighted argument of writeRigidBodyConstraints. Set by the "weightedMean" clause.
  enum class WeightedMeanFlag : uint8_t { Unset = 0,
                                           False,
                                           True };
  std::array<WeightedMeanFlag, RigidBodyDOFSet::NDOF> mWeightedMean{};
  /// per global deformation-mode equation (flat index, see NGlobalEq): this volume's free rigid-body
  /// DOFs take part in it. Set by the "globalModes" clause; see writeGlobalModeConstraints.
  std::array<bool, NGlobalEq> mGlobalModes{};
  /// per mode equation: one equation is written over the whole subtree of this volume, together with
  /// the other volumes of its group mGroupModesName (its own symbolic name if the clause gave none).
  /// Set by the "groupModes" clause; see writeGroupModeConstraints.
  std::array<bool, NGlobalEq> mGroupModes{};
  std::string mGroupModesName{};
  /// set by the "rayModes" clause: the sensors of this branch enter the ray-mode constraints with this
  /// binning (a flagged ancestor takes precedence); see writeRayModeConstraints
  std::shared_ptr<const RayModesConfig> mRayModes{}; //!
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

// Guards the hand-written NGlobalEq against any change to GlobalMode or globalModeNEq()
static_assert(Volume::NGlobalEq == Volume::globalEqFirst(Volume::NGlobalModes),
              "Volume::NGlobalEq is out of step with GlobalMode / globalModeNEq()");

} // namespace o2::alignment

#endif
