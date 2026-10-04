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
#include <fnmatch.h>
#include <cmath>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TMath.h>
#include <nlohmann/json.hpp>

#include "O2Align/Volume.h"
#include "O2Align/Constraint.h"
#include "O2Align/TimeSlotsSet.h"
#include "ITSBase/GeometryTGeo.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"

namespace o2::alignrs
{

Volume::Volume(const char* symName, uint32_t label, uint32_t det, bool sens, bool virt) : mSymName(symName), mLabel(det, label, sens), mVirtual(virt)
{
  init();
}

Volume::Volume(const char* symName, Label label, bool virt) : mSymName(symName), mLabel(label), mVirtual(virt)
{
  init();
}

void Volume::init()
{
  mCalibLabels.assign(1, mLabel.asCalib()); // a single calibration slot unless setCalibSlots says otherwise
  if (mVirtual) { // fictitious volume, has no counterpart in the geometry
    return;
  }
  // check if this sym volume actually exists
  mPNE = gGeoManager->GetAlignableEntry(mSymName.c_str());
  if (mPNE == nullptr) {
    LOGP(fatal, "Symbolic volume '{}' has no corresponding alignable entry!", mSymName);
  }
  mPN = mPNE->GetPhysicalNode();
  if (mPN == nullptr) {
    LOGP(debug, "Adding physical node to {}", mSymName);
    mPN = gGeoManager->MakePhysicalNode(mPNE->GetTitle());
    if (mPN == nullptr) {
      LOGP(fatal, "Failed to make physical node for {}", mSymName);
    }
  }
}

Volume::Ptr Volume::makeRoot(const char* symName)
{
  // the root belongs to no detector and is fictitious: its L2G is the identity of the global frame,
  // left as provided by the default-constructed mL2G, since defineMatrixL2G() is a no-op here
  auto root = std::make_unique<Volume>(symName, Label(Label::DET_GLOBAL, 0, false), true);
  root->setRigidBodyAllowed(false);
  return root;
}

const TGeoHMatrix& Volume::getMatrixL2G() const
{
  // sensors may redefine the L2G matrix and fictitious volumes have no geometry to take it from,
  // in both cases mL2G is filled by defineMatrixL2G()
  if (!isLeaf() && !mVirtual) {
    return *mPN->GetMatrix(); // global matrix (including possible pre-alignment)
  }
  return mL2G;
}

void Volume::computeJacobianL2T(const TGeoHMatrix& t2l, const double* posLoc, Matrix66& jac)
{
  jac.setZero();
  Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> rotT2L(t2l.GetRotationMatrix());
  Eigen::Matrix3d skew, rotL2T = rotT2L.transpose();
  skew << 0, -posLoc[2], posLoc[1], posLoc[2], 0, -posLoc[0], -posLoc[1], posLoc[0], 0;
  jac.topLeftCorner<3, 3>() = rotL2T;
  jac.topRightCorner<3, 3>() = -rotL2T * skew;
  jac.bottomRightCorner<3, 3>() = rotL2T;
}

TGeoHMatrix Volume::computeT2L(double alpha) const
{
  TGeoHMatrix t2l; // TRK -> GLO is the rotation by alpha
  t2l.RotateZ(alpha * TMath::RadToDeg());
  t2l.MultiplyLeft(mG2L); // GLO -> LOC
  return t2l;
}

void Volume::finalise(uint8_t level)
{
  if (level == 0 && !isRoot()) {
    LOGP(fatal, "Finalise should be called only from the root node!");
  }
  mLevel = level;
  if (mRigidBody && !mRBAllowed) {
    LOGP(fatal, "Volume {} is not rigid-body alignable but was assigned rigid-body DOFs", mSymName);
  }
  for (const auto* dofSet : {mRigidBody.get(), mCalib.get()}) {
    if (dofSet && dofSet->nDOFs() > static_cast<int>(Label::DOF_MAX) + 1) {
      LOGP(fatal, "Volume {} has a DOF set with {} DOFs, the Millepede label can encode at most {}", mSymName, dofSet->nDOFs(), Label::DOF_MAX + 1);
    }
  }
  if (isLeaf() && !mLabel.sens() && !mVirtual) {
    // being a leaf is what makes a volume define its own matrices instead of taking them from the
    // geometry, hence a childless volume which is neither a sensor nor fictitious is a dead branch:
    // its L2G would silently become the identity instead of the matrix of its alignable entry
    LOGP(fatal, "Volume {} is childless but neither a sensor nor fictitious: dead branch of the hierarchy", mSymName);
  }
  if (isLeaf()) {
    // for sensors we need also to define the transformation from the measurment (TRK) to the local frame (LOC)
    // need to it with including possible pre-alignment to allow for iterative convergence
    // (TRK) is defined wrt global z-axis
    defineMatrixL2G();
    mG2L = getMatrixL2G().Inverse();
    defineMatrixT2L();
  } else if (mVirtual) {
    // a fictitious intermediate volume must define its L2G before the children derive their L2P from it
    defineMatrixL2G();
  }
  if (!isRoot()) {
    // prepare the transformation matrices, e.g. from child frame to parent frame
    // this is not necessarily just one level transformation
    TGeoHMatrix mat = getMatrixL2G();
    auto inv = mParent->getMatrixL2G().Inverse(); // global (including possible pre-alignment) from the parent to the global frame
    mat.MultiplyLeft(inv);                        // left mult. effectively subtracts the parent transformation which is included in the the childs
    mL2P = mat;                                   // now this is directly the child to the parent transformation (LOC) (including possible pre-alignment)

    // prepare jacobian from child to parent frame
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> rotL2P(mL2P.GetRotationMatrix());
    Eigen::Matrix3d rotInv = rotL2P.transpose(); // parent-to-child rotation
    const double* t = mL2P.GetTranslation();     // child origin in parent frame
    Eigen::Matrix3d skewT;
    skewT << 0, -t[2], t[1], t[2], 0, -t[0], -t[1], t[0], 0;
    mJL2P.setZero();
    mJL2P.topLeftCorner<3, 3>() = rotInv;
    mJL2P.topRightCorner<3, 3>() = -rotInv * skewT;
    mJL2P.bottomRightCorner<3, 3>() = rotInv;
    mJP2L = mJL2P.inverse();
  }
  if (!isLeaf()) {
    // depth first
    for (const auto& c : mChildren) {
      c->finalise(level + 1);
    }
    // A free parent RB DOF with no active children is NOT auto-disabled: it is fully observable
    // from leaf residuals regardless of whether any descendant carries a DOFSet (buildPointGlobals
    // transports the parent-to-child jacobian unconditionally), and it only needs the mean-of-children
    // constraint (writeChildrenMeanConstraints) when a DOF is free at both this level and a child
    // level, which is the actual degeneracy that mechanism breaks. See doc/DOFConfig_rules.md.
    // if (mRigidBody) {
    //   int nActiveChildren = 0;
    //   for (const auto& c : mChildren) {
    //     if (c->isActive()) {
    //       ++nActiveChildren;
    //     }
    //   }
    //   if (!nActiveChildren) {
    //     for (int iDOF = 0; iDOF < mRigidBody->nDOFs(); ++iDOF) {
    //       if (mRigidBody->isFree(iDOF)) {
    //         LOGP(warn, "Auto-disabling DOF {} for {} since no active children",
    //              mRigidBody->dofName(iDOF), mSymName);
    //         mRigidBody->setFree(iDOF, false);
    //       }
    //     }
    //   }
    // }
  }
}

void Volume::writeRigidBodyConstraints(std::ostream& os) const
{
  // a volume which is not rigid-body alignable (the root of the hierarchy, a detector envelope)
  // imposes no constraint on its branches, but these may still constrain their own children
  if (!isLeaf() && mRBAllowed) {
    writeChildrenMeanConstraints(os);
  }
  for (const auto& c : mChildren) {
    c->writeRigidBodyConstraints(os);
  }
}

// Thesis A.83: for every free DOF of this volume, the mean of the child DOFs transported to this
// frame with J^-1 = getJP2L() must vanish. The same equation, which never refers to this volume's
// own label, also serves a second purpose: for a DOF flagged by mPinChildrenMean (the
// "pinChildrenMean" DOF config clause) but fixed here, or even w/o a rigid-body DOF set at all, it
// pins the common mode of the children, there being no parameter of this volume to absorb it.
void Volume::writeChildrenMeanConstraints(std::ostream& os) const
{
  const auto nActiveChildren = std::count_if(mChildren.begin(), mChildren.end(), [](const auto& c) { return c->isActive(); });
  for (int iDOF = 0; iDOF < RigidBodyDOFSet::NDOF; ++iDOF) {
    const bool free = mRigidBody && mRigidBody->isFree(iDOF);
    if (!free && !mPinChildrenMean[iDOF]) {
      continue;
    }
    const char* dofName = RigidBodyDOFSet::RigidBodyDOFNames[iDOF];
    // a free DOF with no active children needs no constraint: there is no degeneracy to break
    // (see finalise()) and the fall-through below (con.getSize() == 0, free == true) stays silent
    Constraint con(std::format("DOF {} for {}{}", dofName, mSymName, free ? "" : " (pinned)"), 0.);
    if (nActiveChildren > 0) {
      addChildrenMeanTerms(con, iDOF, 1.0 / static_cast<double>(nActiveChildren));
    }
    // a single term is a legitimate constraint as well: e.g. for a parent with a single active
    // child it removes the degeneracy of their common rotation
    if (con.getSize() > 0) {
      con.write(os);
    } else if (!free) { // w/o free child DOFs feeding it, a free parent DOF has no degeneracy to break
      LOGP(warn, "Ignoring pinChildrenMean of DOF {} for {}: no free child DOF contributes to it", dofName, mSymName);
    }
  }
}

// add to con the free child DOFs contributing to the DOF iDOF of this volume, scaled by weight
void Volume::addChildrenMeanTerms(Constraint& con, int iDOF, double weight) const
{
  for (const auto& c : mChildren) {
    if (!c->mRigidBody) {
      continue;
    }
    for (int jDOF = 0; jDOF < c->mRigidBody->nDOFs(); ++jDOF) {
      if (!c->mRigidBody->isFree(jDOF)) {
        continue;
      }
      const double coeff = weight * c->getJP2L()(iDOF, jDOF);
      if (std::abs(coeff) > 1e-16) {
        con.add(c->getLabel().raw(jDOF), coeff);
      }
    }
  }
}

void Volume::setCalibSlots(const TimeSlotsSet& slots)
{
  const auto slotIDs = slots.getSlotIDs();
  if (slotIDs.empty()) {
    LOGP(fatal, "Cannot declare an empty set of calibration slots for {}", mSymName);
  }
  mCalibLabels.clear();
  mCalibLabels.reserve(slotIDs.size());
  for (int slotID : slotIDs) {
    // the slot ID plays the role of the volume ID: the parameters of every slot are independent
    mCalibLabels.emplace_back(Label(mLabel.det(), static_cast<uint32_t>(slotID), mLabel.sens(), true));
  }
  mActiveCalibSlot = 0;
  mCalibSlots = &slots;
}

void Volume::setActiveCalibSlot(int slotID)
{
  const auto lbl = Label(mLabel.det(), static_cast<uint32_t>(slotID), mLabel.sens(), true);
  const auto it = std::find(mCalibLabels.begin(), mCalibLabels.end(), lbl);
  if (it == mCalibLabels.end()) {
    LOGP(fatal, "Calibration slot {} was not declared for {}", slotID, mSymName);
  }
  mActiveCalibSlot = static_cast<size_t>(std::distance(mCalibLabels.begin(), it));
}

void Volume::writeParameters(std::ostream& os) const
{
  if (isRoot()) {
    os << "Parameter\n";
  }
  if (!mIsPseudo) {
    if (mRigidBody) {
      for (int iDOF = 0; iDOF < mRigidBody->nDOFs(); ++iDOF) {
        os << std::format("{:<10} {:>+15g} {:>+15g} ! {} {} ",
                          mLabel.raw(iDOF), 0.0, (mRigidBody->isFree(iDOF) ? 0.0 : -1.0),
                          (mRigidBody->isFree(iDOF) ? 'V' : 'F'), mRigidBody->dofName(iDOF))
           << mSymName << '\n';
      }
    }
    if (mCalib) {
      // one independent set of parameters per calibration time slot
      for (const auto& calibLbl : mCalibLabels) {
        for (int iDOF = 0; iDOF < mCalib->nDOFs(); ++iDOF) {
          os << std::format("{:<10} {:>+15g} {:>+15g} ! {} {:<5} ",
                            calibLbl.raw(iDOF), 0.0, (mCalib->isFree(iDOF) ? 0.0 : -1.0),
                            (mCalib->isFree(iDOF) ? 'V' : 'F'), mCalib->dofName(iDOF))
             << mSymName;
          if (mCalibLabels.size() > 1) {
            os << std::format(" slot {}", calibLbl.id());
          }
          os << '\n';
        }
      }
    }
  }
  for (const auto& c : mChildren) {
    c->writeParameters(os);
  }
}

void Volume::writeTree(std::ostream& os, int indent) const
{
  os << std::string(static_cast<size_t>(indent * 2), ' ') << mSymName << (mLabel.sens() ? " (sens)" : " (pasv)");
  if (mIsPseudo) {
    os << " pseudo";
  } else {
    int nFreeDofs{0};
    if (mRigidBody && mRigidBody->nFreeDOFs()) {
      nFreeDofs += mRigidBody->nFreeDOFs();
      os << " RB[";
      for (int i = 0; i < mRigidBody->nDOFs(); ++i) {
        if (mRigidBody->isFree(i)) {
          os << " " << mRigidBody->dofName(i) << "(" << mLabel.raw(i) << ")";
        }
      }
      os << " ]";
    }
    if (mCalib && mCalib->nFreeDOFs()) {
      nFreeDofs += mCalib->nFreeDOFs();
      for (const auto& calibLbl : mCalibLabels) {
        os << " CAL";
        if (mCalibLabels.size() > 1) {
          os << "@" << calibLbl.id();
        }
        os << "[";
        for (int i = 0; i < mCalib->nDOFs(); ++i) {
          if (mCalib->isFree(i)) {
            os << " " << mCalib->dofName(i) << "(" << calibLbl.raw(i) << ")";
          }
        }
        os << " ]";
      }
    }
    if (!nFreeDofs) {
      os << " no DOFs";
    }
  }
  os << '\n';
  for (const auto& c : mChildren) {
    c->writeTree(os, indent + 2);
  }
}

// Apply the free/fixed part of a calibration rule, common to all the calibration DOF set types:
// "fixed" fixes the whole set, "free" fixes everything but the listed DOFs and "fix" fixes the
// listed ones. A DOF is named either by its index or by its name, e.g. "L(1,0)" or "VDRIFT".
void Volume::applyFreeFixConfig(DOFSet& dofSet, const nlohmann::json& cal)
{
  auto setFreeByKey = [&dofSet](const nlohmann::json& item, bool free) {
    if (item.is_number_integer()) {
      const int idx = item.get<int>();
      if (idx < 0 || idx >= dofSet.nDOFs()) {
        LOGP(warn, "Ignoring the DOF index {} of the calib rule: out of the [0, {}) range", idx, dofSet.nDOFs());
        return;
      }
      dofSet.setFree(idx, free);
      return;
    }
    if (!item.is_string()) {
      return;
    }
    const auto name = item.get<std::string>();
    bool found = false;
    for (int k = 0; k < dofSet.nDOFs(); ++k) {
      if (dofSet.dofName(k) == name) {
        dofSet.setFree(k, free);
        found = true;
      }
    }
    if (!found) {
      LOGP(warn, "Ignoring the DOF '{}' of the calib rule: no such DOF in this set", name);
    }
  };

  if (cal.value("fixed", false)) {
    dofSet.setAllFree(false);
  }
  if (cal.contains("free")) { // an explicit list of free DOFs fixes all the others
    dofSet.setAllFree(false);
    for (const auto& item : cal["free"]) {
      setFreeByKey(item, true);
    }
  }
  if (cal.contains("fix")) {
    for (const auto& item : cal["fix"]) {
      setFreeByKey(item, false);
    }
  }
}

// Build the rigid-body DOF set of a rule: "all"/"free" or "fixed", an array of the free DOF names,
// or an object {"dofs": "all" | [names], "fixed": bool}. Returns nullptr for an invalid clause.
std::unique_ptr<DOFSet> Volume::makeRigidBodyDOFSet(const nlohmann::json& rb, const std::string& pattern)
{
  auto dofSet = std::make_unique<RigidBodyDOFSet>();
  // free (or fix, if fixed) the named DOFs, all the others being fixed
  auto setListed = [&](const nlohmann::json& names, bool fixed) {
    dofSet->setAllFree(false);
    for (const auto& name : names) {
      const auto idx = RigidBodyDOFSet::dofIndex(name.get<std::string>());
      if (idx < 0) {
        LOGP(fatal, "Unknown rigid-body DOF '{}' in the rule '{}', allowed are TX,TY,TZ,RX,RY,RZ", name.get<std::string>(), pattern);
      }
      dofSet->setFree(idx, !fixed);
    }
  };
  if (rb.is_string()) {
    const auto s = rb.get<std::string>();
    if (s == "fixed") {
      dofSet->setAllFree(false);
    } else if (s != "all" && s != "free") {
      LOGP(fatal, "Unknown rigidBody value '{}' in the rule '{}', allowed are all, free, fixed or a list of DOFs", s, pattern);
    }
    return dofSet;
  }
  if (rb.is_array()) {
    setListed(rb, false);
    return dofSet;
  }
  if (rb.is_object()) {
    const bool fixed = rb.value("fixed", false);
    const auto dofs = rb.value("dofs", nlohmann::json("all"));
    if (dofs.is_array()) {
      setListed(dofs, false);
      if (fixed) { // the listed DOFs are declared but fixed
        dofSet->setAllFree(false);
      }
    } else if (dofs == "all") {
      dofSet->setAllFree(!fixed);
    } else {
      LOGP(fatal, "Invalid rigidBody.dofs in the rule '{}'", pattern);
    }
    return dofSet;
  }
  LOGP(fatal, "Invalid rigidBody clause in the rule '{}'", pattern);
  return nullptr;
}

// Apply the "pinChildrenMean" clause of a rule: a bool (all DOFs), "all", or an array of DOF names.
// Unlike rigidBody/calib this never creates a DOF set, it only flags DOFs of the volume's own
// (possibly absent or fixed) RigidBodyDOFSet for writeRigidBodyConstraints. A later matching rule
// replaces the flags of an earlier one, like the rigidBody/calib clauses.
void Volume::applyPinChildrenMeanConfig(const nlohmann::json& pin, const std::string& pattern)
{
  if (!isRigidBodyAllowed()) {
    if (pattern.find('*') == std::string::npos) {
      LOGP(warn, "Ignoring the pinChildrenMean rule '{}': {} is not rigid-body alignable", pattern, mSymName);
    }
    return;
  }
  mPinChildrenMean.fill(false);
  if (pin.is_boolean()) {
    mPinChildrenMean.fill(pin.get<bool>());
    return;
  }
  if (pin.is_string() && pin.get<std::string>() == "all") {
    mPinChildrenMean.fill(true);
    return;
  }
  if (!pin.is_array()) {
    LOGP(fatal, "Invalid pinChildrenMean clause in the rule '{}', allowed are a bool, \"all\" or a list of DOF names", pattern);
  }
  for (const auto& name : pin) {
    const auto idx = RigidBodyDOFSet::dofIndex(name.get<std::string>());
    if (idx < 0) {
      LOGP(fatal, "Unknown rigid-body DOF '{}' in the pinChildrenMean rule '{}', allowed are TX,TY,TZ,RX,RY,RZ", name.get<std::string>(), pattern);
    }
    mPinChildrenMean[idx] = true;
  }
}

// Build the calibration DOF set of a rule, with its free/fixed clauses applied
std::unique_ptr<DOFSet> Volume::makeCalibDOFSet(const nlohmann::json& cal, const std::string& pattern, const std::string& sym)
{
  const auto calType = cal.value("type", std::string(""));
  std::unique_ptr<DOFSet> dofSet;
  if (calType == "legendre") {
    dofSet = std::make_unique<LegendreDOFSet>(cal.value("order", 3));
  } else if (calType == "inextensional") {
    dofSet = std::make_unique<InextensionalDOFSet>(cal.value("order", 2));
  } else if (calType == "tpcvdrift") {
    dofSet = std::make_unique<TPCVDriftDOFSet>();
  } else if (calType == "meanvertex") {
    dofSet = std::make_unique<MeanVertexDOFSet>();
  } else {
    LOGP(warn, "Ignoring the calib rule '{}' of {}: unknown calibration type '{}'", pattern, sym, calType);
    return nullptr;
  }
  applyFreeFixConfig(*dofSet, cal);
  return dofSet;
}

void Volume::applyDOFConfig(Volume* root, const std::string& jsonPath)
{
  using json = nlohmann::json;
  std::ifstream f(jsonPath);
  if (!f.is_open()) {
    LOGP(fatal, "Cannot open DOF config file: {}", jsonPath);
  }
  auto data = json::parse(f);
  json rules = data.is_array() ? data : data.value("rules", json::array());
  if (data.is_object() && data.contains("defaults")) {
    json defRule = data["defaults"];
    defRule["match"] = "*";
    rules.insert(rules.begin(), defRule);
  }

  // the pattern is also tried with an implicit leading "*"
  auto matchPattern = [](const std::string& pattern, const std::string& sym) -> bool {
    return fnmatch(pattern.c_str(), sym.c_str(), 0) == 0 || fnmatch(("*" + pattern).c_str(), sym.c_str(), 0) == 0;
  };

  root->traverse([&](Volume* vol) {
    if (vol->isPseudo()) {
      return;
    }
    const std::string& sym = vol->getSymName();
    for (const auto& rule : rules) { // in order: a later matching rule replaces the DOF set of an earlier one
      const auto pattern = rule["match"].get<std::string>();
      if (!matchPattern(pattern, sym)) {
        continue;
      }
      if (rule.contains("rigidBody")) {
        // silently not applicable to the volumes which are not rigid-body alignable (the root of the
        // hierarchy, the envelopes w/o own geometry), for which only the calibration DOFs can be
        // configured. An explicitly named volume is likely a configuration mistake.
        if (vol->isRigidBodyAllowed()) {
          vol->setRigidBody(makeRigidBodyDOFSet(rule["rigidBody"], pattern));
        } else if (pattern.find('*') == std::string::npos) {
          LOGP(warn, "Ignoring the rigidBody rule '{}': {} is not rigid-body alignable", pattern, sym);
        }
      }
      if (rule.contains("calib")) {
        if (auto dofSet = makeCalibDOFSet(rule["calib"], pattern, sym)) {
          vol->setCalib(std::move(dofSet));
        }
      }
      if (rule.contains("pinChildrenMean")) {
        vol->applyPinChildrenMeanConfig(rule["pinChildrenMean"], pattern);
      }
    }
  });
}

Volume::InjectedMisalignment Volume::loadInjectedMisalignment(const std::string& injectedJsonPath)
{
  using json = nlohmann::json;
  InjectedMisalignment inj;
  if (injectedJsonPath.empty()) {
    return inj;
  }
  std::ifstream injFile(injectedJsonPath);
  if (!injFile.is_open()) {
    LOGP(warn, "Cannot open injected misalignment file: {}, writing absolute values", injectedJsonPath);
    return inj;
  }
  json injData = json::parse(injFile);
  for (const auto& item : injData) {
    int id = item["id"].get<int>();
    if (item.contains("rigidBody")) {
      inj.rigidBody[id] = item["rigidBody"].get<std::vector<double>>();
    }
    if (item.contains("matrix")) {
      inj.matrix[id] = item["matrix"].get<std::vector<std::vector<double>>>();
    }
    if (item.contains("inextensional")) {
      InjectedMisalignment::Inextensional ii;
      const auto& inex = item["inextensional"];
      if (inex.contains("modes")) {
        for (auto& [key, val] : inex["modes"].items()) {
          ii.modes[std::stoi(key)] = val.get<std::array<double, 4>>();
        }
      }
      if (inex.contains("alpha")) {
        ii.alpha = inex["alpha"].get<double>();
      }
      if (inex.contains("beta")) {
        ii.beta = inex["beta"].get<double>();
      }
      inj.inextensional[id] = ii;
    }
  }
  LOGP(info, "Loaded injected misalignment for {} sensors", injData.size());
  return inj;
}

// The caller (writeMillepedeResults) only invokes this when getRigidBody() has free DOFs.
bool Volume::MP2JSON_RB(const std::map<uint32_t, double>& labelToValue, const std::vector<double>* inj, nlohmann::json& entry) const
{
  auto rbArr = nlohmann::json::array();
  for (int i = 0; i < getRigidBody()->nDOFs(); ++i) {
    const double fitted = getFittedValue(labelToValue, mLabel.raw(i));
    const double ref = (inj && i < static_cast<int>(inj->size())) ? (*inj)[i] : 0.0;
    rbArr.push_back(fitted - ref);
  }
  entry["rigidBody"] = rbArr;
  return true;
}

// The caller (writeMillepedeResults) only invokes this when getCalib() has free DOFs.
bool Volume::MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const
{
  if (inj) {
    LOGP(warn, "The injected calibration of {} (calib type {}) is not subtracted, its fitted values are absolute", getSymName(), static_cast<int>(getCalib()->type()));
  }
  auto slots = calibSlotsToJSON(labelToValue);
  if (slots.empty()) {
    return false;
  }
  entry["calib"] = std::move(slots);
  return true;
}

nlohmann::json Volume::calibSlotsToJSON(const std::map<uint32_t, double>& labelToValue) const
{
  auto slotArr = nlohmann::json::array();
  for (const auto& lbl : mCalibLabels) {
    const int slotID = static_cast<int>(lbl.id());
    nlohmann::json rec;
    rec["slot"] = slotID;
    if (mCalibSlots) {
      const auto& slot = mCalibSlots->getSlotByID(slotID);
      rec["run"] = slot.runNumber;
      rec["tsS"] = slot.timeStampS;
      rec["tsE"] = slot.timeStampE;
    }
    bool anyFitted = false;
    for (int i = 0; i < mCalib->nDOFs(); ++i) {
      if (!mCalib->isFree(i)) {
        continue;
      }
      const auto it = labelToValue.find(lbl.raw(i));
      rec[mCalib->dofName(i)] = it != labelToValue.end() ? it->second : 0.0;
      anyFitted = anyFitted || it != labelToValue.end();
    }
    if (anyFitted) {
      slotArr.push_back(std::move(rec));
    }
  }
  return slotArr;
}

bool Volume::MP2ROOT_RB(const std::map<uint32_t, double>&, const std::vector<double>*, nlohmann::json&) const
{
  LOGP(warn, "MP2ROOT_RB is not implemented for volume {}", getSymName());
  return false;
}

bool Volume::MP2ROOT_Calib(const std::map<uint32_t, double>&, const InjectedMisalignment*, nlohmann::json&) const
{
  LOGP(warn, "MP2ROOT_Calib is not implemented for volume {} (calib type {})", getSymName(), static_cast<int>(getCalib()->type()));
  return false;
}

std::map<uint32_t, double> Volume::readMillepedeResults(const std::string& milleResPath)
{
  // parse millepede.res: label fittedValue presigma [...]
  std::ifstream fin(milleResPath);
  if (!fin.is_open()) {
    LOGP(fatal, "Cannot open millepede result file: {}", milleResPath);
  }
  std::map<uint32_t, double> labelToValue;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '!' || line[0] == '*') {
      continue;
    }
    if (line.find("Parameter") != std::string::npos) {
      continue;
    }
    std::istringstream iss(line);
    uint32_t label = 0;
    double value = NAN, presigma = NAN;
    if (!(iss >> label >> value >> presigma)) {
      continue;
    }
    if (presigma >= 0.0) { // skip fixed parameters
      labelToValue[label] = value;
    }
  }
  fin.close();
  LOGP(info, "Parsed {} not fixed parameters from {}", labelToValue.size(), milleResPath);
  return labelToValue;
}

void Volume::writeMillepedeResults(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outJsonPath, const std::string& injectedJsonPath)
{
  using json = nlohmann::json;

  // load injected misalignment if provided (same format as closure test input), indexed by sensorID.
  // Its use is optional: an empty injectedJsonPath yields an empty struct, and MP2JSON_RB/Calib then
  // write the absolute fitted values (nullptr inj pointers below).
  const InjectedMisalignment injMisal = loadInjectedMisalignment(injectedJsonPath);
  const bool haveInj = !injectedJsonPath.empty();

  // collect results per volume that has RB or calib DOFs
  json output = json::array();
  root->traverse([&](Volume* vol) {
    const int id = vol->getSensorId();
    json entry;
    entry["symName"] = vol->getSymName();
    entry["id"] = id;

    bool write = false;
    if (vol->getRigidBody() && vol->getRigidBody()->nFreeDOFs()) {
      const auto itRB = injMisal.rigidBody.find(id);
      const std::vector<double>* injRBPtr = (haveInj && itRB != injMisal.rigidBody.end()) ? &itRB->second : nullptr;
      write |= vol->MP2JSON_RB(labelToValue, injRBPtr, entry);
    }
    if (vol->getCalib() && vol->getCalib()->nFreeDOFs()) {
      const InjectedMisalignment* injPtr = haveInj ? &injMisal : nullptr;
      write |= vol->MP2JSON_Calib(labelToValue, injPtr, entry);
    }

    if (write) {
      output.push_back(entry);
    }
  });

  std::ofstream fout(outJsonPath);
  if (!fout.is_open()) {
    LOGP(fatal, "Cannot open output file: {}", outJsonPath);
  }
  fout << output.dump(2) << '\n';
  fout.close();
  LOGP(info, "Wrote millepede results to {}", outJsonPath);
}

} // namespace o2::alignrs
