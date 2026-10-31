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
#include <regex>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoBBox.h>
#include <TMath.h>
#include <nlohmann/json.hpp>

#include "O2Align/Volume.h"
#include "O2Align/AlignGeomUtils.h"
#include "O2Align/Constraint.h"
#include "O2Align/TimeSlotsSet.h"
#include "ITSBase/GeometryTGeo.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

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
    // constraint (writeChildrenMeanConstraints) when a DOF is free at both this level and a descendant
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

bool Volume::hasAnyWeightedMean() const
{
  for (auto f : mWeightedMean) {
    if (f == WeightedMeanFlag::True) {
      return true;
    }
  }
  for (const auto& c : mChildren) {
    if (c->hasAnyWeightedMean()) {
      return true;
    }
  }
  return false;
}

void Volume::writeRigidBodyConstraints(std::ostream& os, bool statWeighted, const ModeEquationOptions& modeOpt) const
{
  if (isRoot()) {
    // a global deformation mode spans the whole tree, so its equation is written once here and not
    // during the per-volume recursion below
    writeGlobalModeConstraints(os, modeOpt);
    writeRayModeConstraints(os);
  }
  // a volume which is not rigid-body alignable (the root of the hierarchy, a detector envelope)
  // imposes no constraint on its branches, but these may still constrain their own children
  if (!isLeaf() && mRBAllowed) {
    writeChildrenMeanConstraints(os, statWeighted);
  }
  writeRigidBodyDOFMeasurements(os); // a sensor may be measured too
  writeCalibDOFMeasurements(os);     // a calibration DOF (e.g. the mean vertex) may be measured too
  for (const auto& c : mChildren) {
    c->writeRigidBodyConstraints(os, statWeighted);
  }
}

// Thesis A.83: for every free DOF of this volume, the mean of the child DOFs transported to this
// frame with J^-1 = getJP2L() must vanish. A child w/o a free RB DOF is looked through, its own
// children (with the jacobians chained) standing in for it, since the degeneracy is then with them.
// Limitation: the walk stops at the first descendant with any free DOF, so the deeper free DOFs
// feeding a DOF fixed at that level are not included. The same equation, which never refers to this volume's
// own label, also serves a second purpose: for a DOF flagged by mPinChildrenMean (the
// "pinChildrenMean" DOF config clause) but fixed here, or even w/o a rigid-body DOF set at all, it
// pins the common mode of the children, there being no parameter of this volume to absorb it.
// A DOF flagged by mPinChildrenMeanNoCross ("pinChildrenMeanNoCross") gets, in addition, the
// cross-term-free variant of the same equation -- see addChildrenMeanTerms for why that matters.
void Volume::writeChildrenMeanConstraints(std::ostream& os, bool statWeighted) const
{
  // the degeneracy is with the nearest free level below, which needs not be the direct children
  std::vector<MeanContributor> contributors;
  collectMeanContributors(Matrix66::Identity(), contributors);
  const size_t n = contributors.size();
  const std::vector<double> equalWeights(n, n ? 1.0 / static_cast<double>(n) : 0.);
  // the per-DOF choice (statWeighted default or a "weightedMean" override) only selects between this
  // and equalWeights below: the underlying statistics do not depend on iDOF, so they are computed once
  std::vector<double> statWeights;
  bool haveStatWeights = false;
  if (n) {
    statWeights.resize(n);
    double sum = 0.;
    for (size_t i = 0; i < n; ++i) {
      sum += (statWeights[i] = static_cast<double>(contributors[i].vol->getDataCounter()));
    }
    if (sum > 0.) {
      for (auto& w : statWeights) {
        w /= sum;
      }
      haveStatWeights = true;
    }
  }
  for (int iDOF = 0; iDOF < RigidBodyDOFSet::NDOF; ++iDOF) {
    const bool free = mRigidBody && mRigidBody->isFree(iDOF);
    const char* dofName = RigidBodyDOFSet::RigidBodyDOFNames[iDOF];
    const bool wantStat = mWeightedMean[iDOF] == WeightedMeanFlag::Unset ? statWeighted : (mWeightedMean[iDOF] == WeightedMeanFlag::True);
    if (wantStat && n && !haveStatWeights) {
      LOGP(warn, "weightedMean of DOF {} for {}: every contributor has zero accumulated statistics, falling back to an equal split", dofName, mSymName);
    }
    const std::vector<double>& weights = (wantStat && haveStatWeights) ? statWeights : equalWeights;
    // the transported (cross-term) equation, for a free or pinChildrenMean-flagged DOF
    if (free || mPinChildrenMean[iDOF]) {
      // w/o a free descendant DOF realising it, a free parent DOF has no degeneracy to break. A merely
      // nonzero coupling is not enough: e.g. for OB chips with only TX free under a module with TX,TY,TZ
      // free, the module TY/TZ equations would be the TX one scaled by the chip/module misrotation,
      // making the constraint matrix singular.
      if (maxChildrenCoupling(iDOF, contributors) < MinChildrenCoupling) {
        if (free) {
          LOGP(debug, "No mean constraint of DOF {} for {}: no free descendant DOF realises it", dofName, mSymName);
        } else {
          LOGP(warn, "Ignoring pinChildrenMean of DOF {} for {}: no free descendant DOF realises it", dofName, mSymName);
        }
      } else {
        // a single term is a legitimate constraint as well: e.g. for a parent with a single active
        // child it removes the degeneracy of their common rotation
        Constraint con(std::format("DOF {} for {}{}", dofName, mSymName, free ? "" : " (pinned)"), 0.);
        addChildrenMeanTerms(con, iDOF, contributors, weights, false);
        con.write(os);
      }
    }
    // the cross-term-free equation, an independent one: the plain mean of the children's own iDOF
    if (mPinChildrenMeanNoCross[iDOF]) {
      Constraint con(std::format("DOF {} for {} (pinned, no cross-terms)", dofName, mSymName), 0.);
      if (!contributors.empty()) {
        addChildrenMeanTerms(con, iDOF, contributors, weights, true);
      }
      if (con.getSize() > 0) {
        con.write(os);
      } else {
        LOGP(warn, "Ignoring pinChildrenMeanNoCross of DOF {} for {}: no descendant has that DOF free", dofName, mSymName);
      }
    }
  }
}

void Volume::writeRigidBodyDOFMeasurements(std::ostream& os) const
{
  if (!mRigidBody) {
    return;
  }
  for (int iDOF = 0; iDOF < mRigidBody->nDOFs(); ++iDOF) {
    const auto& m = mRigidBody->getMeasurement(iDOF);
    if (m.sigma <= 0.) {
      continue;
    }
    const auto dofName = mRigidBody->dofName(iDOF);
    if (!mRigidBody->isFree(iDOF)) { // the label is not a parameter of the fit
      LOGP(warn, "Ignoring the measurement of DOF {} for {}: the DOF is not free", dofName, mSymName);
      continue;
    }
    Constraint con(std::format("Measurement of DOF {} for {}", dofName, mSymName), m.value, m.sigma);
    con.add(mLabel.raw(iDOF), 1.);
    con.write(os);
  }
}

void Volume::writeCalibDOFMeasurements(std::ostream& os) const
{
  if (!mCalib) {
    return;
  }
  for (int iDOF = 0; iDOF < mCalib->nDOFs(); ++iDOF) {
    const auto& m = mCalib->getMeasurement(iDOF);
    if (m.sigma <= 0.) {
      continue;
    }
    const auto dofName = mCalib->dofName(iDOF);
    if (!mCalib->isFree(iDOF)) { // the label is not a parameter of the fit
      LOGP(warn, "Ignoring the measurement of DOF {} for {}: the DOF is not free", dofName, mSymName);
      continue;
    }
    // the same external measurement applies to every calibration slot of this volume
    for (const auto& calibLbl : mCalibLabels) {
      Constraint con(std::format("Measurement of DOF {} for {}{}", dofName, mSymName,
                                  mCalibLabels.size() > 1 ? std::format(" slot {}", calibLbl.id()) : ""),
                      m.value, m.sigma);
      con.add(calibLbl.raw(iDOF), 1.);
      con.write(os);
    }
  }
}

void Volume::collectMeanContributors(const Matrix66& jToThis, std::vector<MeanContributor>& out) const
{
  for (const auto& c : mChildren) {
    if (c->mIsPseudo) { // pseudo volumes (ITS3 tiles) carry no RB DOFs and have no children
      continue;
    }
    const Matrix66 j = jToThis * c->getJP2L();
    if (c->mRigidBody && c->mRigidBody->nFreeDOFs() > 0) {
      out.push_back({c.get(), j});
    } else {
      c->collectMeanContributors(j, out);
    }
  }
}

void Volume::addChildrenMeanTerms(Constraint& con, int iDOF, const std::vector<MeanContributor>& contributors, const std::vector<double>& weights, bool noCross)
{
  for (size_t ic = 0; ic < contributors.size(); ++ic) {
    const auto& [vol, j] = contributors[ic];
    const double weight = weights[ic];
    if (noCross) {
      // only the contributor's own iDOF, with no transport: the equation then says that the plain mean
      // of that one DOF vanishes, which no other DOF of the contributor can pay for
      if (iDOF < vol->mRigidBody->nDOFs() && vol->mRigidBody->isFree(iDOF)) {
        con.add(vol->getLabel().raw(iDOF), weight);
      }
      continue;
    }
    for (int jDOF = 0; jDOF < vol->mRigidBody->nDOFs(); ++jDOF) {
      if (!vol->mRigidBody->isFree(jDOF)) {
        continue;
      }
      const double coeff = weight * j(iDOF, jDOF);
      if (std::abs(coeff) > 1e-16) {
        con.add(vol->getLabel().raw(jDOF), coeff);
      }
    }
  }
}

double Volume::maxChildrenCoupling(int iDOF, const std::vector<MeanContributor>& contributors)
{
  double maxCoupling = 0.;
  for (const auto& [vol, j] : contributors) {
    for (int jDOF = 0; jDOF < vol->mRigidBody->nDOFs(); ++jDOF) {
      if (vol->mRigidBody->isFree(jDOF)) {
        maxCoupling = std::max(maxCoupling, std::abs(j(iDOF, jDOF)));
      }
    }
  }
  return maxCoupling;
}

// Displacement field u(x) of a global deformation mode of unit amplitude, in the GLOBAL frame,
// together with the local rotation omega(x) = 0.5*curl(u) that the field induces at that point.
//
// A volume sitting at x realises the mode by translating by u(x) *and* rotating by omega(x): the
// rotation part is not optional decoration, it is the first-order content of the field's gradient
// and it is what makes e.g. the twist constrain the staves' RX/RY and not only their TX/TY.
// Projecting u and omega on the volume's own local axes gives the coefficients of its 6 rigid-body
// DOFs in the mode equation (see collectGlobalModeTerms).
//
// Cylindrical coordinates, z along the beam: r = sqrt(x^2+y^2), rHat = (cos phi, sin phi, 0),
// phiHat = (-sin phi, cos phi, 0). Each case below states how omega follows from u.
// Returns false if the field is not defined at x (on the beam axis, where rHat and the 1/r rotation
// terms of the non-axisymmetric modes are meaningless); such a volume is then simply left out.
bool Volume::globalModeField(int mode, int phase, const double* pos, double* u, double* om)
{
  constexpr double rMin = 1e-4; // cm
  const double x = pos[0], y = pos[1], z = pos[2];
  const double r = std::sqrt(x * x + y * y);
  u[0] = u[1] = u[2] = om[0] = om[1] = om[2] = 0.;
  switch (mode) {
    case Volume::RadialScaling:
      // u = r*rHat = (x, y, 0) = grad(r^2/2): a gradient field, hence curl-free -- a pure dilation
      // carries no rotation at all, which is exactly why it is so nearly invisible to a helix fit
      u[0] = x;
      u[1] = y;
      break;
    case Volume::ZScaling:
      // u = z*zHat = (0, 0, z) = grad(z^2/2): curl-free as well
      u[2] = z;
      break;
    case Volume::Curl:
      // u = r*phiHat = (-y, x, 0), the rigid rotation about the beam: curl u = (0, 0, 2), om = zHat.
      // Kept available although the top-level pin of RZ normally removes it already
      u[0] = -y;
      u[1] = x;
      om[2] = 1.;
      break;
    case Volume::Twist:
      // u = z*r*phiHat = (-y*z, x*z, 0): curl u = (-x, -y, 2z), so om = (-x/2, -y/2, z).
      // The rotation grows linearly with z, i.e. each end of the barrel is rotated the opposite way
      u[0] = -y * z;
      u[1] = x * z;
      om[0] = -0.5 * x;
      om[1] = -0.5 * y;
      om[2] = z;
      break;
    case Volume::Elliptical: {
      // u_r = cos(2phi) (phase 0) or sin(2phi) (phase 1), u_phi = u_z = 0: each layer is squashed
      // into an ellipse, the two phases covering the two orientations of its axes.
      // In cylindrical coordinates only curl_z survives, curl_z = -(1/r) d(u_r)/dphi, hence
      // om = +sin(2phi)/r * zHat for phase 0 and -cos(2phi)/r * zHat for phase 1
      if (r < rMin) {
        return false;
      }
      const double c = x / r, s = y / r;                 // cos phi, sin phi
      const double c2 = c * c - s * s, s2 = 2. * c * s;  // cos 2phi, sin 2phi
      const double ur = (phase == 0) ? c2 : s2;
      u[0] = ur * c;
      u[1] = ur * s;
      om[2] = ((phase == 0) ? s2 : -c2) / r;
      break;
    }
    case Volume::Bowing: {
      // u = z^2*rHat: the barrel bulges radially towards its ends. Only curl_phi survives,
      // curl_phi = d(u_r)/dz = 2z, hence om = z*phiHat = z*(-sin phi, cos phi, 0)
      if (r < rMin) {
        return false;
      }
      const double c = x / r, s = y / r;
      u[0] = z * z * c;
      u[1] = z * z * s;
      om[0] = -z * s;
      om[1] = z * c;
      break;
    }
    case Volume::RadialHarmonic:
    case Volume::TangentialHarmonic: {
      // The generic azimuthal-harmonic family, amplitude proportional to r:
      //   RadialHarmonic:     u_r   = r*C(phi),  u_phi = 0
      //   TangentialHarmonic: u_phi = r*C(phi),  u_r   = 0
      // with C = cos(l*phi) or sin(l*phi), l = globalHarmonicL(phase) and the cos/sin selected by
      // globalHarmonicIsSin(phase). The amplitude is proportional to r (the field is r*C, not C),
      // the shape in which these modes are coherent across detector layers.
      //
      // Only curl_z survives, curl_z = (1/r)[d(r*u_phi)/dr - d(u_r)/dphi], hence
      //   RadialHarmonic:     om_z = -C'(phi)/2   (C' = dC/dphi)
      //   TangentialHarmonic: om_z = +C(phi)
      // Both are finite on the axis, unlike Elliptical whose r-independent field needs om ~ 1/r.
      // At l = 0 the two families reduce exactly to RadialScaling and Curl, which is why the
      // harmonic orders start at l = 1.
      if (r < rMin) {
        return false;
      }
      const int l = globalHarmonicL(phase);
      const double c = x / r, s = y / r; // cos phi, sin phi
      const double lphi = l * std::atan2(s, c);
      const double C = globalHarmonicIsSin(phase) ? std::sin(lphi) : std::cos(lphi);
      const double dC = l * (globalHarmonicIsSin(phase) ? std::cos(lphi) : -std::sin(lphi));
      if (mode == Volume::RadialHarmonic) {
        u[0] = r * C * c;
        u[1] = r * C * s;
        om[2] = -0.5 * dC;
      } else {
        u[0] = -r * C * s;
        u[1] = r * C * c;
        om[2] = C;
      }
      break;
    }
    case Volume::GlobalTilt: {
      // Rigid rotation of the whole branch about the global X (phase 0) or Y (phase 1) axis:
      // u = e x r, om = e. These are the two transverse partners of Curl, which is the same thing
      // about Z, and they are needed as explicit constraints because no mean-of-children pin can
      // remove them: a tilt about Y shows up as a stave TZ pattern ~ -r*cos(phi) (and a stave RX
      // pattern ~ cos(phi)), whose mean vanishes over a half-shell *and* over the full circle, so
      // the cross-term pin that would see it is instead paid off between the two.
      const double e[3] = {phase == 0 ? 1. : 0., phase == 0 ? 0. : 1., 0.};
      u[0] = e[1] * z - e[2] * y;
      u[1] = e[2] * x - e[0] * z;
      u[2] = e[0] * y - e[1] * x;
      om[0] = e[0];
      om[1] = e[1];
      om[2] = e[2];
      break;
    }
    case Volume::RadialPolar:
    case Volume::TangentialPolar:
    case Volume::LongitudinalPolar: {
      // Axisymmetric fields of amplitude f = r^p * t^k, t = z/r (the tgLambda of a track from the
      // nominal vertex), along rHat, phiHat or zHat; see isPolarDuplicate for the members left out.
      // With dt/dz = 1/r and dt/dr = -t/r, the cylindrical curl of each field is
      //   radial:       curl_phi = d(u_r)/dz                         =  k     r^(p-1) t^(k-1)
      //   tangential:   curl_r   = -d(u_phi)/dz                      = -k     r^(p-1) t^(k-1)
      //                 curl_z   = (1/r) d(r u_phi)/dr               = (p+1-k) r^(p-1) t^k
      //   longitudinal: curl_phi = -d(u_z)/dr                        = (k-p)  r^(p-1) t^k
      // and om = curl/2. Checks: tangential (1,0) is Curl (om_z = 1), longitudinal (1,0) the telescope
      // dz = r (om = -phiHat/2), longitudinal (1,1) = z is curl-free as ZScaling.
      if (r < rMin) {
        return false;
      }
      const auto [p, k] = polarPK(mode, phase);
      if (p < 0) {
        return false;
      }
      const double c = x / r, s = y / r, t = z / r;
      const double rp1 = p == 1 ? 1. : 1. / r; // r^(p-1)
      const double tk = std::pow(t, k), tk1 = k > 0 ? std::pow(t, k - 1) : 0.;
      const double f = (p == 1 ? r : 1.) * tk;
      double omR = 0., omPhi = 0., omZ = 0.;
      if (mode == Volume::RadialPolar) {
        u[0] = f * c;
        u[1] = f * s;
        omPhi = 0.5 * k * rp1 * tk1;
      } else if (mode == Volume::TangentialPolar) {
        u[0] = -f * s;
        u[1] = f * c;
        omR = -0.5 * k * rp1 * tk1;
        omZ = 0.5 * (p + 1 - k) * rp1 * tk;
      } else {
        u[2] = f;
        omPhi = 0.5 * (k - p) * rp1 * tk;
      }
      om[0] = omR * c - omPhi * s;
      om[1] = omR * s + omPhi * c;
      om[2] = omZ;
      break;
    }
    case Volume::GlobalShift:
      // Rigid translation along the global X, Y or Z axis (phase 0, 1, 2): u = e, om = 0. The
      // translation partners of GlobalTilt. Over a whole detector the top-level pins remove it
      // already; it is meant for "groupModes", where it pins the net shift of a part of the
      // detector at every level below it, which no mean-of-children pin can (their contributor
      // walk stops at the nearest free descendants)
      u[phase] = 1.;
      break;
    default:
      return false;
  }
  return true;
}

int Volume::globalEqMode(int eq)
{
  for (int m = 0; m < NGlobalModes; ++m) {
    if (eq < globalEqFirst(m) + globalModeNEq(m)) {
      return m;
    }
  }
  return -1;
}

int Volume::globalEqPhase(int eq)
{
  const int m = globalEqMode(eq);
  return m < 0 ? -1 : eq - globalEqFirst(m);
}

std::string Volume::globalEqName(int eq)
{
  const int m = globalEqMode(eq), ph = globalEqPhase(eq);
  if (m < 0) {
    return "?";
  }
  const std::string base = globalModeName(m);
  if (m == Elliptical) {
    return base + (ph == 0 ? "_cos2phi" : "_sin2phi");
  }
  if (m == GlobalTilt || m == GlobalShift) {
    return base + (ph == 0 ? "_x" : (ph == 1 ? "_y" : "_z"));
  }
  if (m == RadialHarmonic || m == TangentialHarmonic) {
    const int l = globalHarmonicL(ph);
    return std::format("{}{}_{}{}phi", base, l, globalHarmonicIsSin(ph) ? "sin" : "cos", l);
  }
  if (isPolarMode(m)) {
    const auto [p, k] = polarPK(m, ph);
    return std::format("{}_r{}t{}", base, p, k); // amplitude r^p * (z/r)^k
  }
  return base;
}

const char* Volume::globalModeName(int mode)
{
  static constexpr const char* names[NGlobalModes] = {"radialScaling", "zScaling", "curl", "twist", "elliptical", "bowing", "radialHarmonic", "tangentialHarmonic", "globalTilt",
                                                      "radialPolar", "tangentialPolar", "longitudinalPolar", "globalShift"};
  return (mode >= 0 && mode < NGlobalModes) ? names[mode] : "?";
}

void Volume::collectGlobalModeTerms(int eq, std::vector<std::pair<uint32_t, double>>& terms, int& nFlagged, bool allInSubtree) const
{
  const int mode = globalEqMode(eq), phase = globalEqPhase(eq);
  // only the volumes explicitly flagged for this equation contribute (every volume of the subtree for
  // a group equation), and only through DOFs which are actually free: a fixed DOF is not a parameter
  // of the fit and must not appear in a constraint
  if ((allInSubtree || mGlobalModes[eq]) && !mIsPseudo && mRigidBody && mRigidBody->nFreeDOFs() > 0) {
    ++nFlagged;
    const auto& l2g = getMatrixL2G();
    const double* pos = l2g.GetTranslation();    // this volume's origin in the global frame
    const double* rot = l2g.GetRotationMatrix(); // row-major 3x3; column k = local axis k in GLO
    double u[3], om[3];
    if (globalModeField(mode, phase, pos, u, om)) {
      const int nDOF = std::min(mRigidBody->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
      for (int iDOF = 0; iDOF < nDOF; ++iDOF) {
        if (!mRigidBody->isFree(iDOF)) {
          continue;
        }
        // The 6 local rigid-body DOFs are mutually orthogonal unit directions of the parameter
        // space (3 translations along orthonormal axes, 3 rotations about the same axes), so the
        // component of the mode along a DOF is simply the projection of u (for TX,TY,TZ) or of
        // omega (for RX,RY,RZ) on the corresponding local axis, i.e. on that column of L2G.
        // For an ITS sensor local Y is the sensor normal (~ radial), so for radialScaling this
        // reduces to r*(nHat.rHat) ~ r on TY and to ~0 on the in-plane TX, TZ.
        const double* v = (iDOF < 3) ? u : om;
        const int ax = iDOF % 3;
        // kept unconditionally: a term which is only numerically nonzero (the matrices of volumes
        // that are exactly axis-aligned leave ~1e-16 residues) is dropped later, relative to the
        // largest coefficient of the equation, since there is no meaningful absolute scale here
        terms.emplace_back(mLabel.raw(iDOF), v[0] * rot[ax] + v[1] * rot[3 + ax] + v[2] * rot[6 + ax]);
      }
    }
  }
  for (const auto& c : mChildren) {
    c->collectGlobalModeTerms(eq, terms, nFlagged, allInSubtree);
  }
}

Volume::ModeField Volume::modeField(int eq)
{
  const int mode = globalEqMode(eq), phase = globalEqPhase(eq);
  return [mode, phase](const double* pos, double* u) {
    double om[3];
    if (!globalModeField(mode, phase, pos, u, om)) {
      u[0] = u[1] = u[2] = 0.;
    }
  };
}

void Volume::collectSensorSamplePoints(std::vector<std::array<double, 3>>& pts) const
{
  if (!isLeaf()) {
    for (const auto& c : mChildren) {
      c->collectSensorSamplePoints(pts);
    }
    return;
  }
  // the box of the geometry volume, centred on the origin of the leaf's own L2G frame: a sensor may
  // define its L2G differently from its geometry volume (the ITS chip vs its sensitive layer), but
  // only the extent matters here, which serves to weight the rotations of the DOFs by a lever arm
  const auto& l2g = getMatrixL2G();
  const auto* box = mPN ? dynamic_cast<const TGeoBBox*>(mPN->GetShape()) : nullptr;
  if (box == nullptr) {
    const double* o = l2g.GetTranslation();
    pts.push_back({o[0], o[1], o[2]});
    return;
  }
  for (const double sx : {-1., 1.}) {
    for (const double sy : {-1., 1.}) {
      for (const double sz : {-1., 1.}) {
        const double loc[3] = {sx * box->GetDX(), sy * box->GetDY(), sz * box->GetDZ()};
        auto& g = pts.emplace_back();
        l2g.LocalToMaster(loc, g.data());
      }
    }
  }
}

// A unit value of the local DOF k of a volume with origin O and local axes e_k (columns of L2G)
// moves a point p by d(p) = e_k for a translation and by e_k x (p - O) for a rotation (see
// computeJacobianL2T: dX = t + omega x p in the local frame). Its coefficient sum_p u(p).d(p) is then
// U.e_k, resp. (M - O x U).e_k, with U = sum u(p) and M = sum p x u(p), sums which add up over the
// subtrees: one pass over the sensors gives the coefficients of every level.
std::array<double, 6> Volume::collectProjectedModeTerms(const ModeField& field, std::vector<std::pair<uint32_t, double>>& terms, int& nFlagged, bool allInSubtree, int eq) const
{
  std::array<double, 6> sums{}; // U, M
  if (isLeaf()) {
    std::vector<std::array<double, 3>> pts;
    collectSensorSamplePoints(pts);
    for (const auto& p : pts) {
      double u[3];
      field(p.data(), u);
      for (int i = 0; i < 3; ++i) {
        sums[i] += u[i];
        sums[3 + i] += p[(i + 1) % 3] * u[(i + 2) % 3] - p[(i + 2) % 3] * u[(i + 1) % 3];
      }
    }
  }
  for (const auto& c : mChildren) {
    const auto cs = c->collectProjectedModeTerms(field, terms, nFlagged, allInSubtree, eq);
    for (int i = 0; i < 6; ++i) {
      sums[i] += cs[i];
    }
  }
  // same selection as collectGlobalModeTerms
  if ((allInSubtree || mGlobalModes[eq]) && !mIsPseudo && mRigidBody && mRigidBody->nFreeDOFs() > 0) {
    ++nFlagged;
    const auto& l2g = getMatrixL2G();
    const double* o = l2g.GetTranslation();
    const double* rot = l2g.GetRotationMatrix(); // row-major 3x3; column k = local axis k in GLO
    double r[3];                                   // M - O x U
    for (int i = 0; i < 3; ++i) {
      r[i] = sums[3 + i] - (o[(i + 1) % 3] * sums[(i + 2) % 3] - o[(i + 2) % 3] * sums[(i + 1) % 3]);
    }
    const int nDOF = std::min(mRigidBody->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
    for (int iDOF = 0; iDOF < nDOF; ++iDOF) {
      if (mRigidBody->isFree(iDOF)) {
        const double* v = (iDOF < 3) ? sums.data() : r;
        const int ax = iDOF % 3;
        terms.emplace_back(mLabel.raw(iDOF), v[0] * rot[ax] + v[1] * rot[3 + ax] + v[2] * rot[6 + ax]);
      }
    }
  }
  return sums;
}

// u'(p) = u(p) - sum_j beta_j g_j(p), g_j the displacement field of the free DOF j of this volume and
// beta the least-squares solution of G beta = b, G_ij = sum_p g_i.g_j, b_j = sum_p u.g_j. By
// construction sum_p u'.g_j = 0: the free DOFs of this volume drop out of the projected equation.
Volume::ModeField Volume::removeFreeRigidMotion(const ModeField& field, std::array<double, 2>& norm2) const
{
  std::vector<std::array<double, 3>> pts;
  collectSensorSamplePoints(pts);
  std::vector<int> dofs;
  if (mRigidBody && !mIsPseudo) {
    const int nDOF = std::min(mRigidBody->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
    for (int iDOF = 0; iDOF < nDOF; ++iDOF) {
      if (mRigidBody->isFree(iDOF)) {
        dofs.push_back(iDOF);
      }
    }
  }
  const auto& l2g = getMatrixL2G();
  const double* o0 = l2g.GetTranslation();
  const double* rot = l2g.GetRotationMatrix();
  const std::array<double, 3> o{o0[0], o0[1], o0[2]};
  std::array<std::array<double, 3>, 3> axes; // local axes in GLO
  for (int k = 0; k < 3; ++k) {
    axes[k] = {rot[k], rot[3 + k], rot[6 + k]};
  }
  auto dofField = [dofs, axes, o](int j, const double* p, double* g) {
    const auto& e = axes[dofs[j] % 3];
    if (dofs[j] < 3) {
      g[0] = e[0], g[1] = e[1], g[2] = e[2];
    } else {
      const double d[3] = {p[0] - o[0], p[1] - o[1], p[2] - o[2]};
      for (int i = 0; i < 3; ++i) {
        g[i] = e[(i + 1) % 3] * d[(i + 2) % 3] - e[(i + 2) % 3] * d[(i + 1) % 3];
      }
    }
  };
  const int n = static_cast<int>(dofs.size());
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
  std::vector<std::array<double, 3>> g(n);
  for (const auto& p : pts) {
    double u[3];
    field(p.data(), u);
    for (int j = 0; j < n; ++j) {
      dofField(j, p.data(), g[j].data());
      b(j) += u[0] * g[j][0] + u[1] * g[j][1] + u[2] * g[j][2];
      for (int i = 0; i <= j; ++i) {
        G(i, j) += g[i][0] * g[j][0] + g[i][1] * g[j][1] + g[i][2] * g[j][2];
      }
    }
  }
  G = G.selfadjointView<Eigen::Upper>();
  // a complete orthogonal decomposition rather than a plain inverse: a volume with all its sensors on
  // a line (or a single point-like one) leaves some of its rotations without a lever arm
  const Eigen::VectorXd beta = n ? Eigen::VectorXd(G.completeOrthogonalDecomposition().solve(b)) : Eigen::VectorXd();
  ModeField reduced = [field, dofField, beta, n](const double* p, double* u) {
    field(p, u);
    double gj[3];
    for (int j = 0; j < n; ++j) {
      dofField(j, p, gj);
      for (int i = 0; i < 3; ++i) {
        u[i] -= beta(j) * gj[i];
      }
    }
  };
  norm2 = {0., 0.};
  for (const auto& p : pts) {
    double u[3], ur[3];
    field(p.data(), u);
    reduced(p.data(), ur);
    norm2[0] += u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    norm2[1] += ur[0] * ur[0] + ur[1] * ur[1] + ur[2] * ur[2];
  }
  return reduced;
}

void Volume::collectGroupModeRoots(std::map<std::string, std::vector<const Volume*>>& groups) const
{
  if (std::any_of(mGroupModes.begin(), mGroupModes.end(), [](bool b) { return b; })) {
    groups[mGroupModesName.empty() ? mSymName : mGroupModesName].push_back(this);
  }
  for (const auto& c : mChildren) {
    c->collectGroupModeRoots(groups);
  }
}

uint32_t Volume::minFreeRigidBodyLabel() const
{
  uint32_t lmin = std::numeric_limits<uint32_t>::max();
  if (mRigidBody) {
    const int nDOF = std::min(mRigidBody->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
    for (int iDOF = 0; iDOF < nDOF; ++iDOF) {
      if (mRigidBody->isFree(iDOF)) {
        lmin = std::min(lmin, mLabel.raw(iDOF));
      }
    }
  }
  for (const auto& c : mChildren) {
    lmin = std::min(lmin, c->minFreeRigidBodyLabel());
  }
  return lmin;
}

void Volume::addPedeGroupingAnchor(Constraint& con, uint32_t anchorLabel)
{
  // against a largest coefficient of 1 this is negligible for the equation and its rank, but nonzero:
  // pede's text reader (addItem) silently drops a zero factor, and the anchor with it
  constexpr double AnchorCoeff = 1e-12;
  const auto& labs = con.getLabels();
  if (anchorLabel != std::numeric_limits<uint32_t>::max() && std::find(labs.begin(), labs.end(), anchorLabel) == labs.end()) {
    con.add(anchorLabel, AnchorCoeff);
  }
}

void Volume::collectSensors(std::vector<const Volume*>& out) const
{
  if (isLeaf()) {
    if (!mVirtual && !mIsPseudo && !isRoot()) {
      out.push_back(this);
    }
    return;
  }
  for (const auto& c : mChildren) {
    c->collectSensors(out);
  }
}

const char* Volume::rayCompName(int comp)
{
  static constexpr const char* names[NRayComps] = {"rdphi", "dz", "dr"};
  return (comp >= 0 && comp < NRayComps) ? names[comp] : "?";
}

// Apply the "rayModes" clause of a rule, e.g.
//   { "match": "ITS", "rayModes": { "nPhi": 24, "nT": 16, "tMax": 4.0, "components": ["rdphi", "dz", "dr"] } }
// true takes the RayModesConfig defaults, false removes the flag, an array of objects gives several
// independent sets of equations (e.g. dz, dr over the OB only and rdphi over IB + OB with combineRanges).
// Like globalModes it frees nothing; it flags a branch whose sensors are binned, see writeRayModeConstraints.
void Volume::applyRayModesConfig(const nlohmann::json& cfg, const std::string& pattern)
{
  if (cfg.is_boolean()) {
    mRayModes = cfg.get<bool>() ? std::make_shared<const std::vector<RayModesConfig>>(1) : nullptr;
    return;
  }
  std::vector<RayModesConfig> sets;
  if (cfg.is_object()) {
    sets.push_back(parseRayModesConfig(cfg, pattern));
  } else if (cfg.is_array() && !cfg.empty()) {
    for (const auto& item : cfg) {
      if (!item.is_object()) {
        LOGP(fatal, "Invalid rayModes clause in the rule '{}': the array elements must be objects", pattern);
      }
      sets.push_back(parseRayModesConfig(item, pattern));
    }
  } else {
    LOGP(fatal, "Invalid rayModes clause in the rule '{}': a bool, an object or a non-empty array of objects is expected", pattern);
  }
  mRayModes = std::make_shared<const std::vector<RayModesConfig>>(std::move(sets));
}

Volume::RayModesConfig Volume::parseRayModesConfig(const nlohmann::json& cfg, const std::string& pattern)
{
  RayModesConfig c;
  for (const auto& [key, val] : cfg.items()) {
    if (key == "nPhi") {
      c.nPhi = val.get<int>();
    } else if (key == "phi0") {
      c.phi0 = val.get<double>();
    } else if (key == "nT") {
      c.nT = val.get<int>();
    } else if (key == "tMax") {
      c.tMax = val.get<double>();
    } else if (key == "rSpanMin") {
      c.rSpanMin = val.get<double>();
    } else if (key == "rMin") {
      c.rMin = val.get<double>();
    } else if (key == "rMax") {
      c.rMax = val.get<double>();
    } else if (key == "rEdges") {
      c.rEdges = val.get<std::vector<double>>();
    } else if (key == "components") {
      c.comps.fill(false);
      for (const auto& item : val) {
        const auto name = item.get<std::string>();
        int comp = 0;
        while (comp < NRayComps && name != rayCompName(comp)) {
          ++comp;
        }
        if (comp == NRayComps) {
          LOGP(fatal, "Unknown rayModes component '{}' in the rule '{}', allowed are rdphi, dz, dr", name, pattern);
        }
        c.comps[comp] = true;
      }
    } else if (key == "combineRanges") {
      c.combineRanges = val.get<bool>();
    } else {
      LOGP(fatal, "Unknown rayModes key '{}' in the rule '{}', allowed are nPhi, phi0, nT, tMax, rSpanMin, rEdges, rMin, rMax, components, combineRanges", key, pattern);
    }
  }
  if (c.nPhi < 1 || c.nT < 1 || c.tMax <= 0. || !std::is_sorted(c.rEdges.begin(), c.rEdges.end())) {
    LOGP(fatal, "Invalid rayModes binning in the rule '{}': nPhi, nT >= 1, tMax > 0 and increasing rEdges are required", pattern);
  }
  if (c.combineRanges && c.rEdges.empty()) {
    LOGP(fatal, "Invalid rayModes clause in the rule '{}': combineRanges needs rEdges", pattern);
  }
  return c;
}

namespace
{
// Linearised displacement of a point x (global frame) per unit value of each rigid-body DOF of a
// volume: d(x) = A*x + b, stored row-major as A (9) then b (3). Built from the very matrix the fitted
// corrections are applied with (AlignParam in the local frame, conjugated by L2G, as in
// getFittedLocalDelta and the result writing), by central differences, so the rotation convention
// cannot drift from the one of the fit.
using DOFDisp = std::array<double, 12>;
std::array<DOFDisp, RigidBodyDOFSet::NDOF> rigidBodyDisplacements(const TGeoHMatrix& l2g)
{
  std::array<DOFDisp, RigidBodyDOFSet::NDOF> out{};
  const TGeoHMatrix g2l = l2g.Inverse();
  for (int i = 0; i < RigidBodyDOFSet::NDOF; ++i) {
    const double eps = i < 3 ? 1e-4 : 1e-6; // cm, rad
    TGeoHMatrix gpm[2];
    for (int sgn = 0; sgn < 2; ++sgn) {
      std::array<double, RigidBodyDOFSet::NDOF> par{};
      par[i] = sgn == 0 ? eps : -eps;
      const auto delta = o2::detectors::AlignParam("", -1, par[0], par[1], par[2], par[3], par[4], par[5], true, false).createMatrix();
      gpm[sgn] = l2g * delta * g2l;
    }
    const double *rp = gpm[0].GetRotationMatrix(), *rm = gpm[1].GetRotationMatrix();
    const double *tp = gpm[0].GetTranslation(), *tm = gpm[1].GetTranslation();
    for (int k = 0; k < 9; ++k) {
      out[i][k] = (rp[k] - rm[k]) / (2. * eps);
    }
    for (int k = 0; k < 3; ++k) {
      out[i][9 + k] = (tp[k] - tm[k]) / (2. * eps);
    }
  }
  return out;
}
} // namespace

// For a track from the vertex at azimuth phi and t = z/r, a displacement r*F(phi, t)*e of every
// sensor it crosses (e = phiHat: its phi changes by F; e = zHat: its tgLambda by F; e = rHat: its
// points slide along it, i.e. dz = -t*dr = -r*t*F absorbed by tgLambda, and with B != 0 the curvature)
// costs no chi2 whatever F is. The global modes constrain a few smooth members of this family; with
// DOFs distributed in z the fit finds the others. Discretised: per cell of (phi, t) and component,
//   sum over the sensors s of the cell of (e_s . d_s) / r_s = 0,
// d_s = sum over the free rigid-body DOFs p on the path from s to the root of (dd_s/dp) * p, the
// linearised displacement of the sensor centre, i.e. the mean angular displacement of the cell seen
// from the vertex must vanish. A common-ray motion eps*r*e of the cell gives eps*nSensors != 0 and is
// removed. Any positive weighting would remove it; the angular one with unit weights is the natural
// measure of a ray and does not favour the outer layers by r^2. What else is removed is one direction
// per cell, the mean angular motion of its sensors, which for a single layer is a measurable motion:
// hence a cell whose sensors do not span rSpanMin in radius is skipped.
// The cells use the nominal t = z/r, i.e. a vertex at z = 0. A track from z_v crosses the layer at
// radius r at t + z_v/r, so the cell indicator is a null direction only where z_v/r is small compared
// to the cell width: at the ITS IB (r ~ 2-4 cm, z_v ~ +-5 cm) the shift exceeds the cell, the IB
// sensors of a cell are not on the rays of its OB sensors, and the equation constrains measurable
// IB-vs-OB motion instead (seen as worse vertex and IB dz residuals). rMin excludes such sensors.
void Volume::writeRayModeConstraints(std::ostream& os) const
{
  std::vector<const Volume*> tops; // the topmost flagged volumes: each one is a branch of its own
  std::function<void(const Volume*)> findTops = [&](const Volume* v) {
    if (v->mRayModes) {
      tops.push_back(v);
      return;
    }
    for (const auto& c : v->mChildren) {
      findTops(c.get());
    }
  };
  findTops(this);
  if (tops.empty()) {
    return;
  }
  const uint32_t anchor = minFreeRigidBodyLabel(); // called on the root: see writeGlobalModeConstraints
  std::unordered_map<const Volume*, std::array<DOFDisp, RigidBodyDOFSet::NDOF>> dispCache;
  for (const auto* top : tops) {
    for (const auto& cfg : *top->mRayModes) {
      writeRayModeSet(os, top, cfg, anchor, dispCache);
    }
  }
}

// With combineRanges the cells of the radial ranges sharing (phi, t) give a single equation,
//   sum over the ranges k of (1/n_k) * sum over the sensors s of range k of (e_s . d_s) / r_s = 0,
// i.e. the mean angular displacement of the cell is averaged over the ranges with equal weights, so
// a range with many sensors (the OB) does not pin the ray motion alone: only the common ray motion of
// all ranges is removed, their relative motion (e.g. an IB/OB kink, which the tracks measure) is left
// to the data. A cell populated in one range only reduces to the equation of that range (counted in
// the log as "incomplete").
int Volume::writeRayModeSet(std::ostream& os, const Volume* top, const RayModesConfig& cfg, uint32_t anchorLabel,
                            std::unordered_map<const Volume*, std::array<DOFDisp, RigidBodyDOFSet::NDOF>>& dispCache) const
{
  constexpr double rMin = 1e-4; // cm
  const int nR = static_cast<int>(cfg.rEdges.size()) + 1;
  struct Cell {
    std::array<std::map<uint32_t, double>, NRayComps> terms;
    double rLow{std::numeric_limits<double>::max()}, rHigh{0.};
    int nSensors{0};
  };
  std::vector<Cell> cells(nR * cfg.nPhi * cfg.nT);
  auto cellAt = [&](int iR, int iPhi, int iT) -> Cell& { return cells[(iR * cfg.nPhi + iPhi) * cfg.nT + iT]; };
  std::vector<const Volume*> sensors;
  top->collectSensors(sensors);
  const double dPhi = 2. * M_PI / cfg.nPhi, dT = 2. * cfg.tMax / cfg.nT;
  for (const auto* s : sensors) {
    const double* x = s->getMatrixL2G().GetTranslation();
    const double r = std::hypot(x[0], x[1]);
    if (r < rMin || r < cfg.rMin || r > cfg.rMax) {
      continue;
    }
    const int iR = static_cast<int>(std::upper_bound(cfg.rEdges.begin(), cfg.rEdges.end(), r) - cfg.rEdges.begin());
    double phi = std::fmod(std::atan2(x[1], x[0]) - cfg.phi0, 2. * M_PI);
    phi += phi < 0. ? 2. * M_PI : 0.;
    const int iPhi = std::min(cfg.nPhi - 1, static_cast<int>(phi / dPhi));
    const int iT = std::clamp(static_cast<int>(std::floor((x[2] / r + cfg.tMax) / dT)), 0, cfg.nT - 1);
    auto& cell = cellAt(iR, iPhi, iT);
    cell.rLow = std::min(cell.rLow, r);
    cell.rHigh = std::max(cell.rHigh, r);
    ++cell.nSensors;
    // angular displacement of the sensor seen from the vertex: e/r with e = phiHat, zHat, rHat
    const double r2 = r * r;
    const double field[NRayComps][3] = {{-x[1] / r2, x[0] / r2, 0.}, {0., 0., 1. / r}, {x[0] / r2, x[1] / r2, 0.}};
    for (const Volume* a = s; a; a = a->mParent) {
      const auto* rb = a->mRigidBody.get();
      if (!rb || rb->nFreeDOFs() == 0) {
        continue;
      }
      auto it = dispCache.find(a);
      if (it == dispCache.end()) {
        it = dispCache.emplace(a, rigidBodyDisplacements(a->getMatrixL2G())).first;
      }
      const int nDOF = std::min(rb->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
      for (int i = 0; i < nDOF; ++i) {
        if (!rb->isFree(i)) {
          continue;
        }
        const auto& D = it->second[i];
        double d[3];
        for (int k = 0; k < 3; ++k) {
          d[k] = D[3 * k] * x[0] + D[3 * k + 1] * x[1] + D[3 * k + 2] * x[2] + D[9 + k];
        }
        for (int comp = 0; comp < NRayComps; ++comp) {
          if (cfg.comps[comp]) {
            cell.terms[comp][a->mLabel.raw(i)] += field[comp][0] * d[0] + field[comp][1] * d[1] + field[comp][2] * d[2];
          }
        }
      }
    }
  }

  // the cells entering one equation: a single radial range, or all of them with combineRanges
  int nWritten = 0, nSingleLayer = 0, nEmpty = 0, nIncomplete = 0;
  auto writeEquation = [&](int iRFirst, int iRLast, int iPhi, int iT) {
    double rLow = std::numeric_limits<double>::max(), rHigh = 0.;
    int nSensors = 0, nRanges = 0;
    for (int iR = iRFirst; iR <= iRLast; ++iR) {
      const auto& cell = cellAt(iR, iPhi, iT);
      if (cell.nSensors > 0) {
        rLow = std::min(rLow, cell.rLow);
        rHigh = std::max(rHigh, cell.rHigh);
        nSensors += cell.nSensors;
        ++nRanges;
      }
    }
    if (nSensors == 0) {
      ++nEmpty;
      return;
    }
    if (rHigh - rLow < cfg.rSpanMin) {
      ++nSingleLayer;
      return;
    }
    nIncomplete += nRanges < iRLast - iRFirst + 1;
    for (int comp = 0; comp < NRayComps; ++comp) {
      if (!cfg.comps[comp]) {
        continue;
      }
      std::map<uint32_t, double> terms;
      for (int iR = iRFirst; iR <= iRLast; ++iR) {
        const auto& cell = cellAt(iR, iPhi, iT);
        const double w = iRFirst == iRLast ? 1. : 1. / std::max(1, cell.nSensors); // equal weight per range
        for (const auto& [lab, coeff] : cell.terms[comp]) {
          terms[lab] += w * coeff;
        }
      }
      if (terms.size() < 2) {
        continue;
      }
      double norm = 0.;
      for (const auto& [lab, coeff] : terms) {
        norm = std::max(norm, std::abs(coeff));
      }
      const auto rangeName = iRFirst == iRLast ? std::to_string(iRFirst) : std::format("{}-{} combined", iRFirst, iRLast);
      Constraint con(std::format("Ray mode {} over {}: r-range {} phi-cell {} t-cell {} ({} sensors)", rayCompName(comp), top->mSymName, rangeName, iPhi, iT, nSensors), 0.);
      for (const auto& [lab, coeff] : terms) {
        if (std::abs(coeff) > 1e-9 * norm) { // same relative cut as the global modes
          con.add(lab, coeff / norm);
        }
      }
      if (con.getSize() < 2) {
        continue;
      }
      addPedeGroupingAnchor(con, anchorLabel);
      con.write(os);
      ++nWritten;
    }
  };
  for (int iPhi = 0; iPhi < cfg.nPhi; ++iPhi) {
    for (int iT = 0; iT < cfg.nT; ++iT) {
      if (cfg.combineRanges) {
        writeEquation(0, nR - 1, iPhi, iT);
      } else {
        for (int iR = 0; iR < nR; ++iR) {
          writeEquation(iR, iR, iPhi, iT);
        }
      }
    }
  }
  LOGP(info, "Wrote {} ray-mode constraints over {}: {} sensors in {} x {} x {} cells{}, {} cells empty, {} skipped as spanning less than {} cm in radius{}",
       nWritten, top->mSymName, sensors.size(), nR, cfg.nPhi, cfg.nT, cfg.combineRanges ? " (radial ranges combined)" : "", nEmpty, nSingleLayer, cfg.rSpanMin,
       cfg.combineRanges ? std::format(", {} missing some radial range", nIncomplete) : std::string{});
  return nWritten;
}

void Volume::writeGlobalModeConstraints(std::ostream& os, const ModeEquationOptions& opt) const
{
  // Every mode equation gets the lowest label of any free rigid-body DOF of the tree as an anchor term
  // (coefficient 1e-12, not 0, which pede's reader drops), so that pede's constraint grouping (which
  // processes constraints by their lowest parameter index) meets it first and forms one group; see
  // addPedeGroupingAnchor. Any other long constraint (ray modes) needs the same anchor.
  const uint32_t firstLabel = minFreeRigidBodyLabel();
  std::vector<std::pair<uint32_t, double>> terms;
  for (int eq = 0; eq < NGlobalEq; ++eq) {
    terms.clear();
    int nFlagged = 0;
    if (opt.sensorProjection) {
      collectProjectedModeTerms(modeField(eq), terms, nFlagged, false, eq);
    } else {
      collectGlobalModeTerms(eq, terms, nFlagged);
    }
    if (nFlagged == 0) { // the mode was not requested, or nothing is free where it was requested
      continue;
    }
    writeModeEquation(os, std::format("Global mode {} over {}", globalEqName(eq), mSymName), terms, nFlagged, firstLabel);
  }
  writeGroupModeConstraints(os, firstLabel, opt);
}

// The "groupModes" equations: for each group (the volumes sharing a group name, or a single volume
// without one) and each requested equation, the mode summed over the free rigid-body DOFs of the
// *whole* subtrees of the group's volumes, at every level. Summing every level is the point: a
// mean-of-children pin only sees the nearest free descendants, so a motion it pins can reappear one
// level further down (the OB half-barrel "scissor" did, through the module TX + RY), while a sum over
// all levels constrains the total, whichever level realises it.
// With ModeEquationOptions::relativeToGroupRoots the mode is first made orthogonal, over the sensors
// of each group volume, to the free rigid motion of that volume: its own free DOFs then drop out and
// the equation constrains the deformation of the volume's content only, e.g. a per-layer stave bow
// z^2 rHat acts as (z^2 - <z^2>) rHat and leaves the radial offset of each stave to the data.
void Volume::writeGroupModeConstraints(std::ostream& os, uint32_t anchorLabel, const ModeEquationOptions& opt) const
{
  if (opt.relativeToGroupRoots && !opt.sensorProjection) {
    LOGP(warn, "relativeToGroupRoots needs sensorProjection for the coefficients of the group modes: ignored");
  }
  const bool relToRoots = opt.relativeToGroupRoots && opt.sensorProjection;
  std::map<std::string, std::vector<const Volume*>> groups;
  collectGroupModeRoots(groups);
  std::vector<std::pair<uint32_t, double>> terms;
  for (const auto& [group, roots] : groups) {
    // the members of a named group must request the same equations: it is one set of equations
    for (const auto* v : roots) {
      if (v->mGroupModes != roots.front()->mGroupModes) {
        LOGP(fatal, "The groupModes group '{}' is requested with different mode lists for {} and {}", group, roots.front()->mSymName, v->mSymName);
      }
    }
    for (int eq = 0; eq < NGlobalEq; ++eq) {
      if (!roots.front()->mGroupModes[eq]) {
        continue;
      }
      terms.clear();
      int nFlagged = 0;
      std::array<double, 2> norm2{}; // of the mode before and after removing the free motion of the roots
      for (const auto* v : roots) {
        if (!opt.sensorProjection) {
          v->collectGlobalModeTerms(eq, terms, nFlagged, true);
          continue;
        }
        auto field = modeField(eq);
        if (relToRoots) {
          std::array<double, 2> n2{};
          field = v->removeFreeRigidMotion(field, n2);
          norm2[0] += n2[0];
          norm2[1] += n2[1];
        }
        v->collectProjectedModeTerms(field, terms, nFlagged, true, eq);
      }
      if (nFlagged == 0) {
        LOGP(warn, "Ignoring the group mode {} over {}: no free rigid-body DOF in its {} subtree(s)", globalEqName(eq), group, roots.size());
        continue;
      }
      // the mode is (numerically) a rigid motion of each group volume, which its own free DOFs carry:
      // nothing is left to constrain, and what writeModeEquation would normalise is rounding noise
      if (relToRoots && norm2[1] <= 1e-12 * norm2[0]) {
        LOGP(warn, "Ignoring the group mode {} over {}: it is a free rigid motion of the group volumes", globalEqName(eq), group);
        continue;
      }
      writeModeEquation(os, std::format("Group mode {} over {}", globalEqName(eq), group), terms, nFlagged, anchorLabel);
    }
  }
}

bool Volume::writeModeEquation(std::ostream& os, const std::string& name, const std::vector<std::pair<uint32_t, double>>& terms, int nFlagged, uint32_t anchorLabel)
{
  // the equation is homogeneous (= 0), so its overall scale is arbitrary; normalise the largest
  // coefficient to 1 to keep the constraint matrix well conditioned, since the raw coefficients
  // carry the physical lever arms and span orders of magnitude between modes (r ~ 40 cm for
  // radialScaling against z*r ~ 10^3 cm^2 for the twist)
  const double norm = terms.empty() ? 0. : std::abs(std::max_element(terms.begin(), terms.end(), [](const auto& a, const auto& b) {
                                                      return std::abs(a.second) < std::abs(b.second);
                                                    })->second);
  Constraint con(name, 0.);
  if (norm > 0.) {
    for (const auto& [lab, coeff] : terms) {
      // relative cut: everything above it is a real projection of the mode, everything below is
      // a rounding residue of the L2G matrices (the two are separated by many orders of magnitude)
      if (std::abs(coeff) > 1e-9 * norm) {
        con.add(lab, coeff / norm);
      }
    }
  }
  // A mode can legitimately have no handle on the flagged level: z scaling projects on the local Z
  // (beam) DOF with a coefficient ~ z of the volume (origin evaluation) or ~ sum z of its sensors
  // (sensor projection), which vanishes for an ITS stave centred at z = 0 -- such modes are only
  // observable with free DOFs that are *distributed* in z (module/chip level). Say so rather than silently writing nothing. One term
  // alone is not a mode either, it is a measurement of a single DOF.
  if (con.getSize() < 2) {
    LOGP(warn, "Ignoring the constraint '{}': only {} of the {} flagged volumes' free DOFs project on this mode",
         name, con.getSize(), nFlagged);
    return false;
  }
  addPedeGroupingAnchor(con, anchorLabel);
  con.write(os);
  LOGP(info, "Wrote the constraint '{}': {} terms over {} flagged volumes", name, con.getSize(), nFlagged);
  return true;
}

// Local delta d^n of the volume from its fitted rigid-body parameters (step 1). Returns false if
// the volume has no free rigid-body DOF or all its fitted values are 0 (labels absent from
// labelToValue count as 0).
bool Volume::getFittedLocalDelta(const std::map<uint32_t, double>& labelToValue, TGeoHMatrix& delta) const
{
  // RSTODO In principle, if the volume has no rigid-body DOF, it should not have any fitted parameters either. 
  // So, in the alignment results writing session one could skip initializing the same RB DOFs just for this check,
  // instead, one can check isRigidBodyAllowed() and then loop as for (int i = 0; i < RigidBodyDOFSet::NDOF; ++i)
  const auto* rb = getRigidBody();
  if (rb == nullptr || !rb->nFreeDOFs()) {
    return false;
  }
  std::array<double, RigidBodyDOFSet::NDOF> par{};
  bool nonZero = false;
  for (int i = 0; i < rb->nDOFs(); ++i) {
    if (auto it = labelToValue.find(getLabel().raw(i)); it != labelToValue.end()) {
      par[i] = it->second;
      nonZero |= par[i] != 0.;
    }
  }
  if (!nonZero) {
    return false;
  }
  // the RigidBodyDOFSet parameters coincide with the AlignParam (x,y,z,psi,theta,phi), which builds
  // the matrix with the same angle convention used for the derivatives (A.47)
  using RB = RigidBodyDOFSet;
  delta = o2::detectors::AlignParam("", -1, par[RB::TX], par[RB::TY], par[RB::TZ], par[RB::RX], par[RB::RY], par[RB::RZ], true, false).createMatrix();
  return true;
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
          os << std::format("{:<10} {:>+15g} {:>+15g} ! {} {:<6} ",
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

namespace
{
// A rule pattern naming one volume explicitly (no regex metacharacters). Only for such a pattern an
// inapplicable rule is worth a warning: a regex is typically meant to be broad, matching volumes
// of several kinds.
bool isLiteralPattern(const std::string& pattern)
{
  return pattern.find_first_of("\\^$.|?*+()[]{}") == std::string::npos;
}
} // namespace

// Apply the "pinChildrenMean" clause of a rule: a bool (all DOFs), "all", or an array of DOF names.
// Unlike rigidBody/calib this never creates a DOF set, it only flags DOFs of the volume's own
// (possibly absent or fixed) RigidBodyDOFSet for writeRigidBodyConstraints. A later matching rule
// replaces the flags of an earlier one, like the rigidBody/calib clauses.
void Volume::applyPinChildrenMeanConfig(const nlohmann::json& pin, const std::string& pattern, bool noCross)
{
  const char* clause = noCross ? "pinChildrenMeanNoCross" : "pinChildrenMean";
  auto& flags = noCross ? mPinChildrenMeanNoCross : mPinChildrenMean;
  if (!isRigidBodyAllowed()) {
    if (isLiteralPattern(pattern)) {
      LOGP(warn, "Ignoring the {} rule '{}': {} is not rigid-body alignable", clause, pattern, mSymName);
    }
    return;
  }
  flags.fill(false);
  if (pin.is_boolean()) {
    flags.fill(pin.get<bool>());
    return;
  }
  if (pin.is_string() && pin.get<std::string>() == "all") {
    flags.fill(true);
    return;
  }
  if (!pin.is_array()) {
    LOGP(fatal, "Invalid {} clause in the rule '{}', allowed are a bool, \"all\" or a list of DOF names", clause, pattern);
  }
  for (const auto& name : pin) {
    const auto idx = RigidBodyDOFSet::dofIndex(name.get<std::string>());
    if (idx < 0) {
      LOGP(fatal, "Unknown rigid-body DOF '{}' in the {} rule '{}', allowed are TX,TY,TZ,RX,RY,RZ", name.get<std::string>(), clause, pattern);
    }
    flags[idx] = true;
  }
}

// Apply the "weightedMean" clause of a rule: a bool (all DOFs), "all", or an array of DOF names.
// Unlike pinChildrenMean's flags (which are meaningless without a separate "on" switch), each DOF
// here is tri-state: Unset (no matching rule yet) falls back to the statWeighted default passed to
// writeRigidBodyConstraints, while a matching rule always sets every DOF explicitly (True or False),
// "all" or a bool touching every DOF and an array touching only the named ones -- false for the rest,
// mirroring pinChildrenMean's "replace, don't merge" reset of the whole array on every match.
void Volume::applyWeightedMeanConfig(const nlohmann::json& cfg, const std::string& pattern)
{
  if (!isRigidBodyAllowed()) {
    if (isLiteralPattern(pattern)) {
      LOGP(warn, "Ignoring the weightedMean rule '{}': {} is not rigid-body alignable", pattern, mSymName);
    }
    return;
  }
  if (cfg.is_boolean()) {
    mWeightedMean.fill(cfg.get<bool>() ? WeightedMeanFlag::True : WeightedMeanFlag::False);
    return;
  }
  if (cfg.is_string() && cfg.get<std::string>() == "all") {
    mWeightedMean.fill(WeightedMeanFlag::True);
    return;
  }
  if (!cfg.is_array()) {
    LOGP(fatal, "Invalid weightedMean clause in the rule '{}', allowed are a bool, \"all\" or a list of DOF names", pattern);
  }
  mWeightedMean.fill(WeightedMeanFlag::False);
  for (const auto& name : cfg) {
    const auto idx = RigidBodyDOFSet::dofIndex(name.get<std::string>());
    if (idx < 0) {
      LOGP(fatal, "Unknown rigid-body DOF '{}' in the weightedMean rule '{}', allowed are TX,TY,TZ,RX,RY,RZ", name.get<std::string>(), pattern);
    }
    mWeightedMean[idx] = WeightedMeanFlag::True;
  }
}

// Apply the "globalModes" clause of a rule: a bool (all equations), "all", or an array of names, each
// of which may be a mode name (radialScaling, zScaling, curl, twist, elliptical, bowing,
// radialHarmonic, tangentialHarmonic, globalTilt, radialPolar, tangentialPolar, longitudinalPolar,
// globalShift), one order of a harmonic family (radialHarmonic3), or one equation (elliptical_cos2phi,
// tangentialHarmonic2_sin2phi, radialPolar_r0t1, globalShift_x).
// Like pinChildrenMean this never creates a DOF set and never frees anything: it only flags this
// volume as a participant of the global deformation-mode equations written by
// writeGlobalModeConstraints, over whichever of its rigid-body DOFs happen to be free at the time of
// writing. A later matching rule replaces the flags of an earlier one, like the other clauses.
// The clause is meant for a level which is *free*, typically the stave level, e.g.
//   { "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+",
//     "globalModes": ["radialScaling", "curl", "elliptical", "twist",
//                     "radialHarmonic", "tangentialHarmonic"] }
// Flagging several nested levels at once is legitimate and not double counting: the equation then
// constrains the total amplitude of the mode summed over all of them, which is the quantity that is
// degenerate, regardless of the level that realises it.
void Volume::applyGlobalModesConfig(const nlohmann::json& modes, const std::string& pattern)
{
  if (!isRigidBodyAllowed()) {
    if (isLiteralPattern(pattern)) {
      LOGP(warn, "Ignoring the globalModes rule '{}': {} is not rigid-body alignable", pattern, mSymName);
    }
    return;
  }
  mGlobalModes = parseGlobalModeList(modes, "globalModes", pattern);
}

// Apply the "groupModes" clause of a rule. It takes the same mode list as "globalModes", for which
// this volume forms a group of its own, or {"modes": <list>, "group": "<name>"}, for which all the
// volumes given the same name form one group. Either way the equations are written over the free
// rigid-body DOFs of the whole subtree of every volume of the group, at all levels, by
// writeGroupModeConstraints; nothing has to be flagged below. The volume itself needs no free DOF
// (it can be the fixed half barrel whose children are to be constrained), and like the other clauses
// a later matching rule replaces the modes and the group name of an earlier one, e.g.
//   { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]", "groupModes": ["globalShift_x", "globalShift_y"] },
//   { "match": "ITS/ITSULayer3/ITSUHalfBarrel[0-1]", "groupModes": {"group": "OB-Lr3", "modes": ["bowing"]} }
void Volume::applyGroupModesConfig(const nlohmann::json& cfg, const std::string& pattern)
{
  mGroupModesName.clear();
  if (cfg.is_object()) {
    if (!cfg.contains("modes")) {
      LOGP(fatal, "The groupModes clause of the rule '{}' has no \"modes\" list", pattern);
    }
    for (const auto& [key, val] : cfg.items()) {
      if (key != "modes" && key != "group") {
        LOGP(fatal, "Unknown key '{}' in the groupModes clause of the rule '{}', allowed are \"modes\" and \"group\"", key, pattern);
      }
    }
    mGroupModes = parseGlobalModeList(cfg["modes"], "groupModes", pattern);
    if (cfg.contains("group")) {
      mGroupModesName = cfg["group"].get<std::string>();
    }
  } else {
    mGroupModes = parseGlobalModeList(cfg, "groupModes", pattern);
  }
}

std::array<bool, Volume::NGlobalEq> Volume::parseGlobalModeList(const nlohmann::json& modes, const std::string& clause, const std::string& pattern)
{
  std::array<bool, NGlobalEq> flags{};
  if (modes.is_boolean()) {
    flags.fill(modes.get<bool>());
    return flags;
  }
  if (modes.is_string()) {
    if (modes.get<std::string>() != "all") {
      LOGP(fatal, "Unknown {} value '{}' in the rule '{}', allowed are a bool, \"all\" or a list of mode names", clause, modes.get<std::string>(), pattern);
    }
    flags.fill(true);
    return flags;
  }
  if (!modes.is_array()) {
    LOGP(fatal, "Invalid {} clause in the rule '{}', allowed are a bool, \"all\" or a list of mode names", clause, pattern);
  }
  for (const auto& item : modes) {
    const auto name = item.get<std::string>();
    int nSet = 0;
    for (int eq = 0; eq < NGlobalEq; ++eq) {
      const int m = globalEqMode(eq);
      const std::string mode = globalModeName(m);
      // a name selects either one equation ("radialHarmonic3_cos3phi", "elliptical_sin2phi"), or a
      // whole mode ("elliptical" = both phases, "radialHarmonic" = every order), or one order of a
      // harmonic family ("radialHarmonic3" = its cos and sin)
      const bool hit = name == globalEqName(eq) || name == mode ||
                       ((m == RadialHarmonic || m == TangentialHarmonic) &&
                        name == mode + std::to_string(globalHarmonicL(globalEqPhase(eq))));
      if (hit) {
        flags[eq] = true;
        ++nSet;
      }
    }
    if (nSet == 0) {
      LOGP(fatal, "Unknown global mode '{}' in the {} rule '{}'. Allowed are the mode names radialScaling, zScaling, curl, twist, "
                  "elliptical, bowing, radialHarmonic, tangentialHarmonic, globalTilt, radialPolar, tangentialPolar, "
                  "longitudinalPolar, globalShift, a single harmonic order (e.g. radialHarmonic3) or a single equation "
                  "(e.g. elliptical_cos2phi, globalTilt_y, globalShift_x, radialPolar_r0t1)",
           name, clause, pattern);
    }
  }
  return flags;
}

// Apply the "measurement" clause of a rule: an object DOF name -> [value, sigma] or {"value", "sigma"}.
// Each name is resolved against whichever of this volume's rigid-body or calibration DOF set has a
// matching dofName() (the two never share a name, e.g. TX..RZ vs X,Y,Z,SlopeX,SlopeY for the mean
// vertex). Like pinChildrenMean it creates no DOF set: the measurement is written only if the DOF is
// free at the time of writing. A later matching rule replaces all the measurements of both DOF sets.
void Volume::applyMeasurementConfig(const nlohmann::json& meas, const std::string& pattern)
{
  if (!mRigidBody && !mCalib) { // nothing to measure on this volume, e.g. under a broadly matching pattern
    if (isLiteralPattern(pattern)) {
      LOGP(warn, "Ignoring the measurement rule '{}': {} has neither a rigid-body nor a calibration DOF set", pattern, mSymName);
    }
    return;
  }
  if (!meas.is_object()) {
    LOGP(fatal, "Invalid measurement clause in the rule '{}', an object DOF name -> [value, sigma] is expected", pattern);
  }
  if (mRigidBody) {
    mRigidBody->clearMeasurements();
  }
  if (mCalib) {
    mCalib->clearMeasurements();
  }
  for (const auto& [name, item] : meas.items()) {
    DOFSet* target = nullptr;
    int idx = -1;
    if (mRigidBody && (idx = mRigidBody->dofIndex(name)) >= 0) {
      target = mRigidBody.get();
    } else if (mCalib && (idx = mCalib->dofIndex(name)) >= 0) {
      target = mCalib.get();
    }
    if (!target) {
      LOGP(fatal, "Unknown DOF '{}' in the measurement rule '{}' for {}: no such DOF in its rigid-body or calibration set", name, pattern, mSymName);
    }
    double value{}, sigma{};
    if (item.is_array() && item.size() == 2) {
      value = item[0].get<double>();
      sigma = item[1].get<double>();
    } else if (item.is_object() && item.contains("value") && item.contains("sigma")) {
      value = item["value"].get<double>();
      sigma = item["sigma"].get<double>();
    } else {
      LOGP(fatal, "Invalid measurement of {} in the rule '{}', allowed are [value, sigma] or {{\"value\": v, \"sigma\": s}}", name, pattern);
    }
    if (!(sigma > 0.)) {
      LOGP(fatal, "Measurement of {} in the rule '{}' needs sigma > 0, got {}", name, pattern, sigma);
    }
    target->setMeasurement(idx, value, sigma);
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
  const json rules = geomutils::readJsonRules(jsonPath, "DOF config", true);

  // a rigidBody/calib clause defines a DOF set (and replaces an earlier one), a measurement only adds
  // a prior to DOFs which are already free: in one rule the DOF list would be read as the definition
  // and silently fix all the unlisted DOFs, hence the measurement must come in a rule of its own,
  // after the one defining the DOFs
  for (const auto& rule : rules) {
    if (rule.contains("measurement") && (rule.contains("rigidBody") || rule.contains("calib"))) {
      LOGP(fatal, "The rule '{}' has both a measurement clause and a rigidBody or calib clause: put the measurement in a separate rule following the DOF definition",
           rule.value("match", std::string("?")));
    }
  }

  // Each pattern is also tried with an implicit leading ".*", so e.g. a
  // bare "ITSUStave3" still matches regardless of its ancestors. Compiled once per rule here,
  // rather than once per (rule, volume) pair in the traversal below.
  std::vector<std::string> patterns;
  std::vector<geomutils::VolumePattern> compiledPatterns;
  patterns.reserve(rules.size());
  compiledPatterns.reserve(rules.size());
  for (const auto& rule : rules) {
    const auto pattern = rule["match"].get<std::string>();
    try {
      compiledPatterns.emplace_back(pattern);
    } catch (const std::regex_error& e) {
      LOGP(fatal, "Invalid regex in the rule match '{}': {}", pattern, e.what());
    }
    patterns.push_back(pattern);
  }

  root->traverse([&](Volume* vol) {
    if (vol->isPseudo()) {
      return;
    }
    const std::string& sym = vol->getSymName();
    for (size_t i = 0; i < rules.size(); ++i) { // in order: a later matching rule replaces the DOF set of an earlier one
      if (!compiledPatterns[i].matches(sym)) {
        continue;
      }
      const auto& rule = rules[i];
      const auto& pattern = patterns[i];
      if (rule.contains("rigidBody")) {
        // silently not applicable to the volumes which are not rigid-body alignable (the root of the
        // hierarchy, the envelopes w/o own geometry), for which only the calibration DOFs can be
        // configured. An explicitly named volume is likely a configuration mistake.
        if (vol->isRigidBodyAllowed()) {
          vol->setRigidBody(makeRigidBodyDOFSet(rule["rigidBody"], pattern));
        } else if (isLiteralPattern(pattern)) {
          LOGP(warn, "Ignoring the rigidBody rule '{}': {} is not rigid-body alignable", pattern, sym);
        }
      }
      if (rule.contains("calib")) {
        if (auto dofSet = makeCalibDOFSet(rule["calib"], pattern, sym)) {
          vol->setCalib(std::move(dofSet));
        }
      }
      if (rule.contains("pinChildrenMean")) {
        vol->applyPinChildrenMeanConfig(rule["pinChildrenMean"], pattern, false);
      }
      if (rule.contains("pinChildrenMeanNoCross")) {
        vol->applyPinChildrenMeanConfig(rule["pinChildrenMeanNoCross"], pattern, true);
      }
      if (rule.contains("weightedMean")) {
        vol->applyWeightedMeanConfig(rule["weightedMean"], pattern);
      }
      if (rule.contains("globalModes")) {
        vol->applyGlobalModesConfig(rule["globalModes"], pattern);
      }
      if (rule.contains("groupModes")) {
        vol->applyGroupModesConfig(rule["groupModes"], pattern);
      }
      if (rule.contains("rayModes")) {
        vol->applyRayModesConfig(rule["rayModes"], pattern);
      }
      if (rule.contains("measurement")) {
        vol->applyMeasurementConfig(rule["measurement"], pattern);
      }
    }
  });
}

// The caller (writeMillepedeResults) only invokes this when getRigidBody() has free DOFs.
bool Volume::MP2JSON_RB(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const
{
  auto rbArr = nlohmann::json::array();
  for (int i = 0; i < getRigidBody()->nDOFs(); ++i) {
    rbArr.push_back(getFittedValue(labelToValue, mLabel.raw(i)));
  }
  entry["rigidBody"] = rbArr;
  return true;
}

// The caller (writeMillepedeResults) only invokes this when getCalib() has free DOFs.
bool Volume::MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const
{
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

bool Volume::MP2ROOT_RB(const std::map<uint32_t, double>&, nlohmann::json&) const
{
  LOGP(warn, "MP2ROOT_RB is not implemented for volume {}", getSymName());
  return false;
}

bool Volume::MP2ROOT_Calib(const std::map<uint32_t, double>&, nlohmann::json&) const
{
  LOGP(warn, "MP2ROOT_Calib is not implemented for volume {} (calib type {})", getSymName(), static_cast<int>(getCalib()->type()));
  return false;
}

std::map<uint32_t, double> Volume::readMillepedeResults(const std::string& milleResPath)
{
  // parse millepede.res: label fittedValue presigma [...]. Fixed parameters (presigma < 0) are kept
  // too: pede reports them at their preset value, which is a real correction when the steering fixed
  // them at a non-zero one (e.g. a fit restarted on top of a previous solution); otherwise it is 0.
  std::ifstream fin(milleResPath);
  if (!fin.is_open()) {
    LOGP(fatal, "Cannot open millepede result file: {}", milleResPath);
  }
  std::map<uint32_t, double> labelToValue;
  size_t nFixed = 0, nFixedNonZero = 0;
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
    labelToValue[label] = value;
    nFixed += presigma < 0.0;
    nFixedNonZero += presigma < 0.0 && value != 0.0;
  }
  fin.close();
  LOGP(info, "Parsed {} parameters from {}, of which {} fixed ({} at a non-zero value)", labelToValue.size(), milleResPath, nFixed, nFixedNonZero);
  return labelToValue;
}

void Volume::writeMillepedeResults(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outJsonPath)
{
  using json = nlohmann::json;

  // collect results per volume that has RB or calib DOFs
  json output = json::array();
  root->traverse([&](Volume* vol) {
    const int id = vol->getSensorId();
    json entry;
    entry["symName"] = vol->getSymName();
    entry["id"] = id;

    bool write = false;
    if (vol->getRigidBody() && vol->getRigidBody()->nFreeDOFs()) {
      write |= vol->MP2JSON_RB(labelToValue, entry);
    }
    if (vol->getCalib() && vol->getCalib()->nFreeDOFs()) {
      write |= vol->MP2JSON_Calib(labelToValue, entry);
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
