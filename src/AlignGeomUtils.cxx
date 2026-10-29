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

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>

#include <Eigen/Dense>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>

#include "Framework/Logger.h"
#include "DetectorsBase/GeometryManager.h"
#include "O2Align/AlignGeomUtils.h"

namespace o2::alignrs::geomutils
{

nlohmann::json readJsonRules(const std::string& path, const char* what, bool allowDefaults)
{
  using json = nlohmann::json;
  std::ifstream f(path);
  if (!f.is_open()) {
    LOGP(fatal, "Cannot open {} file: {}", what, path);
  }
  const auto data = parseJsonWithComments(f);
  json rules = data.is_array() ? data : data.value("rules", json::array());
  if (data.is_object() && data.contains("defaults")) {
    if (!allowDefaults) {
      LOGP(fatal, "{} file {}: the 'defaults' member is not supported, give the rules explicitly", what, path);
    }
    json defRule = data["defaults"];
    defRule["match"] = ".*";
    rules.insert(rules.begin(), defRule);
  }
  // elements which had nothing but comments
  rules.erase(std::remove_if(rules.begin(), rules.end(), [](const json& r) { return r.is_object() && r.empty(); }), rules.end());
  return rules;
}

VolumePattern::VolumePattern(const std::string& pattern)
  : mDirect(pattern), mAnyPrefix(".*(?:" + pattern + ")")
{
}

bool VolumePattern::matches(const std::string& symName) const
{
  return std::regex_match(symName, mDirect) || std::regex_match(symName, mAnyPrefix);
}

std::string getAlignablePath(const std::string& symName)
{
  const auto* pne = gGeoManager->GetAlignableEntry(symName.c_str());
  if (pne == nullptr) {
    LOGP(fatal, "Symbolic volume '{}' has no corresponding alignable entry!", symName);
  }
  return pne->GetTitle();
}

int getPathDepth(const std::string& path)
{
  return static_cast<int>(std::count(path.begin(), path.end(), '/'));
}

const TGeoHMatrix& getGeometryMatrix(const std::string& symName)
{
  const auto* m = o2::base::GeometryManager::getMatrix(symName.c_str());
  if (m == nullptr) {
    LOGP(fatal, "Cannot get the matrix of {}", symName);
  }
  return *m;
}

TGeoHMatrix getOnPath(const PathMatrices& matrices, std::string path, bool includeSelf)
{
  auto stripLast = [&path]() {
    const auto pos = path.find_last_of('/');
    path.resize(pos == std::string::npos ? 0 : pos);
  };
  if (!includeSelf) {
    stripLast();
  }
  while (!path.empty()) {
    if (auto it = matrices.find(path); it != matrices.end()) {
      return it->second;
    }
    stripLast();
  }
  return {};
}

TGeoHMatrix getInitialGlobalDelta(const o2::detectors::AlignParam& par)
{
  auto delta = par.createMatrix();
  if (par.isGlobal()) {
    return delta;
  }
  const auto& gp = getGeometryMatrix(par.getSymName());
  return gp * delta * gp.Inverse();
}

PathMatrices collectInitialCumulativeDeltas(const std::vector<o2::detectors::AlignParam>* initial)
{
  std::vector<std::pair<std::string, const o2::detectors::AlignParam*>> byDepth;
  if (!initial) {
    return {};
  }
  for (const auto& par : *initial) {
    byDepth.emplace_back(getAlignablePath(par.getSymName()), &par);
  }
  std::stable_sort(byDepth.begin(), byDepth.end(), [](const auto& a, const auto& b) { return getPathDepth(a.first) < getPathDepth(b.first); });
  PathMatrices byPath;
  for (const auto& [path, par] : byDepth) {
    if (byPath.contains(path)) {
      LOGP(fatal, "Initial alignment contains {} more than once", par->getSymName());
    }
    byPath[path] = getInitialGlobalDelta(*par) * getOnPath(byPath, path, false);
  }
  return byPath;
}

std::vector<Alignable> collectAlignablesBelow(const std::unordered_set<std::string>& topPaths)
{
  std::vector<Alignable> result;
  for (int i = 0; i < gGeoManager->GetNAlignable(); ++i) {
    const auto* pne = gGeoManager->GetAlignableEntry(i);
    const std::string path = pne->GetTitle();
    for (std::string p = path; !p.empty();) {
      if (topPaths.contains(p)) {
        result.emplace_back(pne->GetName(), path);
        break;
      }
      const auto pos = p.find_last_of('/');
      p.resize(pos == std::string::npos ? 0 : pos);
    }
  }
  std::stable_sort(result.begin(), result.end(), [](const Alignable& a, const Alignable& b) { return getPathDepth(a.path) < getPathDepth(b.path); });
  return result;
}

std::unordered_map<std::string, std::vector<size_t>> collectAlignableChildren(const std::vector<Alignable>& alignables)
{
  std::unordered_set<std::string> paths;
  for (const auto& a : alignables) {
    paths.insert(a.path);
  }
  std::unordered_map<std::string, std::vector<size_t>> children;
  for (size_t i = 0; i < alignables.size(); ++i) {
    std::string p = alignables[i].path;
    while (true) {
      const auto pos = p.find_last_of('/');
      if (pos == std::string::npos || pos == 0) {
        break;
      }
      p.resize(pos);
      if (paths.contains(p)) {
        children[p].push_back(i);
        break;
      }
    }
  }
  return children;
}

namespace
{
using Vec3 = std::array<double, 3>;
using Vector6 = Eigen::Matrix<double, 6, 1>;

/// small-angle rotation vector of the rotation part of m (its antisymmetric part)
Vec3 getRotationVector(const TGeoHMatrix& m)
{
  const double* r = m.GetRotationMatrix();
  return {0.5 * (r[7] - r[5]), 0.5 * (r[2] - r[6]), 0.5 * (r[3] - r[1])};
}

/// rigid motion x -> R(w) * x + t, R(w) being the rotation by |w| around w (Rodrigues formula)
TGeoHMatrix makeRigidMotion(const Vec3& t, const Vec3& w)
{
  double rot[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  if (const double angle = std::hypot(w[0], w[1], w[2]); angle > 0.) {
    const Vec3 k{w[0] / angle, w[1] / angle, w[2] / angle};
    const double c = std::cos(angle), s = std::sin(angle);
    const double kx[9] = {0., -k[2], k[1], k[2], 0., -k[0], -k[1], k[0], 0.};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rot[3 * i + j] = (i == j ? c : 0.) + s * kx[3 * i + j] + (1. - c) * k[i] * k[j];
      }
    }
  }
  TGeoHMatrix m;
  m.SetRotation(rot);
  m.SetTranslation(t.data());
  return m;
}

/// corners of the bounding box of the TGeo volume at path, in its local frame
std::vector<Vec3> getBoundingBoxCorners(const std::string& path)
{
  if (!gGeoManager->cd(path.c_str())) {
    LOGP(fatal, "Cannot navigate to {}", path);
  }
  auto* volume = gGeoManager->GetCurrentVolume();
  auto* box = static_cast<TGeoBBox*>(volume->GetShape()); // every shape is a TGeoBBox
  if (volume->IsAssembly()) {
    box->ComputeBBox(); // computed only on demand for an assembly (e.g. an ITS stave), 0 until then
  }
  const double* o = box->GetOrigin();
  std::vector<Vec3> corners;
  for (const double sx : {-1., 1.}) {
    for (const double sy : {-1., 1.}) {
      for (const double sz : {-1., 1.}) {
        corners.push_back({o[0] + sx * box->GetDX(), o[1] + sy * box->GetDY(), o[2] + sz * box->GetDZ()});
      }
    }
  }
  return corners;
}

/// jacobian of the displacement T + w x x of the point x wrt (T, w)
Eigen::Matrix<double, 3, 6> getRigidMotionJacobian(const Vec3& x)
{
  Eigen::Matrix<double, 3, 6> a;
  a << 1., 0., 0., 0., x[2], -x[1],
    0., 1., 0., -x[2], 0., x[0],
    0., 0., 1., x[1], -x[0], 0.;
  return a;
}

/// Small rigid motion (T, w), x -> x + T + w x x, best describing the displacements u of the points x
/// in the L1 sense, i.e. minimizing the sum of the Euclidean norms of the residuals, by iteratively
/// reweighted least squares. This is the robust (median-like) counterpart of the least squares fit:
/// the points moving differently from the majority do not bias it. A combination of (T, w) not
/// constrained by the points (e.g. a rotation around the line they all lie on) is left at 0.
/// Returns the sum of the norms of the residuals.
double fitRobustRigidMotion(const std::vector<Vec3>& x, const std::vector<Vec3>& u, Vec3& t, Vec3& w)
{
  constexpr int MaxIterations = 200;
  constexpr double MinResidual = 1e-7; // cm, regularizes the weight of the points fitted exactly
  constexpr double Tolerance = 1e-11;
  std::vector<double> weights(x.size(), 1.); // the 1st iteration is the least squares fit
  Vector6 theta = Vector6::Zero();
  double cost = 0.;
  for (int it = 0; it < MaxIterations; ++it) {
    Eigen::Matrix<double, 6, 6> normal = Eigen::Matrix<double, 6, 6>::Zero();
    Vector6 rhs = Vector6::Zero();
    for (size_t k = 0; k < x.size(); ++k) {
      const auto a = getRigidMotionJacobian(x[k]);
      normal += weights[k] * a.transpose() * a;
      rhs += weights[k] * a.transpose() * Eigen::Vector3d(u[k].data());
    }
    const Vector6 next = normal.completeOrthogonalDecomposition().solve(rhs);
    const bool converged = (next - theta).cwiseAbs().maxCoeff() < Tolerance;
    theta = next;
    cost = 0.;
    for (size_t k = 0; k < x.size(); ++k) {
      const double res = (Eigen::Vector3d(u[k].data()) - getRigidMotionJacobian(x[k]) * theta).norm();
      cost += res;
      weights[k] = 1. / std::max(res, MinResidual);
    }
    if (converged) {
      break;
    }
  }
  t = {theta[0], theta[1], theta[2]};
  w = {theta[3], theta[4], theta[5]};
  return cost;
}

/// Displacements u, in the frame of the parent whose current L2G is gPar, of the points x of the
/// children: with D_i = P_i * P_p^-1 the global delta of child i relative to its parent, its motion in
/// the parent frame is E_i = G_p^-1 * D_i * G_p and u = E_i(x) - x.
std::vector<Vec3> getChildrenDisplacements(const TGeoHMatrix& gPar, const TGeoHMatrix& pPar, const std::vector<TGeoHMatrix>& pChildren, const std::vector<std::vector<Vec3>>& points, std::vector<Vec3>& allPoints)
{
  const auto gParInv = gPar.Inverse();
  const auto pParInv = pPar.Inverse();
  std::vector<Vec3> displacements;
  allPoints.clear();
  for (size_t i = 0; i < pChildren.size(); ++i) {
    const TGeoHMatrix e = gParInv * pChildren[i] * pParInv * gPar;
    for (const auto& x : points[i]) {
      Vec3 moved;
      e.LocalToMaster(x.data(), moved.data());
      displacements.push_back({moved[0] - x[0], moved[1] - x[1], moved[2] - x[2]});
      allPoints.push_back(x);
    }
  }
  return displacements;
}

/// L1 sizes of the children motions relative to their parent, before and after the transfer
struct TransferCost {
  double before = 0.;
  double after = 0.;
};

/// Move the parent by the rigid motion robustly describing the displacements of the bounding box
/// corners of its children (fitRobustRigidMotion), keeping the cumulative deltas of the children,
/// i.e. their placement in the lab. Iterated, the displacements being linearized in the fit.
/// Returns the transferred motion in the parent frame.
TGeoHMatrix moveParentByChildrenMotion(PathMatrices& cumulative, const std::string& parentPath, const std::vector<std::string>& childPaths, const PathMatrices& ideal, TransferCost& cost)
{
  constexpr int MaxIterations = 5;
  constexpr double RotTolerance = 1e-14, ShiftTolerance = 1e-12;
  const auto& gIdPar = ideal.at(parentPath);
  std::vector<TGeoHMatrix> pChildren;
  std::vector<std::vector<Vec3>> points; // bounding box corners of the children at their ideal place in the parent frame
  for (const auto& path : childPaths) {
    pChildren.push_back(cumulative.at(path));
    const TGeoHMatrix rel = gIdPar.Inverse() * ideal.at(path);
    auto& corners = points.emplace_back(getBoundingBoxCorners(path));
    for (auto& c : corners) {
      Vec3 x;
      rel.LocalToMaster(c.data(), x.data());
      c = x;
    }
  }
  TGeoHMatrix pPar = cumulative.at(parentPath), total;
  std::vector<Vec3> allPoints;
  for (int it = 0; it < MaxIterations; ++it) {
    const TGeoHMatrix gPar = pPar * gIdPar;
    const auto u = getChildrenDisplacements(gPar, pPar, pChildren, points, allPoints);
    if (it == 0) {
      for (const auto& d : u) {
        cost.before += std::hypot(d[0], d[1], d[2]);
      }
    }
    Vec3 t, w;
    cost.after = fitRobustRigidMotion(allPoints, u, t, w);
    const auto step = makeRigidMotion(t, w);
    pPar = gPar * step * gPar.Inverse() * pPar; // G_p -> G_p * step
    total = total * step;
    if (std::max({std::abs(w[0]), std::abs(w[1]), std::abs(w[2])}) < RotTolerance && std::max({std::abs(t[0]), std::abs(t[1]), std::abs(t[2])}) < ShiftTolerance) {
      break;
    }
  }
  cumulative[parentPath] = pPar;
  return total;
}

} // namespace

void transferChildrenCommonMotion(PathMatrices& cumulative, const std::vector<Alignable>& alignables, const PathMatrices& ideal, const ChildrenTransferConfig& cfg)
{
  const auto children = collectAlignableChildren(alignables);
  std::vector<size_t> parents;
  for (size_t i = 0; i < alignables.size(); ++i) {
    const auto it = children.find(alignables[i].path);
    if (it == children.end() || static_cast<int>(it->second.size()) < cfg.minChildren) {
      continue;
    }
    if (cfg.parents.empty() || std::any_of(cfg.parents.begin(), cfg.parents.end(), [&](const auto& pat) { return pat.matches(alignables[i].symName); })) {
      parents.push_back(i);
    }
  }
  // deepest first, for the common motion of the children of a parent to be passed further up
  std::stable_sort(parents.begin(), parents.end(), [&](size_t a, size_t b) { return getPathDepth(alignables[a].path) > getPathDepth(alignables[b].path); });

  // an explicit cumulative delta for a volume moving so far with its nearest ancestor having one
  auto makeExplicit = [&cumulative](const std::string& path) {
    if (!cumulative.contains(path)) {
      cumulative[path] = getOnPath(cumulative, path, true);
    }
  };
  double maxShift = 0., maxRot = 0.;
  TransferCost totalCost;
  size_t nPoints = 0;
  for (const auto ip : parents) {
    const auto& parent = alignables[ip];
    std::vector<std::string> childPaths;
    for (const auto ic : children.at(parent.path)) {
      childPaths.push_back(alignables[ic].path);
      makeExplicit(childPaths.back()); // before the parent is moved
    }
    makeExplicit(parent.path);
    TransferCost cost;
    const auto moved = moveParentByChildrenMotion(cumulative, parent.path, childPaths, ideal, cost);
    const auto w = getRotationVector(moved);
    const double* t = moved.GetTranslation();
    const size_t np = 8 * childPaths.size();
    LOGP(debug, "{}: transferred common motion of {} children: T = ({:+.3e}, {:+.3e}, {:+.3e}) cm, R = ({:+.3e}, {:+.3e}, {:+.3e}) rad, mean child displacement {:.3e} -> {:.3e} cm",
         parent.symName, childPaths.size(), t[0], t[1], t[2], w[0], w[1], w[2], cost.before / np, cost.after / np);
    maxShift = std::max({maxShift, std::abs(t[0]), std::abs(t[1]), std::abs(t[2])});
    maxRot = std::max({maxRot, std::abs(w[0]), std::abs(w[1]), std::abs(w[2])});
    totalCost.before += cost.before;
    totalCost.after += cost.after;
    nPoints += np;
  }
  LOGP(info, "Transferred the common motion of the children to {} parents, max transferred shift {:.3e} cm, rotation {:.3e} rad; mean displacement of the children wrt their parents {:.3e} -> {:.3e} cm",
       parents.size(), maxShift, maxRot, nPoints ? totalCost.before / nPoints : 0., nPoints ? totalCost.after / nPoints : 0.);
}

} // namespace o2::alignrs::geomutils
