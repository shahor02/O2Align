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

#include "O2Align/MisalignmentUtils.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <vector>
#include <array>

#include <TMatrixD.h>
#include <nlohmann/json.hpp>
#include "O2Align/AlignGeomUtils.h"

#include "Framework/Logger.h"
#include "CommonConstants/MathConstants.h"
#include "ITS3Base/SpecsV2.h"
#include "O2Align/DOFSet.h"
#include "O2Align/Volume.h"

namespace o2::alignrs
{

namespace
{
using json = nlohmann::json;
using DOFValues = std::vector<std::pair<int, double>>; // DOF index, value

std::string describe(const std::string& pattern, const Volume* vol)
{
  return std::format("rule '{}' on volume {}", pattern, vol->getSymName());
}

// rigidBody clause: an object DOF name -> value, or an array of the 6 values (TX..RZ)
DOFValues readRigidBodyClause(const json& rb, const Volume* vol, const std::string& pattern)
{
  const auto where = describe(pattern, vol);
  if (!vol->isRigidBodyAllowed()) {
    LOGP(fatal, "Misalignment {}: the volume is not rigid-body alignable", where);
  }
  std::vector<std::pair<std::string, json>> items;
  if (rb.is_array()) {
    if (rb.size() != RigidBodyDOFSet::NDOF) {
      LOGP(fatal, "Misalignment {}: the rigidBody array has {} values instead of {}", where, rb.size(), static_cast<int>(RigidBodyDOFSet::NDOF));
    }
    for (int i = 0; i < RigidBodyDOFSet::NDOF; ++i) {
      items.emplace_back(RigidBodyDOFSet::RigidBodyDOFNames[i], rb[i]);
    }
  } else if (rb.is_object()) {
    for (const auto& [name, val] : rb.items()) {
      items.emplace_back(name, val);
    }
  } else {
    LOGP(fatal, "Misalignment {}: the rigidBody clause must be an object {{DOF: value}} or an array of 6 values", where);
  }
  const auto* dofSet = vol->getRigidBody();
  DOFValues out;
  for (const auto& [name, val] : items) {
    const int idx = RigidBodyDOFSet::dofIndex(name);
    if (idx < 0) {
      LOGP(fatal, "Misalignment {}: unknown rigid-body DOF '{}', expected TX,TY,TZ,RX,RY,RZ", where, name);
    }
    if (!val.is_number() || !std::isfinite(val.get<double>())) {
      LOGP(fatal, "Misalignment {}: the value of {} is not a finite number", where, name);
    }
    if (val.get<double>() == 0. && rb.is_array()) {
      continue; // an array lists all 6 values, a 0 is no request to inject
    }
    if (!dofSet || !dofSet->isFree(idx)) {
      LOGP(fatal, "Misalignment {}: the rigid-body DOF {} is {} in the DOF configuration, the fit cannot recover it: free it in the DOF config or drop it from the misalignment",
           where, name, dofSet ? "fixed" : "not defined (the volume has no rigid-body DOF set)");
    }
    out.emplace_back(idx, val.get<double>());
  }
  return out;
}

struct CalibClause {
  DOFSet::Type type{DOFSet::Type::Legendre};
  bool exact{false};
  DOFValues values;
};

DOFSet::Type calibTypeFromName(const std::string& name, const std::string& where)
{
  if (name == "legendre") {
    return DOFSet::Type::Legendre;
  }
  if (name == "inextensional") {
    return DOFSet::Type::Inextensional;
  }
  if (name == "tpcvdrift") {
    return DOFSet::Type::TPCVDrift;
  }
  if (name == "meanvertex") {
    LOGP(fatal, "Misalignment {}: calib type meanvertex is not supported, the vertex frames do not go through the hierarchy derivatives", where);
  }
  LOGP(fatal, "Misalignment {}: unknown calib type '{}', expected legendre, inextensional or tpcvdrift", where, name);
  return DOFSet::Type::Legendre;
}

// an ITS3 half-barrel: the sensor of the exact model, the parent of the pseudo tiles
bool isIT3HalfBarrel(const Volume* vol)
{
  if (vol->getNChildren() == 0 || vol->getSensorId() < 0 || vol->getSensorId() >= 6) {
    return false;
  }
  bool allPseudo = true;
  vol->forEachChild([&](const Volume* c) { allPseudo = allPseudo && c->isPseudo(); });
  return allPseudo;
}

bool hasMeasurementLeafChild(const Volume* vol)
{
  bool has = false;
  vol->forEachChild([&](const Volume* c) { has = has || c->isLeaf(); });
  return has;
}

// calib clause: {"type": ..., "model": "linear" (default) | "exact", "values": {DOF name: value}}
CalibClause readCalibClause(const json& cal, const Volume* vol, const std::string& pattern, bool exactLegendreAllowed)
{
  const auto where = describe(pattern, vol);
  if (!cal.is_object()) {
    LOGP(fatal, "Misalignment {}: the calib clause must be an object with type, values and optionally model", where);
  }
  for (const auto& [key, val] : cal.items()) {
    if (key != "type" && key != "model" && key != "values") {
      LOGP(fatal, "Misalignment {}: unknown key '{}' in the calib clause, allowed: type, model, values", where, key);
    }
  }
  if (!cal.contains("type") || !cal["type"].is_string() || !cal.contains("values") || !cal["values"].is_object()) {
    LOGP(fatal, "Misalignment {}: the calib clause needs a string 'type' and an object 'values'", where);
  }
  CalibClause out;
  out.type = calibTypeFromName(cal["type"].get<std::string>(), where);
  const auto model = cal.value("model", std::string("linear"));
  if (model != "linear" && model != "exact") {
    LOGP(fatal, "Misalignment {}: unknown model '{}', expected linear or exact", where, model);
  }
  out.exact = model == "exact";

  const auto* dofSet = vol->getCalib();
  if (!dofSet) {
    LOGP(fatal, "Misalignment {}: the volume has no calibration DOF set in the DOF configuration, the fit cannot recover a calib misalignment", where);
  }
  if (dofSet->type() != out.type) {
    LOGP(fatal, "Misalignment {}: the calib type {} differs from the configured one of the volume", where, cal["type"].get<std::string>());
  }
  if (!hasMeasurementLeafChild(vol)) {
    LOGP(fatal, "Misalignment {}: none of the children of the volume is a measurement leaf, its calibration DOFs do not enter the fit (see buildPointGlobals)", where);
  }
  if (out.exact) {
    if (out.type != DOFSet::Type::Legendre && out.type != DOFSet::Type::Inextensional) {
      LOGP(fatal, "Misalignment {}: the exact model exists only for the legendre and inextensional types", where);
    }
    if (!isIT3HalfBarrel(vol)) {
      LOGP(fatal, "Misalignment {}: the exact model is defined only for the ITS3 half-barrel sensors", where);
    }
    if (out.type == DOFSet::Type::Legendre && !exactLegendreAllowed) {
      LOGP(fatal, "Misalignment {}: the exact Legendre model needs the MC truth of the tracks, i.e. MC input (do not use --disable-mc); use the linear model for data", where);
    }
  }
  for (const auto& [name, val] : cal["values"].items()) {
    const int idx = dofSet->dofIndex(name);
    if (idx < 0) {
      LOGP(fatal, "Misalignment {}: unknown DOF '{}' of the configured calibration set, see the names in mp2param.txt (e.g. L(1,0), a_2, alpha, VDRIFT)", where, name);
    }
    if (!val.is_number() || !std::isfinite(val.get<double>())) {
      LOGP(fatal, "Misalignment {}: the value of {} is not a finite number", where, name);
    }
    if (!dofSet->isFree(idx)) {
      LOGP(fatal, "Misalignment {}: the calibration DOF {} is fixed in the DOF configuration, the fit cannot recover it: free it in the DOF config or drop it from the misalignment", where, name);
    }
    out.values.emplace_back(idx, val.get<double>());
  }
  return out;
}

// exact model of an ITS3 half-barrel from the values of the DOFs of its configured calibration set
SensorMisalignment buildExactModel(const DOFSet& dofSet, const DOFValues& values)
{
  SensorMisalignment model;
  if (const auto* leg = dynamic_cast<const LegendreDOFSet*>(&dofSet)) {
    const int n = leg->order() + 1;
    TMatrixD m(n, n);
    for (const auto& [idx, val] : values) {
      int i = 0; // the DOF is L(i,j), see LegendreDOFSet::dofName
      while ((i + 1) * (i + 2) / 2 <= idx) {
        ++i;
      }
      m(i, idx - (i * (i + 1) / 2)) = val;
    }
    model.legendre = o2::math_utils::Legendre2DPolynominal(m);
    model.hasLegendre = true;
  } else if (const auto* inex = dynamic_cast<const InextensionalDOFSet*>(&dofSet)) {
    for (const auto& [idx, val] : values) {
      if (idx == inex->alphaIdx()) {
        model.inextensional.alpha = val;
      } else if (idx == inex->betaIdx()) {
        model.inextensional.beta = val;
      } else {
        model.inextensional.modes[(idx / 4) + 2][idx % 4] = val; // a_n,b_n,c_n,d_n
      }
    }
    model.hasInextensional = true;
  }
  return model;
}

struct CompiledRule {
  std::string pattern;
  geomutils::VolumePattern matcher;
  const json* rule;
  bool matched{false};
};

std::vector<CompiledRule> compileRules(const json& rules, const std::string& jsonPath)
{
  std::vector<CompiledRule> out;
  out.reserve(rules.size());
  for (const auto& rule : rules) {
    if (!rule.is_object()) {
      LOGP(fatal, "Misalignment file {}: every rule must be an object, got {}", jsonPath, rule.dump());
    }
    for (const auto& [key, val] : rule.items()) {
      if (key != "match" && key != "rigidBody" && key != "calib") {
        LOGP(fatal, "Misalignment file {}: unknown key '{}' in the rule {}, allowed: match, rigidBody, calib", jsonPath, key, rule.dump());
      }
    }
    if (!rule.contains("match") || !rule["match"].is_string()) {
      LOGP(fatal, "Misalignment file {}: the rule {} has no string 'match'", jsonPath, rule.dump());
    }
    if (!rule.contains("rigidBody") && !rule.contains("calib")) {
      LOGP(fatal, "Misalignment file {}: the rule '{}' has neither rigidBody nor calib clause, nothing to inject", jsonPath, rule["match"].get<std::string>());
    }
    const auto pattern = rule["match"].get<std::string>();
    try {
      out.push_back({pattern, geomutils::VolumePattern(pattern), &rule});
    } catch (const std::regex_error& e) {
      LOGP(fatal, "Misalignment file {}: invalid regex in the rule match '{}': {}", jsonPath, pattern, e.what());
    }
  }
  return out;
}
} // namespace

MisalignmentConfig MisalignmentConfig::load(Volume* root, const std::string& jsonPath, bool exactLegendreAllowed)
{
  MisalignmentConfig cfg;
  const json rules = geomutils::readJsonRules(jsonPath, "misalignment", false);
  auto compiled = compileRules(rules, jsonPath);

  root->traverse([&](Volume* vol) {
    if (vol->isPseudo() || vol->isRoot()) {
      return;
    }
    // in order: a later matching rule replaces the rigid-body (resp. calib) misalignment set by an earlier one
    const CompiledRule* rbRule = nullptr;
    const CompiledRule* calRule = nullptr;
    for (auto& cr : compiled) {
      if (!cr.matcher.matches(vol->getSymName())) {
        continue;
      }
      cr.matched = true;
      rbRule = cr.rule->contains("rigidBody") ? &cr : rbRule;
      calRule = cr.rule->contains("calib") ? &cr : calRule;
    }
    if (rbRule) {
      for (const auto& [idx, val] : readRigidBodyClause((*rbRule->rule)["rigidBody"], vol, rbRule->pattern)) {
        cfg.addLinear(vol->getLabel().raw(idx), vol->getSymName(), vol->getRigidBody()->dofName(idx), val);
      }
    }
    if (calRule) {
      const auto cal = readCalibClause((*calRule->rule)["calib"], vol, calRule->pattern, exactLegendreAllowed);
      for (const auto& [idx, val] : cal.values) {
        for (const auto& lbl : vol->getCalibLabels()) { // the same value in every calibration time slot
          cfg.addCalib(lbl.raw(idx), vol->getSymName(), vol->getCalib()->dofName(idx), val, cal.exact);
        }
      }
      if (cal.exact) {
        cfg.mExact.emplace(vol, buildExactModel(*vol->getCalib(), cal.values));
      }
    }
  });
  for (const auto& cr : compiled) {
    if (!cr.matched) {
      LOGP(fatal, "Misalignment file {}: the rule '{}' matches no volume of the hierarchy", jsonPath, cr.pattern);
    }
  }
  if (cfg.empty()) {
    LOGP(fatal, "Misalignment file {} injects nothing", jsonPath);
  }
  cfg.print();
  return cfg;
}

void MisalignmentConfig::addLinear(uint32_t label, const std::string& vol, const std::string& dof, double value)
{
  mLinear[label] = value;
  mAllLabels[label] = value;
  mEntries.push_back({vol, dof, value, false});
}

void MisalignmentConfig::addCalib(uint32_t label, const std::string& vol, const std::string& dof, double value, bool exact)
{
  if (!exact) {
    mLinear[label] = value;
  }
  mAllLabels[label] = value;
  mEntries.push_back({vol, dof, value, exact});
}

void MisalignmentConfig::print() const
{
  LOGP(info, "Injected misalignment: {} DOF values, {} of them with the exact ITS3 model", mAllLabels.size(), std::count_if(mEntries.begin(), mEntries.end(), [](const Entry& e) { return e.exact; }));
  for (const auto& e : mEntries) {
    LOGP(info, "  {:<48s} {:<10s} {:+.6e}{}", e.volume, e.dof, e.value, e.exact ? "  (exact model)" : "");
  }
}

MisalignmentShift evaluateLegendreShift(const SensorMisalignment& sensor, const MisalignmentFrame& frame, const TrackSlopes& slopes)
{
  MisalignmentShift shift;
  if (!sensor.hasLegendre) {
    return shift;
  }

  const double gloX = frame.x * std::cos(frame.alpha);
  const double gloY = frame.x * std::sin(frame.alpha);
  const double gloZ = frame.z;
  auto [u, v] = LegendreDOFSet::computeUV(gloX, gloY, gloZ, frame.sensorID, o2::its3::constants::radii[frame.layerID]);
  const double h = sensor.legendre(u, v);

  // this is the shift due to back-projection of the track on the ideal surface
  shift.dy = slopes.dydx * h;
  shift.dz = slopes.dzdx * h;

  if (std::abs(u) > o2::constants::math::Almost0) {
    // account for additional tangential movement due to radial shift
    // we have to approximate the difference in arc-length from the reference pnt on the deformed surface
    // this is done by integrating the height function via Gauss-Legendre quadrature (from Numerical recipes 4.6 [1])
    const auto& x = LegendreDOFSet::GaussNodes;
    const auto& w = LegendreDOFSet::GaussWeights;
    const double mid = 0.5 * u;
    const double half = 0.5 * u;
    double integral = 0.;
    for (int i = 0; i < 8; ++i) {
      const double up = mid + (half * x[i]);
      integral += w[i] * sensor.legendre(up, v);
    }
    integral *= half;
    shift.dy += 0.5 * LegendreDOFSet::getSensorPhiWidth(frame.sensorID, o2::its3::constants::radii[frame.layerID]) * integral;
  }

  const double newGloY = gloY + (shift.dy * std::cos(frame.alpha));
  const double newGloX = gloX - (shift.dy * std::sin(frame.alpha));
  const double newGloZ = gloZ + shift.dz;
  auto [uNew, vNew] = LegendreDOFSet::computeUV(newGloX, newGloY, newGloZ, frame.sensorID, o2::its3::constants::radii[frame.layerID]);
  shift.accepted = std::abs(uNew) <= 1. && std::abs(vNew) <= 1.;
  return shift;
}

MisalignmentShift evaluateInextensionalShift(const SensorMisalignment& sensor, const MisalignmentFrame& frame, const TrackSlopes& slopes)
{
  MisalignmentShift shift;
  if (!sensor.hasInextensional) {
    return shift;
  }

  const double r = o2::its3::constants::radii[frame.layerID];
  const double phi = std::atan2(r * std::sin(frame.alpha), r * std::cos(frame.alpha));
  const double z = frame.z;
  const auto& inex = sensor.inextensional;

  double uz = 0., uphi = 0., ur = 0.;
  for (const auto& [n, coeffs] : inex.modes) {
    const double a_n = coeffs[0], b_n = coeffs[1], c_n = coeffs[2], d_n = coeffs[3];
    const double sn = std::sin(n * phi);
    const double cn = std::cos(n * phi);
    const int n2 = n * n;

    const double fn = (a_n * cn) + (b_n * sn);
    const double fpn = (-n * a_n * sn) + (n * b_n * cn);
    const double fppn = (-n2 * a_n * cn) - (n2 * b_n * sn);
    const double gn = (c_n * cn) + (d_n * sn);
    const double gpn = (-n * c_n * sn) + (n * d_n * cn);

    uz += fn;
    uphi += -(z / r) * fpn + gn;
    ur += (z / r) * fppn - gpn;
  }

  uz += inex.alpha * phi;
  uphi += -(z / r) * inex.alpha + inex.beta * phi;
  ur += -inex.beta;

  shift.dy = -uphi + (slopes.dydx * ur);
  shift.dz = -uz + (slopes.dzdx * ur);
  return shift;
}

} // namespace o2::alignrs
