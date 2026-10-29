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
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <TGeoManager.h>
#include <TGeoMatrix.h>

#include "Framework/Logger.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "O2Align/AlignGeomUtils.h"
#include "O2Align/Detector.h"
#include "O2Align/Params.h"
#include "DetectorsBase/Propagator.h"

namespace o2::alignrs
{

/*
  Conversion of the fitted rigid-body corrections to o2::detectors::AlignParam (Felix Schlepper thesis A.7, A.8)
  ---------------------------------------------------------------------------------------------

  Notation, for a TGeo alignable volume j:
    G^id_j      ideal L2G matrix
    G^p_j       L2G matrix of the geometry the fit was done on (ideal + initial alignment)
    a(j)        nearest TGeo ancestor of j carrying an alignment object
    D_j         global delta of j: AlignParam::applyToGeometry applies it on top of the already
                aligned ancestors (GeometryManager::applyAlignment sorts by level), G <- D_j * G
    d_j         local delta of j: applyToGeometry sets the local matrix to L^id_j * d_j
    P_j         cumulative global delta, P_j = D_j * P_{a(j)}, such that G_j = P_j * G^id_j
  Superscripts p, n, c stand for the initial, the new (fitted) and the combined corrections.

  1) New corrections. The fitted parameters (TX,TY,TZ,RX,RY,RZ) of volume j are the displacement
     of the volume in the frame of the L2G matrix G^p_j in which its derivatives were computed
     (Volume::getMatrixL2G; for SensorITS it is offset from the chip frame by the epitaxial layer
     correction, which the conjugation below absorbs). The derivatives being those of the
     prediction (A.57), the fitted values need no sign flip: the corrected volume is G^p_j * d^n_j,
     d^n_j being built with the angle convention R = Rx(RX)*Ry(RY)*Rz(RZ) of A.47, the one of
     AlignParam. The volume moves with its hierarchy parents first, then by d^n_j in its own frame
     (A.75-A.77), so that the cumulative global delta relative to the pre-aligned geometry is
        P^n_j = P^n_{par(j)} * G^p_j * d^n_j * (G^p_j)^-1
     A TGeo alignable which is not part of the hierarchy follows rigidly its nearest TGeo ancestor
     which is. Virtual volumes have no TGeo counterpart: their correction is only transported to
     their non-virtual descendants.
  2) Initial alignment. P^p_j = D^p_j * P^p_{a(j)}. An initial local delta d^p_j was applied as
     G^p_j = Z_j * d^p_j with Z_j = P^p_{a(j)} * G^id_j, its global equivalent is (A.77)
        D^p_j = Z_j * d^p_j * Z_j^-1 = G^p_j * d^p_j * (G^p_j)^-1        (Z_j = G^p_j * (d^p_j)^-1)
     with G^p_j the TGeo matrix of the loaded geometry.
  3) Combination. Since G^c_j = P^n_j * G^p_j = P^n_j * P^p_j * G^id_j:
        P^c_j = P^n_j * P^p_j                                            (A.79)
        D^c_j = P^c_j * (P^c_{a(j)})^-1                                  (A.80)
     The output convention is chosen by the caller (Params::writeLocalAlignParams) for all objects,
     independently of the convention of the initial ones. In the local one, G^c_j = Z^c_j * d^c_j
     with Z^c_j = P^c_{a(j)} * G^id_j = P^c_{a(j)} * (P^p_j)^-1 * G^p_j, hence
        d^c_j = (Z^c_j)^-1 * D^c_j * Z^c_j
     With no new correction this gives back the initial objects, converted to the requested convention.

  The ideal geometry is never needed: G^id_j = (P^p_j)^-1 * G^p_j wherever it enters.

  4) Optional regularization (OutputOpt::Regularize), geomutils::transferChildrenCommonMotion: before the
     conversion to objects, the P^c of the parents are modified to absorb the common (robust) motion of their
     TGeo alignable children, while the P^c of the children, hence the placement of every sensor, are
     kept. The volumes made explicit by this (a parent or a child without object) are added to the
     output, and the sensor placement given by the final objects is checked against the one before.
*/
namespace
{
using namespace o2::alignrs::geomutils;
using AlgParVec = std::vector<o2::detectors::AlignParam>;

bool isSameMatrix(const TGeoHMatrix& a, const TGeoHMatrix& b, double eps = 1e-12)
{
  for (int i = 0; i < 9; ++i) {
    if (std::abs(a.GetRotationMatrix()[i] - b.GetRotationMatrix()[i]) > eps) {
      return false;
    }
  }
  for (int i = 0; i < 3; ++i) {
    if (std::abs(a.GetTranslation()[i] - b.GetTranslation()[i]) > eps) {
      return false;
    }
  }
  return true;
}

/// Step 1: cumulative global delta P^n of the new corrections for every non-virtual volume of the
/// branch below top (included). The branch is traversed parents first, the parent of top being
/// the common root, which has no DOFs: the corrections of different detectors never mix.
PathMatrices collectNewCumulativeDeltas(Volume* top, const std::map<uint32_t, double>& labelToValue)
{
  std::unordered_map<const Volume*, TGeoHMatrix> cumulative;
  PathMatrices byPath;
  top->traverse([&](Volume* vol) {
    TGeoHMatrix pn = vol == top ? TGeoHMatrix() : cumulative.at(vol->getParent());
    TGeoHMatrix delta;
    if (vol->getFittedLocalDelta(labelToValue, delta)) {
      const auto& l2g = vol->getMatrixL2G();
      pn = pn * l2g * delta * l2g.Inverse();
      if (vol->isVirtual() && !vol->isLeaf()) {
        LOGP(warn, "Virtual volume {} has no geometry counterpart, its rigid-body correction is transported to its descendants", vol->getSymName());
      }
    }
    if (!vol->isVirtual()) {
      byPath[getAlignablePath(vol->getSymName())] = pn;
    }
    cumulative[vol] = pn;
  });
  return byPath;
}

/// volume for which an alignment object is to be written
struct OutputEntry {
  std::string symName;
  std::string path;
  int algID{-1};
};

/// Volumes to write: all those of the initial alignment, then the volumes of the branch which the
/// new corrections move differently from their TGeo parent. Ordered parents first, as needed by step 3.
std::vector<OutputEntry> selectOutputVolumes(Volume* top, const AlgParVec* initial, const PathMatrices& newCum)
{
  std::vector<OutputEntry> entries;
  std::unordered_set<std::string> selected;
  if (initial) {
    for (const auto& par : *initial) {
      auto& e = entries.emplace_back(par.getSymName(), getAlignablePath(par.getSymName()), par.getAlignableID());
      selected.insert(e.path);
    }
  }
  top->traverse([&](Volume* vol) {
    if (vol->isVirtual()) {
      return;
    }
    auto path = getAlignablePath(vol->getSymName());
    if (!selected.contains(path) && !isSameMatrix(newCum.at(path), getOnPath(newCum, path, false))) {
      selected.insert(path);
      // the alignable ID of a sensor is not known here: the DET field of the Label is the index of
      // the Detector plugin, not the DetID. GeometryManager::applyAlignment does not use it.
      entries.emplace_back(vol->getSymName(), std::move(path), -1);
    }
  });
  std::stable_sort(entries.begin(), entries.end(), [](const OutputEntry& a, const OutputEntry& b) { return getPathDepth(a.path) < getPathDepth(b.path); });
  return entries;
}

/// Step 3: combined cumulative delta P^c of every output volume
PathMatrices combineCumulativeDeltas(const std::vector<OutputEntry>& entries, const PathMatrices& iniCum, const PathMatrices& newCum)
{
  PathMatrices combinedCum;
  for (const auto& e : entries) {
    combinedCum[e.path] = getOnPath(newCum, e.path, true) * getOnPath(iniCum, e.path, true);
  }
  return combinedCum;
}

/// Step 3: the alignment objects of the output volumes (ordered parents first), D^c_j = P^c_j * (P^c_{a(j)})^-1,
/// converted to the local convention if requested
AlgParVec makeAlignParams(const std::vector<OutputEntry>& entries, const PathMatrices& combinedCum, const PathMatrices& iniCum, bool writeLocal)
{
  AlgParVec result;
  result.reserve(entries.size());
  for (const auto& e : entries) {
    const auto& pc = combinedCum.at(e.path);
    const auto pcParent = getOnPath(combinedCum, e.path, false);
    TGeoHMatrix delta = pc * pcParent.Inverse();
    if (writeLocal) {
      // d^c = (Z^c)^-1 * D^c * Z^c, with Z^c = P^c_{a(j)} * G^id_j the matrix of j before its own delta
      const TGeoHMatrix z = pcParent * getOnPath(iniCum, e.path, true).Inverse() * getGeometryMatrix(e.symName);
      delta = z.Inverse() * delta * z;
    }
    result.emplace_back(e.symName.c_str(), e.algID, delta, !writeLocal, false);
  }
  return result;
}

/// paths of the non-virtual volumes of the branch below top (included)
std::unordered_set<std::string> collectBranchPaths(Volume* top)
{
  std::unordered_set<std::string> paths;
  top->traverse([&](Volume* vol) {
    if (!vol->isVirtual()) {
      paths.insert(getAlignablePath(vol->getSymName()));
    }
  });
  return paths;
}

/// ideal L2G of every alignable, G^id_j = (P^p_j)^-1 * G^p_j
PathMatrices collectIdealMatrices(const std::vector<Alignable>& alignables, const PathMatrices& iniCum)
{
  PathMatrices ideal;
  for (const auto& a : alignables) {
    ideal[a.path] = getOnPath(iniCum, a.path, true).Inverse() * getGeometryMatrix(a.symName);
  }
  return ideal;
}

/// Step 4: add to the output the volumes which got an explicit cumulative delta, keeping the order parents first
void addRegularizedEntries(std::vector<OutputEntry>& entries, const std::vector<Alignable>& alignables, const PathMatrices& combinedCum)
{
  std::unordered_set<std::string> selected;
  for (const auto& e : entries) {
    selected.insert(e.path);
  }
  for (const auto& a : alignables) {
    if (combinedCum.contains(a.path) && !selected.contains(a.path)) {
      entries.emplace_back(a.symName, a.path, -1);
    }
  }
  std::stable_sort(entries.begin(), entries.end(), [](const OutputEntry& a, const OutputEntry& b) { return getPathDepth(a.path) < getPathDepth(b.path); });
}

/// Step 4 closure: the cumulative deltas obtained by applying the output objects (parents first) to the
/// ideal geometry must place every sensor (alignable without alignable children) as before the regularization
void checkSensorPlacement(const AlgParVec& result, const PathMatrices& before, const std::vector<Alignable>& alignables, const PathMatrices& ideal, const std::string& detName)
{
  PathMatrices applied;
  for (const auto& par : result) {
    const auto path = getAlignablePath(par.getSymName());
    const auto pa = getOnPath(applied, path, false);
    auto delta = par.createMatrix();
    if (!par.isGlobal()) {
      const TGeoHMatrix z = pa * ideal.at(path);
      delta = z * delta * z.Inverse();
    }
    applied[path] = delta * pa;
  }
  const auto children = collectAlignableChildren(alignables);
  double maxShift = 0., maxRot = 0.;
  size_t nSensors = 0;
  for (const auto& a : alignables) {
    if (children.contains(a.path)) {
      continue;
    }
    const TGeoHMatrix diff = getOnPath(applied, a.path, true) * ideal.at(a.path) * (getOnPath(before, a.path, true) * ideal.at(a.path)).Inverse();
    for (int i = 0; i < 3; ++i) {
      maxShift = std::max(maxShift, std::abs(diff.GetTranslation()[i]));
    }
    for (int i = 0; i < 9; ++i) {
      maxRot = std::max(maxRot, std::abs(diff.GetRotationMatrix()[i] - (i % 4 ? 0. : 1.)));
    }
    ++nSensors;
  }
  if (maxShift > 1e-9 || maxRot > 1e-11) {
    LOGP(error, "Regularization of {} changed the sensor placement: max deviation over {} sensors {:.3e} cm, {:.3e} in rotation", detName, nSensors, maxShift, maxRot);
  } else {
    LOGP(info, "Regularization of {} keeps the sensor placement: max deviation over {} sensors {:.3e} cm, {:.3e} in rotation", detName, nSensors, maxShift, maxRot);
  }
}

} // namespace

std::vector<o2::detectors::AlignParam> Detector::MP2AlignParams(const std::map<uint32_t, double>& labelToValue, bool writeLocal, const geomutils::ChildrenTransferConfig* regularize) const
{
  if (getO2DetID() < 0) {
    LOGP(fatal, "MP2AlignParams called for {} which has no DetID", getDetName());
  }
  if (!mTopVolume) {
    LOGP(fatal, "MP2AlignParams called for {} before attachTo or for a discarded branch", getDetName());
  }
  auto iniAlg = o2::base::GRPGeomHelper::instance().getAlignment(getO2DetID());
  if (!iniAlg) {
    LOGP(warn, "Detector {} has no initial alignment, the MP2AlignParams output will contain only the fitted corrections", getDetName());
  } else {
    LOGP(info, "Building MP2AlignParams for detector {} on top of initial alignment", getDetName());
  }
  // a branch made only of virtual volumes (TPC, mean vertex) produces no new correction, its
  // initial objects, if any, are still rewritten in the requested convention
  const auto newCum = collectNewCumulativeDeltas(mTopVolume, labelToValue);
  const auto iniCum = collectInitialCumulativeDeltas(iniAlg);
  auto entries = selectOutputVolumes(mTopVolume, iniAlg, newCum);
  auto combinedCum = combineCumulativeDeltas(entries, iniCum, newCum);
  AlgParVec result;
  if (regularize && !entries.empty()) {
    const auto alignables = collectAlignablesBelow(collectBranchPaths(mTopVolume));
    const auto ideal = collectIdealMatrices(alignables, iniCum);
    const auto before = combinedCum;
    transferChildrenCommonMotion(combinedCum, alignables, ideal, *regularize);
    addRegularizedEntries(entries, alignables, combinedCum);
    result = makeAlignParams(entries, combinedCum, iniCum, writeLocal);
    checkSensorPlacement(result, before, alignables, ideal, getDetName());
  } else {
    result = makeAlignParams(entries, combinedCum, iniCum, writeLocal);
  }
  LOGP(info, "Converted rigid-body corrections of {} to {} {} AlignParam objects ({} in the initial alignment){}", getDetName(), result.size(), writeLocal ? "local" : "global", iniAlg ? iniAlg->size() : 0, regularize ? ", regularized" : "");
  return result;
}

Volume* Detector::attachTo(Volume* root, Volume::SensorMapping& sensorMap)
{
  if (root == nullptr) {
    LOGP(fatal, "Cannot attach the hierarchy of {} to a null root", getDetName());
  }
  auto top = buildHierarchy(sensorMap);
  if (!top) {
    LOGP(fatal, "Detector {} provided no hierarchy", getDetName());
  }
  // a branch which provides no measurement at all is not attached: this happens when the geometry
  // has no volume of this detector (all its supermodules disabled), leaving an envelope alone.
  // A detector consisting of a single sensor (the mean vertex) is perfectly legitimate: it needs no
  // envelope above it, its only volume being directly the top of its branch.
  bool hasSensor = false;
  top->traverse([&hasSensor](Volume* vol) { hasSensor = hasSensor || vol->getLabel().sens(); });
  if (!hasSensor) {
    LOGP(warn, "Detector {} has no sensor volume, its branch {} is discarded", getDetName(), top->getSymName());
    return nullptr;
  }
  mTopVolume = root->adoptChild(std::move(top));
  LOGP(info, "Attached the branch {} of {} to {}", mTopVolume->getSymName(), getDetName(), root->getSymName());
  return mTopVolume;
}

const std::string& Detector::getTimeSlotsJson() const
{
  static const std::string noSlots{}; // a detector without time-sliced calibration
  return noSlots;
}

void Detector::loadTimeSlots()
{
  const auto& jsonPath = getTimeSlotsJson();
  if (jsonPath.empty()) {
    return;
  }
  mTimeSlots = std::make_unique<TimeSlotsSet>();
  if (mTimeSlots->readSlotsFromFile(jsonPath) < 0) {
    LOGP(fatal, "Failed to load the calibration time slots of {} from {}", getDetName(), jsonPath);
  }
  LOGP(info, "{} will be calibrated in {} time slots of {}", getDetName(), mTimeSlots->slots.size(), jsonPath);
}

bool Detector::setTimeStamp(long tsMS, bool ignoreMismatch)
{
  if (!mTimeSlots) { // no slots: everything belongs to the single slot 0
    return false;
  }
  const int slotID = mTimeSlots->getSlotID(tsMS);
  if (slotID < 0) {
    if (!ignoreMismatch) {
      LOGP(fatal, "Timestamp {} is not covered by any calibration slot of {} from {}", tsMS, getDetName(), getTimeSlotsJson());
    } else {
      LOGP(warn, "Timestamp {} is not covered by any calibration slot of {} from {}, ignoring has requested", tsMS, getDetName(), getTimeSlotsJson());
      return false;
    }
  }
  if (slotID == mSlotID) {
    return false;
  }
  mSlotID = slotID;
  onSlotChange(slotID);
  return true;
}

std::string Detector::reportDOFSummary() const
{
  if (mTopVolume == nullptr) {
    return {};
  }
  struct Counts {
    int nRB{0}, nCal{0}, nTot{0};
    void add(const Volume* v)
    {
      ++nTot;
      nRB += v->getRigidBody() && v->getRigidBody()->nFreeDOFs() > 0;
      nCal += v->getCalib() && v->getCalib()->nFreeDOFs() > 0;
    }
  };
  const int topLevel = mTopVolume->getLevel();
  std::map<int, Counts> byLevel;
  mTopVolume->traverse([&byLevel, topLevel](Volume* v) { byLevel[v->getLevel() - topLevel].add(v); });
  std::ostringstream oss;
  oss << getDetName() << ":";
  for (const auto& [lvl, c] : byLevel) {
    oss << " level" << lvl << ": " << c.nRB << '/' << c.nCal << '/' << c.nTot << ',';
  }
  return oss.str();
}

// The KF (re)fit of the prepared track moved to Track::fitTrack / Track::continueFitOutward

bool Detector::propagateCosmicRef(const o2::track::TrackParCovD& ref, double alpha, double x, o2::track::TrackParCovD& trk)
{
  const auto& params = Params::Instance();
  trk = ref;
  return o2::base::PropagatorD::Instance()->propagateToAlphaX(trk, nullptr, alpha, x, false, params.maxSnp, params.maxStep, 1, o2::base::PropagatorD::MatCorrType::USEMatCorrNONE);
}

} // namespace o2::alignrs
