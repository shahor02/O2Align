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
#include "O2Align/Detector.h"

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
*/
namespace
{
using AlgParVec = std::vector<o2::detectors::AlignParam>;
/// matrices keyed by the TGeo path of the alignable volume they belong to
using PathMatrices = std::unordered_map<std::string, TGeoHMatrix>;

/// TGeo path of an alignable volume, used to navigate the TGeo hierarchy by path prefixes
std::string getAlignablePath(const std::string& symName)
{
  const auto* pne = gGeoManager->GetAlignableEntry(symName.c_str());
  if (pne == nullptr) {
    LOGP(fatal, "Symbolic volume '{}' has no corresponding alignable entry!", symName);
  }
  return pne->GetTitle();
}

/// TGeo level of a path, used to order the volumes parents first
int getPathDepth(const std::string& path)
{
  return static_cast<int>(std::count(path.begin(), path.end(), '/'));
}

/// L2G matrix of an alignable volume in the loaded (pre-aligned) geometry
const TGeoHMatrix& getGeometryMatrix(const std::string& symName)
{
  const auto* m = o2::base::GeometryManager::getMatrix(symName.c_str());
  if (m == nullptr) {
    LOGP(fatal, "Cannot get the matrix of {}", symName);
  }
  return *m;
}

/// Matrix of the TGeo volume `path` (if includeSelf) or else of its nearest TGeo ancestor present in
/// the map, the identity if there is none. Walks up the path component by component, which is
/// O(depth) per lookup rather than a scan over all entries (the ITS alone has ~25k alignables).
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

/// Local delta d^n of the volume from its fitted rigid-body parameters (step 1). Returns false if
/// the volume has no free rigid-body DOF or all its fitted values are 0 (fixed parameters are
/// absent from labelToValue, hence count as 0).
bool getFittedLocalDelta(const Volume* vol, const std::map<uint32_t, double>& labelToValue, TGeoHMatrix& delta)
{
  // RSTODO In principle, if the volume has no rigid-body DOF, it should not have any fitted parameters either. 
  // So, in the alignment results writing session one could skip initializing the same RB DOFs just for this check,
  // instead, one can check isRigidBodyAllowed() and then loop as for (int i = 0; i < RigidBodyDOFSet::NDOF; ++i)
  const auto* rb = vol->getRigidBody();
  if (rb == nullptr || !rb->nFreeDOFs()) {
    return false;
  }
  std::array<double, RigidBodyDOFSet::NDOF> par{};
  bool nonZero = false;
  for (int i = 0; i < rb->nDOFs(); ++i) {
    if (auto it = labelToValue.find(vol->getLabel().raw(i)); it != labelToValue.end()) {
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
    if (getFittedLocalDelta(vol, labelToValue, delta)) {
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

/// global delta D^p equivalent to an initial object of either convention (step 2)
TGeoHMatrix getInitialGlobalDelta(const o2::detectors::AlignParam& par)
{
  auto delta = par.createMatrix();
  if (par.isGlobal()) {
    return delta;
  }
  const auto& gp = getGeometryMatrix(par.getSymName());
  return gp * delta * gp.Inverse();
}

/// Step 2: cumulative global delta P^p of the initial alignment for every volume it contains, the
/// objects being accumulated parents first, like GeometryManager::applyAlignment applies them
PathMatrices collectInitialCumulativeDeltas(const AlgParVec* initial)
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

} // namespace

std::vector<o2::detectors::AlignParam> Detector::MP2AlignParams(const std::map<uint32_t, double>& labelToValue, bool writeLocal) const
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
  const auto entries = selectOutputVolumes(mTopVolume, iniAlg, newCum);

  // step 3, parents being processed first so that P^c_{a(j)} is known when j is reached
  PathMatrices combinedCum;
  AlgParVec result;
  result.reserve(entries.size());
  for (const auto& e : entries) {
    const auto pp = getOnPath(iniCum, e.path, true);
    const auto pc = getOnPath(newCum, e.path, true) * pp;
    const auto pcParent = getOnPath(combinedCum, e.path, false);
    combinedCum[e.path] = pc;
    TGeoHMatrix delta = pc * pcParent.Inverse();
    if (writeLocal) {
      // d^c = (Z^c)^-1 * D^c * Z^c, with Z^c = P^c_{a(j)} * G^id_j the matrix of j before its own delta
      const TGeoHMatrix z = pcParent * pp.Inverse() * getGeometryMatrix(e.symName);
      delta = z.Inverse() * delta * z;
    }
    result.emplace_back(e.symName.c_str(), e.algID, delta, !writeLocal, false);
  }
  LOGP(info, "Converted rigid-body corrections of {} to {} {} AlignParam objects ({} in the initial alignment)", getDetName(), result.size(), writeLocal ? "local" : "global", iniAlg ? iniAlg->size() : 0);
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

} // namespace o2::alignrs
