#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TGeoMatrix.h>
#include <TH1F.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <regex>
#include <string>
#include <unordered_map>

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSBase/GeometryTGeo.h"
#endif

/*
  Alignment formalism:

  Vector l in the local frame of the volume_j (assuming hierarchy of nested
  volumes 0...J from most coarse to the end volume) is transformed to master
  frame vector g = G_j*l_j Matrix G_j is Local2Global matrix (L2G in the code).
  If the volume has a parent volume j-1, the global vector g can be transformed
  to the local volume of j-1 as l_{j-1} = G^-1_{j-1}* g Hence, the
  transformation from volume j to j-1 is l_{j-1} = G^-1_{j-1}*G_j l_j = R_j*l_j

  The alignment corrections in general can be defined either as a

  1) local delta:   l'_j = delta_j * l_j
  hence g'  = G_j * delta_j = G'_j*l_j
  or as
  2) global Delta:  g' = Delta_j * G_j * l_j = G'_j*l_j

  Hence Delta and delta are linked as
  Delta_j = G_j delta_j G^-1_j
  delta_j = G^-1_j Delta_j G_j

  In case the whole chain of nested volumes is aligned, the corrections pile-up
  as:

  G_0*delta_0 ... G^-1_{j-2}*G_{j-1}*delta_{j-1}*G^-1_{j-1}*G_j*delta_j =
  Delta_0 * Delta_{1} ... Delta_{j-1}*Delta_{j}... * G_j

  From this by induction one gets relation between local and global deltas:

  Delta_j = Z_j * delta_j * Z^-1_j

  where Z_j = [ Prod_{k=0}^{j-1} (G_k * delta_k * G^-1_k) ] * G_j

  By convention, aliroot alignment framework stores global Deltas !

  In case the geometry was already prealigned by PDelta_j matrices, the result
  of the new incremental alignment Delta_j must be combined with PDelta_j to
  resulting matrix TDelta_j before writing new alignment object.

  Derivation: if G_j and IG_j are final and ideal L2G matrices for level j, then

  G_j = TDelta_j * TDelta_{j-1} ... TDelta_0 * IG_j
  =     (Delta_j * Delta_{j-1} ... Delta_0)  * (PDelta_j * PDelta_{j-1} ...
  PDelta_0) * IG_j

  Hence:
  TDelta_j = [Prod_{i=j}^0 Delta_i ] * [Prod_{k=j}^0 PDelta_k ] *
  [Prod_{l=0}^{j-1} TDelta_l]

  By induction we get combination rule:

  TDelta_j = Delta_j * X_{j-1} * PDelta_j * X^-1_{j-1}

  where X_i = Delta_i * Delta_{i-1} ... Delta_0

  ---------------------------------

  This alignment framework internally allows to find geometry corrections either
  in the volume LOCAL frame or in its TRACKING frame. The latter is defined for
  sensors as lab frame, rotated by the angle alpha in such a way that the X axis
  is normal to the sensor plane (note, that for ITS the rotated X axis origin is
  also moved to the sensor) For the non-sensor volumes the TRACKING frame is
  defined by rotation of the lab frame with the alpha angle = average angle of
  centers of its children, seen from the origin.

  The TRACKING and IDEAL LOCAL (before misalignment) frames are related by the
  tracking-to-local matrix (T2L in the code), i.e. the vectors in local and
  tracking frames are related as l = T2L * t

  The alignment can be done using both frames for different volumes of the same
  geometry branch. The alignments deltas in local and tracking frames are
  related as:

  l' = T2L * delta_t * t
  l' = delta_l * T2L * t
  -> delta_l = T2L * delta_t * T2L^-1

 */

using AlgParams = std::vector<o2::detectors::AlignParam>;

void loadGeom(std::string const& name) {
  TFile f(name.c_str());
  f.Get("ccdb_object");
}

void matrixToAngles(const double* rot, double& psi, double& theta,
                    double& phi) {
  if (std::abs(rot[0]) < 1e-7 || std::abs(rot[8]) < 1e-7) {
    std::cerr << "Failed to extract roll-pitch-yall angles!";
    exit(41);
  }
  psi = std::atan2(-rot[5], rot[8]);
  theta = std::asin(rot[2]);
  phi = std::atan2(-rot[1], rot[0]);
}

void anglesToMatrix(double psi, double theta, double phi, double* rot) {
  double sinpsi = std::sin(psi);
  double cospsi = std::cos(psi);
  double sinthe = std::sin(theta);
  double costhe = std::cos(theta);
  double sinphi = std::sin(phi);
  double cosphi = std::cos(phi);

  rot[0] = costhe * cosphi;
  rot[1] = -costhe * sinphi;
  rot[2] = sinthe;
  rot[3] = sinpsi * sinthe * cosphi + cospsi * sinphi;
  rot[4] = -sinpsi * sinthe * sinphi + cospsi * cosphi;
  rot[5] = -costhe * sinpsi;
  rot[6] = -cospsi * sinthe * cosphi + sinpsi * sinphi;
  rot[7] = cospsi * sinthe * sinphi + sinpsi * cosphi;
  rot[8] = costhe * cospsi;
}

void convertZXZEulerToXYZCardano(double rotA, double rotB, double rotC,
                                 double& roll, double& pitch, double& yaw) {
  static TRotation aRot;  // this how Alex minimizes his parameter
  aRot.SetXEulerAngles(rotA, rotB, rotC);
  static double rota[9];
  for (int i{0}; i < 3; ++i) {
    for (int j{0}; j < 3; ++j) {
      rota[i * 3 + j] = aRot(i, j);
    }
  }
  matrixToAngles(rota, roll, pitch, yaw);
}

double rectify(double v, double zero = 1e-12) {
  if (std::abs(v) < zero) {
    return 0.;
  }
  return v;
}

using DeltaTree = boost::property_tree::ptree;
using Delta = std::array<double, 6>;
void store(DeltaTree& node, const std::string& path, const Delta& par) {
  DeltaTree anode;
  for (const auto& v : par) {
    DeltaTree item{};
    item.put("", rectify(v));
    anode.push_back(std::make_pair("", item));
  }
  node.put_child(path, anode);
}

Delta load(const DeltaTree& node, const std::string& path) {
  Delta a;
  const auto& array_node = node.get_child(path);
  size_t i = 0;
  for (const auto& item : array_node) {
    a[i++] = rectify(item.second.get_value<double>());
    if (i > 5) {
      break;
    }
  }
  return a;
}

Delta toArray(const TGeoHMatrix& m) {
  Delta del;
  const Double_t *tra, *rot;
  tra = m.GetTranslation();
  std::memcpy(&del[0], &tra[0], 3 * sizeof(Double_t));
  rot = m.GetRotationMatrix();
  matrixToAngles(rot, del[3], del[4], del[5]);
  return del;
}

TGeoHMatrix toMatrix(const Delta& del) {
  TGeoHMatrix m;
  Double_t rot[9];
  anglesToMatrix(del[3], del[4], del[5], rot);
  m.SetRotation(rot);
  m.SetTranslation(del.data());
  return m;
}

static std::unordered_map<std::string, std::string> cacheSym2Path;
void cacheBuildSym2Path(const char* geom) {
  LOGP(info, "Building sym <-> path cache!");
  loadGeom(geom);
  auto gman = o2::its::GeometryTGeo::Instance();  // just for navigation!

  auto add = [&](const char* sym) {
    TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(sym);
    cacheSym2Path[sym] = pne->GetPath();
    LOGP(debug, "From '{}' to '{}'", sym, pne->GetPath());
  };

  add(gman->composeSymNameITS());
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      add(gman->composeSymNameHalfBarrel(ilr, ihb));
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        add(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          add(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            add(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }
  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    add(o2::base::GeometryManager::getSymbolicName(0, ich));
  }
}

static std::unordered_map<std::string, TGeoHMatrix> cacheLocalIdeal;
void cacheBuildLocalIdeal(const char* geom) {
  LOGP(info, "Building local ideal matrix cache!");
  if (gSystem->AccessPathName(geom)) {  // if needed, create geometry
    std::cout << geom << " does not exist. Will create it on the fly\n";
    std::stringstream str;
    // constructing an **unaligned** geom (Geant3 used since faster
    // initialization) --> can be avoided by passing an existing geometry
    str << "${O2_ROOT}/bin/o2-sim-serial -n 0 -e TGeant3 --field 0";
    gSystem->Exec(str.str().c_str());
  }
  loadGeom(geom);
  auto gman = o2::its::GeometryTGeo::Instance();  // just for navigation!

  auto add = [&](const char* sym) {
    if (!gGeoManager->cd(cacheSym2Path[sym].c_str())) {
      LOGP(warn, "From '{}' to '{}' is not valid!", sym, cacheSym2Path[sym]);
    }
    const auto node = gGeoManager->GetCurrentNode();
    cacheLocalIdeal[sym] = node->GetMatrix();
  };

  add(gman->composeSymNameITS());
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      add(gman->composeSymNameHalfBarrel(ilr, ihb));
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        add(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          add(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            add(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }
  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    add(o2::base::GeometryManager::getSymbolicName(0, ich));
  }
}

static std::unordered_map<std::string, Delta> gLocalDelta;
static std::unordered_map<std::string, Delta> gLocalDeltaDelta;

// Extracts local deltas of ITS geometry.
// Returns n-ary tree of local deltas
DeltaTree getLocalDeltas(const std::string& geom) {
  LOGP(info, "Getting local deltas '{}'", geom);
  loadGeom(geom);
  DeltaTree deltas;
  std::unordered_map<std::string, TGeoHMatrix> cache;  // local delta cache
  auto gman = o2::its::GeometryTGeo::Instance();       // just for navigation!

  auto addDelta = [&](const char* sym) {
    TGeoHMatrix orig, *cur, delta, z;
    if (!(cur = o2::base::GeometryManager::getMatrix(sym))) {
      LOGP(fatal, "cannot get aligned matrix");  // this should not happen
    }
    if (!o2::base::GeometryManager::getOriginalMatrix(sym, orig)) {
      LOGP(fatal, "cannot get original matrix");  // this should not happen
    }
    // local delta is Del_loc,j = G^-1_orig,j * G_cur,j (including parent
    // cont.!)
    delta = orig.Inverse();
    delta.Multiply(cur);

    // subtract with del_loc,j = Del_loc,j * del_loc,j-1^-1 * ... * del_loc,0^-1
    //               del_loc,j = Del_loc,j * z^-1
    // -> z = del_loc,0 * ... * del_loc,j-1
    const std::regex patLayer{"ITSULayer\\d+$"};
    std::string spath{sym};
    size_t pos;
    while ((pos = spath.find_last_of('/')) != std::string::npos) {
      spath = spath.substr(0, pos);
      if (std::regex_search(spath, patLayer)) {
        continue;
      }
      if (cache.find(spath) == cache.end()) {
        LOGP(fatal, "No cached for {}", spath);
      }
      // Get parent's original L2G matrix
      TGeoHMatrix G_parent_orig;
      if (!o2::base::GeometryManager::getOriginalMatrix(spath.c_str(),
                                                        G_parent_orig)) {
        LOGP(fatal, "Original matrix missing for parent: {}", spath);
      }

      // Transform parent delta to child's frame: Δ_parent_child = G_child⁻¹ *
      // G_parent * Δ_parent * G_parent⁻¹ * G_child
      // Works but prone to float-imp
      TGeoHMatrix G_child_orig = orig,
                  G_child_orig_inv = G_child_orig.Inverse();
      TGeoHMatrix deltaParentInChildFrame =
          G_child_orig_inv * G_parent_orig * cache[spath] *
          G_parent_orig.Inverse() * G_child_orig;
      // -> Since G_child_orig_inv * G_parent_orig = L_j^-1 TODO
      // TGeoHMatrix deltaParentInChildFrame = cacheLocalI[]
      // Apply inverse of parent delta in child's frame
      delta *= deltaParentInChildFrame.Inverse();
    }
    cache[sym] = delta;
    const auto adelta = toArray(delta);
    std::string path{sym};
    gLocalDelta[path] = adelta;
    gLocalDeltaDelta[path] = adelta;
    std::replace(path.begin(), path.end(), '/', '.');
    store(deltas, path, adelta);
  };

  addDelta(gman->composeSymNameITS());
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      addDelta(gman->composeSymNameHalfBarrel(ilr, ihb));
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        addDelta(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          addDelta(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            addDelta(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }

  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    addDelta(o2::base::GeometryManager::getSymbolicName(0, ich));
  }

  boost::property_tree::write_json("local.json", deltas);
  return deltas;
}

// build updated local delta tree
DeltaTree getLocalDeltas(
    const DeltaTree& locDeltas,
    const std::unordered_map<std::string, TGeoHMatrix>& updatedLocDeltas = {}) {
  DeltaTree deltas;
  std::unordered_map<std::string, TGeoHMatrix> cache;  // global delta cache
  auto gman = o2::its::GeometryTGeo::Instance();       // just for navigation

  auto getLocalDelta = [&](const char* path) {
    std::string spath{path}, dpath{path};
    if (auto it = updatedLocDeltas.find(spath); it != updatedLocDeltas.end()) {
      return it->second;
    } else {
      // get local delta from tree
      std::replace(dpath.begin(), dpath.end(), '/', '.');
      const auto ldelta = load(locDeltas, dpath);
      return toMatrix(ldelta);
    }
  };

  auto addDelta = [&](const char* p) {
    auto delta = getLocalDelta(p);
    std::string path{p};
    std::replace(path.begin(), path.end(), '/', '.');
    store(deltas, path, toArray(delta));
  };

  addDelta(gman->composeSymNameITS());
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      addDelta(gman->composeSymNameHalfBarrel(ilr, ihb));
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        addDelta(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          addDelta(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            addDelta(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }

  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    addDelta(o2::base::GeometryManager::getSymbolicName(0, ich));
  }

  boost::property_tree::write_json("local.json", deltas);
  return deltas;
}

// Build global deltas from local delta tree
DeltaTree getGlobalDeltas(
    const DeltaTree& locDeltas,
    const std::unordered_map<std::string, TGeoHMatrix>& updatedLocDeltas = {}) {
  DeltaTree deltas;
  std::unordered_map<std::string, TGeoHMatrix> cache;  // global delta cache
  auto gman = o2::its::GeometryTGeo::Instance();       // just for navigation

  auto getLocalDelta = [&](const char* path) {
    std::string spath{path}, dpath{path};
    if (auto it = updatedLocDeltas.find(spath); it != updatedLocDeltas.end()) {
      return it->second;
    } else {
      // get local delta from tree
      std::replace(dpath.begin(), dpath.end(), '/', '.');
      const auto ldelta = load(locDeltas, dpath);
      return toMatrix(ldelta);
    }
  };

  auto getGlobalDelta = [&](const char* path) {
    std::string spath{path}, dpath{path};
    auto delta = getLocalDelta(path);

    size_t pos;
    TGeoHMatrix z, orig;
    const std::regex patLayer{"ITSULayer\\d+$"};
    while ((pos = spath.find_last_of('/')) != std::string::npos) {
      spath = spath.substr(0, pos);
      if (std::regex_search(spath, patLayer)) {
        continue;
      }
      if (cache.find(spath) != cache.end()) {
        z.MultiplyLeft(cache[spath]);
      } else {
        LOGP(fatal, "No cached for {}", spath);
      }
    }
    const auto zi = z.Inverse();
    if (!o2::base::GeometryManager::getOriginalMatrix(path, orig)) {
      LOGP(fatal, "cannot get original matrix");  // this should not happen
    }
    orig.MultiplyLeft(z);  // Z = parent_Z · G_j_ideal
    const TGeoHMatrix origi = orig.Inverse();
    delta.MultiplyLeft(orig);  // Z · Δ_local
    delta.Multiply(origi);     // Z · Δ_local · Z⁻¹
    cache[path] = delta;
    return delta;
  };

  auto addDelta = [&](const char* p) {
    auto delta = getGlobalDelta(p);
    std::string path{p};
    std::replace(path.begin(), path.end(), '/', '.');
    store(deltas, path, toArray(delta));
  };

  addDelta(gman->composeSymNameITS());
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      addDelta(gman->composeSymNameHalfBarrel(ilr, ihb));
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        addDelta(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          addDelta(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            addDelta(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }

  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    addDelta(o2::base::GeometryManager::getSymbolicName(0, ich));
  }

  boost::property_tree::write_json("global.json", deltas);
  return deltas;
}

// Alex provides the corrections as a global correction since they include the
// deltas of the parents
static std::unordered_map<std::string, TGeoHMatrix> alexParamGlobalDelta;
static std::unordered_map<std::string, TGeoHMatrix> alexParamLocalDelta;
static DeltaTree alexLocalDeltas;
void alexParamsRead(const char* file, const char* geom) {
  Info("", "Reading params from %s with geom %s", file, geom);
  auto localDeltas = getLocalDeltas(geom);
  auto gman = o2::its::GeometryTGeo::Instance();  // just for navigation

  // Prints out all alignable objects
  // const TObjArray *arr = gGeoManager->GetListOfPhysicalNodes();
  // const int nev = arr->GetEntries();
  // for (int i{0}; i < nev; ++i) {
  //   ((TGeoPhysicalNode *)arr->At(i))->Print();
  // }

  auto alexFile = TFile::Open(file);
  auto alexParams = alexFile->Get<TTree>("ITSalignParams");
  float fLayer, fStave, fHalfStave, fDx, fDy, fDz, fDRotA, fDRotB, fDRotC;
  alexParams->SetBranchAddress("ITSAlayer", &fLayer);
  alexParams->SetBranchAddress("ITSAstave", &fStave);
  bool hasHS = !alexParams->SetBranchAddress("ITSAhalfstave", &fHalfStave);
  if (hasHS) {
    LOGP(info, "Applying HS alignment");
  }
  alexParams->SetBranchAddress("ITSADX", &fDx);
  alexParams->SetBranchAddress("ITSADY", &fDy);
  alexParams->SetBranchAddress("ITSADZ", &fDz);
  alexParams->SetBranchAddress("ITSARotA", &fDRotA);
  alexParams->SetBranchAddress("ITSARotB", &fDRotB);
  alexParams->SetBranchAddress("ITSARotC", &fDRotC);
  double psi, theta, phi;
  TGeoHMatrix* cur;
  for (int iEntry{0}; alexParams->LoadTree(iEntry) >= 0; ++iEntry) {
    if (alexParams->GetEntry(iEntry) < 0) {
      continue;
    }
    convertZXZEulerToXYZCardano(fDRotA, fDRotB, fDRotC, psi, theta, phi);

    const int nstaves = gman->getNumberOfStaves((int)fLayer) / 2;
    const int lay = (int)fLayer;
    const int sta = ((int)fStave) % nstaves;
    const int hba = (int)fStave >= nstaves;
    int hs = (int)fHalfStave;

    std::string name;
    name = gman->composeSymNameStave(lay, hba, sta);
    if (hasHS && lay > 2) {
      name = gman->composeSymNameHalfStave(lay, hba, sta, hs);
    } else {
      name = gman->composeSymNameStave(lay, hba, sta);
    }

    // global delta
    auto& g = alexParamGlobalDelta[name] = TGeoHMatrix();
    double rot[9] = {};
    double tr[3] = {fDx, fDy, fDz};
    g.SetTranslation(tr);
    anglesToMatrix(psi, theta, phi, rot);
    g.SetRotation(rot);

    // local delta
    auto& l = alexParamLocalDelta[name] = g;
    if (!(cur = o2::base::GeometryManager::getMatrix(name.c_str()))) {
      LOGP(fatal, "cannot get aligned matrix");  // this should not happen
    }
    l.Multiply(cur);
    l.MultiplyLeft(cur->Inverse());
  }

  // This is a global shift of the envelope volume
  // these should be absorbed into the TPC distortion maps
  // const Double_t itsXYZ[3] = {-0.02763, 0.1835, 0.22586}; // actual TPC drift
  // const Double_t itsXYZ[3] = {0.02763, -0.1835, -0.22586};  // reversed!
  // std::string itsTopSym = gman->composeSymNameITS();
  // TGeoHMatrix mat;
  // mat.SetTranslation(itsXYZ);
  // alexParamLocalDelta[itsTopSym] = mat;

  Info("", "Recalculated local deltas");

  // now we have local deltas from the refGeom and new to be applied local
  // deltas they need to be combined del_loc_comb = del_loc_ref * del_loc_new
  // then we build the global deltas and save them as AlignParam vector
  std::unordered_map<std::string, TGeoHMatrix> updatedLocDeltas;
  for (const auto& [sym, dloc_new] : alexParamLocalDelta) {
    std::string dpath{sym};
    std::replace(dpath.begin(), dpath.end(), '/', '.');
    const auto dlocp = load(localDeltas, dpath);
    auto dloc = toMatrix(dlocp);
    dloc.Multiply(dloc_new);
    updatedLocDeltas[sym] = dloc;
  }
  Info("", "Updated local deltas with newpar");

  AlgParams algNewParams;
  auto locDeltas = getLocalDeltas(localDeltas, updatedLocDeltas);
  auto addPar = [&](const char* p, int id = -1, bool print = false) {
    std::string path{p};
    std::replace(path.begin(), path.end(), '/', '.');
    auto t = load(locDeltas, path);
    auto& par = algNewParams.emplace_back(p, id, t[0], t[1], t[2], t[3], t[4],
                                          t[5], false, false);
    if (print) {
      par.print();
    }
  };
  // auto globalDeltas = getGlobalDeltas(localDeltas, updatedLocDeltas);
  // auto addPar = [&](const char *p, int id = -1) {
  //   std::string path{p};
  //   std::replace(path.begin(), path.end(), '/', '.');
  //   auto t = load(globalDeltas, path);
  //   algNewParams.emplace_back(p, id, t[0], t[1], t[2], t[3], t[4], t[5]);
  // };

  addPar(gman->composeSymNameITS(), -1, true);
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      addPar(gman->composeSymNameHalfBarrel(ilr, ihb), -1, true);
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        addPar(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          addPar(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            addPar(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }

  for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
    addPar(o2::base::GeometryManager::getSymbolicName(0, ich),
           o2::base::GeometryManager::getSensID(0, ich));
  }

  LOGP(info, "saving new parameters in snapshot.root");
  auto algNewFile = TFile::Open("snapshot.root", "RECREATE");
  algNewFile->WriteObjectAny(
      &algNewParams, "std::vector<o2::detectors::AlignParam>", "ccdb_object");
  algNewFile->Close();
}

void createDummy() {
  auto gman = o2::its::GeometryTGeo::Instance();  // just for navigation
  AlgParams algNewParams;
  using Pars = std::array<float, 6>;
  auto addPar = [&](const char* c, Pars p, int id = -1, bool print = false) {
    auto& par = algNewParams.emplace_back(c, id, p[0], p[1], p[2], p[3], p[4],
                                          p[5], false, false);
    if (print) {
      par.print();
    }
  };

  const Pars itsXYZ = {-0.02763, 0.1835, 0.22586, 0, 0, 0.1};
  const Pars itsHB = {-0.02763, 0.1835, 0.22586, 0, 0, 0.1};

  addPar(gman->composeSymNameITS(), itsXYZ, -1, true);
  for (int ilr{0}; ilr < gman->getNumberOfLayers(); ++ilr) {
    for (int ihb{0}; ihb < gman->getNumberOfHalfBarrels(); ++ihb) {
      addPar(gman->composeSymNameHalfBarrel(ilr, ihb), itsHB, -1, true);
      const int nst = gman->getNumberOfStaves(ilr) / 2;
      for (int ist{0}; ist < nst; ++ist) {
        // addPar(gman->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < gman->getNumberOfHalfStaves(ilr); ++ihst) {
          // addPar(gman->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < gman->getNumberOfModules(ilr); imd++) {
            // addPar(gman->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }

  // for (int ich = 0; ich < gman->getNumberOfChips(); ich++) {
  //   addPar(o2::base::GeometryManager::getSymbolicName(0, ich),
  //          o2::base::GeometryManager::getSensID(0, ich));
  // }

  LOGP(info, "saving new parameters in snapshot.root");
  auto algNewFile = TFile::Open("snapshot.root", "RECREATE");
  algNewFile->WriteObjectAny(
      &algNewParams, "std::vector<o2::detectors::AlignParam>", "ccdb_object");
  algNewFile->Close();
}

void closureLocalTest() {
  Info("", "running closure test");
  auto algCCDBFile = TFile::Open("ccdb_its_alignment.root");
  if (!algCCDBFile) {
    Error("", "did not open file");
  }
  AlgParams* ccdbPars;
  algCCDBFile->GetObject("ccdb_object", ccdbPars);
  AlgParams algNewParams;
  for (const auto& p : *ccdbPars) {
    TGeoHMatrix loc;
    if (!p.createLocalMatrix(loc)) {
      Error("", "!!!");
      return;
    }
    const Double_t *tra, *rot;
    tra = loc.GetTranslation();
    rot = loc.GetRotationMatrix();
    double psi, theta, phi;
    matrixToAngles(rot, psi, theta, phi);
    algNewParams.emplace_back(p.getSymName().c_str(), p.getAlignableID(),
                              tra[0], tra[1], tra[2], psi, theta, phi, false,
                              false);
  }

  LOGP(info, "saving new parameters in snapshot.root");
  auto algNewFile = TFile::Open("snapshot.root", "RECREATE");
  algNewFile->WriteObjectAny(
      &algNewParams, "std::vector<o2::detectors::AlignParam>", "ccdb_object");
  algNewFile->Close();
}

void applyAlignmentRed(
    const char* alexRefGeom = nullptr,     // reference geometry upon which
                                           // Alex’s changes are based
    const char* alexParamsFile = nullptr,  // Alex’s file containing corrections
    const char* refGeom =
        "o2sim_geometry-aligned_CCDB.root",  // just use to built correct matrix
                                             // paths!
    const char* idealGeom = "o2sim_geometry.root", bool dummy = false,
    bool closure = false) {
  if (alexRefGeom) {
    cacheBuildSym2Path(alexRefGeom);
    cacheBuildLocalIdeal(idealGeom);
    alexParamsRead(alexParamsFile, alexRefGeom);
  } else if (dummy) {
    loadGeom(idealGeom);
    createDummy();
  } else if (closure) {
    loadGeom(refGeom);
    closureLocalTest();
  }
}
