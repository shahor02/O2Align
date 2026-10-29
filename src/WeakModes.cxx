#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <fstream>
#include <regex>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TGeoMatrix.h>

#include "Framework/Logger.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "O2Align/AlignGeomUtils.h"
#include "O2Align/Detector.h"
#include "O2Align/WeakModes.h"
#include "O2Align/Volume.h"

namespace o2::alignrs
{
namespace
{
constexpr double CmToUm = 1e4;
constexpr int NRigid = 6;                   // global shift (3) and rotation (3) absorbed before the modes
constexpr int NCols = Volume::NGlobalEq;    // one column per mode *equation*, see Volume::NGlobalEq
constexpr int NPar = NRigid + NCols;

// A column of the mode part of the design matrix is one global-mode equation, identified by its flat
// index; the mode and the phase are decoded from it, so adding a mode to Volume::GlobalMode
// propagates here with no change.
std::string colName(int c) { return Volume::globalEqName(c); }

// Displacement of the centre of a sensor between the geometry the fit started from and the corrected one
struct SensorDisp {
  const Volume* vol{};
  std::array<double, 3> x{};  // position, cm
  std::array<double, 3> d{};  // displacement, cm
  std::array<double, 3> w{};  // rotation vector of the correction, rad
  int layer{-1};
};

struct Branch {
  Volume* top{};
  std::vector<SensorDisp> sensors;
};

int layerOf(const std::string& symName)
{
  static const std::regex re("Layer([0-9]+)");
  std::smatch m;
  return std::regex_search(symName, m, re) ? std::stoi(m[1]) : -1;
}

SensorDisp makeSensorDisp(const Volume* vol, const TGeoHMatrix& cumulative)
{
  SensorDisp s;
  s.vol = vol;
  s.layer = layerOf(vol->getSymName());
  const double* pos = vol->getMatrixL2G().GetTranslation();
  double moved[3];
  cumulative.LocalToMaster(pos, moved);
  const double* R = cumulative.GetRotationMatrix();
  for (int i = 0; i < 3; ++i) {
    s.x[i] = pos[i];
    s.d[i] = moved[i] - pos[i];
  }
  s.w = {0.5 * (R[7] - R[5]), 0.5 * (R[2] - R[6]), 0.5 * (R[3] - R[1])};
  return s;
}

// Compose the fitted corrections down the tree like Detector::MP2AlignParams does (cumulative global
// delta, parents first) and collect the displacement of every sensor, grouped by detector branch.
// withHalfBarrels: also the volumes whose children are all pseudo (ITS3 half barrels) are sensors, which
// the weak-modes diagnostics do not use.
bool isSensorVolume(const Volume* vol, bool withHalfBarrels)
{
  if (vol->isRoot() || vol->isVirtual() || vol->isPseudo()) {
    return false;
  }
  if (vol->isLeaf()) {
    return true;
  }
  bool allPseudo = true;
  vol->forEachChild([&](const Volume* c) { allPseudo = allPseudo && c->isPseudo(); });
  return withHalfBarrels && allPseudo;
}

std::vector<Branch> collectBranches(Volume* root, const std::map<uint32_t, double>& labelToValue, bool withHalfBarrels = false)
{
  std::vector<Branch> branches;
  std::unordered_map<const Volume*, TGeoHMatrix> cumulative;
  std::unordered_map<const Volume*, int> branchOf;
  root->traverse([&](Volume* vol) {
    TGeoHMatrix pn;
    if (!vol->isRoot()) {
      pn = cumulative.at(vol->getParent());
      branchOf[vol] = vol->getParent()->isRoot() ? static_cast<int>(branches.size()) : branchOf.at(vol->getParent());
      if (vol->getParent()->isRoot()) {
        branches.push_back({vol, {}});
      }
    }
    TGeoHMatrix delta;
    if (vol->getFittedLocalDelta(labelToValue, delta)) {
      const auto& l2g = vol->getMatrixL2G();
      pn = pn * l2g * delta * l2g.Inverse();
    }
    cumulative[vol] = pn;
    if (isSensorVolume(vol, withHalfBarrels)) {
      branches[branchOf.at(vol)].sensors.push_back(makeSensorDisp(vol, pn));
    }
  });
  return branches;
}

// Displacement caused by the initial (CCDB) alignment alone, relative to the ideal geometry. Pp
// maps the ideal sensor-centre position to the loaded (pre-aligned) one, in the global frame (same
// convention as Detector's P^p, see Detector.cxx).
SensorDisp makeInitialSensorDisp(const Volume* vol, const TGeoHMatrix& Pp)
{
  SensorDisp s;
  s.vol = vol;
  s.layer = layerOf(vol->getSymName());
  const double* posLoaded = vol->getMatrixL2G().GetTranslation();
  double posIdeal[3];
  Pp.Inverse().LocalToMaster(posLoaded, posIdeal);
  const double* R = Pp.GetRotationMatrix();
  for (int i = 0; i < 3; ++i) {
    s.x[i] = posIdeal[i];
    s.d[i] = posLoaded[i] - posIdeal[i];
  }
  s.w = {0.5 * (R[7] - R[5]), 0.5 * (R[2] - R[6]), 0.5 * (R[3] - R[1])};
  return s;
}

// Same branch grouping as collectBranches, but Pp comes from each detector's CCDB initial
// alignment (geomutils::collectInitialCumulativeDeltas) instead of from a millepede fit.
std::vector<Branch> collectInitialBranches(Volume* root, const std::vector<Detector*>& detectors)
{
  std::vector<Branch> branches;
  std::vector<geomutils::PathMatrices> iniCumPerBranch;
  std::unordered_map<const Volume*, int> branchOf;
  root->traverse([&](Volume* vol) {
    if (vol->isRoot()) {
      return;
    }
    if (vol->getParent()->isRoot()) {
      branchOf[vol] = static_cast<int>(branches.size());
      branches.push_back({vol, {}});
      const auto it = std::find_if(detectors.begin(), detectors.end(), [vol](const Detector* d) { return d->getTopVolume() == vol; });
      const auto* iniAlg = (it != detectors.end() && (*it)->getO2DetID() >= 0) ? o2::base::GRPGeomHelper::instance().getAlignment((*it)->getO2DetID()) : nullptr;
      iniCumPerBranch.push_back(geomutils::collectInitialCumulativeDeltas(iniAlg));
    } else {
      branchOf[vol] = branchOf.at(vol->getParent());
    }
    if (vol->isLeaf() && !vol->isRoot() && !vol->isVirtual() && !vol->isPseudo()) {
      const int ib = branchOf.at(vol);
      const auto Pp = geomutils::getOnPath(iniCumPerBranch[ib], geomutils::getAlignablePath(vol->getSymName()), true);
      branches[ib].sensors.push_back(makeInitialSensorDisp(vol, Pp));
    }
  });
  return branches;
}

// Metric A of one mode column over the subtree of top: sum(c*p)/sum(|c|) in um, and the number of terms
struct MetricA {
  double value{0.};
  int nTerms{0};
};

MetricA metricA(Volume* top, const std::map<uint32_t, double>& labelToValue, int eq)
{
  const int mode = Volume::globalEqMode(eq), phase = Volume::globalEqPhase(eq);
  double sumCP = 0., sumC = 0.;
  int nTerms = 0;
  top->traverse([&](Volume* vol) {
    const auto* rb = vol->getRigidBody();
    if (!rb || rb->nFreeDOFs() == 0 || vol->isPseudo()) {
      return;
    }
    const auto* l2g = &vol->getMatrixL2G();
    double u[3], om[3];
    if (!Volume::globalModeField(mode, phase, l2g->GetTranslation(), u, om)) {
      return;
    }
    const double* rot = l2g->GetRotationMatrix();
    const int nDOF = std::min(rb->nDOFs(), static_cast<int>(RigidBodyDOFSet::NDOF));
    for (int i = 0; i < nDOF; ++i) {
      if (!rb->isFree(i)) {
        continue;
      }
      const double* v = (i < 3) ? u : om;
      const int ax = i % 3;
      const double c = v[0] * rot[ax] + v[1] * rot[3 + ax] + v[2] * rot[6 + ax];
      sumCP += c * Volume::getFittedValue(labelToValue, vol->getLabel().raw(i));
      sumC += std::abs(c);
      ++nTerms;
    }
  });
  return {sumC > 0. ? CmToUm * sumCP / sumC : 0., nTerms};
}

// Result of the fit of the displacements of one branch
struct ModeFit {
  std::array<double, NCols> amp{};   // displacement of the unit-scale mode field, cm: the amplitude at the largest lever arm
  std::array<double, NCols> err{};
  std::array<double, NCols> scale{}; // lever arm of the mode, in units of cm^power, to convert amp to the native amplitude
  std::array<double, NCols> explained{};
  std::array<double, NRigid> rigid{};
  Eigen::Matrix<double, NCols, NCols> corr = Eigen::Matrix<double, NCols, NCols>::Zero();
  double rmsRaw{0.}, rmsRigid{0.}, rmsModes{0.}; // cm
  double condition{0.};                          // smallest/largest singular value of the (column-normalised) design matrix
  int rank{0};
  // per-sensor model: rigid part and every mode contribution, cm
  std::vector<std::array<double, 3>> dRigid;
  std::vector<std::array<std::array<double, 3>, NCols>> dMode;
};

std::array<double, NCols> modeScales(const std::vector<SensorDisp>& sensors)
{
  double rMax = 0., zMax = 0.;
  for (const auto& s : sensors) {
    rMax = std::max(rMax, std::hypot(s.x[0], s.x[1]));
    zMax = std::max(zMax, std::abs(s.x[2]));
  }
  rMax = rMax > 0. ? rMax : 1.;
  zMax = zMax > 0. ? zMax : 1.;
  std::array<double, NCols> sc{};
  for (int c = 0; c < NCols; ++c) {
    switch (Volume::globalEqMode(c)) {
      case Volume::RadialScaling:
      case Volume::Curl:
      case Volume::RadialHarmonic:     // the harmonic families carry an explicit factor r, so their
      case Volume::TangentialHarmonic: // native amplitude is the dimensionless eps as well
        sc[c] = rMax;
        break;
      case Volume::ZScaling:
        sc[c] = zMax;
        break;
      case Volume::Twist:
        sc[c] = rMax * zMax;
        break;
      case Volume::Bowing:
        sc[c] = zMax * zMax;
        break;
      case Volume::GlobalTilt: // a rigid rotation: the displacement is the angle times the lever arm
        sc[c] = std::max(rMax, zMax);
        break;
      case Volume::RadialPolar:
      case Volume::TangentialPolar:
      case Volume::LongitudinalPolar: // r^p*(z/r)^k: no closed form, take the largest |u| over the branch
        sc[c] = 0.;
        for (const auto& s : sensors) {
          double u[3], om[3];
          if (Volume::globalModeField(Volume::globalEqMode(c), Volume::globalEqPhase(c), s.x.data(), u, om)) {
            sc[c] = std::max(sc[c], std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
          }
        }
        sc[c] = sc[c] > 0. ? sc[c] : 1.;
        break;
      default: // Elliptical: the field is already a length
        sc[c] = 1.;
    }
  }
  return sc;
}

// d_i = t + r x x_i + sum_m a_m u_m(x_i)/scale_m, solved by SVD with the columns normalised; the modes
// which have no handle on the branch (e.g. z scaling when all sensors sit at z=0) come out as 0.
ModeFit fitModes(const std::vector<SensorDisp>& sensors)
{
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  ModeFit fit;
  fit.scale = modeScales(sensors);
  const int nS = static_cast<int>(sensors.size());
  MatrixXd A = MatrixXd::Zero(3 * nS, NPar);
  VectorXd b(3 * nS);
  for (int i = 0; i < nS; ++i) {
    const auto& s = sensors[i];
    for (int k = 0; k < 3; ++k) {
      b(3 * i + k) = s.d[k];
      A(3 * i + k, k) = 1.; // shift
      // rotation about axis k: e_k x x_i
      const double* x = s.x.data();
      const double ek[3] = {double(k == 0), double(k == 1), double(k == 2)};
      for (int j = 0; j < 3; ++j) {
        A(3 * i + j, 3 + k) = ek[(j + 1) % 3] * x[(j + 2) % 3] - ek[(j + 2) % 3] * x[(j + 1) % 3];
      }
    }
    for (int c = 0; c < NCols; ++c) {
      double u[3], om[3];
      if (Volume::globalModeField(Volume::globalEqMode(c), Volume::globalEqPhase(c), s.x.data(), u, om)) {
        for (int k = 0; k < 3; ++k) {
          A(3 * i + k, NRigid + c) = u[k] / fit.scale[c];
        }
      }
    }
  }
  VectorXd colNorm(NPar);
  for (int j = 0; j < NPar; ++j) {
    colNorm(j) = A.col(j).norm() > 0. ? A.col(j).norm() : 1.;
    A.col(j) /= colNorm(j);
  }
  // JacobiSVD, not BDCSVD: the design matrix is exactly rank deficient by construction (curl and
  // globalTilt duplicate the rigid rotation columns, globalShift the rigid shift ones), and BDCSVD
  // (Eigen 3.4) returned NaN singular values on it, which left rank 0 and the whole fit void. ~50
  // columns, so the cost is irrelevant.
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const auto& sv = svd.singularValues();
  const double tol = 1e-8 * sv(0);
  Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(NPar, NPar); // of the normalised parameters, unit sigma
  VectorXd y = VectorXd::Zero(NPar);
  for (int j = 0; j < NPar; ++j) {
    if (sv(j) > tol) {
      ++fit.rank;
      y += svd.matrixV().col(j) * (svd.matrixU().col(j).dot(b) / sv(j));
      cov += svd.matrixV().col(j) * svd.matrixV().col(j).transpose() / (sv(j) * sv(j));
    }
  }
  fit.condition = sv(NPar - 1) / sv(0);
  VectorXd theta = y.cwiseQuotient(colNorm);
  const VectorXd model = (A * y);
  const double rss = (b - model).squaredNorm();
  const double sigma2 = rss / std::max(1, 3 * nS - fit.rank);
  for (int j = 0; j < NRigid; ++j) {
    fit.rigid[j] = theta(j);
  }
  for (int c = 0; c < NCols; ++c) {
    const int j = NRigid + c;
    fit.amp[c] = theta(j);
    fit.err[c] = std::sqrt(sigma2 * cov(j, j)) / colNorm(j);
  }
  for (int c = 0; c < NCols; ++c) {
    for (int d = 0; d < NCols; ++d) {
      const double vc = cov(NRigid + c, NRigid + c), vd = cov(NRigid + d, NRigid + d);
      fit.corr(c, d) = (vc > 0. && vd > 0.) ? cov(NRigid + c, NRigid + d) / std::sqrt(vc * vd) : 0.;
    }
  }
  // per-sensor contributions and the power explained by every mode after the rigid part is removed
  fit.dRigid.resize(nS);
  fit.dMode.resize(nS);
  double pRaw = 0., pRigid = 0., pModes = 0.;
  std::array<double, NCols> pMode{};
  for (int i = 0; i < nS; ++i) {
    for (int k = 0; k < 3; ++k) {
      double rig = 0.;
      for (int j = 0; j < NRigid; ++j) {
        rig += A(3 * i + k, j) * y(j);
      }
      fit.dRigid[i][k] = rig;
      double mod = 0.;
      for (int c = 0; c < NCols; ++c) {
        const double v = A(3 * i + k, NRigid + c) * y(NRigid + c);
        fit.dMode[i][c][k] = v;
        pMode[c] += v * v;
        mod += v;
      }
      pRaw += b(3 * i + k) * b(3 * i + k);
      pRigid += std::pow(b(3 * i + k) - rig, 2);
      pModes += std::pow(b(3 * i + k) - rig - mod, 2);
    }
  }
  const double n = std::max(1, nS);
  fit.rmsRaw = std::sqrt(pRaw / n);
  fit.rmsRigid = std::sqrt(pRigid / n);
  fit.rmsModes = std::sqrt(pModes / n);
  for (int c = 0; c < NCols; ++c) {
    fit.explained[c] = pRigid > 0. ? pMode[c] / pRigid : 0.;
  }
  return fit;
}

nlohmann::json branchToJSON(const Branch& br, const ModeFit& fit, const std::array<MetricA, NCols>* mA)
{
  nlohmann::json j;
  j["detector"] = br.top->getSymName();
  j["nSensors"] = br.sensors.size();
  j["rmsDisplacement_um"] = CmToUm * fit.rmsRaw;
  j["rmsAfterRigid_um"] = CmToUm * fit.rmsRigid;
  j["rmsAfterModes_um"] = CmToUm * fit.rmsModes;
  j["rank"] = fit.rank;
  j["conditioning"] = fit.condition;
  j["rigid"] = {{"TX_um", CmToUm * fit.rigid[0]}, {"TY_um", CmToUm * fit.rigid[1]}, {"TZ_um", CmToUm * fit.rigid[2]},
                {"RX_urad", 1e6 * fit.rigid[3]}, {"RY_urad", 1e6 * fit.rigid[4]}, {"RZ_urad", 1e6 * fit.rigid[5]}};
  auto modes = nlohmann::json::array();
  auto corr = nlohmann::json::array();
  for (int c = 0; c < NCols; ++c) {
    nlohmann::json m{{"name", colName(c)},
                     {"metricB_um", CmToUm * fit.amp[c]},
                     {"metricBError_um", CmToUm * fit.err[c]},
                     {"explainedFraction", fit.explained[c]},
                     {"leverArm", fit.scale[c]},
                     {"nativeAmplitude", fit.amp[c] / fit.scale[c]}};
    if (mA) {
      m["metricA_um"] = (*mA)[c].value;
      m["nTermsA"] = (*mA)[c].nTerms;
    }
    modes.push_back(std::move(m));
    std::vector<double> row(NCols);
    for (int d = 0; d < NCols; ++d) {
      row[d] = fit.corr(c, d);
    }
    corr.push_back(std::move(row));
  }
  j["modes"] = std::move(modes);
  j["correlation"] = std::move(corr);
  return j;
}

void logBranch(const Branch& br, const ModeFit& fit, const std::array<MetricA, NCols>* mA)
{
  LOGP(info, "Weak modes of {}: {} sensors, displacement rms {:.1f} um -> {:.1f} um after rigid shift/rotation -> {:.1f} um after modes (rank {}/{}, cond {:.1e})",
       br.top->getSymName(), br.sensors.size(), CmToUm * fit.rmsRaw, CmToUm * fit.rmsRigid, CmToUm * fit.rmsModes, fit.rank, NPar, fit.condition);
  if (mA) {
    LOGP(info, "  {:<20} {:>12} {:>10} {:>10} {:>14}", "mode", "B [um@lever]", "+-", "explained", "A [um] (nDOF)");
    for (int c = 0; c < NCols; ++c) {
      LOGP(info, "  {:<20} {:>12.2f} {:>10.2f} {:>9.1f}% {:>9.2f} ({})", colName(c), CmToUm * fit.amp[c], CmToUm * fit.err[c], 100. * fit.explained[c], (*mA)[c].value, (*mA)[c].nTerms);
    }
  } else {
    LOGP(info, "  {:<20} {:>12} {:>10} {:>10}", "mode", "B [um@lever]", "+-", "explained");
    for (int c = 0; c < NCols; ++c) {
      LOGP(info, "  {:<20} {:>12.2f} {:>10.2f} {:>9.1f}%", colName(c), CmToUm * fit.amp[c], CmToUm * fit.err[c], 100. * fit.explained[c]);
    }
  }
}

// ROOT output: one row per sensor (all branches), and one per branch and mode
void writeRoot(const std::string& fileName, const std::vector<Branch>& branches, const std::vector<ModeFit>& fits, const std::vector<std::array<MetricA, NCols>>* mA)
{
  TFile fout(fileName.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    LOGP(fatal, "Cannot open the weak modes output {}", fileName);
  }
  int iBranch = 0, layer = 0, sensorId = 0;
  double x[3], d[3], w[3], r, phi, dr, rdphi, dRigid[3], uMode[NCols * 3];
  char symName[256];
  TTree tSens("sensors", "per-sensor displacement of the fitted corrections [cm, rad]");
  tSens.Branch("branch", &iBranch);
  tSens.Branch("layer", &layer);
  tSens.Branch("sensorId", &sensorId);
  tSens.Branch("symName", symName, "symName/C");
  tSens.Branch("x", &x[0], "x/D");
  tSens.Branch("y", &x[1], "y/D");
  tSens.Branch("z", &x[2], "z/D");
  tSens.Branch("r", &r);
  tSens.Branch("phi", &phi);
  tSens.Branch("dx", &d[0], "dx/D");
  tSens.Branch("dy", &d[1], "dy/D");
  tSens.Branch("dz", &d[2], "dz/D");
  tSens.Branch("dr", &dr);
  tSens.Branch("rdphi", &rdphi);
  tSens.Branch("wx", &w[0], "wx/D");
  tSens.Branch("wy", &w[1], "wy/D");
  tSens.Branch("wz", &w[2], "wz/D");
  tSens.Branch("dRigid", dRigid, "dRigid[3]/D");
  tSens.Branch("uMode", uMode, std::format("uMode[{}]/D", NCols * 3).c_str());
  for (int ib = 0; ib < static_cast<int>(branches.size()); ++ib) {
    iBranch = ib;
    for (int is = 0; is < static_cast<int>(branches[ib].sensors.size()); ++is) {
      const auto& s = branches[ib].sensors[is];
      layer = s.layer;
      sensorId = s.vol->getSensorId();
      std::snprintf(symName, sizeof(symName), "%s", s.vol->getSymName().c_str());
      std::copy(s.x.begin(), s.x.end(), x);
      std::copy(s.d.begin(), s.d.end(), d);
      std::copy(s.w.begin(), s.w.end(), w);
      r = std::hypot(x[0], x[1]);
      phi = std::atan2(x[1], x[0]);
      const double c = x[0] / (r > 0. ? r : 1.), sn = x[1] / (r > 0. ? r : 1.);
      dr = d[0] * c + d[1] * sn;
      rdphi = -d[0] * sn + d[1] * c;
      std::copy(fits[ib].dRigid[is].begin(), fits[ib].dRigid[is].end(), dRigid);
      for (int m = 0; m < NCols; ++m) {
        for (int k = 0; k < 3; ++k) {
          uMode[3 * m + k] = fits[ib].dMode[is][m][k];
        }
      }
      tSens.Fill();
    }
  }
  tSens.Write();

  std::string detName, modeName;
  double ampB, errB, explained, ampA, lever;
  int nTermsA;
  TTree tMode("modes", "fit results per branch and mode, um");
  tMode.Branch("branch", &iBranch);
  tMode.Branch("detector", &detName);
  tMode.Branch("mode", &modeName);
  tMode.Branch("metricB", &ampB);
  tMode.Branch("metricBErr", &errB);
  tMode.Branch("explained", &explained);
  tMode.Branch("metricA", &ampA);
  tMode.Branch("nTermsA", &nTermsA);
  tMode.Branch("leverArm", &lever);
  for (int ib = 0; ib < static_cast<int>(branches.size()); ++ib) {
    iBranch = ib;
    detName = branches[ib].top->getSymName();
    for (int c = 0; c < NCols; ++c) {
      modeName = colName(c);
      ampB = CmToUm * fits[ib].amp[c];
      errB = CmToUm * fits[ib].err[c];
      explained = fits[ib].explained[c];
      ampA = mA ? (*mA)[ib][c].value : 0.;
      nTermsA = mA ? (*mA)[ib][c].nTerms : 0;
      lever = fits[ib].scale[c];
      tMode.Fill();
    }
  }
  tMode.Write();
}
} // namespace

std::vector<SensorDisplacement> computeSensorDisplacements(Volume* root, const std::map<uint32_t, double>& labelToValue)
{
  std::vector<SensorDisplacement> out;
  for (const auto& br : collectBranches(root, labelToValue, true)) {
    for (const auto& s : br.sensors) {
      out.push_back({br.top, s.vol, s.layer, s.x, s.d, s.w});
    }
  }
  return out;
}

void diagnoseWeakModes(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outBase)
{
  auto branches = collectBranches(root, labelToValue);
  // keep the branches with enough sensors for a fit of all the parameters
  branches.erase(std::remove_if(branches.begin(), branches.end(), [](const Branch& b) { return b.sensors.size() < NPar; }), branches.end());
  std::vector<ModeFit> fits;
  std::vector<std::array<MetricA, NCols>> mA;
  auto jsonOut = nlohmann::json::array();
  for (const auto& br : branches) {
    fits.push_back(fitModes(br.sensors));
    auto& a = mA.emplace_back();
    for (int c = 0; c < NCols; ++c) {
      a[c] = metricA(br.top, labelToValue, c);
    }
    logBranch(br, fits.back(), &a);
    jsonOut.push_back(branchToJSON(br, fits.back(), &a));
  }
  if (branches.empty()) {
    LOGP(warn, "No detector branch with enough sensors for the weak modes diagnostics");
    return;
  }
  std::ofstream fout(outBase + ".json");
  fout << jsonOut.dump(2) << '\n';
  writeRoot(outBase + ".root", branches, fits, &mA);
  LOGP(info, "Wrote the weak modes diagnostics to {}.json and {}.root", outBase, outBase);
}

void diagnoseInitialAlignment(Volume* root, const std::vector<Detector*>& detectors, const std::string& outBase)
{
  auto branches = collectInitialBranches(root, detectors);
  branches.erase(std::remove_if(branches.begin(), branches.end(), [](const Branch& b) { return b.sensors.size() < NPar; }), branches.end());
  std::vector<ModeFit> fits;
  auto jsonOut = nlohmann::json::array();
  for (const auto& br : branches) {
    fits.push_back(fitModes(br.sensors));
    logBranch(br, fits.back(), nullptr);
    jsonOut.push_back(branchToJSON(br, fits.back(), nullptr));
  }
  if (branches.empty()) {
    LOGP(warn, "No detector branch with enough sensors for the initial alignment weak modes diagnostics");
    return;
  }
  std::ofstream fout(outBase + ".json");
  fout << jsonOut.dump(2) << '\n';
  writeRoot(outBase + ".root", branches, fits, nullptr);
  LOGP(info, "Wrote the initial alignment weak modes diagnostics to {}.json and {}.root", outBase, outBase);
}

} // namespace o2::alignrs
