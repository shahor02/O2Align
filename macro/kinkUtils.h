// Common helpers of the IB/OB kink macros (kinkMap.C, kinkResponse.C,
// makeKinkCure.C), see doc/handovers/kink_cure.md.
//
// The input is the `res` tree of o2-check-resid-workflow (checkResid): per
// track, the IB sub-track (vertex + Lr0..Lr2, trIBOut) and the OB sub-track
// (Lr6..Lr3, trOBInw), both at R = 12 cm and in the same frame, plus the
// clusters (points.{dy,dz,phi,z,sens,lr}). The kink observables are the IB - OB
// differences dPar0 (y, um) and dPar2 (snp, urad). Sensor frames are read
// either from an aligned geometry file (o2sim_geometry-aligned.root) or from an
// algDump output (trees gm/gmx/gmy/gmz, branch "pos"), both giving the ITS chip
// alignable matrix.
#ifndef O2ALIGN_KINKUTILS_H
#define O2ALIGN_KINKUTILS_H

#include "DetectorsBase/GeometryManager.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include <TComplex.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace kink {
constexpr int NP = 24; // phi bins (15 deg)
constexpr int NT = 8;  // tgl bins
constexpr double TE[NT + 1] = {-2,  -1.2, -0.8, -0.4, 0,
                               0.4, 0.8,  1.2,  2}; // tgl bin edges
constexpr double RComp =
    12.; // radius of the IB/OB comparison, cm (checkresid.rCompIBOB)
constexpr double RLayer[7] = {2.24, 3.01, 3.78, 19.6,
                              24.5, 34.4, 39.3}; // nominal layer radii, cm

// global azimuth of the IB sub-track position at R = 12 cm, in [0, 2pi)
const char *PhiIB =
    "TMath::ATan2(TMath::Sin(trIBOut.mAlpha+TMath::ATan2(trIBOut.mP[0],trIBOut."
    "mX)),TMath::Cos(trIBOut.mAlpha+TMath::ATan2(trIBOut.mP[0],trIBOut.mX)))+"
    "TMath::Pi()*(1-TMath::Sign(1.,TMath::Sin(trIBOut.mAlpha+TMath::ATan2("
    "trIBOut.mP[0],trIBOut.mX))))";
const char *SelIOB =
    "trIBOut.mX>1&&trOBInw.mX>1"; // both sub-tracks defined (7-cluster tracks)

int tglBin(double t) {
  for (int i = 0; i < NT; i++) {
    if (t < TE[i + 1]) {
      return i;
    }
  }
  return NT - 1;
}

double tglCentre(int i) { return 0.5 * (TE[i] + TE[i + 1]); }

/// Gaussian-fitted mean (and its error) of `expr` per (tgl bin, phi bin) cell,
/// multiplied by `scale`. `range` is the half-range of the histogrammed
/// quantity in its own units; `extraSel` is ANDed with SelIOB.
void cellMeans(TTree *t, const char *expr, double range, double scale,
               double m[NT][NP], double e[NT][NP], const char *extraSel = "",
               long nmax = TTree::kMaxEntries) {
  TH3F h("kinkCells", "", NP, 0, 2 * M_PI, NT, -2, 2, 100, -range, range);
  h.GetYaxis()->Set(NT, TE);
  std::string sel = SelIOB;
  if (extraSel && *extraSel) {
    sel += std::string("&&(") + extraSel + ")";
  }
  t->Draw(Form("%s:track.mP[3]:%s>>kinkCells", expr, PhiIB), sel.c_str(),
          "goff", nmax);
  for (int it = 0; it < NT; it++) {
    TH2D h2("kinkCells2", "", NP, 0, 2 * M_PI, 100, -range, range);
    for (int i = 1; i <= NP; i++) {
      for (int k = 1; k <= 100; k++) {
        h2.SetBinContent(i, k, h.GetBinContent(i, it + 1, k));
      }
    }
    TObjArray arr;
    TF1 g("g", "gaus");
    h2.FitSlicesY(&g, 0, -1, 0, "QN", &arr);
    auto mm = (TH1 *)arr.At(1);
    for (int i = 0; i < NP; i++) {
      m[it][i] = mm->GetBinContent(i + 1) * scale;
      e[it][i] = mm->GetBinError(i + 1) * scale;
    }
    arr.Delete();
  }
}

/// dPar0 [um] / dPar2 [urad] cell means
void dParMeans(TTree *t, int ipar, double m[NT][NP], double e[NT][NP],
               const char *extraSel = "", long nmax = TTree::kMaxEntries) {
  cellMeans(t, Form("trIBOut.mP[%d]-trOBInw.mP[%d]", ipar, ipar),
            ipar == 0 ? 0.02 : 4e-3, ipar == 0 ? 1e4 : 1e6, m, e, extraSel,
            nmax);
}

/// harmonic n of a phi profile sampled at the NP bin centres: n=0 -> mean, n>0
/// -> (cos, sin) amplitudes
TComplex harmonic(const double *v, int n) {
  TComplex s(0, 0);
  for (int i = 0; i < NP; i++) {
    double ph = (i + 0.5) * 2 * M_PI / NP;
    s += TComplex(v[i] * cos(n * ph), v[i] * sin(n * ph));
  }
  return s * (n ? 2. / NP : 1. / NP);
}

/// chip frame: centre and local axes (unit vectors) in the global frame
struct Frame {
  TVector3 c, ax, ay, az;
};

/// ITS chip frames by chip (sensor) ID, from an aligned geometry or from an
/// algDump file. Fills `symNames` if given (geometry input only).
std::map<int, Frame>
loadFrames(const char *fname, std::map<int, std::string> *symNames = nullptr) {
  std::map<int, Frame> frames;
  TFile f(fname);
  if (auto gm = (TTree *)f.Get("gm")) { // algDump output
    const char *tn[4] = {"gm", "gmx", "gmy", "gmz"};
    TTree *t[4];
    TVector3 *p[4] = {};
    int det, sid;
    for (int k = 0; k < 4; k++) {
      t[k] = (TTree *)f.Get(tn[k]);
      t[k]->SetBranchAddress("pos", &p[k]);
    }
    gm->SetBranchAddress("det", &det);
    gm->SetBranchAddress("sid", &sid);
    for (Long64_t i = 0; i < gm->GetEntries(); i++) {
      for (int k = 0; k < 4; k++) {
        t[k]->GetEntry(i);
      }
      if (det == o2::detectors::DetID::ITS) {
        frames[sid] = {*p[0], *p[1] - *p[0], *p[2] - *p[0], *p[3] - *p[0]};
      }
    }
    if (symNames) {
      printf("loadFrames: %s is an algDump file, no symbolic names\n", fname);
    }
    return frames;
  }
  f.Close();
  delete gGeoManager;
  TGeoManager::Import(fname);
  o2::itsmft::ChipMappingITS mp;
  o2::detectors::DetID its("ITS");
  for (int ic = 0; ic < mp.getNChips(); ic++) {
    auto m = o2::base::GeometryManager::getMatrix(its, ic);
    const double *r = m->GetRotationMatrix();
    const double *tr = m->GetTranslation();
    frames[ic] = {TVector3(tr[0], tr[1], tr[2]), TVector3(r[0], r[3], r[6]),
                  TVector3(r[1], r[4], r[7]), TVector3(r[2], r[5], r[8])};
    if (symNames) {
      (*symNames)[ic] =
          o2::base::GeometryManager::getPNEntry(its, ic)->GetName();
    }
  }
  return frames;
}

/// apparent shift [um] along +phi of the hit when chip B (reference geometry)
/// is replaced by chip A (other geometry), for the straight ray from the origin
/// through the hit point P on B: local x of the ray crossing in A vs in B.
double apparentShift(const Frame &A, const Frame &B, const TVector3 &P) {
  TVector3 u = P.Unit();
  double t = A.c.Dot(A.ay) / u.Dot(A.ay);
  double xB = (P - B.c).Dot(B.ax) / B.ax.Mag2(),
         xA = (t * u - A.c).Dot(A.ax) / A.ax.Mag2();
  TVector3 ph(-P.Y(), P.X(), 0);
  return (xB - xA) * B.ax.Dot(ph.Unit()) * 1e4;
}

/// Per (tgl, phi) cell and layer, the mean apparent r*dphi shift [um] of the
/// hits of the tracks of `t` (the hit pattern) for each geometry of `others`
/// with respect to `ref`. Cells are those of the IB sub-track (phi at 12 cm,
/// track tgl). D[g][it][ip][lr]; layers without hits in a cell get 0.
using CellShifts =
    std::vector<std::vector<std::array<double, 7>>>; // [it][ip][lr]
std::vector<CellShifts>
cellHitShifts(TTree *t, const std::map<int, Frame> &ref,
              const std::vector<const std::map<int, Frame> *> &others,
              long nmax) {
  long nt = std::min<long>(nmax, t->GetEntries());
  t->SetEstimate(nt * 9 + 1);
  long np = t->Draw("Entry$:points.lr:points.sens:points.phi:points.z",
                    "points.lr>=0", "goff", nt);
  std::vector<double> E(t->GetVal(0), t->GetVal(0) + np),
      L(t->GetVal(1), t->GetVal(1) + np), S(t->GetVal(2), t->GetVal(2) + np),
      P(t->GetVal(3), t->GetVal(3) + np), Z(t->GetVal(4), t->GetVal(4) + np);
  long ntr = t->Draw(Form("Entry$:%s:track.mP[3]", PhiIB), SelIOB, "goff", nt);
  std::map<long, std::pair<int, int>> cell; // track entry -> (tgl bin, phi bin)
  for (long i = 0; i < ntr; i++) {
    int ip = int(t->GetVal(1)[i] / (2 * M_PI) * NP);
    if (ip >= 0 && ip < NP) {
      cell[long(t->GetVal(0)[i])] = {tglBin(t->GetVal(2)[i]), ip};
    }
  }
  std::vector<CellShifts> D(others.size(),
                            CellShifts(NT, std::vector<std::array<double, 7>>(
                                               NP, std::array<double, 7>{})));
  std::vector<std::vector<std::array<int, 7>>> N(
      NT, std::vector<std::array<int, 7>>(NP, std::array<int, 7>{}));
  for (long i = 0; i < np; i++) {
    auto it = cell.find(long(E[i]));
    int lr = int(L[i]), sid = int(S[i]);
    auto fb = ref.find(sid);
    if (it == cell.end() || fb == ref.end()) {
      continue;
    }
    const Frame &B = fb->second; // hit point on chip B from the measured phi, z
                                 // (radial position from the chip plane)
    TVector3 dir(cos(P[i]), sin(P[i]), 0);
    double rho = (B.c.Dot(B.ay) - Z[i] * B.ay.Z()) / dir.Dot(B.ay);
    TVector3 hit(rho * cos(P[i]), rho * sin(P[i]), Z[i]);
    auto [a, b] = it->second;
    for (size_t g = 0; g < others.size(); g++) {
      D[g][a][b][lr] += apparentShift(others[g]->at(sid), B, hit);
    }
    N[a][b][lr]++;
  }
  for (auto &Dg : D) {
    for (int a = 0; a < NT; a++) {
      for (int b = 0; b < NP; b++) {
        for (int l = 0; l < 7; l++) {
          if (N[a][b][l]) {
            Dg[a][b][l] /= N[a][b][l];
          }
        }
      }
    }
  }
  return D;
}

/// split "a,b,c" into tokens
std::vector<std::string> tokenize(const std::string &s, char sep = ',') {
  std::vector<std::string> v;
  size_t p = 0;
  while (p <= s.size()) {
    size_t q = s.find(sep, p);
    if (q == std::string::npos) {
      q = s.size();
    }
    if (q > p) {
      v.push_back(s.substr(p, q - p));
    }
    p = q + 1;
  }
  return v;
}
} // namespace kink

#endif
