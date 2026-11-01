// Linear response of the IB/OB kink dPar2 to the apparent r*dphi shifts of the
// hits, see doc/handovers/kink_cure.md.
//
// Model: for a geometry G replacing a reference geometry R, the change of the
// cell mean of dPar2 (per tgl bin x 24 phi bins, as in kinkMap.C) is  delta
// dPar2 = sum_L c_L d_L,  d_L being the mean apparent shift [um] along +phi of
// the layer-L hits of that cell (kink::cellHitShifts: straight ray from the
// origin through the hit, local x crossing in G vs in R). The c_L [urad/um]
// (and the analogous coefficients of dPar0 [um/um]) are fitted over the cells
// with |tgl| < 1.2, using several geometries G_i against the same R, each with
// its own checkResid tree on the same data (the hit pattern is taken from the
// first one). The fit uses only the c_L, so it does not depend on what the
// alignments changed.
//
// kinkResponse(): fit the coefficients.
//   root -b -q 'kinkResponse.C+("/data/devAlg/v10n/algdump.root",
//   "/data/devAlg/v10n/checkResid_LHC25ab_562862_AlgV10n_PV1.root",
//     "/data/devAlg/solA2v18/algdump.root,/data/devAlg/solA2v19/algdump.root,/data/devAlg/solA2v19rmc/algdump.root",
//     "/data/devAlg/solA2v18/checkResid_LHC25ab_562862_AlgRS_PV1.root,/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root,/data/devAlg/solA2v19rmc/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
//   (pp B0 run 562862, giving c_L = -33.9 -42.3 75.7 20.7 -21.6 1.2 0.04 for
//   Lr0..Lr6)
// kinkPredict(): predicted dPar2 harmonics of a modified geometry (e.g. a cure)
// from the measured ones of a solution:
//   root -b -q -e '.L kinkResponse.C+' -e
//   'kinkPredict("/data/devAlg/solA2v19/o2sim_geometry-aligned.root",
//   "/data/devAlg/solA2v19cure/cure012/o2sim_geometry-aligned.root",
//     "/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
// Geometry arguments are aligned geometry files (o2sim_geometry-aligned.root)
// or algDump outputs.
#include "kinkUtils.h"
#include <TDecompSVD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

// default c_L [urad/um], fitted on pp B0 562862 with V10n as reference and v18,
// v19, v19rmc
const std::array<double, 7> CLpp{-33.86, -42.28, 75.73, 20.72,
                                 -21.64, 1.18,   0.04};

void kinkResponse(const char *refGeom, const char *refRes, const char *geoms,
                  const char *ress, long nmax = 600000) {
  using namespace kink;
  auto gl = tokenize(geoms), rl = tokenize(ress);
  if (gl.size() != rl.size() || gl.empty()) {
    printf("need the same number (>0) of geometries and res trees\n");
    return;
  }
  int ng = gl.size();
  auto R = loadFrames(refGeom);
  std::vector<std::map<int, Frame>> G(ng);
  std::vector<const std::map<int, Frame> *> pg;
  for (int g = 0; g < ng; g++) {
    G[g] = loadFrames(gl[g].c_str());
    pg.push_back(&G[g]);
  }
  std::vector<CellShifts> D;
  {
    TFile f(rl[0].c_str());
    D = cellHitShifts((TTree *)f.Get("res"), R, pg, nmax);
  }
  // observed cell means: reference and each geometry
  std::vector<std::array<std::array<std::array<double, NP>, NT>, 2>> M(ng + 1),
      E(ng + 1); // [set][dPar0/2][it][ip]
  for (int s = 0; s <= ng; s++) {
    TFile f(s ? rl[s - 1].c_str() : refRes);
    auto t = (TTree *)f.Get("res");
    double m[NT][NP], e[NT][NP];
    for (int o = 0; o < 2; o++) {
      dParMeans(t, o ? 2 : 0, m, e);
      for (int a = 0; a < NT; a++) {
        for (int b = 0; b < NP; b++) {
          M[s][o][a][b] = m[a][b], E[s][o][a][b] = e[a][b];
        }
      }
    }
  }
  for (int o = 0; o < 2; o++) {
    TMatrixDSym A(7);
    TVectorD bv(7);
    double chi0 = 0;
    int nc = 0;
    for (int s = 0; s < ng; s++) {
      for (int a = 1; a < NT - 1; a++) { // |tgl| < 1.2
        for (int b = 0; b < NP; b++) {
          double y = M[s + 1][o][a][b] - M[0][o][a][b],
                 err = hypot(E[s + 1][o][a][b], E[0][o][a][b]);
          if (!std::isfinite(y) || err <= 0) {
            continue;
          }
          double w = 1. / (err * err);
          for (int i = 0; i < 7; i++) {
            bv[i] += w * D[s][a][b][i] * y;
            for (int j = 0; j < 7; j++) {
              A(i, j) += w * D[s][a][b][i] * D[s][a][b][j];
            }
          }
          chi0 += w * y * y;
          nc++;
        }
      }
    }
    TDecompSVD svd(A);
    bool ok;
    TVectorD c = svd.Solve(bv, ok);
    TVectorD Ac = A * c;
    double chi = chi0 - 2. * (c * bv) + (c * Ac);
    printf("%s change vs reference, %d cells (|tgl|<1.2): chi2 without model "
           "%.0f, with model %.0f\n  c_L [%s/um] Lr0..Lr6:",
           o ? "dPar2" : "dPar0", nc, chi0, chi, o ? "urad" : "um");
    for (int i = 0; i < 7; i++) {
      printf(" %7.2f", c[i]);
    }
    printf("\n");
    if (o) { // response to an OB ray rotation omega [urad] about R = 12 cm: d_L
             // = omega*(r_L-12)*1e-2 um
      double gain = 0;
      for (int l = 3; l < 7; l++) {
        gain += c[l] * (RLayer[l] - RComp) * 1e-2;
      }
      printf("  delta dPar2 = %.3f * omega for an OB rphi shift "
             "omega*(r-12cm): cure gain (makeKinkCure.C) = %.3f\n",
             gain, -gain);
    }
  }
}

void kinkPredict(const char *solGeom, const char *newGeom, const char *solRes,
                 long nmax = 600000) {
  using namespace kink;
  auto B = loadFrames(solGeom), A = loadFrames(newGeom);
  TFile f(solRes);
  auto t = (TTree *)f.Get("res");
  auto D = cellHitShifts(t, B, {&A}, nmax)[0];
  double m[NT][NP], e[NT][NP], p[NT][NP];
  dParMeans(t, 2, m, e);
  for (int a = 0; a < NT; a++) {
    for (int b = 0; b < NP; b++) {
      p[a][b] = m[a][b];
      for (int l = 0; l < 7; l++) {
        p[a][b] += CLpp[l] * D[a][b][l];
      }
    }
  }
  printf("dPar2 harmonics [urad] per tgl bin, measured on %s -> predicted for "
         "%s\n   tgl bin    |  n0          | |n1|         | |n2|\n",
         solGeom, newGeom);
  for (int a = 0; a < NT; a++) {
    printf("  [%+4.1f,%+4.1f] | %5.0f -> %5.0f | %4.0f -> %4.0f | %4.0f -> "
           "%4.0f\n",
           TE[a], TE[a + 1], harmonic(m[a], 0).Re(), harmonic(p[a], 0).Re(),
           TComplex::Abs(harmonic(m[a], 1)), TComplex::Abs(harmonic(p[a], 1)),
           TComplex::Abs(harmonic(m[a], 2)), TComplex::Abs(harmonic(p[a], 2)));
  }
}
