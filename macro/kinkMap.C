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

// IB/OB kink detection from checkResid `res` trees, see doc/handovers/kink_cure.md.
//
// For each input tree, per tgl bin (kink::TE) and in 24 phi bins (phi of the IB sub-track at R = 12 cm), the gaussian-fitted
// mean of dPar0 = y_IB - y_OB [um] and dPar2 = snp_IB - snp_OB [urad] at R = 12 cm, and their phi harmonics n = 0, |n=1|,
// n = 2 (cos, sin). R0 = 12 - dPar0/dPar2 of the n = 2 term is the radius where the two sub-tracks cross: R0 ~ 12 cm means a
// purely directional kink, i.e. a relative IB/OB rotation about R = 12 cm.
// With curvCorr = true (B != 0 data) dPar2 is also shown with both sub-tracks constrained to their common q/pT
// (kink::DPar2QC): without it dPar2 is dominated by the IB sub-track curvature error, the cell errors being ~3x larger.
// The statistical errors of the harmonics are printed as well (from the cell errors, n0: e/sqrt(NP), n > 0: e*sqrt(2/NP)).
//
// Usage (O2 environment, e.g. after `o2d`):
//   root -b -q 'kinkMap.C+("v19=/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root,cure=/data/devAlg/solA2v19cure/cure012/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
//   root -b -q 'kinkMap.C+("v19Bp=/data/devAlg/solA2v19/checkResid_LHC25ab_562837_AlgRS_PV1.root", true)'
// nmax limits the number of tracks read per tree.
#include "kinkUtils.h"

void kinkMapOne(const char* fname, const char* tag, bool curvCorr, long nmax)
{
  using namespace kink;
  TFile f(fname);
  auto t = (TTree*)f.Get("res");
  if (!t) {
    printf("no res tree in %s\n", fname);
    return;
  }
  double m0[NT][NP], e0[NT][NP], m2[NT][NP], e2[NT][NP];
  dParMeans(t, 0, m0, e0, "", nmax);
  dParMeans(t, 2, m2, e2, "", nmax);
  printf("%s (%s)\n   tgl bin    | dPar0[um] n0  |n1|  n2(cos,sin)  | dPar2[urad] n0 |n1|  n2(cos,sin)  | R0(n2) [cm]\n", tag, fname);
  for (int it = 0; it < NT; it++) {
    TComplex p0 = harmonic(m0[it], 2), p2 = harmonic(m2[it], 2), r = p0 / p2; // um/urad -> x100 cm
    printf("  [%+4.1f,%+4.1f] | %5.1f %5.1f (%5.1f,%5.1f) | %5.0f %5.0f (%5.0f,%5.0f) | %6.1f\n", TE[it], TE[it + 1],
           harmonic(m0[it], 0).Re(), TComplex::Abs(harmonic(m0[it], 1)), p0.Re(), p0.Im(),
           harmonic(m2[it], 0).Re(), TComplex::Abs(harmonic(m2[it], 1)), p2.Re(), p2.Im(), RComp - 100 * r.Re());
  }
  auto printErr = [&](const char* title) {
    printf("  %s [urad], harmonics with statistical errors\n   tgl bin    | n0         | n2 cos     | n2 sin\n", title);
    for (int it = 0; it < NT; it++) {
      double se = 0;
      int nc = 0;
      for (int i = 0; i < NP; i++) {
        if (std::isfinite(e2[it][i]) && e2[it][i] > 0) {
          se += e2[it][i] * e2[it][i], nc++;
        }
      }
      double ec = nc ? sqrt(se / nc) : 0;
      TComplex p2 = harmonic(m2[it], 2);
      printf("  [%+4.1f,%+4.1f] | %5.0f±%4.0f | %5.0f±%4.0f | %5.0f±%4.0f\n", TE[it], TE[it + 1], harmonic(m2[it], 0).Re(), ec / sqrt(NP),
             p2.Re(), ec * sqrt(2. / NP), p2.Im(), ec * sqrt(2. / NP));
    }
  };
  printErr("dPar2");
  if (!curvCorr) {
    return;
  }
  dParMeans(t, 2, m2, e2, "", nmax, true);
  printErr("dPar2, sub-tracks constrained to the common q/pT");
}

/// inputs: comma-separated "label=file" list (a bare file name is labelled by itself)
void kinkMap(const char* inputs, bool curvCorr = false, long nmax = TTree::kMaxEntries)
{
  for (const auto& tok : kink::tokenize(inputs)) {
    auto eq = tok.find('=');
    std::string tag = eq == std::string::npos ? tok : tok.substr(0, eq), fn = eq == std::string::npos ? tok : tok.substr(eq + 1);
    kinkMapOne(fn.c_str(), tag.c_str(), curvCorr, nmax);
  }
}
