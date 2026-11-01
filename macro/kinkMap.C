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
// With curvCorr = true (B != 0 data only) dPar2 is also shown after subtracting k * dPar4, k being the regression slope of
// dPar2 on dPar4 over all tracks: in B+ OO data dPar2 is dominated by the IB/OB curvature difference.
//
// Usage (O2 environment, e.g. after `o2d`):
//   root -b -q 'kinkMap.C+("v19=/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root,cure=/data/devAlg/solA2v19cure/cure012/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
//   root -b -q 'kinkMap.C+("v19OO=/data/devAlg/solA2v19/checkResid_LHC25ae_AlgRS_PV1.root", true)'
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
  if (!curvCorr) {
    return;
  }
  // regression slope of dPar2 on dPar4 over the core
  const char* core = "trIBOut.mX>1&&trOBInw.mX>1&&abs(trIBOut.mP[2]-trOBInw.mP[2])<3e-3&&abs(trIBOut.mP[4]-trOBInw.mP[4])<0.3";
  long nt = std::min<long>(nmax, t->GetEntries());
  t->SetEstimate(nt + 1);
  long n = t->Draw("trIBOut.mP[2]-trOBInw.mP[2]:trIBOut.mP[4]-trOBInw.mP[4]", core, "goff", nt);
  double sx = 0, sy = 0, sxy = 0, syy = 0;
  for (long i = 0; i < n; i++) {
    double x = t->GetV1()[i], y = t->GetV2()[i];
    sx += x, sy += y, sxy += x * y, syy += y * y;
  }
  double k = (sxy / n - sx / n * sy / n) / (syy / n - sy / n * sy / n);
  cellMeans(t, Form("trIBOut.mP[2]-trOBInw.mP[2]-(%g)*(trIBOut.mP[4]-trOBInw.mP[4])", k), 4e-3, 1e6, m2, e2, "", nmax);
  printf("  dPar2 - k*dPar4, k = %.1f urad per 1e-3 (GeV/c)^-1\n   tgl bin    | n0  |n1|  n2(cos,sin)\n", k * 1e3);
  for (int it = 0; it < NT; it++) {
    TComplex p2 = harmonic(m2[it], 2);
    printf("  [%+4.1f,%+4.1f] | %5.0f %5.0f (%5.0f,%5.0f)\n", TE[it], TE[it + 1], harmonic(m2[it], 0).Re(), TComplex::Abs(harmonic(m2[it], 1)), p2.Re(), p2.Im());
  }
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
