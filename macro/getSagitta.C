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

// Estimate of the curvature (sagitta) bias of the ITS-only tracks from the charge-mirror symmetry of their q/pT spectrum.
//
//   root -l -b -q 'macro/getSagitta.C+("checkResid_LHC25as_568664_AlgRS_PV1.root", 6)'
//
// Input: the `res` tree of the CheckResid study (o2-check-resid-workflow). The q/pT used is that of `tr.track`, the
// ITS-only inward fit at the vertex, before the vertex update; the azimuth is that of the track direction at the vertex,
// phi = alpha + asin(snp), in [0, 2pi).
//
// Output: a header line with the phi sectors and units, then one line with the offset delta (and N+/N-) for all phi and
// for each of the nPhi equal sectors. The q/pT vs phi histogram is saved to `outFile` (if not empty).
//
// Method
// ------
// A geometric curvature bias shifts every measured q/pT by the same delta: x_meas = x_true + delta. The q/pT spectrum
// has a deep hole around 0 (it falls steeply with pT), mirror symmetric up to the charge ratio. For a trial delta,
// the counts in mirrored |x| windows around delta, N+(k) = N(delta + [a_k, b_k]) and N-(k) = N(delta - [b_k, a_k]),
// are compared: at the true delta their ratio N+/N- is the same in every window k. The ratio is fitted as a constant
// per trial delta, and delta is the value minimising the chi2 of this constant fit (scanned in steps of 1e-4 /GeV,
// counts in windows obtained by linear interpolation of the cumulative spectrum, so delta is not quantised by the
// histogram binning). The quoted error is the half width of the chi2_min + 1 interval.
//
// Sensitivity: with a local spectrum slope n = -dlnN/dln|x| (~2.6 at |q/pT| ~ 0.4), delta changes the mirror ratio by
// ~2*n*delta/|x|, i.e. most at high pT (small |x|). 1 unit of the curl/sagitta mode k (r*dphi = k*1e-6/cm * r^2)
// corresponds to ~13e-4 /GeV at 0.5 T.
//
// Charge abundances: the charge ratio is a free parameter, so a pT-independent charge asymmetry is absorbed. Only a
// change of N+/N- with pT inside the windows (|q/pT| 0.06-0.85, pT ~1.2-17 GeV/c) biases delta: a 1% drift across the
// window fakes ~1.7e-4 /GeV. Checks:
//  - a geometric bias flips its sign with the field polarity, the physics charge asymmetry does not:
//    (delta(B+) + delta(B-))/2 bounds the non-geometric part;
//  - restricting the windows to low or high pT must give compatible values;
//  - the difference between two geometries on the same tracks is free of it.
// The per-sector values are less safe as absolute numbers: charge-dependent acceptance holes (dead chips, the cluster
// requirement, TPC sector edges) remove positive and negative tracks at different phi (they bend in opposite
// directions); this cancels when integrating over phi, but not per sector, and it flips with the field like a geometric
// bias. Differences between geometries on the same data are fine per sector as well.

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

namespace sagitta
{
// edges of the |q/pT - delta| windows, (GeV/c)^-1
const std::vector<double> Edges = {0.06, 0.10, 0.14, 0.18, 0.22, 0.27, 0.33, 0.40, 0.50, 0.65, 0.85};
constexpr double DeltaMax = 0.025, DeltaStep = 1e-4; // scanned delta range and step, (GeV/c)^-1
constexpr double MinCounts = 20.;                    // a window pair is used only if both sides have this many tracks
constexpr int NPhiBins = 36;                         // phi binning of the stored histogram (nPhi must divide it)

struct Result {
  double delta = 0., error = 0., chargeRatio = 0.;
};

/// cumulative spectrum: cum[i] = number of entries in the bins 1..i
std::vector<double> cumulative(const TH1D& h)
{
  std::vector<double> cum(h.GetNbinsX() + 1, 0.);
  for (int i = 1; i <= h.GetNbinsX(); i++) {
    cum[i] = cum[i - 1] + h.GetBinContent(i);
  }
  return cum;
}

/// number of entries in [lo, hi], with linear interpolation of the cumulative spectrum inside the bins
double countIn(const std::vector<double>& cum, double xmin, double binW, double lo, double hi)
{
  auto atX = [&](double x) {
    double u = (x - xmin) / binW;
    int i = std::floor(u);
    if (i < 0) {
      return 0.;
    }
    if (i >= (int)cum.size() - 1) {
      return cum.back();
    }
    return cum[i] + (u - i) * (cum[i + 1] - cum[i]);
  };
  return atX(hi) - atX(lo);
}

/// chi2 of the constant fit of the mirror ratios N+(k)/N-(k) for the trial offset delta; ratio = fitted constant
double mirrorChi2(const std::vector<double>& cum, double xmin, double binW, double delta, double& ratio)
{
  std::vector<double> r, er;
  for (size_t k = 0; k + 1 < Edges.size(); k++) {
    double np = countIn(cum, xmin, binW, delta + Edges[k], delta + Edges[k + 1]);
    double nm = countIn(cum, xmin, binW, delta - Edges[k + 1], delta - Edges[k]);
    if (np < MinCounts || nm < MinCounts) {
      continue;
    }
    r.push_back(np / nm);
    er.push_back(np / nm * std::sqrt(1. / np + 1. / nm));
  }
  double sw = 0., swr = 0.;
  for (size_t k = 0; k < r.size(); k++) {
    sw += 1. / (er[k] * er[k]);
    swr += r[k] / (er[k] * er[k]);
  }
  ratio = sw > 0. ? swr / sw : 0.;
  double chi2 = 0.;
  for (size_t k = 0; k < r.size(); k++) {
    chi2 += std::pow((r[k] - ratio) / er[k], 2);
  }
  return r.size() > 1 ? chi2 : 1e30;
}

/// scan delta, return the chi2 minimum and the half width of its chi2_min+1 interval
Result estimate(const TH1D& h)
{
  auto cum = cumulative(h);
  double xmin = h.GetXaxis()->GetXmin(), binW = h.GetXaxis()->GetBinWidth(1);
  std::vector<double> deltas, chi2s;
  Result res;
  double chi2Min = 1e30;
  for (double d = -DeltaMax; d <= DeltaMax + 1e-12; d += DeltaStep) {
    double ratio = 0.;
    double c2 = mirrorChi2(cum, xmin, binW, d, ratio);
    deltas.push_back(d);
    chi2s.push_back(c2);
    if (c2 < chi2Min) {
      chi2Min = c2;
      res.delta = d;
      res.chargeRatio = ratio;
    }
  }
  double lo = 1., hi = -1.;
  for (size_t i = 0; i < chi2s.size(); i++) {
    if (chi2s[i] < chi2Min + 1.) {
      lo = std::min(lo, deltas[i]);
      hi = std::max(hi, deltas[i]);
    }
  }
  res.error = (hi - lo) / 2.;
  return res;
}
} // namespace sagitta

void getSagitta(const char* fname, int nPhi = 6, const char* outFile = "hq2pt.root")
{
  using namespace sagitta;
  if (nPhi < 1 || NPhiBins % nPhi) {
    printf("nPhi=%d must divide %d\n", nPhi, NPhiBins);
    return;
  }
  auto fin = TFile::Open(fname);
  if (!fin || fin->IsZombie()) {
    printf("cannot open %s\n", fname);
    return;
  }
  auto tree = (TTree*)fin->Get("res");
  if (!tree) {
    printf("no res tree in %s\n", fname);
    return;
  }
  auto h = new TH2F("hq", "q/pT vs phi_dir;#phi;q/p_{T}", NPhiBins, 0, 2 * TMath::Pi(), 2000, -2, 2);
  tree->Draw("tr.track.mP[4]:fmod(tr.track.mAlpha+asin(tr.track.mP[2])+4*TMath::Pi(),2*TMath::Pi())>>hq", "", "goff");
  h->SetDirectory(nullptr); // detach only after Draw, which looks the histogram up in gDirectory

  // header: sectors in degrees, units
  printf("delta(q/pT) [1e-4 (GeV/c)^-1] (+- chi2+1 half width), N+/N- for all phi; file %s\n", fname);
  printf("%-27s", "all phi");
  int group = NPhiBins / nPhi;
  for (int i = 0; i < nPhi; i++) {
    printf(" %14s", Form("[%3.0f,%3.0f)", 360. * i / nPhi, 360. * (i + 1) / nPhi));
  }
  printf("\n");

  std::unique_ptr<TH1D> pAll(h->ProjectionY("pall"));
  auto all = estimate(*pAll);
  printf("%6.1f +- %4.1f  N+/N- %.3f", all.delta * 1e4, all.error * 1e4, all.chargeRatio);
  for (int i = 0; i < nPhi; i++) {
    std::unique_ptr<TH1D> p(h->ProjectionY(Form("p%d", i), i * group + 1, (i + 1) * group));
    auto sec = estimate(*p);
    printf(" %6.1f +- %4.1f", sec.delta * 1e4, sec.error * 1e4);
  }
  printf("\n");

  if (outFile && outFile[0]) {
    TFile fo(outFile, "recreate");
    h->Write();
    fo.Close();
  }
  delete h;
}
