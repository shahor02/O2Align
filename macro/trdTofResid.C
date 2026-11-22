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

// TRD / TOF residuals of the `gbltrack` tree of debug_alg.root (see doc/DebugInfo.md), in an o2::HistoManager.
//
// trdTofResid(): per TRD chamber (sector, layer, stack) dY vs snp (nSnp bins in |snp| < maxSnp) and dZ vs tgl of the track
// (nTglBins bins between Zmin/X and Zmax/X of the chamber, X = nominal tracking X of the layer), per TOF sector and
// Z slice (nTOFSlices equal slices of the strip centre Z, i.e. groups of whole strips) the same, the tgl range being
// the Z limits of the slice divided by the strip radius. All the residuals are in cm. Only points with hasMeas() are used.
// The TH2F are fitted with FitSlicesY (gaus in [-resRange, resRange]), the mean and sigma TH1 are kept in the manager. A slice
// whose fit has its mean outside the fit window, or a sigma not within (0, resRange], takes the profile mean and spread instead.
// Per TRD layer and stack / TOF Z slice, the mean and sigma vs sector of the whole residual distribution of each chamber / sector
// slice (gaussian fit with the same validity check, otherwise the mean and the standard deviation) are added as well.
//
// o2::HistoManager IDs:   Sect*10000 + Lr*1000 + St*100 + iYZ*10  (+1: mean, +2: sigma of the fitted slices)
//   Sect 0:17, Lr 0:5 = TRD layer, 6 = TOF, St = TRD stack or TOF Z slice, iYZ 0 (Y) or 1 (Z).
//   per-sector summaries: 200000 + Lr*1000 + St*100 + iYZ*10 (+1: mean, +2: sigma), 18 bins in sector.
//
// drawTRDTOFResid(): superimposes several such managers, the same way as `o2-check-resid-workfow --draw-external-only`
// does it with the ext_hm_list, ext_leg_list, ext_color_list, ext_marker_list keys of CheckResidConfig (see
// doc/reResidAll.sh). First one page of per-sector summaries per TRD layer and one for TOF, then one page per sector and layer;
// each page has 4 rows (mean dY, sigma dY, mean dZ, sigma dZ) and one column per TRD stack or TOF Z slice.
//
// Usage (O2 environment, libO2AlignRS and the O2Align headers findable, ROOT_INCLUDE_PATH):
//   root -b -q 'trdTofResid.C+("debug_alg.root", "resTRDTOF_Unb_hman.root", "o2sim_geometry-aligned.root", 1)'
//   root -b -q 'trdTofResid.C+("files.txt", "resTRDTOF_Pre_hman.root", "", 0)'    // text file: one root file per line
//   root -l 'trdTofResid.C+' -e 'drawTRDTOFResid("checkresid.ext_hm_list=a_hman.root,b_hman.root;checkresid.ext_leg_list=A,B;checkresid.ext_color_list=kBlack,kRed;checkresid.ext_marker_list=20,24;", "cmp.pdf")'

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TH1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fmt/format.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "CommonUtils/StringUtils.h"
#include "DetectorsBase/GeometryManager.h"
#include "GlobalTrackingStudy/CheckResidConfig.h"
#include "GlobalTrackingStudy/HistoManager.h"
#include "TOFBase/Geo.h"
#include "TRDBase/Geometry.h"
#include "O2Align/AlignmentDbg.h"

namespace trdtof
{

constexpr int NSect = 18, NTRDLr = 6, NTRDSt = 5, LrTOF = 6, DetTRD = 3, DetTOF = 4;
constexpr int NStripSM = o2::tof::Geo::NSTRIPXSECTOR;
constexpr int MaxTOFSlices = 10;
constexpr int SummaryBase = 200000; // IDs of the per-sector summaries, above those of the per-chamber histograms
constexpr std::array<const char*, 3> ResTypeName{"Pre", "Unb", "Fit"};

inline int histoID(int sect, int lr, int st, int iyz) { return sect * 10000 + lr * 1000 + st * 100 + iyz * 10; }
inline int summaryID(int lr, int col, int iyz) { return SummaryBase + lr * 1000 + col * 100 + iyz * 10; }

/// strips grouped in Z slices of the TOF sector, from the ideal geometry
struct TOFSlices {
  int nSlices = 0;
  std::vector<int> sliceOfStrip;             // by continuous strip number, -1 for missing strips
  std::vector<std::array<float, 3>> zLoHiR;  // per slice: Z min, Z max of the slice, mean radius of its strips
  std::vector<std::array<int, 2>> stripRange; // per slice and sector (index slice*NSect + sect): first, last strip number in the sector
};

/// chain of the `gbltrack` trees from a root file or from a text file with one root file per line (not recursive)
TChain* loadDbgChain(const char* inpData, const char* chName)
{
  auto* chain = new TChain(chName);
  TString inp = inpData;
  if (inp.EndsWith(".root")) {
    chain->AddFile(inpData);
  } else {
    std::ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n", inpData);
      delete chain;
      return nullptr;
    }
    TString flName;
    while (flName.ReadLine(inpf)) {
      flName = flName.Strip(TString::kBoth, ' ');
      if (flName.IsNull() || flName.BeginsWith("//") || flName.BeginsWith("#")) {
        continue;
      }
      flName = flName.Strip(TString::kBoth, ',');
      flName = flName.Strip(TString::kBoth, '"');
      printf("Adding %s\n", flName.Data());
      chain->AddFile(flName.Data());
    }
  }
  if (chain->GetEntries() < 1) {
    printf("Obtained chain is empty\n");
    delete chain;
    return nullptr;
  }
  printf("Opened %s chain with %lld entries\n", chName, chain->GetEntries());
  return chain;
}

TOFSlices makeTOFSlices(int nSlices)
{
  TOFSlices ts;
  ts.nSlices = nSlices;
  ts.sliceOfStrip.assign(NSect * NStripSM, -1);
  std::vector<float> zc(NSect * NStripSM, 0.f), rc(NSect * NStripSM, 0.f);
  float zMin = 1e9, zMax = -1e9;
  for (int isect = 0; isect < NSect; ++isect) {
    for (int istr = 0; istr < NStripSM; ++istr) {
      auto* pne = gGeoManager->GetAlignableEntry(Form("TOF/sm%02d/strip%02d", isect, istr + 1));
      if (!pne) {
        continue;
      }
      auto* pn = pne->GetPhysicalNode();
      if (!pn) {
        pn = gGeoManager->MakeAlignablePN(pne);
      }
      const double loc[3] = {0, 0, 0};
      double glo[3];
      pn->GetMatrix()->LocalToMaster(loc, glo);
      const int sid = isect * NStripSM + istr;
      zc[sid] = glo[2];
      rc[sid] = std::hypot(glo[0], glo[1]);
      zMin = std::min(zMin, zc[sid]);
      zMax = std::max(zMax, zc[sid]);
    }
  }
  const float halfPitch = 0.5f * (zMax - zMin) / (NStripSM - 1);
  zMin -= halfPitch;
  zMax += halfPitch;
  const float dz = (zMax - zMin) / nSlices;
  ts.zLoHiR.assign(nSlices, {0.f, 0.f, 0.f});
  std::vector<int> nStr(nSlices, 0);
  ts.stripRange.assign(nSlices * NSect, {NStripSM, -1});
  for (int is = 0; is < nSlices; ++is) {
    ts.zLoHiR[is][0] = zMin + is * dz;
    ts.zLoHiR[is][1] = zMin + (is + 1) * dz;
  }
  for (int sid = 0; sid < NSect * NStripSM; ++sid) {
    if (ts.sliceOfStrip[sid] < 0 && rc[sid] > 0.f) {
      const int is = std::min(nSlices - 1, std::max(0, int((zc[sid] - zMin) / dz)));
      ts.sliceOfStrip[sid] = is;
      ts.zLoHiR[is][2] += rc[sid];
      nStr[is]++;
      auto& rng = ts.stripRange[is * NSect + sid / NStripSM];
      rng[0] = std::min(rng[0], sid % NStripSM);
      rng[1] = std::max(rng[1], sid % NStripSM);
    }
  }
  for (int is = 0; is < nSlices; ++is) {
    ts.zLoHiR[is][2] = nStr[is] ? ts.zLoHiR[is][2] / nStr[is] : 380.f;
  }
  return ts;
}

/// books dY vs snp and dZ vs tgl for one chamber / slice, tgl in [tglMin, tglMax]
void bookPair(o2::HistoManager& hm, int id, const std::string& name, const std::string& title, int resType, int nSnp, float maxSnp,
              int nTgl, float tglMin, float tglMax, int nBinRes, float resRangeY, float resRangeZ)
{
  for (int iyz = 0; iyz < 2; ++iyz) {
    const std::string rname = fmt::format("d{}{}", iyz ? 'Z' : 'Y', ResTypeName[resType]);
    const std::string nm = fmt::format("{}_{}", name, rname);
    const std::string tit = fmt::format("{} {}", title, rname);
    const std::string ytit = fmt::format("{} (cm)", iyz ? "#DeltaZ" : "#DeltaY");
    TH2F* h2 = iyz == 0 ? new TH2F(nm.c_str(), fmt::format("{};snp;{}", tit, ytit).c_str(), nSnp, -maxSnp, maxSnp, nBinRes, -resRangeY, resRangeY)
                        : new TH2F(nm.c_str(), fmt::format("{};tg#lambda;{}", tit, ytit).c_str(), nTgl, tglMin, tglMax, nBinRes, -resRangeZ, resRangeZ);
    hm.addHisto(h2, id + iyz * 10);
  }
}

void bookHistos(o2::HistoManager& hm, const TOFSlices& ts, int resType, int nSnp, float maxSnp, int nTglBins, int nBinRes, float resRangeY, float resRangeZ)
{
  auto* trdGeo = o2::trd::Geometry::instance();
  for (int isect = 0; isect < NSect; ++isect) {
    for (int ilr = 0; ilr < NTRDLr; ++ilr) {
      for (int ist = 0; ist < NTRDSt; ++ist) {
        const auto symChamber = fmt::format("TRD/sm{:02d}/st{}/pl{}", isect, ist, ilr);
        if (!gGeoManager->GetAlignableEntry(symChamber.c_str())) { // missing chamber
          continue;
        }
        const int det = o2::trd::Geometry::getDetector(ilr, ist, isect);
        const float x = o2::trd::Geometry::getTime0(ilr);
        const float z0 = trdGeo->getRowEnd(ilr, ist), z1 = trdGeo->getRow0(ilr, ist);
        bookPair(hm, histoID(isect, ilr, ist, 0), fmt::format("trd_s{:02d}_l{}_st{}", isect, ilr, ist), fmt::format("{} ch{}", symChamber, det),
                 resType, nSnp, maxSnp, nTglBins, std::min(z0, z1) / x, std::max(z0, z1) / x, nBinRes, resRangeY, resRangeZ);
      }
    }
    for (int is = 0; is < ts.nSlices; ++is) {
      const auto& rng = ts.stripRange[is * NSect + isect];
      if (rng[1] < 0) {
        continue;
      }
      const auto& lim = ts.zLoHiR[is];
      bookPair(hm, histoID(isect, LrTOF, is, 0), fmt::format("tof_s{:02d}_z{}", isect, is), fmt::format("TOF/sm{:02d}/strip{}:{}", isect, rng[0], rng[1]),
               resType, nSnp, maxSnp, nTglBins, lim[0] / lim[2], lim[1] / lim[2], nBinRes, resRangeY, resRangeZ);
    }
  }
}

/// the ID of the TH2F receiving the point, -1 if the point is not TRD/TOF or has no booked histogram
int pointToHistoID(const o2::alignrs::DbgPoint& p, const TOFSlices& ts)
{
  if (p.det == DetTRD && p.sens >= 0) {
    return histoID(p.sens / (NTRDLr * NTRDSt), p.sens % NTRDLr, (p.sens % (NTRDLr * NTRDSt)) / NTRDLr, 0);
  }
  if (p.det == DetTOF && p.sens >= 0 && p.sens < (int)ts.sliceOfStrip.size() && ts.sliceOfStrip[p.sens] >= 0) {
    return histoID(p.sens / NStripSM, LrTOF, ts.sliceOfStrip[p.sens], 0);
  }
  return -1;
}

void fillHistos(TChain* chain, o2::HistoManager& hm, const TOFSlices& ts, int resType, bool onlyAccepted, long maxTracks)
{
  auto* dt = new o2::alignrs::DbgTrack;
  chain->SetBranchAddress("t", &dt);
  const long nEnt = maxTracks > 0 ? std::min<long>(maxTracks, chain->GetEntries()) : chain->GetEntries();
  for (long i = 0; i < nEnt; ++i) {
    chain->GetEntry(i);
    if (onlyAccepted && !dt->accepted) {
      continue;
    }
    const float tgl = dt->getTgl();
    for (const auto& p : dt->points) {
      if (!p.hasMeas()) {
        continue;
      }
      const int id = pointToHistoID(p, ts);
      if (id < 0) {
        continue;
      }
      float dy = p.dyPre, dz = p.dzPre;
      if (resType == 1) {
        if (p.sigYUnb <= 0.f) { // the unbiased residuals were not computed
          continue;
        }
        dy = p.dyUnb;
        dz = p.dzUnb;
      } else if (resType == 2) {
        dy = p.dyFit;
        dz = p.dzFit;
      }
      if (auto* hy = hm.getHisto2F(id)) {
        hy->Fill(p.snp, dy);
      }
      if (auto* hz = hm.getHisto2F(id + 10)) {
        hz->Fill(tgl, dz);
      }
    }
  }
  delete dt;
  chain->ResetBranchAddresses();
}

/// profile of the 2D histogram converted to TH1D: mean (error of the mean) or sigma (spread, error sigma/sqrt(2n)) per bin; empty bins have zero error
TH1* profileToTH1(TH2F* h2, bool sigma, const char* name)
{
  auto* prf = h2->ProfileX(Form("%s_pfx%d", h2->GetName(), (int)sigma), 1, -1, sigma ? "s" : "");
  auto* h1 = new TH1D(name, h2->GetTitle(), h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  h1->SetDirectory(nullptr);
  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    const double n = prf->GetBinEntries(i);
    if (n < 2) { // the spread of a single entry is undefined
      continue;
    }
    h1->SetBinContent(i, sigma ? prf->GetBinError(i) : prf->GetBinContent(i));
    h1->SetBinError(i, sigma ? prf->GetBinError(i) / std::sqrt(2. * n) : prf->GetBinError(i));
  }
  delete prf;
  return h1;
}

/// a gaussian fit result is accepted if its mean is within the fit window [-rng, rng] and its sigma within (0, rng]
inline bool isGoodGausFit(double mean, double meanErr, double sigma, double rng)
{
  return std::isfinite(mean) && std::isfinite(meanErr) && std::isfinite(sigma) && meanErr > 0. && std::abs(mean) <= rng && sigma > 0. && sigma <= rng;
}

/// FitSlicesY of h2 (gaus within the Y range of the histogram) into hmean / hsig; a fitted slice failing isGoodGausFit is replaced
/// by its profile mean and spread (or cleared if the slice has fewer than 2 entries)
void fitSlicesWithFallback(TH2F* h2, int minSliceStat, TH1*& hmean, TH1*& hsig)
{
  TObjArray arr;
  arr.SetOwner(true);
  const double rng = h2->GetYaxis()->GetXmax();
  TF1 gs("gs", "gaus", -rng, rng);
  h2->FitSlicesY(&gs, 0, -1, minSliceStat, "QNR", &arr);
  hmean = (TH1*)arr.RemoveAt(1);
  hsig = (TH1*)arr.RemoveAt(2);
  arr.Delete();
  if (!hmean || !hsig) {
    return;
  }
  std::unique_ptr<TH1> pmean(profileToTH1(h2, false, Form("%s_pm", h2->GetName()))), psig(profileToTH1(h2, true, Form("%s_ps", h2->GetName())));
  int nRepl = 0;
  for (int ib = 1; ib <= hmean->GetNbinsX(); ++ib) {
    const bool fitted = hmean->GetBinContent(ib) != 0. || hmean->GetBinError(ib) != 0.;
    const double sig = std::abs(hsig->GetBinContent(ib));
    if (!fitted || isGoodGausFit(hmean->GetBinContent(ib), hmean->GetBinError(ib), sig, rng)) {
      hsig->SetBinContent(ib, sig);
      continue;
    }
    hmean->SetBinContent(ib, pmean->GetBinContent(ib));
    hmean->SetBinError(ib, pmean->GetBinError(ib));
    hsig->SetBinContent(ib, psig->GetBinContent(ib));
    hsig->SetBinError(ib, psig->GetBinError(ib));
    nRepl++;
  }
  if (nRepl) {
    printf("%s: %d failed slice fit(s) replaced by the profile\n", h2->GetName(), nRepl);
  }
}

/// mean, sigma and their errors of the whole Y distribution of h2: gaussian fit within the Y range if it has at least minEntriesFit
/// entries and the fit is good (isGoodGausFit), otherwise the mean and the standard deviation. False if fewer than 2 entries.
bool fitWholeY(TH2F* h2, int minEntriesFit, std::array<double, 4>& res)
{
  std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("%s_py", h2->GetName()), 1, h2->GetNbinsX()));
  proj->SetDirectory(nullptr);
  const double n = proj->GetEntries();
  if (n < 2) {
    return false;
  }
  const double rng = h2->GetYaxis()->GetXmax();
  if (n >= minEntriesFit) {
    TF1 gs("gsw", "gaus", -rng, rng);
    gs.SetParameters(proj->GetMaximum(), proj->GetMean(), proj->GetStdDev());
    if (proj->Fit(&gs, "QNR0") == 0) {
      const double sig = std::abs(gs.GetParameter(2));
      if (isGoodGausFit(gs.GetParameter(1), gs.GetParError(1), sig, rng)) {
        res = {gs.GetParameter(1), gs.GetParError(1), sig, gs.GetParError(2)};
        return true;
      }
    }
  }
  const double sd = proj->GetStdDev();
  res = {proj->GetMean(), sd / std::sqrt(n), sd, sd / std::sqrt(2. * n)};
  return true;
}

/// per TRD layer and stack / TOF Z slice and per Y/Z: mean and sigma vs sector of the whole residual distribution of each chamber /
/// sector slice (fitWholeY), IDs summaryID(lr, col, iyz) + 1 (mean) / + 2 (sigma). Booked also if all sectors are empty.
void makeSectorSummaries(o2::HistoManager& hm, int minEntriesFit)
{
  for (int ilr = 0; ilr <= LrTOF; ++ilr) {
    const int nCol = ilr == LrTOF ? MaxTOFSlices : NTRDSt;
    for (int ic = 0; ic < nCol; ++ic) {
      // the TOF slices exist up to nTOFSlices of the filling; the tgl range of the slice is that of its dZ histograms
      TH2F* h2ref[2] = {nullptr, nullptr};
      for (int isect = 0; isect < NSect && !(h2ref[0] && h2ref[1]); ++isect) {
        for (int iyz = 0; iyz < 2; ++iyz) {
          if (!h2ref[iyz]) {
            h2ref[iyz] = hm.getHisto2F(histoID(isect, ilr, ic, iyz));
          }
        }
      }
      if (ilr == LrTOF && !h2ref[0] && !h2ref[1]) {
        continue;
      }
      for (int iyz = 0; iyz < 2; ++iyz) {
        const std::string res = h2ref[iyz] ? std::string(h2ref[iyz]->GetName()).substr(std::string(h2ref[iyz]->GetName()).rfind('_') + 1) : (iyz ? "dZ" : "dY");
        const std::string where = ilr == LrTOF ? (h2ref[1] ? fmt::format("TOF tg#lambda {:.2f}:{:.2f}", h2ref[1]->GetXaxis()->GetXmin(), h2ref[1]->GetXaxis()->GetXmax()) : fmt::format("TOF z{}", ic))
                                               : fmt::format("TRD L{} st{}", ilr, ic);
        const std::string base = ilr == LrTOF ? fmt::format("sumSect_tof_z{}_{}", ic, res) : fmt::format("sumSect_trd_l{}_st{}_{}", ilr, ic, res);
        auto* hmean = new TH1D(Form("%s_1", base.c_str()), Form("%s;sector;<%s> (cm)", where.c_str(), res.c_str()), NSect, -0.5, NSect - 0.5);
        auto* hsig = new TH1D(Form("%s_2", base.c_str()), Form("%s;sector;#sigma(%s) (cm)", where.c_str(), res.c_str()), NSect, -0.5, NSect - 0.5);
        for (int isect = 0; isect < NSect; ++isect) {
          auto* h2 = hm.getHisto2F(histoID(isect, ilr, ic, iyz));
          std::array<double, 4> r{};
          if (!h2 || !fitWholeY(h2, minEntriesFit, r)) {
            continue;
          }
          hmean->SetBinContent(isect + 1, r[0]);
          hmean->SetBinError(isect + 1, r[1]);
          hsig->SetBinContent(isect + 1, r[2]);
          hsig->SetBinError(isect + 1, r[3]);
        }
        const int id = summaryID(ilr, ic, iyz);
        hm.addHisto(hmean, id + 1);
        hm.addHisto(hsig, id + 2);
      }
    }
  }
}

/// mean/sigma of all the base 2D histograms of the manager, added with IDs +1/+2: FitSlicesY with the profile fallback for the failed
/// slices (fitSlicesWithFallback) if it has at least minEntriesFit entries, otherwise the profile converted to TH1D (mean and spread);
/// then the per-sector summaries (makeSectorSummaries)
void postProcessHistos(o2::HistoManager& hm, int minSliceStat, int minEntriesFit)
{
  printf("Fitting histos\n");
  for (int isect = 0; isect < NSect; ++isect) {
    for (int ilr = 0; ilr <= LrTOF; ++ilr) {
      for (int ist = 0; ist < (ilr == LrTOF ? MaxTOFSlices : NTRDSt); ++ist) {
        for (int iyz = 0; iyz < 2; ++iyz) {
          const int id = histoID(isect, ilr, ist, iyz);
          auto* h2 = hm.getHisto2F(id);
          if (!h2 || h2->GetEntries() < 1) {
            continue;
          }
          TH1 *hmean = nullptr, *hsig = nullptr;
          if (h2->GetEntries() < minEntriesFit) {
            hmean = profileToTH1(h2, false, Form("%s_1", h2->GetName()));
            hsig = profileToTH1(h2, true, Form("%s_2", h2->GetName()));
          } else {
            fitSlicesWithFallback(h2, minSliceStat, hmean, hsig);
          }
          if (hmean) {
            hmean->SetTitle(Form("<%s>", h2->GetTitle()));
            hmean->SetXTitle(h2->GetXaxis()->GetTitle());
            hmean->SetYTitle(Form("<%s>", h2->GetYaxis()->GetTitle()));
            hm.addHisto(hmean, id + 1);
          }
          if (hsig) {
            hsig->SetTitle(Form("#sigma(%s)", h2->GetTitle()));
            hsig->SetXTitle(h2->GetXaxis()->GetTitle());
            hsig->SetYTitle(Form("#sigma(%s)", h2->GetYaxis()->GetTitle()));
            hm.addHisto(hsig, id + 2);
          }
        }
      }
    }
  }
  makeSectorSummaries(hm, minEntriesFit);
}

} // namespace trdtof

/// Fill and fit the TRD/TOF residual histograms.
/// inpName   : debug_alg.root, or a text file with one root file per line
/// outName   : output file of the o2::HistoManager
/// geomName  : geometry file (alignable entries select the existing chambers/strips and give the TOF strip positions), "" = default
/// resType   : 0: dPre, 1: dUnb, 2: dFit
/// maxSnp, nSnp          : |snp| range and number of snp bins
/// nTglBins              : number of tgl bins (limits from the Z extent of the chamber / TOF slice divided by X)
/// nTOFSlices            : number of TOF Z slices
/// nBinRes, resRangeY/Z   : bins and half-ranges (cm) of the residual axes of dY and dZ, also the ranges of the gaussian fits
/// minEntriesFit         : histograms with fewer entries use the profile (mean and spread) instead of FitSlicesY
/// minSliceStat          : minimal number of non-empty bins in the slice to be fitted
/// onlyAccepted          : skip the tracks rejected by the chi2 cut of the GBL fit
/// maxTracks             : process at most that many tracks (<=0: all)
void trdTofResid(const char* inpName, const char* outName = "trdTofResid_hman.root", const char* geomName = "", int resType = 0,
                 float maxSnp = 0.4, int nSnp = 7, int nTglBins = 5, int nTOFSlices = 5, int nBinRes = 100, float resRangeY = 3,
                 float resRangeZ = 5, int minSliceStat = 3, int minEntriesFit = 100, bool onlyAccepted = true, long maxTracks = -1)
{
  using namespace trdtof;
  if (nTOFSlices < 1 || nTOFSlices > MaxTOFSlices) {
    printf("nTOFSlices must be within 1:10\n");
    return;
  }
  if (resType < 0 || resType > 2) {
    printf("resType must be 0 (dPre), 1 (dUnb) or 2 (dFit)\n");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geomName);
  auto* chain = loadDbgChain(inpName, "gbltrack");
  if (!chain) {
    return;
  }
  const auto ts = makeTOFSlices(nTOFSlices);
  o2::HistoManager hm("", outName);
  hm.SetName(gSystem->BaseName(outName));
  bookHistos(hm, ts, resType, nSnp, maxSnp, nTglBins, nBinRes, resRangeY, resRangeZ);
  fillHistos(chain, hm, ts, resType, onlyAccepted, maxTracks);
  delete chain;
  postProcessHistos(hm, minSliceStat, minEntriesFit);
  hm.write(outName);
  printf("Saved %d histograms to %s\n", hm.getNHistos(), outName);
}

//____________________________________________________________________________________________
// Drawing

namespace trdtof
{
// the same EColor name parsing as CheckResidSpec::drawHistos
std::optional<int> parseColor(const std::string& token)
{
  static constexpr std::array<std::pair<std::string_view, int>, 15> rootColors{{{"kWhite", kWhite}, {"kBlack", kBlack}, {"kGray", kGray}, {"kRed", kRed}, {"kGreen", kGreen}, {"kBlue", kBlue}, {"kYellow", kYellow}, {"kMagenta", kMagenta}, {"kCyan", kCyan}, {"kOrange", kOrange}, {"kSpring", kSpring}, {"kTeal", kTeal}, {"kAzure", kAzure}, {"kViolet", kViolet}, {"kPink", kPink}}};
  std::string normalized;
  for (unsigned char ch : token) {
    if (!std::isspace(ch)) {
      normalized.push_back(ch);
    }
  }
  for (const auto& [name, value] : rootColors) {
    if (normalized.compare(0, name.size(), name) != 0) {
      continue;
    }
    const auto increment = std::string_view(normalized).substr(name.size());
    if (increment.empty()) {
      return value;
    }
    if (increment.size() < 2 || (increment.front() != '+' && increment.front() != '-')) {
      return std::nullopt;
    }
    int offset = 0;
    const auto* begin = increment.data() + 1;
    const auto* end = increment.data() + increment.size();
    auto [parsedEnd, error] = std::from_chars(begin, end, offset);
    if (error != std::errc{} || parsedEnd != end) {
      return std::nullopt;
    }
    return increment.front() == '-' ? value - offset : value + offset;
  }
  return std::nullopt;
}

/// appends the contents of the bins with an error (i.e. a successful fit) of the histograms with the given ID of all the managers
void collectFitted(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int id, std::vector<double>& vals)
{
  for (const auto& hm : hmv) {
    auto* h = hm ? hm->getHisto(id) : nullptr;
    if (!h) {
      continue;
    }
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      if (h->GetBinError(ib) > 0.) {
        vals.push_back(h->GetBinContent(ib));
      }
    }
  }
}

/// median of v (reordered)
double median(std::vector<double>& v)
{
  const size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + n, v.end());
  double m = v[n];
  if (v.size() % 2 == 0) {
    m = 0.5 * (m + *std::max_element(v.begin(), v.begin() + n));
  }
  return m;
}

/// median and MAD of vals (reordered and overwritten), {0, 0} if empty
std::array<double, 2> medianMAD(std::vector<double>& vals)
{
  if (vals.empty()) {
    return {0., 0.};
  }
  const double med = median(vals);
  for (auto& v : vals) {
    v = std::abs(v - med);
  }
  return {med, median(vals)};
}

/// range of a row of pads (one sector and layer, all stacks/slices and managers) from the bins with a successful fit:
/// mean: +-(|median| + 4 MAD); sigma: [max(0, median - 4 MAD), median + 4 MAD]. Returns {0, 0} if there are no bins.
std::array<float, 2> rowRange(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int isect, int ilr, int nCol, int iyz, bool isSigma)
{
  std::vector<double> vals;
  for (int ic = 0; ic < nCol; ++ic) {
    collectFitted(hmv, histoID(isect, ilr, ic, iyz) + 1 + isSigma, vals);
  }
  if (vals.empty()) {
    return {0.f, 0.f};
  }
  const auto [med, mad] = medianMAD(vals);
  if (isSigma) {
    return {(float)std::max(0., med - 4. * mad), (float)(med + 4. * mad)};
  }
  const float lim = std::abs(med) + 4. * mad;
  return {-lim, lim};
}

/// common range of one type (iyz, isSigma) of the per-sector summaries of the layers [lrMin, lrMax], all columns and managers:
/// mean: +-(|median| + 5 MAD), sigma: [0, |median| + 5 MAD]. Returns {0, 0} if there are no bins.
std::array<float, 2> summaryRange(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int lrMin, int lrMax, int iyz, bool isSigma)
{
  std::vector<double> vals;
  for (int ilr = lrMin; ilr <= lrMax; ++ilr) {
    for (int ic = 0; ic < MaxTOFSlices; ++ic) {
      collectFitted(hmv, summaryID(ilr, ic, iyz) + 1 + isSigma, vals);
    }
  }
  if (vals.empty()) {
    return {0.f, 0.f};
  }
  const auto [med, mad] = medianMAD(vals);
  const float lim = std::abs(med) + 5. * mad;
  return {isSigma ? 0.f : -lim, lim};
}

/// graph of the bins of h with an error (i.e. a successful fit), shifted by dx
TGraphErrors* fittedBinsGraph(const TH1* h, int color, int marker, double dx = 0.)
{
  auto* gr = new TGraphErrors();
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    if (h->GetBinError(ib) > 0.) {
      const int n = gr->GetN();
      gr->SetPoint(n, h->GetBinCenter(ib) + dx, h->GetBinContent(ib));
      gr->SetPointError(n, 0., h->GetBinError(ib));
    }
  }
  gr->SetMarkerStyle(marker);
  gr->SetMarkerSize(1.0);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  return gr;
}

/// frame histogram (drawn empty, its line clipped away) with the binning, the titles and the vertical range [lo, hi]
TH1F* makeFrame(const char* name, const char* title, int nb, double xmin, double xmax, const char* xtit, const char* ytit, float lo, float hi)
{
  auto* frame = new TH1F(name, title, nb, xmin, xmax);
  frame->SetDirectory(nullptr);
  frame->SetStats(false);
  frame->SetXTitle(xtit);
  frame->SetYTitle(ytit);
  frame->SetMinimum(lo);
  frame->SetMaximum(hi);
  for (int ib = 1; ib <= nb; ++ib) {
    frame->SetBinContent(ib, lo - 10. * (hi - lo));
  }
  frame->SetLineColorAlpha(kWhite, 0.);
  return frame;
}

/// legend with one entry per loaded manager (owned by the current pad)
TLegend* makeLegend(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, const std::vector<int>& colors, const std::vector<int>& markers,
                    double x1, double y1, double x2, double y2, int nColumns)
{
  auto* lg = new TLegend(x1, y1, x2, y2);
  lg->SetNColumns(std::max(1, nColumns));
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  for (size_t i = 0; i < hmv.size(); ++i) {
    if (hmv[i]) {
      auto* mk = new TGraph(); // legend marker
      mk->SetMarkerStyle(markers[i]);
      mk->SetMarkerColor(colors[i]);
      mk->SetLineColor(colors[i]);
      lg->AddEntry(mk, hmv[i]->GetName(), "lp")->SetTextColor(colors[i]);
    }
  }
  return lg;
}

/// one page of per-sector summaries of layer ilr (TRD layer or TOF): 4 rows (mean dY, sigma dY, mean dZ, sigma dZ), one column per
/// TRD stack / TOF Z slice, all drawn even if empty. The frames of a row abut and have the same width: only the leftmost pad has a
/// left margin and the Y axis labels and title. ranges: vertical range of each row.
void drawSummaryPage(TCanvas& canv, const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int ilr, int nCol, const std::array<std::array<float, 2>, 4>& ranges,
                     const std::vector<int>& colors, const std::vector<int>& markers)
{
  constexpr double XL = 0.08, XR = 0.01, YT = 0.04, YB = 0.005; // canvas fractions: left (Y axis), right, top (header), bottom
  constexpr double PadTop = 0.09, PadBot = 0.18;                // pad margins (fractions of the pad height)
  constexpr int LabSize = 26, TitSize = 30;                     // pixels, font precision 3: the same on pads of different widths
  canv.Clear();
  canv.cd();
  const double w = (1. - XL - XR) / nCol, h = (1. - YT - YB) / 4.;
  auto* hdr = new TLatex(0.5, 1. - 0.45 * YT, ilr == LrTOF ? "TOF, per sector" : Form("TRD layer %d, per sector", ilr));
  hdr->SetNDC();
  hdr->SetTextFont(63);
  hdr->SetTextSize(40);
  hdr->SetTextAlign(22);
  hdr->Draw();
  canv.cd();
  makeLegend(hmv, colors, markers, 0.70, 1. - 0.9 * YT, 1. - XR, 1. - 0.1 * YT, hmv.size())->Draw();
  for (int ir = 0; ir < 4; ++ir) {
    const int iyz = ir / 2, isSigma = ir % 2;
    const auto [lo, hi] = ranges[ir];
    for (int ic = 0; ic < nCol; ++ic) {
      const double x1 = ic == 0 ? 0. : XL + ic * w, x2 = XL + (ic + 1) * w + (ic == nCol - 1 ? XR : 0.);
      const double y2 = 1. - YT - ir * h, y1 = y2 - h;
      canv.cd();
      auto* pad = new TPad(Form("psum_%d_%d", ir, ic), "", x1, y1, x2, y2);
      pad->SetLeftMargin(ic == 0 ? XL / (x2 - x1) : 0.);
      pad->SetRightMargin(ic == nCol - 1 ? XR / (x2 - x1) : 0.);
      pad->SetTopMargin(PadTop);
      pad->SetBottomMargin(PadBot);
      pad->SetFillStyle(4000);
      pad->Draw();
      pad->cd();
      const int id = summaryID(ilr, ic, iyz) + 1 + isSigma;
      const TH1* href = nullptr;
      for (const auto& hm : hmv) {
        if (!href && hm) {
          href = hm->getHisto(id);
        }
      }
      const std::string title = href ? href->GetTitle() : (ilr == LrTOF ? Form("TOF z%d", ic) : Form("TRD L%d st%d", ilr, ic));
      const std::string ytit = href ? href->GetYaxis()->GetTitle() : (isSigma ? (iyz ? "#sigma(dZ) (cm)" : "#sigma(dY) (cm)") : (iyz ? "<dZ> (cm)" : "<dY> (cm)"));
      auto* frame = makeFrame(Form("fsum_%d_%d_%d", ilr, ir, ic), "", NSect, -0.5, NSect - 0.5, "sector", ytit.c_str(), lo, hi);
      for (auto* ax : {frame->GetXaxis(), frame->GetYaxis()}) {
        ax->SetLabelFont(43);
        ax->SetLabelSize(LabSize);
        ax->SetTitleFont(43);
        ax->SetTitleSize(TitSize);
        ax->CenterTitle();
      }
      frame->GetXaxis()->SetTitleOffset(1.8);
      frame->GetYaxis()->SetTitleOffset(2.8);
      for (int isect = 0; isect < NSect; ++isect) { // integer sector labels at the bin centres, every 3rd one
        frame->GetXaxis()->SetBinLabel(isect + 1, isect % 3 ? "" : Form("%d", isect));
      }
      frame->GetXaxis()->LabelsOption("h");
      if (ic) { // the Y axis is labelled only on the leftmost pad
        frame->GetYaxis()->SetLabelSize(0);
        frame->GetYaxis()->SetTitleSize(0);
      }
      frame->Draw();
      for (size_t i = 0; i < hmv.size(); ++i) {
        auto* hs = hmv[i] ? hmv[i]->getHisto(id) : nullptr;
        if (hs) {
          auto* gr = fittedBinsGraph(hs, colors[i], markers[i], 0.25 * (i - 0.5 * (hmv.size() - 1)) / std::max<size_t>(1, hmv.size() - 1)); // managers side by side
          if (gr->GetN()) {
            gr->Draw("P");
          }
        }
      }
      pad->SetGrid();
      const double lm = pad->GetLeftMargin(), rm = pad->GetRightMargin();
      auto* tl = new TLatex(0.5 * (lm + 1. - rm), 1. - 0.5 * PadTop, title.c_str());
      tl->SetNDC();
      tl->SetTextFont(43);
      tl->SetTextSize(TitSize);
      tl->SetTextAlign(22);
      tl->Draw();
    }
  }
  canv.cd();
}

/// draws, superimposed, the histograms with the given ID of all the managers on the current pad, skipping the bins w/o a successful
/// fit; yRange (lo < hi) is the vertical range, otherwise it is set from the content.
/// Returns the graph of each manager actually drawn (the pad owns them).
void drawOnPad(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int id, bool isSigma, std::array<float, 2> yRange, const std::vector<int>& colors, const std::vector<int>& markers, std::vector<TGraphErrors*>& drawn)
{
  drawn.assign(hmv.size(), nullptr);
  float vmin = 1e9, vmax = -1e9;
  const TH1* href = nullptr;
  for (size_t i = 0; i < hmv.size(); ++i) {
    auto* h = hmv[i] ? hmv[i]->getHisto(id) : nullptr;
    if (!h) {
      continue;
    }
    drawn[i] = fittedBinsGraph(h, colors[i], markers[i]);
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      if (h->GetBinError(ib) > 0.) {
        vmin = std::min(vmin, (float)(h->GetBinContent(ib) - h->GetBinError(ib)));
        vmax = std::max(vmax, (float)(h->GetBinContent(ib) + h->GetBinError(ib)));
      }
    }
    if (!href) {
      href = h;
    }
  }
  if (!href) {
    return;
  }
  float lo = vmin, hi = vmax;
  if (yRange[1] > yRange[0]) {
    lo = yRange[0];
    hi = yRange[1];
  } else if (vmax > vmin) {
    const float mrg = 0.15f * (vmax - vmin);
    lo = isSigma ? 0.f : vmin - mrg;
    hi = vmax + mrg;
  } else {
    lo = isSigma ? 0.f : -1.f;
    hi = 1.f;
  }
  auto* frame = makeFrame(Form("%s_frame", href->GetName()), href->GetTitle(), href->GetNbinsX(), href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax(),
                          href->GetXaxis()->GetTitle(), href->GetYaxis()->GetTitle(), lo, hi);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.03);
  for (auto* ax : {frame->GetXaxis(), frame->GetYaxis()}) {
    ax->SetTitleSize(0.065);
    ax->SetLabelSize(0.055);
  }
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->Draw();
  for (auto* gr : drawn) {
    if (gr && gr->GetN()) {
      gr->Draw("P");
    }
  }
  gPad->SetGrid();
}
} // namespace trdtof

/// Superimpose the HistoManagers listed in the CheckResidConfig keys ext_hm_list, ext_leg_list, ext_color_list, ext_marker_list
/// (same format as in the --configKeyValues of o2-check-resid-workfow --draw-external-only, e.g. in doc/reResidAll.sh), written
/// by trdTofResid, into pdfName: one page per sector and TRD layer / TOF, 4 rows (mean dY, sigma dY, mean dZ, sigma dZ) and
/// one column per TRD stack / TOF Z slice, preceded by one page of per-sector summaries per TRD layer and one for TOF, in the same
/// layout, with abutting pads and all the pads drawn even if empty. Managers written without the summaries get them computed here.
/// cfgKeyValues : "checkresid.ext_hm_list=a.root,b.root;checkresid.ext_leg_list=...;..." (empty: the keys are taken as they are)
/// rerunPostProcess : refit/reconvert the 2D histograms of the loaded managers (minSliceStat, minEntriesFit as in trdTofResid)
/// maxMeanY, maxSigY, maxMeanZ, maxSigZ : fixed vertical ranges (cm), 0 = automatic: per row of a page (all stacks and managers)
///   median and MAD of the plotted bins, mean: +-(|median| + 4 MAD), sigma: [max(0, median - 4 MAD), median + 4 MAD];
///   for the summary pages the range is common to all the summaries of one type of TRD (resp. TOF), from all its layers, columns
///   and managers: mean: +-(|median| + 5 MAD), sigma: [0, |median| + 5 MAD]
void drawTRDTOFResid(const char* cfgKeyValues = "", const char* pdfName = "trdTofResid.pdf", bool rerunPostProcess = false, float maxMeanY = 0.f,
                     float maxSigY = 0.f, float maxMeanZ = 0.f, float maxSigZ = 0.f, int sectMin = 0, int sectMax = 17, int minSliceStat = 3,
                     int minEntriesFit = 100)
{
  using namespace trdtof;
  if (cfgKeyValues && cfgKeyValues[0]) {
    o2::conf::ConfigurableParam::updateFromString(cfgKeyValues);
  }
  const auto& params = o2::checkresid::CheckResidConfig::Instance();
  const auto hmNames = o2::utils::Str::tokenize(params.ext_hm_list, ",;");
  const auto vecLegends = o2::utils::Str::tokenize(params.ext_leg_list, ',');
  const auto colorNames = o2::utils::Str::tokenize(params.ext_color_list, ",;");
  const auto markerIDs = o2::utils::Str::tokenize(params.ext_marker_list, ",;");
  if (hmNames.empty()) {
    printf("checkresid.ext_hm_list is empty, nothing to draw\n");
    return;
  }
  const int nhm = hmNames.size();
  bool useLeg = vecLegends.size() == hmNames.size();
  bool useColors = colorNames.size() == hmNames.size() && !colorNames.empty();
  bool useMarkers = markerIDs.size() == hmNames.size() && !markerIDs.empty();
  if (!useLeg) {
    printf("Warning: %zu legend names provided for %d external histomanagers, will use file names as legends\n", vecLegends.size(), nhm);
  }
  if (!useColors) {
    printf("Warning: %zu color names provided for %d external histomanagers, will use automatic colors\n", colorNames.size(), nhm);
  }
  if (!useMarkers) {
    printf("Warning: %zu marker IDs provided for %d external histomanagers, will use automatic markers\n", markerIDs.size(), nhm);
  }
  std::vector<int> colors, markers;
  if (useColors) { // an unrecognized name falls back to the automatic color of that entry only
    for (size_t i = 0; i < colorNames.size(); ++i) {
      auto c = parseColor(colorNames[i]);
      if (!c) {
        printf("Warning: unrecognized EColor name '%s', will use the automatic color for entry %zu\n", colorNames[i].c_str(), i);
      }
      colors.push_back(c ? *c : -1);
    }
  }
  if (useMarkers) {
    for (const auto& mk : markerIDs) {
      int marker = 0;
      auto [parsedEnd, error] = std::from_chars(mk.data(), mk.data() + mk.size(), marker);
      if (error != std::errc{} || parsedEnd != mk.data() + mk.size()) {
        printf("Warning: unrecognized marker ID '%s', will use automatic markers\n", mk.c_str());
        useMarkers = false;
        break;
      }
      markers.push_back(marker);
    }
  }
  constexpr std::array<unsigned int, 3> hcol{kRed, kBlue, kGreen + 2};
  std::vector<std::unique_ptr<o2::HistoManager>> hmv;
  std::vector<int> hmColor(nhm), hmMarker(nhm);
  for (int i = 0; i < nhm; ++i) {
    printf("Loading external HistoManager %s\n", hmNames[i].c_str());
    auto& hm = hmv.emplace_back(std::make_unique<o2::HistoManager>("", hmNames[i], true));
    if (hm->GetLast() < 1) {
      printf("Failed to load histograms from %s\n", hmNames[i].c_str());
      hm.reset();
      continue;
    }
    if (rerunPostProcess) {
      postProcessHistos(*hm, minSliceStat, minEntriesFit);
    } else if (!hm->getHisto(summaryID(0, 0, 0) + 1)) { // written before the per-sector summaries existed
      makeSectorSummaries(*hm, minEntriesFit);
    }
    hm->SetName(useLeg ? vecLegends[i].c_str() : hmNames[i].c_str());
    hmColor[i] = useColors && colors[i] >= 0 ? colors[i] : hcol[i % hcol.size()];
    hmMarker[i] = useMarkers ? markers[i] : 20 + i + (i % 2) * 4;
    hm->setMarkerStyle(hmMarker[i], 0.6);
    hm->setColor(hmColor[i]);
  }

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  // pages: sector and layer with at least one histogram, nCol = number of columns (TRD stacks or TOF Z slices) to draw
  struct Page {
    int sect, lr, nCol;
  };
  std::vector<Page> pages;
  for (int isect = sectMin; isect <= sectMax; ++isect) {
    for (int ilr = 0; ilr <= LrTOF; ++ilr) {
      int nCol = 0;
      for (int ist = 0; ist < (ilr == LrTOF ? MaxTOFSlices : NTRDSt); ++ist) {
        for (const auto& hm : hmv) {
          if (hm && hm->getHisto(histoID(isect, ilr, ist, 0) + 1)) {
            nCol = ist + 1;
          }
        }
      }
      if (nCol) {
        nCol = ilr == LrTOF ? nCol : NTRDSt; // all the stacks of a TRD layer
        pages.push_back({isect, ilr, nCol});
      }
    }
  }
  // per-sector summary pages: every TRD layer, TOF if any manager has its summaries; nCol as above
  std::vector<std::array<int, 2>> sumPages;
  for (int ilr = 0; ilr <= LrTOF; ++ilr) {
    int nCol = 0;
    for (int ic = 0; ic < (ilr == LrTOF ? MaxTOFSlices : NTRDSt); ++ic) {
      for (const auto& hm : hmv) {
        if (hm && hm->getHisto(summaryID(ilr, ic, 0) + 1)) {
          nCol = ic + 1;
        }
      }
    }
    if (nCol) {
      sumPages.push_back({ilr, ilr == LrTOF ? nCol : NTRDSt});
    }
  }
  const int nPages = sumPages.size() + pages.size();
  if (!nPages) {
    printf("No histograms to draw\n");
    return;
  }
  TCanvas canv("ctrdtof", "", 2000, 2600); // same aspect ratio as the ROOT pdf page (20x26), else the content is misplaced
  canv.SetCanvasSize(2000, 2600);              // the exact size, w/o it the pdf content is shifted
  int iPage = 0;
  auto printPage = [&]() {
    const char* flag = iPage == 0 ? "(" : (iPage + 1 == nPages ? ")" : "");
    canv.Print(nPages == 1 ? pdfName : Form("%s%s", pdfName, flag));
    iPage++;
  };
  const float rowMax[4] = {maxMeanY, maxSigY, maxMeanZ, maxSigZ};
  // the same vertical range for all the summaries of one type of a detector (TRD layers / TOF)
  std::array<std::array<float, 2>, 4> sumRangeTRD{}, sumRangeTOF{};
  for (int ir = 0; ir < 4; ++ir) {
    const int iyz = ir / 2, isSigma = ir % 2;
    sumRangeTRD[ir] = summaryRange(hmv, 0, NTRDLr - 1, iyz, isSigma);
    sumRangeTOF[ir] = summaryRange(hmv, LrTOF, LrTOF, iyz, isSigma);
    for (auto* rng : {&sumRangeTRD[ir], &sumRangeTOF[ir]}) {
      if (rowMax[ir] > 0.f) { // explicit fixed range overrides the median/MAD one
        *rng = {isSigma ? 0.f : -rowMax[ir], rowMax[ir]};
      } else if ((*rng)[1] <= (*rng)[0]) { // no fitted bin at all
        *rng = {isSigma ? 0.f : -1.f, 1.f};
      }
    }
  }
  for (const auto& [ilr, nCol] : sumPages) {
    drawSummaryPage(canv, hmv, ilr, nCol, ilr == LrTOF ? sumRangeTOF : sumRangeTRD, hmColor, hmMarker);
    printPage();
  }
  for (size_t ip = 0; ip < pages.size(); ++ip) {
    const auto [isect, ilr, nCol] = pages[ip];
    canv.Clear();
    canv.Divide(nCol, 4, 0.001, 0.001);
    bool legDrawn = false;
    for (int ir = 0; ir < 4; ++ir) {
      const int iyz = ir / 2, isSigma = ir % 2;
      std::array<float, 2> yRange = rowRange(hmv, isect, ilr, nCol, iyz, isSigma);
      if (rowMax[ir] > 0.f) { // explicit fixed range overrides the median/MAD one
        yRange = {isSigma ? 0.f : -rowMax[ir], rowMax[ir]};
      }
      for (int ic = 0; ic < nCol; ++ic) {
        canv.cd(ir * nCol + ic + 1);
        std::vector<TGraphErrors*> drawn;
        drawOnPad(hmv, histoID(isect, ilr, ic, iyz) + 1 + isSigma, isSigma, yRange, hmColor, hmMarker, drawn);
        if (legDrawn || std::none_of(drawn.begin(), drawn.end(), [](const TGraphErrors* g) { return g != nullptr; })) {
          continue;
        }
        auto* lg = makeLegend(hmv, hmColor, hmMarker, 0.12, 0.13, 0.9, 0.13 + std::min(0.5f, nhm * 0.2f / 3.f), std::round(nhm / 2));
        lg->Draw();
        legDrawn = true;
      }
    }
    canv.cd(0);
    printPage();
  }
  printf("Wrote %d pages to %s\n", nPages, pdfName);
}

/// Called by ROOT when the macro is given as `root trdTofResid.C+ -e '...'` (the function named as the file is executed): only loads the macro.
void trdTofResid()
{
  printf("trdTofResid(inp, out, geom, resType, ...) fills and fits, drawTRDTOFResid(cfgKeyValues, pdf, ...) draws, see the header of trdTofResid.C\n");
}
