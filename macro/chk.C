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

#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <fstream>
#include <sstream>
#include <map>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include <functional>

using namespace o2::base;
using M = TGeoHMatrix;
static const char* DD = "/data0/devAlg/tests/sol1/savITSBarFix/";

struct Diff {
  double dt{0}, dr{0};
};
Diff cmp(const M& a, const M& b)
{
  Diff d;
  for (int i = 0; i < 3; i++)
    d.dt = std::max(d.dt, std::abs(a.GetTranslation()[i] - b.GetTranslation()[i]));
  for (int i = 0; i < 9; i++)
    d.dr = std::max(d.dr, std::abs(a.GetRotationMatrix()[i] - b.GetRotationMatrix()[i]));
  return d;
}
std::string parentOf(const std::string& s)
{
  auto p = s.find_last_of('/');
  return p == std::string::npos ? "" : s.substr(0, p);
}

void chk(int nChipsMax = 24120)
{
  // 1. parameters: label -> value from millepede.res, label -> (name, dof) from the labeled copy
  std::map<uint32_t, double> val;
  {
    std::ifstream f(std::string(DD) + "millepede.res");
    std::string l;
    std::getline(f, l);
    while (std::getline(f, l)) {
      std::istringstream is(l);
      uint32_t lb;
      double v;
      if (is >> lb >> v)
        val[lb] = v;
    }
  }
  const char* dofn[6] = {"TX", "TY", "TZ", "RX", "RY", "RZ"};
  std::map<std::string, std::array<double, 6>> par; // volume -> TX..RZ
  int nMis = 0, nLab = 0;
  {
    std::ifstream f(std::string(DD) + "millepede_labeled.res");
    std::string l;
    std::getline(f, l);
    while (std::getline(f, l)) {
      auto p = l.find('!');
      if (p == std::string::npos)
        continue;
      std::istringstream is(l.substr(0, p));
      uint32_t lb;
      double v;
      is >> lb >> v;
      std::istringstream ic(l.substr(p + 1));
      std::string flag, dof, name;
      ic >> flag >> dof >> name;
      int idof = -1;
      for (int i = 0; i < 6; i++)
        if (dof == dofn[i])
          idof = i;
      if (idof < 0 || name.rfind("ITS", 0) != 0)
        continue;
      nLab++;
      if (!val.count(lb) || val[lb] != v)
        nMis++;
      par[name][idof] = val.count(lb) ? val[lb] : 0.; // trust millepede.res
    }
  }
  printf("ITS rigid-body labels in labeled file: %d, mismatching millepede.res: %d, volumes with params: %zu\n", nLab, nMis, par.size());

  // 2. initial geometry
  GeometryManager::loadGeometry("/data0/ALG25/Pede_v10n/GLO/Config/GeometryAligned/snapshot.root");
  std::vector<std::string> chipName;
  std::vector<M> chipIni;
  for (int i = 0; i < nChipsMax; i++) {
    std::string nm = GeometryManager::getSymbolicName(o2::detectors::DetID::ITS, i);
    if (nm.empty())
      break;
    chipName.push_back(nm);
    chipIni.push_back(*GeometryManager::getMatrix(o2::detectors::DetID::ITS, i));
  }
  printf("chips: %zu\n", chipName.size());
  std::map<std::string, M> volIni; // all volumes with params and all their ancestors
  std::map<std::string, bool> isChip;
  for (auto& n : chipName)
    isChip[n] = true;
  auto fetch = [&](const std::string& n) {
    if (!volIni.count(n)) {
      auto* m = GeometryManager::getMatrix(n.c_str());
      volIni[n] = m ? *m : M();
      if (!m)
        printf("no matrix for %s\n", n.c_str());
    }
  };
  for (auto& [n, p] : par) {
    for (std::string s = n; !s.empty(); s = parentOf(s))
      fetch(s);
  }
  // chip matrix of the volume vs fitting frame
  const double eps = o2::itsmft::SegmentationAlpide::SensorLayerThickness - o2::itsmft::SegmentationAlpide::SensorLayerThicknessEff;
  TGeoTranslation tra(0., 0.5 * eps, 0.);

  // 3. cumulative predicted global delta P(vol) = P(parent) * L2G * d * L2G^-1
  std::map<std::string, M> cum;
  std::function<M(const std::string&)> getCum = [&](const std::string& n) -> M {
    if (n.empty())
      return M();
    if (auto it = cum.find(n); it != cum.end())
      return it->second;
    M pn = getCum(parentOf(n));
    if (auto it = par.find(n); it != par.end()) {
      auto& v = it->second;
      bool nz = false;
      for (double x : v)
        nz |= x != 0.;
      if (nz) {
        M d = o2::detectors::AlignParam("", -1, v[0], v[1], v[2], v[3], v[4], v[5], true, false).createMatrix();
        M l2g;
        if (isChip.count(n)) {
          l2g = volIni.count(n) ? volIni[n] : M();
          l2g *= tra;
        } else
          l2g = volIni[n];
        pn = pn * l2g * d * l2g.Inverse();
      }
    }
    return cum[n] = pn;
  };
  // pre-fetch chips needing params
  std::vector<M> chipPred(chipName.size());
  for (size_t i = 0; i < chipName.size(); i++) {
    if (par.count(chipName[i]))
      volIni[chipName[i]] = chipIni[i];
  }
  for (size_t i = 0; i < chipName.size(); i++) {
    // ancestors of chips without own params may not be in volIni (only those with params were fetched), fine: they have no params
    chipPred[i] = getCum(chipName[i]) * chipIni[i];
  }
  std::vector<M> volPred;
  std::vector<std::string> volNames;
  for (auto& [n, p] : par)
    if (!isChip.count(n)) {
      volNames.push_back(n);
      volPred.push_back(getCum(n) * volIni[n]);
    }
  std::vector<M> volIniV;
  for (auto& n : volNames)
    volIniV.push_back(volIni[n]);

  // 4. corrected geometry
  GeometryManager::loadGeometry((std::string(DD) + "o2sim_geometry-aligned.root").c_str());
  double mxT = 0, mxR = 0, mxT0 = 0, mxR0 = 0, sT = 0;
  int nBad = 0, nMoved = 0;
  std::string worst;
  for (size_t i = 0; i < chipName.size(); i++) {
    M nw = *GeometryManager::getMatrix(o2::detectors::DetID::ITS, i);
    auto d = cmp(nw, chipPred[i]), d0 = cmp(nw, chipIni[i]);
    if (d.dt > mxT) {
      mxT = d.dt;
      worst = chipName[i];
    }
    mxR = std::max(mxR, d.dr);
    mxT0 = std::max(mxT0, d0.dt);
    mxR0 = std::max(mxR0, d0.dr);
    sT += d.dt;
    if (d0.dt > 1e-9 || d0.dr > 1e-9)
      nMoved++;
    if (d.dt > 1e-6 || d.dr > 1e-8)
      nBad++;
  }
  printf("CHIPS: moved by correction: %d / %zu\n", nMoved, chipName.size());
  printf("  |new - ini|  max: dT=%.3e cm  dR=%.3e   (size of the correction)\n", mxT0, mxR0);
  printf("  |new - pred| max: dT=%.3e cm  dR=%.3e   mean dT=%.3e, worst %s, #bad(>1e-6cm/1e-8)=%d\n", mxT, mxR, sT / chipName.size(), worst.c_str(), nBad);
  double vT = 0, vR = 0, vT0 = 0;
  std::string vw;
  for (size_t i = 0; i < volNames.size(); i++) {
    M nw = *GeometryManager::getMatrix(volNames[i].c_str());
    auto d = cmp(nw, volPred[i]), d0 = cmp(nw, volIniV[i]);
    if (d.dt > vT) {
      vT = d.dt;
      vw = volNames[i];
    }
    vR = std::max(vR, d.dr);
    vT0 = std::max(vT0, d0.dt);
  }
  printf("VOLUMES with params (staves/halfstaves/...): %zu, max|new-ini| dT=%.3e, max|new-pred| dT=%.3e dR=%.3e (worst %s)\n", volNames.size(), vT0, vT, vR, vw.c_str());
}
