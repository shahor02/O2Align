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

// Visualisation of the shift of the ITS half-barrel shells (ITS/ITSULayer<l>/ITSUHalfBarrel<h>) between two geometries.
//
//   root -l -b -q 'macro/drawHBShift.C+("/data/devAlg/best/IBEdge/o2sim_geometry-aligned.root", \
//                                       "/data/devAlg/solPbHBA0/o2sim_geometry-aligned.root", 500, "hbShift.png")'
//
// One pad per layer-half-barrel (row = half-barrel 0/1, column = layer). In every pad:
//  - solid, shaded half-cylinder: the half-barrel shell in the reference geometry (file 1);
//  - outline: the same shell in the corrected geometry (file 2), where the displacement of every point of the shell
//    (new - old, lab frame) is multiplied by `magnification`;
//  - lab x,y,z axes from the origin, ticked in cm (the scale of the shell), with the tick step given in the pad.
// The shell (radius, azimuth span, Z extent) is derived from the chips of the half-barrel in the reference geometry; the
// shift is the exact difference of the half-barrel TGeo matrices, i.e. all the rigid-body motion of the HB is included.
//
// View: by default lab X to the right, lab Y up, lab Z towards the viewer (half-barrels in their true lab orientation).
// rotX, rotY, rotZ (degrees) rotate the scene about these fixed screen axes, in this order.

#include <TArrow.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TPolyLine.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace hbshift
{
using V3 = std::array<double, 3>;
constexpr int NLayers = 7, NHB = 2;

struct Shell {
  TGeoHMatrix mOld, mNew; // HB local-to-lab in the reference / corrected geometry
  double radius = 0., phiMid = 0., phiHalf = 0., zMin = 1e9, zMax = -1e9; // HB local frame
  int nChips = 0;
  double sumR = 0.;
  double sumSin = 0., sumCos = 0.;
  std::vector<double> chipPhi; // corners' phi, to find the span
};

std::string hbName(int l, int h) { return "ITS/ITSULayer" + std::to_string(l) + "/ITSUHalfBarrel" + std::to_string(h); }

bool cdTo(const std::string& name, TGeoHMatrix& m)
{
  auto* e = gGeoManager->GetAlignableEntry(name.c_str());
  if (!e || !gGeoManager->cd(e->GetTitle())) {
    return false;
  }
  m = *gGeoManager->GetCurrentMatrix();
  return true;
}

// HB matrices, plus (if collectShape) the shell shape from the chips of the geometry
bool readGeometry(const char* file, std::map<int, Shell>& shells, bool second)
{
  if (!TGeoManager::Import(file)) {
    std::cerr << "Cannot load geometry " << file << '\n';
    return false;
  }
  for (int l = 0; l < NLayers; l++) {
    for (int h = 0; h < NHB; h++) {
      auto& s = shells[l * NHB + h];
      if (!cdTo(hbName(l, h), second ? s.mNew : s.mOld)) {
        std::cerr << "No alignable " << hbName(l, h) << " in " << file << '\n';
        return false;
      }
    }
  }
  if (!second) { // shell shape from the chips
    for (int i = 0; i < gGeoManager->GetNAlignable(); i++) {
      auto* e = gGeoManager->GetAlignableEntry(i);
      const std::string name = e ? e->GetName() : "";
      if (name.find("ITSUChip") == std::string::npos || name.compare(0, 13, "ITS/ITSULayer") != 0) {
        continue;
      }
      const int l = name[13] - '0';
      const auto hp = name.find("ITSUHalfBarrel");
      if (l < 0 || l >= NLayers || hp == std::string::npos || !gGeoManager->cd(e->GetTitle())) {
        continue;
      }
      auto& s = shells[l * NHB + (name[hp + 14] - '0')];
      auto* box = dynamic_cast<TGeoBBox*>(gGeoManager->GetCurrentVolume()->GetShape());
      const TGeoHMatrix* c2l = gGeoManager->GetCurrentMatrix();
      const double zero[3] = {0, 0, 0};
      double g[3], loc[3];
      c2l->LocalToMaster(zero, g);
      s.mOld.MasterToLocal(g, loc);
      s.sumR += std::hypot(loc[0], loc[1]);
      s.nChips++;
      for (int c = 0; c < 8; c++) { // corners, for the extents
        const double d[3] = {(c & 1 ? 1. : -1.) * (box ? box->GetDX() : 0.), (c & 2 ? 1. : -1.) * (box ? box->GetDY() : 0.), (c & 4 ? 1. : -1.) * (box ? box->GetDZ() : 0.)};
        c2l->LocalToMaster(d, g);
        s.mOld.MasterToLocal(g, loc);
        s.zMin = std::min(s.zMin, loc[2]);
        s.zMax = std::max(s.zMax, loc[2]);
        s.chipPhi.push_back(std::atan2(loc[1], loc[0]));
        s.sumSin += std::sin(s.chipPhi.back());
        s.sumCos += std::cos(s.chipPhi.back());
      }
    }
    for (auto& [k, s] : shells) {
      if (!s.nChips) {
        std::cerr << "No chips for shell " << k << '\n';
        return false;
      }
      s.radius = s.sumR / s.nChips;
      s.phiMid = std::atan2(s.sumSin, s.sumCos);
      for (double p : s.chipPhi) {
        s.phiHalf = std::max(s.phiHalf, std::abs(TVector2::Phi_mpi_pi(p - s.phiMid)));
      }
    }
  }
  delete gGeoManager;
  gGeoManager = nullptr;
  return true;
}

// ---- projection
struct View {
  V3 d, r, u; // screen axes in the lab frame: to the viewer (depth), right, up
  // The scene is rotated about the fixed screen axes: first by rx about screen X, then ry about Y, then rz about Z (degrees).
  // Zero angles: lab X to the right, lab Y up, lab Z towards the viewer.
  void init(double rxDeg, double ryDeg, double rzDeg)
  {
    const double ax = rxDeg * TMath::DegToRad(), ay = ryDeg * TMath::DegToRad(), az = rzDeg * TMath::DegToRad();
    const double cx = std::cos(ax), sx = std::sin(ax), cy = std::cos(ay), sy = std::sin(ay), cz = std::cos(az), sz = std::sin(az);
    // R = Rz Ry Rx; rows of R are the screen axes expressed in the lab frame
    const double R[3][3] = {{cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx},
                            {sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx},
                            {-sy, cy * sx, cy * cx}};
    r = {R[0][0], R[0][1], R[0][2]};
    u = {R[1][0], R[1][1], R[1][2]};
    d = {R[2][0], R[2][1], R[2][2]};
  }
  static V3 cross(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
  static double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  static V3 norm(V3 a)
  {
    const double n = std::sqrt(dot(a, a));
    for (auto& x : a) {
      x /= n;
    }
    return a;
  }
  void xy(const V3& p, double& x, double& y) const
  {
    x = dot(p, r);
    y = dot(p, u);
  }
};

V3 toLab(const TGeoHMatrix& m, double x, double y, double z)
{
  const double l[3] = {x, y, z};
  double g[3];
  m.LocalToMaster(l, g);
  return {g[0], g[1], g[2]};
}

double niceStep(double range)
{
  const double raw = range / 4., p = std::pow(10., std::floor(std::log10(raw))), f = raw / p;
  return p * (f < 1.5 ? 1. : f < 3.5 ? 2. : f < 7.5 ? 5. : 10.);
}

int shadeColor(int idx, double frac, int base) // base 0: solid shell colour family
{
  static const double rgb[2][3] = {{0.20, 0.45, 0.80}, {0.85, 0.50, 0.15}}; // HB0 blue, HB1 orange
  const double f = 0.65 + 0.35 * frac, white = 0.40; // lighter fill
  auto* col = new TColor(idx, white + (1 - white) * rgb[base][0] * f, white + (1 - white) * rgb[base][1] * f, white + (1 - white) * rgb[base][2] * f);
  return idx;
}

struct Quad {
  double depth, shade;
  V3 p[4];
};

void drawPad(TPad* pad, int l, int h, const Shell& s, double mag, double rotX, double rotY, double rotZ, int nPhi)
{
  pad->cd();
  pad->SetMargin(0., 0., 0., 0.); // the drawing range covers the whole pad: nothing is clipped by the frame
  View v;
  v.init(rotX, rotY, rotZ);

  auto point = [&](const TGeoHMatrix& m, double phi, double z) { return toLab(m, s.radius * std::cos(phi), s.radius * std::sin(phi), z); };
  // displaced point: old + mag * (new - old)
  auto shifted = [&](double phi, double z) {
    const V3 o = point(s.mOld, phi, z), n = point(s.mNew, phi, z);
    return V3{o[0] + mag * (n[0] - o[0]), o[1] + mag * (n[1] - o[1]), o[2] + mag * (n[2] - o[2])};
  };
  const double ph0 = s.phiMid - s.phiHalf, dph = 2. * s.phiHalf / nPhi;

  // shift statistics (true, not magnified), um
  V3 mean{0, 0, 0};
  double maxd = 0.;
  int np = 0;
  for (int i = 0; i <= nPhi; i++) {
    for (double z : {s.zMin, s.zMax}) {
      const V3 o = point(s.mOld, ph0 + i * dph, z), n = point(s.mNew, ph0 + i * dph, z);
      double d2 = 0;
      for (int k = 0; k < 3; k++) {
        mean[k] += n[k] - o[k];
        d2 += (n[k] - o[k]) * (n[k] - o[k]);
      }
      maxd = std::max(maxd, std::sqrt(d2));
      np++;
    }
  }

  // solid reference shell: shaded quads, painter's order
  std::vector<Quad> quads;
  for (int i = 0; i < nPhi; i++) {
    Quad q;
    q.p[0] = point(s.mOld, ph0 + i * dph, s.zMin);
    q.p[1] = point(s.mOld, ph0 + (i + 1) * dph, s.zMin);
    q.p[2] = point(s.mOld, ph0 + (i + 1) * dph, s.zMax);
    q.p[3] = point(s.mOld, ph0 + i * dph, s.zMax);
    V3 c{0, 0, 0};
    for (auto& p : q.p) {
      for (int k = 0; k < 3; k++) {
        c[k] += 0.25 * p[k];
      }
    }
    const V3 nrm = View::norm(View::cross({q.p[1][0] - q.p[0][0], q.p[1][1] - q.p[0][1], q.p[1][2] - q.p[0][2]}, {q.p[3][0] - q.p[0][0], q.p[3][1] - q.p[0][1], q.p[3][2] - q.p[0][2]}));
    q.depth = View::dot(c, v.d);
    q.shade = std::abs(View::dot(nrm, View::norm({v.d[0] + 0.5 * v.r[0] + 0.5 * v.u[0], v.d[1] + 0.5 * v.r[1] + 0.5 * v.u[1], v.d[2] + 0.5 * v.r[2] + 0.5 * v.u[2]})));
    quads.push_back(q);
  }
  std::sort(quads.begin(), quads.end(), [](const Quad& a, const Quad& b) { return a.depth < b.depth; });

  // extents of everything drawn, in screen coordinates
  double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
  auto grow = [&](const V3& p) {
    double x, y;
    v.xy(p, x, y);
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  };
  for (auto& q : quads) {
    for (auto& p : q.p) {
      grow(p);
    }
  }
  for (int i = 0; i <= nPhi; i++) {
    for (double z : {s.zMin, s.zMax}) {
      grow(shifted(ph0 + i * dph, z));
    }
  }
  const double axLenXY = 1.15 * s.radius, zlo = std::min(0., s.zMin), zhi = std::max(0., s.zMax);
  const double axLenZ = 1.10 * std::max(std::abs(zlo), std::abs(zhi));
  const V3 axEnd[3] = {{axLenXY, 0, 0}, {0, axLenXY, 0}, {0, 0, axLenZ}};
  for (auto& a : axEnd) {
    grow(a);
  }
  const double cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax);
  double hw = 0.5 * (xmax - xmin), hh = 0.5 * (ymax - ymin);
  const double aspect = double(pad->GetWw() * pad->GetAbsWNDC()) / std::max(1., double(pad->GetWh() * pad->GetAbsHNDC()));
  hw = std::max(hw, hh * aspect) * 1.12;
  hh = hw / aspect;
  pad->Range(cx - hw, cy - 1.0 * hh, cx + hw, cy + 1.35 * hh); // headroom for the text
  const double span = hw;                                         // screen length unit for ticks

  for (auto& q : quads) {
    double x[5], y[5];
    for (int k = 0; k < 4; k++) {
      v.xy(q.p[k], x[k], y[k]);
    }
    x[4] = x[0];
    y[4] = y[0];
    static int cidx = 1600;
    const int col = shadeColor(cidx++, q.shade, h);
    auto* pl = new TPolyLine(5, x, y);
    pl->SetFillColor(col);
    pl->SetLineColor(col);
    pl->Draw("f");
  }

  auto polyline = [&](const std::vector<V3>& pts, int col, int width, int style) {
    std::vector<double> x(pts.size()), y(pts.size());
    for (size_t k = 0; k < pts.size(); k++) {
      v.xy(pts[k], x[k], y[k]);
    }
    auto* pl = new TPolyLine(pts.size(), x.data(), y.data());
    pl->SetLineColor(col);
    pl->SetLineWidth(width);
    pl->SetLineStyle(style);
    pl->Draw();
  };
  // reference shell edges, dashed, on top of the fill
  for (double z : {s.zMin, s.zMax}) {
    std::vector<V3> arc;
    for (int i = 0; i <= nPhi; i++) {
      arc.push_back(point(s.mOld, ph0 + i * dph, z));
    }
    polyline(arc, kGray + 2, 1, 1);
  }
  for (int i : {0, nPhi}) {
    polyline({point(s.mOld, ph0 + i * dph, s.zMin), point(s.mOld, ph0 + i * dph, s.zMax)}, kGray + 2, 1, 1);
  }
  // corrected shell outline: end arcs, side edges and a few generators
  for (double z : {s.zMin, s.zMax}) {
    std::vector<V3> arc;
    for (int i = 0; i <= nPhi; i++) {
      arc.push_back(shifted(ph0 + i * dph, z));
    }
    polyline(arc, kBlack, 2, 1);
  }
  for (int i = 0; i <= nPhi; i += std::max(1, nPhi / 6)) { // generators
    polyline({shifted(ph0 + i * dph, s.zMin), shifted(ph0 + i * dph, s.zMax)}, kBlack, 1, 1);
  }
  for (int i : {0, nPhi}) { // both side edges, drawn last
    polyline({shifted(ph0 + i * dph, s.zMin), shifted(ph0 + i * dph, s.zMax)}, kBlack, 2, 1);
  }

  // lab axes with ticks (cm)
  const char* axName[3] = {"x", "y", "z"};
  const int axCol[3] = {kRed + 1, kGreen + 2, kMagenta + 1};
  const double step = niceStep(std::max(axLenXY, axLenZ)), tickLen = 0.012 * span;
  for (int a = 0; a < 3; a++) {
    double x0, y0, x1, y1;
    v.xy({0, 0, 0}, x0, y0);
    v.xy(axEnd[a], x1, y1);
    auto* ar = new TArrow(x0, y0, x1, y1, 0.005, "|>");
    ar->SetLineColor(axCol[a]);
    ar->SetFillColor(axCol[a]);
    ar->SetLineWidth(2);
    ar->Draw();
    auto* t = new TLatex(x1, y1, Form(" %s", axName[a]));
    t->SetTextColor(axCol[a]);
    t->SetTextSize(0.07);
    t->SetTextFont(62);
    t->Draw();
    const double len = a == 2 ? axLenZ : axLenXY;
    for (double pos = step; pos < len; pos += step) {
      V3 p{0, 0, 0};
      p[a] = pos;
      double px, py;
      v.xy(p, px, py);
      const double dx = x1 - x0, dy = y1 - y0, n = std::hypot(dx, dy);
      auto* tk = new TLine(px - tickLen * dy / n, py + tickLen * dx / n, px + tickLen * dy / n, py - tickLen * dx / n);
      tk->SetLineColor(axCol[a]);
      tk->SetLineWidth(2);
      tk->Draw();
      auto* tl = new TLatex(px + 2.2 * tickLen * dy / n, py - 2.2 * tickLen * dx / n, Form("%g", pos));
      tl->SetTextSize(0.04);
      tl->SetTextAlign(22);
      tl->SetTextColor(axCol[a]);
      tl->Draw();
    }
  }

  // annotation
  const double us = 1.e4;
  auto* title = new TLatex(0.03, 0.935, Form("L%d  HB%d", l, h));
  title->SetNDC();
  title->SetTextSize(0.062);
  title->SetTextFont(62);
  title->Draw();
  auto* info = new TLatex(0.03, 0.880, Form("<#Delta>=(%.0f,%.0f,%.0f) #mum  max|#Delta|=%.0f #mum", us * mean[0] / np, us * mean[1] / np, us * mean[2] / np, us * maxd));
  info->SetNDC();
  info->SetTextSize(0.050);
  info->Draw();
  auto* info2 = new TLatex(0.03, 0.830, Form("shift #times%g, ticks %g cm", mag, step));
  info2->SetNDC();
  info2->SetTextSize(0.050);
  info2->SetTextColor(kGray + 2);
  info2->Draw();
}
} // namespace hbshift

// geomOld: reference geometry (solid), geomNew: corrected one (outline); magnification: factor applied to the shift of every
// point of the shell; layer < 0: all layers; outFile: optional image/pdf to save.
TCanvas* drawHBShift(const char* geomOld = "/data/devAlg/best/IBEdge/o2sim_geometry-aligned.root",
                     const char* geomNew = "/data/devAlg/solPbHBA0/o2sim_geometry-aligned.root",
                     double magnification = 500., const char* outFile = "3dview.png",
		     int layer = -1, double rotX = 10., double rotY = -5., double rotZ = 0.,
                     int nPhi = 48)
{
  using namespace hbshift;
  std::map<int, Shell> shells;
  if (!readGeometry(geomOld, shells, false) || !readGeometry(geomNew, shells, true)) {
    return nullptr;
  }
  gStyle->SetOptStat(0);
  const int l0 = layer < 0 ? 0 : layer, l1 = layer < 0 ? NLayers : layer + 1, nl = l1 - l0;
  auto* c = new TCanvas("cHBShift", Form("HB shifts x%g", magnification), 420 * nl, 800);
  c->Divide(nl, NHB, 0.002, 0.002);
  for (int l = l0; l < l1; l++) {
    for (int h = 0; h < NHB; h++) {
      auto* pad = (TPad*)c->cd(h * nl + (l - l0) + 1);
      pad->SetFillColor(kWhite);
      pad->SetBorderMode(0);
      pad->Modified();
      pad->Update();
      drawPad(pad, l, h, shells[l * NHB + h], magnification, rotX, rotY, rotZ, nPhi);
    }
  }
  c->Update();
  if (outFile && outFile[0]) {
    c->SaveAs(outFile);
  }
  return c;
}
