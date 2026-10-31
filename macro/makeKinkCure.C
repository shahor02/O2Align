// Build a millepede.res that removes the IB/OB kink measured by checkResid on a solution, see doc/handovers/kink_cure.md.
//
// The cure is a relative IB/OB rotation about R = 12 cm applied to the OB only, per ray: every OB chip (Lr3..Lr6) gets the
// rphi shift  d(phi, t) = omega(phi, t) * (r - 12 cm),  t = z/r of the chip centre, with
//   omega = b(phi, t) / gain,
// b being the measured dPar2 [urad] of the solution: its phi harmonics n = nmin..nmax per tgl bin, linear in t between the
// bin centres and constant beyond the outermost ones; gain = 0.855 is the response of dPar2 to such a rotation from the pp
// B0 fit of kinkResponse.C. The shift is written as a change of the chip TX label (local x), TX += d / (x_hat . phi_hat),
// all other lines of the input millepede.res being copied unchanged. The labels are matched to the chips through the
// "! V TX <symname>" comments of the mp2param.txt of the same solution, so only free chip TX labels are modified.
//
// The output is then turned into an alignment like any pede result: --output MilleRes with the same DOF config, and
// o2-create-aligned-geometry-workflow for checkResid (doc/handovers/kink_cure.md, section "Producing the cure").
//
//   root -b -q 'makeKinkCure.C+("/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root", "/data/devAlg/solA2v19/o2sim_geometry-aligned.root",
//     "/data/devAlg/solA2v19/mp2param.txt", "/data/devAlg/solA2v19/millepede.res", "millepede_cure012.res", 0, 2)'
#include "kinkUtils.h"
#include <fstream>
#include <sstream>
#include <regex>

/// b [urad] at (phi, t) from the per-tgl-bin harmonics H[it][n]
double kinkAt(const TComplex H[kink::NT][4], double phi, double t, int nmin, int nmax)
{
  using namespace kink;
  int i0 = 0;
  double w = 0;
  if (t >= tglCentre(NT - 1)) {
    i0 = NT - 2, w = 1;
  } else if (t > tglCentre(0)) {
    while (t > tglCentre(i0 + 1)) {
      i0++;
    }
    w = (t - tglCentre(i0)) / (tglCentre(i0 + 1) - tglCentre(i0));
  }
  double v = 0;
  for (int n = nmin; n <= nmax; n++) {
    TComplex c = H[i0][n] * (1 - w) + H[i0 + 1][n] * w;
    v += n ? c.Re() * cos(n * phi) + c.Im() * sin(n * phi) : c.Re();
  }
  return v;
}

/// resTree: checkResid tree of the solution (pp B0 recommended); geom: its aligned geometry (or algDump file of it);
/// param: mp2param.txt of the solution; resIn/resOut: millepede.res to modify / to write; harmonics nmin..nmax (<= 3) of dPar2.
void makeKinkCure(const char* resTree, const char* geom, const char* param, const char* resIn, const char* resOut, int nmin = 0, int nmax = 2, double gain = 0.855)
{
  using namespace kink;
  TComplex H[NT][4];
  {
    TFile f(resTree);
    double m[NT][NP], e[NT][NP];
    dParMeans((TTree*)f.Get("res"), 2, m, e);
    printf("dPar2 harmonics of %s [urad]: n0 | n1(cos,sin) | n2(cos,sin)\n", resTree);
    for (int it = 0; it < NT; it++) {
      for (int n = 0; n < 4; n++) {
        H[it][n] = harmonic(m[it], n);
      }
      printf("  [%+4.1f,%+4.1f] %5.0f | (%5.0f,%5.0f) | (%5.0f,%5.0f)\n", TE[it], TE[it + 1], H[it][0].Re(), H[it][1].Re(), H[it][1].Im(), H[it][2].Re(), H[it][2].Im());
    }
  }
  std::map<int, std::string> sym;
  auto frames = loadFrames(geom, &sym);
  std::map<std::string, int> sym2sid;
  for (const auto& [sid, nm] : sym) {
    sym2sid[nm] = sid;
  }
  // free OB chip TX labels
  std::map<uint32_t, int> lab2sid;
  {
    std::ifstream in(param);
    std::string l;
    std::regex rl("^\\s*(\\d+)\\s+\\S+\\s+\\S+\\s+! V TX (ITS/ITSULayer[3-6]/\\S+ITSUChip\\d+)\\s*$");
    std::smatch m;
    while (std::getline(in, l)) {
      if (std::regex_match(l, m, rl)) {
        auto it = sym2sid.find(m[2]);
        if (it != sym2sid.end()) {
          lab2sid[std::stoul(m[1])] = it->second;
        }
      }
    }
  }
  printf("free OB chip TX labels matched to chips: %zu\n", lab2sid.size());
  std::map<uint32_t, double> dtx;
  double s2[7] = {0}, mx = 0;
  int nn[7] = {0};
  for (auto [lab, sid] : lab2sid) {
    const auto& fr = frames[sid];
    double r = fr.c.Perp(), phi = fr.c.Phi(), t = fr.c.Z() / r;
    TVector3 phiHat(-fr.c.Y() / r, fr.c.X() / r, 0);
    double dphi = kinkAt(H, phi, t, nmin, nmax) / gain * (r - RComp) * 1e-6; // cm along +phi
    dtx[lab] = dphi / fr.ax.Dot(phiHat);
    int lr = r < 22 ? 3 : r < 29 ? 4 : r < 37 ? 5 : 6;
    s2[lr] += pow(dphi * 1e4, 2);
    nn[lr]++;
    mx = std::max(mx, fabs(dphi * 1e4));
  }
  printf("rphi shift rms per layer [um]:");
  for (int l = 3; l < 7; l++) {
    printf(" Lr%d %.1f", l, nn[l] ? sqrt(s2[l] / nn[l]) : 0.);
  }
  printf(", max %.1f\n", mx);
  std::ifstream in(resIn);
  std::ofstream out(resOut);
  std::string l;
  int nmod = 0;
  while (std::getline(in, l)) {
    std::istringstream is(l);
    uint32_t lab;
    double v;
    if (is >> lab >> v) {
      auto it = dtx.find(lab);
      if (it != dtx.end()) {
        std::string rest;
        std::getline(is, rest);
        out << Form("%10u  %14.7E", lab, v + it->second) << rest << "\n";
        nmod++;
        continue;
      }
    }
    out << l << "\n";
  }
  printf("wrote %s: %d labels modified\n", resOut, nmod);
}
