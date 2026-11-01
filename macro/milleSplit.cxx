// Accumulate, from a Mille binary, the normal equations of the ITS IB chips
// with separate DOF copies for their y rows and z rows, after profiling out the
// local (track) parameters of each record. A real displacement gives the same y
// and z copies; milleSplitModels.py compares the rigid, the split and the rigid
// + y-only-incidence models. Needs the Mille data recorded with all 6 chip
// DOFs.
//   g++ -O2 -std=c++17 -o milleSplit milleSplit.cxx
//   python3 -I milleSplitModels.py labels mp2tree.txt iblabels.txt   # chip TX
//   label -> sensor ID
//   ./milleSplit mp2data_0.bin iblabels.txt ms_B0.bin [maxrec]          # ~20
//   min per 6 GB, ~215 MB output
// Records with > 600 locals or chi2/ndf > 6 after the local fit are skipped.
// Parameter index: sens*12 + (isZrow ? 6 : 0) + dof   (dof =
// TX,TY,TZ,RX,RY,RZ), then the mean-vertex calibration (X,Y,Z,SlopeX,SlopeY of
// up to 8 time slots) at NPIB + slot*5 + dof. The mean vertex must be fitted as
// well: the records measure it relative to a reference at (0,0), while the beam
// is ~1.8 mm away; fixing it fakes rphi-only terms. All the other non-IB
// globals are kept at their recorded reference values.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <unordered_map>
#include <vector>
static const int NS = 432, NPS = 12, NPIB = NS * NPS, NMVSLOT = 8,
                 NMV = NMVSLOT * 5, NP = NPIB + NMV;
static const int MVBASE =
    134218112; // label of PVT/meanVertex X, slot 0; slots are 256 apart
struct Row {
  double meas, sig;
  std::vector<std::pair<int, double>> loc, glo;
};
bool cholInv(std::vector<double> &A,
             int n) { // in-place inverse of SPD matrix (row-major)
  std::vector<double> L(n * n, 0.);
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      double s = A[i * n + j];
      for (int k = 0; k < j; k++)
        s -= L[i * n + k] * L[j * n + k];
      if (i == j) {
        if (s <= 1e-30)
          return false;
        L[i * n + i] = std::sqrt(s);
      } else
        L[i * n + j] = s / L[j * n + j];
    }
  std::vector<double> Li(n * n, 0.);
  for (int i = 0; i < n; i++) {
    Li[i * n + i] = 1. / L[i * n + i];
    for (int j = 0; j < i; j++) {
      double s = 0;
      for (int k = j; k < i; k++)
        s -= L[i * n + k] * Li[k * n + j];
      Li[i * n + j] = s / L[i * n + i];
    }
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      double s = 0;
      for (int k = i; k < n; k++)
        s += Li[k * n + i] * Li[k * n + j];
      A[i * n + j] = A[j * n + i] = s;
    }
  return true;
}
int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: milleSplit file.bin labels.txt out.bin [maxrec]\n");
    return 1;
  }
  long maxrec = argc > 4 ? atol(argv[4]) : -1;
  std::unordered_map<int, int> base2sens;
  {
    std::ifstream lf(argv[2]);
    int l, s;
    while (lf >> l >> s)
      base2sens[l] = s;
  }
  std::vector<double> N((size_t)NP * NP, 0.), B(NP, 0.);
  std::vector<long> cnt(NP, 0);
  FILE *f = fopen(argv[1], "rb");
  std::vector<float> vf;
  std::vector<double> vd;
  std::vector<int> vi;
  long nrec = 0, nused = 0, nbad = 0;
  double chi2sum = 0;
  long ndfsum = 0;
  while (maxrec < 0 || nrec < maxrec) {
    int n2;
    if (fread(&n2, 4, 1, f) != 1)
      break;
    bool dbl = n2 < 0;
    int n = std::abs(n2) / 2;
    std::vector<double> val(n);
    if (dbl) {
      vd.resize(n);
      fread(vd.data(), 8, n, f);
      for (int i = 0; i < n; i++)
        val[i] = vd[i];
    } else {
      vf.resize(n);
      fread(vf.data(), 4, n, f);
      for (int i = 0; i < n; i++)
        val[i] = vf[i];
    }
    vi.resize(n);
    fread(vi.data(), 4, n, f);
    nrec++;
    std::vector<Row> rows;
    int state = 0;
    for (int i = 1; i < n; i++) {
      double v = val[i];
      int k = vi[i];
      if (k == 0) {
        if (state == 0 || state == 2) {
          rows.push_back({v, 0, {}, {}});
          state = 1;
        } else {
          rows.back().sig = v;
          state = 2;
        }
      } else if (state == 1)
        rows.back().loc.push_back({k, v});
      else
        rows.back().glo.push_back({k, v});
    }
    std::map<int, int> lidx;
    for (auto &r : rows)
      for (auto &p : r.loc)
        lidx.emplace(p.first, 0);
    int nl = 0;
    for (auto &p : lidx)
      p.second = nl++;
    if (nl == 0 || nl > 600)
      continue;
    std::vector<double> A(nl * nl, 0.), br(nl, 0.);
    for (auto &r : rows) {
      if (r.sig <= 0)
        continue;
      double w = 1. / (r.sig * r.sig);
      for (auto &a : r.loc) {
        int i = lidx[a.first];
        br[i] += w * a.second * r.meas;
        for (auto &b : r.loc)
          A[i * nl + lidx[b.first]] += w * a.second * b.second;
      }
    }
    std::vector<double> Ai = A;
    if (!cholInv(Ai, nl)) {
      nbad++;
      continue;
    }
    std::vector<double> xr(nl, 0.);
    for (int i = 0; i < nl; i++)
      for (int j = 0; j < nl; j++)
        xr[i] += Ai[i * nl + j] * br[j];
    double chi2 = 0;
    int nm = 0;
    for (auto &r : rows) {
      if (r.sig <= 0)
        continue;
      double p = 0;
      for (auto &a : r.loc)
        p += a.second * xr[lidx[a.first]];
      chi2 += std::pow((r.meas - p) / r.sig, 2);
      nm++;
    }
    int ndf = nm - nl;
    if (ndf <= 0 || chi2 / ndf > 6.) {
      nbad++;
      continue;
    }
    chi2sum += chi2;
    ndfsum += ndf;
    // collect IB-chip parameter derivatives per row
    std::map<int, std::vector<double>>
        u; // param -> local cross vector (sum w*g*dloc)
    std::vector<std::vector<std::pair<int, double>>> rg(rows.size());
    for (size_t ir = 0; ir < rows.size(); ir++) {
      auto &r = rows[ir];
      if (r.sig <= 0)
        continue;
      for (auto &g : r.glo) {
        int o = g.first - MVBASE;
        if (o >= 0 && o < NMVSLOT * 256 && o % 256 < 5)
          rg[ir].push_back({NPIB + (o / 256) * 5 + o % 256, g.second});
      }
      std::unordered_map<int, double[6]> per;
      per.reserve(4);
      for (auto &g : r.glo) {
        int base = g.first & ~7, d = g.first & 7;
        auto it = base2sens.find(base);
        if (it == base2sens.end() || d > 5)
          continue;
        per[it->second][d] = g.second;
      }
      for (auto &[s, d] : per) {
        bool zrow = std::fabs(d[2]) > std::fabs(d[0]);
        for (int k = 0; k < 6; k++)
          if (d[k] != 0)
            rg[ir].push_back({s * NPS + (zrow ? 6 : 0) + k, d[k]});
      }
    }
    for (size_t ir = 0; ir < rows.size(); ir++) {
      auto &r = rows[ir];
      if (rg[ir].empty())
        continue;
      double w = 1. / (r.sig * r.sig);
      for (auto &[p, g] : rg[ir]) {
        auto &uv = u[p];
        if (uv.empty())
          uv.assign(nl, 0.);
        for (auto &a : r.loc)
          uv[lidx[a.first]] += w * g * a.second;
        B[p] += w * g * r.meas;
        cnt[p]++;
        for (auto &[q, h] : rg[ir])
          N[(size_t)p * NP + q] += w * g * h;
      }
    }
    if (u.empty())
      continue;
    nused++;
    std::vector<int> ps;
    std::vector<std::vector<double>> xs;
    for (auto &[p, uv] : u) {
      std::vector<double> x(nl, 0.);
      for (int i = 0; i < nl; i++)
        for (int j = 0; j < nl; j++)
          x[i] += Ai[i * nl + j] * uv[j];
      ps.push_back(p);
      xs.push_back(x);
    }
    for (size_t a = 0; a < ps.size(); a++) {
      auto &ua = u[ps[a]];
      double sb = 0;
      for (int i = 0; i < nl; i++)
        sb += ua[i] * xr[i];
      B[ps[a]] -= sb;
      for (size_t b = 0; b < ps.size(); b++) {
        double s = 0;
        for (int i = 0; i < nl; i++)
          s += ua[i] * xs[b][i];
        N[(size_t)ps[a] * NP + ps[b]] -= s;
      }
    }
    if (nrec % 500000 == 0)
      fprintf(stderr, "%ld records, %ld used\n", nrec, nused);
  }
  fprintf(stderr, "records %ld used %ld rejected %ld  <chi2/ndf> %.3f\n", nrec,
          nused, nbad, chi2sum / std::max(1L, ndfsum));
  FILE *o = fopen(argv[3], "wb");
  int np = NP;
  fwrite(&np, 4, 1, o);
  fwrite(N.data(), 8, N.size(), o);
  fwrite(B.data(), 8, NP, o);
  fwrite(cnt.data(), 8, NP, o);
  fclose(o);
}
