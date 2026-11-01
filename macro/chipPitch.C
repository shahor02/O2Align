// Pitch of the chip planes: d(x_tracking)/d(z_local) of every IB chip of an
// aligned geometry (GeometryAligned snapshot). root -l -b -q
// 'chipPitch.C("GLO/Config/GeometryAligned/snapshot.root")'   (function name
// clx2 inside)
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <cmath>
#include <cstdio>
void chipPitch(const char *fn) {
  TFile f(fn);
  auto *g = (TGeoManager *)f.Get("ccdb_object");
  if (!g) {
    printf("no geometry\n");
    return;
  }
  const int nSt[3] = {12, 16, 20};
  for (int lay = 0; lay < 3; lay++) {
    double s0 = 0, s1 = 0, mx = 0;
    int n = 0;
    for (int hb = 0; hb < 2; hb++)
      for (int st = 0; st < nSt[lay] / 2; st++)
        for (int ch = 0; ch < 9; ch++) {
          auto *e = g->GetAlignableEntry(
              Form("ITS/ITSULayer%d/ITSUHalfBarrel%d/ITSUStave%d/ITSUChip%d",
                   lay, hb, st, ch));
          if (!e)
            continue;
          auto *pn = e->GetPhysicalNode();
          if (!pn)
            continue;
          const TGeoHMatrix *m = pn->GetMatrix();
          double lz[3] = {0, 0, 1}, gz[3], lo[3] = {0, 0, 0}, go[3];
          m->LocalToMasterVect(lz, gz);
          m->LocalToMaster(lo, go);
          // alpha as getFlatChipAlpha: foot of the normal from the beam axis
          // along local X
          double A[3] = {-100, 0, 0}, B[3] = {100, 0, 0}, ga[3], gb[3];
          m->LocalToMaster(A, ga);
          m->LocalToMaster(B, gb);
          double dx = gb[0] - ga[0], dy = gb[1] - ga[1];
          double t = (gb[0] * dx + gb[1] * dy) / (dx * dx + dy * dy);
          double xp = gb[0] - dx * t, yp = gb[1] - dy * t, al = atan2(yp, xp);
          double slope = gz[0] * cos(al) + gz[1] * sin(al); // d x_trk / d z_loc
          s0 += slope;
          s1 += slope * slope;
          n++;
          if (fabs(slope) > fabs(mx))
            mx = slope;
          if (n <= 2)
            printf("   L%d hb%d st%d ch%d  r=%.4f  dxtrk/dzloc=%.3e\n", lay, hb,
                   st, ch, hypot(go[0], go[1]), slope);
        }
    if (n)
      printf("Layer %d: <dx_trk/dz_loc> = %.3e rms %.3e max %.3e (n=%d)\n", lay,
             s0 / n, sqrt(fmax(0, s1 / n - s0 * s0 / n / n)), mx, n);
  }
}
