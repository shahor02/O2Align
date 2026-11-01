#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TGeoMatrix.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <cstdio>
#include <regex>
#include <string>
#include <vector>

#include "DetectorsCommonDataFormats/AlignParam.h"
#endif

// Closure test of the ITS IB chip rotation about the sensor normal (local RY).
//
// Reads an ITS alignment object (vector<AlignParam> in the local convention, as
// written by the MilleRes stage to ITS_alignment.root), composes an extra local
// RY rotation on top of the IB chips (layers 0-2),
//   eps(chip) = eps[layer] + gRandom->Gaus(0, sigma)     [rad],
// and writes the result as <outDir>/ITS/Calib/Align/snapshot.root, the layout
// used with
//   --condition-remap file://<outDir>=ITS/Calib/Align
// The injected angles are listed in <outDir>/injectChipRY.txt (symName eps) for
// chipRYSummary.py.
//
// The geometry is changed, not the fit: the data then see a chip misaligned by
// -eps, so a consistent chain (derivatives -> pede -> AlignParam) must return a
// chip RY correction of -eps.
//
// root -l -b -q 'injectChipRY.C("ITS_alignment.root", "inj", 6e-4, -6e-4,
// 3e-4, 1.5e-4)' (to be run in the O2 environment)

void injectChipRY(const char *inFile = "ITS_alignment.root",
                  const char *outDir = "inj", double eps0 = 6e-4,
                  double eps1 = -6e-4, double eps2 = 3e-4,
                  double sigma = 1.5e-4) {
  using o2::detectors::AlignParam;
  const double epsL[3] = {eps0, eps1, eps2};

  TFile fin(inFile);
  auto *vec = fin.Get<std::vector<AlignParam>>("ccdb_object");
  if (!vec) {
    printf("No ccdb_object in %s\n", inFile);
    return;
  }

  gSystem->mkdir(Form("%s/ITS/Calib/Align", outDir), true);
  FILE *fl = fopen(Form("%s/injectChipRY.txt", outDir), "w");
  const std::regex chipRE("ITS/ITSULayer([0-2])/ITSUHalfBarrel[0-9]+/"
                          "ITSUStave[0-9]+/ITSUChip[0-9]+$");
  TRandom3 rnd(12345);
  int nInj = 0;
  for (auto &a : *vec) {
    std::smatch m;
    if (!std::regex_search(a.getSymName(), m, chipRE)) {
      continue;
    }
    if (a.isGlobal()) {
      printf("Entry %s is a global delta: this macro composes local deltas "
             "only, rerun the MilleRes with writeLocalAlignParams=true\n",
             a.getSymName().c_str());
      return;
    }
    const double eps =
        epsL[std::stoi(m[1])] + (sigma > 0. ? rnd.Gaus(0., sigma) : 0.);
    // the fit's convention: AlignParam(psi,theta,phi) = (RX,RY,RZ), local delta
    const TGeoHMatrix dInj =
        AlignParam("", -1, 0., 0., 0., 0., eps, 0., false, false)
            .createMatrix();
    TGeoHMatrix comb = a.createMatrix(); // stored local delta of the chip
    comb *= dInj;                        // l' = d_old * d_inj * l
    a.setGlobalParams(comb); // stores the parameters, the entry stays "local"
    fprintf(fl, "%s %.8e\n", a.getSymName().c_str(), eps);
    nInj++;
  }
  fclose(fl);

  TFile fout(Form("%s/ITS/Calib/Align/snapshot.root", outDir), "RECREATE");
  fout.WriteObjectAny(vec, "std::vector<o2::detectors::AlignParam>",
                      "ccdb_object");
  fout.Close();
  printf("Injected RY into %d IB chips (eps = %g, %g, %g rad + gaus %g) -> "
         "%s/ITS/Calib/Align/snapshot.root\n",
         nInj, eps0, eps1, eps2, sigma, outDir);
  if (nInj == 0) {
    printf("WARNING: no IB chip entry in %s: add them (a chip absent from the "
           "initial alignment is not injected)\n",
           inFile);
  }
}
