// Per-bin MEDIAN and its error (um-scale robust residuals) of every dy/dz vs
// phi,Z,tgl,Pt TH2F of a checkResid _hman file, to a text file. The "_1"
// fitted-mean histograms of some hman files were found to be inflated ~7x
// w.r.t. the projection medians. root -l -b -q
// 'checkResidMedians.C("checkResid_LHC25ae_hman_AlgRS_PV1.root","med.txt")'
void checkResidMedians(const char *f, const char *out) {
  TFile *F = TFile::Open(f);
  FILE *o = fopen(out, "w");
  TIter it(F->GetListOfKeys());
  TKey *k;
  while ((k = (TKey *)it())) {
    TString s = k->GetName();
    if (!(s.BeginsWith("dy_") || s.BeginsWith("dz_")) || s.Contains("pull") ||
        !(s.EndsWith("_phi") || s.EndsWith("_Z") || s.EndsWith("_tgl") ||
          s.EndsWith("_Pt")))
      continue;
    if (k->GetClassName() != TString("TH2F"))
      continue;
    auto *h = (TH2F *)k->ReadObj();
    fprintf(o, "%s", s.Data());
    for (int i = 1; i <= h->GetNbinsX(); i++) {
      auto *p = h->ProjectionY("p", i, i);
      double n = p->GetEntries();
      if (n < 300) {
        fprintf(o, " %g:0:0", h->GetXaxis()->GetBinCenter(i));
        continue;
      }
      double q[3], pr[3] = {0.25, 0.5, 0.75};
      p->GetQuantiles(3, q, pr);
      double sig = (q[2] - q[0]) / 1.349;
      fprintf(o, " %g:%g:%g", h->GetXaxis()->GetBinCenter(i), q[1],
              1.25 * sig / sqrt(n));
    }
    fprintf(o, "\n");
    delete h;
  }
  fclose(o);
}
