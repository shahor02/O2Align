// Run in an O2 environment:
// root -l -b -q 'drawChipResidualMaps.C+("solC1v19/checkResid_LHC25ah_564644_hman_AlgRS_PV1.root")'
// Arguments: input, output prefix, geometry, first layer, last layer, ID offset,
//            float nMads = 4.f.
// geometry: optional local geometry file; null/empty loads GLO/Config/Geometry
// from CCDB for run 566190. An already loaded gGeoManager is reused.
// Color scale, independently for each map:
//   nMads == 0: automatic ROOT scale, without explicit limits.
//   nMads < 0:  minimum/maximum of the non-empty fitted values.
//   nMads > 0:  [max(minValue, median - nMads * MAD),
//                min(maxValue, median + nMads * MAD)], using raw MAD
//                (median absolute deviation) and non-empty data extrema.
// Scale calculations ignore missing fits (zero content AND error) and non-finite
// values. Degenerate ranges get a tiny padding so ROOT can display them.
// COLZ0 saturates colors at both limits; bin contents/errors are unchanged.
// ID offset 0 writes standard residuals and a second pair <output>_iob.{root,pdf}.
// A nonzero ID offset selects only that residual type (100000 selects IOB).
// Each PDF page contains the four maps for one layer / half-barrel.
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "GlobalTrackingStudy/HistoManager.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSBase/GeometryTGeo.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#else
R__LOAD_LIBRARY(libO2GlobalTrackingStudy)
R__LOAD_LIBRARY(libO2ITSBase)
R__LOAD_LIBRARY(libO2ITSMFTReconstruction)
R__LOAD_LIBRARY(libO2DetectorsBase)
R__LOAD_LIBRARY(libO2CCDB)
#endif


void drawChipResidualMaps(const char* input="checkResid_LHC25ab_562862_hman_AlgRS_PV1.root",  // input histomanager file from CheckResid workflow
			  const char* output = "chipResidualMaps_B0", // output name (pdf and root will be added)
			  const char* geometry = nullptr, // if empty, load from CCDB for run 566190
			  int firstLayer = 0, int lastLayer = 6, int idOffset = 0,
			  float nMads = 4.f)
{
  if (firstLayer < 0 || lastLayer > 6 || firstLayer > lastLayer) {
    throw std::runtime_error("Layer range must be within 0..6");
  }
  if (gSystem->AccessPathName(input)) {
    throw std::runtime_error(std::string("Cannot read ") + input);
  }
  // An existing geometry can also be supplied by the caller before invoking us.
  if (!gGeoManager) {
    if (geometry && geometry[0]) {
      o2::base::GeometryManager::loadGeometry(geometry);
    } else {
      gGeoManager = o2::ccdb::BasicCCDBManager::instance().getForRun<TGeoManager>("GLO/Config/GeometryAligned", 566190);
      if (!gGeoManager) {
        throw std::runtime_error("Cannot load GLO/Config/Geometry from CCDB for run 566190");
      }
    }
  }
  auto* geom = o2::its::GeometryTGeo::Instance();
  o2::itsmft::ChipMappingITS mapping;
  o2::HistoManager hm("", input, true);
  // Check both names before opening either output, to avoid overwriting input.
  const std::vector<int> offsets = idOffset == 0 ? std::vector<int>{0, 100000}
                                               : std::vector<int>{idOffset};
  for (const int offset : offsets) {
    const TString prefix = TString(output) + (offset == 100000 ? "_iob" : "");
    if (prefix + ".root" == input) {
      throw std::runtime_error("Output must differ from input");
    }
  }
  for (const int offset : offsets) {
    const TString prefix = TString(output) + (offset == 100000 ? "_iob" : "");
    const TString rootName = prefix + ".root";
    const TString pdfName = prefix + ".pdf";
    TFile out(rootName, "RECREATE");
    if (out.IsZombie()) {
      throw std::runtime_error("Cannot create output ROOT file");
    }
    TCanvas canvas("chipResidualMapsCanvas", "ITS per-chip residual maps", 1600, 1100);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFont(42, "t");
    gStyle->SetTitleFontSize(0.045);
    gStyle->SetTitleH(0.07);
    gStyle->SetTitleW(0.8);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(0.98);
    gStyle->SetNumberContours(100);
    gStyle->SetPaperSize(26, 18);
    canvas.Print(pdfName + "[");
    int nMaps = 0;
    for (int layer = firstLayer; layer <= lastLayer; ++layer) {
      const int nStaves = mapping.getNStavesOnLr(layer);
      const int nHalf = mapping.getNChipsOnLayer(layer) / 2;
      const int nRows = layer < 3 ? 1 : 4;
      const int nZ = mapping.getNChipsOnLayer(layer) / nStaves / nRows;
      const int nPhi = nStaves / 2 * nRows;
      for (int hb = 0; hb < 2; ++hb) {
        const int firstChip = mapping.getFirstChipsOnLayer(layer) + hb * nHalf;
        struct Chip { int id; double z, phi; };
        // Topology identifies the four rows of an OB/MB stave. Geometry orders
        // these rows in increasing phi and the chips within them in increasing Z.
        // This avoids assuming that numerical chip IDs increase along Z.
        std::vector<std::vector<Chip>> rows(nPhi);
        for (int id = firstChip; id < firstChip + nHalf; ++id) {
          int lr, stave, halfStave, module, chip;
          mapping.expandChipInfoSW(id, lr, stave, halfStave, module, chip);
          const int localStave = (id - firstChip) / (nRows * nZ);
          const int row = layer < 3 ? 0 : 2 * halfStave + chip / 7;
          auto* mat = geom->getMatrix(id);
          if (!mat) {
            throw std::runtime_error("Missing chip geometry matrix");
          }
          const auto* xyz = mat->GetTranslation();
          // getMatrix() may return a shared temporary in an ideal geometry.
          // Copy coordinates before the reference-chip query overwrites it.
          const double z = xyz[2];
          double phi = std::atan2(xyz[1], xyz[0]);
          // Unwrap about the first chip in each stave, including across 2*pi.
          auto* ref = geom->getMatrix(firstChip + localStave * nRows * nZ);
          const auto* r = ref->GetTranslation();
          double phi0 = std::atan2(r[1], r[0]);
          if (phi0 < 0) phi0 += 2 * std::acos(-1.);
          phi = phi0 + std::remainder(phi - phi0, 2 * std::acos(-1.));
          rows[localStave * nRows + row].push_back({id, z, phi});
        }
        auto meanPhi = [](const std::vector<Chip>& row) {
          double sum = 0;
          for (const auto& chip : row) sum += chip.phi;
          return sum / row.size();
        };
        for (auto& row : rows) {
          if (int(row.size()) != nZ) throw std::runtime_error("Unexpected chip topology");
          std::sort(row.begin(), row.end(), [](const Chip& a, const Chip& b) { return a.z < b.z; });
        }
        std::sort(rows.begin(), rows.end(), [&](const auto& a, const auto& b) { return meanPhi(a) < meanPhi(b); });
        canvas.Clear();
        canvas.Divide(2, 2);
        std::vector<std::unique_ptr<TH2D>> plots;
        std::vector<std::unique_ptr<TLine>> lines;
        for (int yz = 0; yz < 2; ++yz) {
          for (int stat = 1; stat <= 2; ++stat) {
            const int pad = 2 * yz + stat;
            canvas.cd(pad)->SetRightMargin(0.16);
            gPad->SetLeftMargin(0.18);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.10);
            const int id = offset + 15000 + (2 * (layer + 1) + hb) * 100 + yz * 10 + stat;
            // FitSlicesY produces TH1D, so use the generic HistoManager accessor.
            auto* source = hm.getHisto(id);
            if (!source) {
              Warning("drawChipResidualMaps", "Missing fitted histogram ID %d (L%d HB%d)", id, layer, hb);
              continue;
            }
            if (source->GetNbinsX() != nHalf ||
                std::abs(source->GetXaxis()->GetBinCenter(1) - firstChip) > 0.01) {
              throw std::runtime_error("Source histogram chip axis disagrees with ITS mapping");
            }
            const TString quantity = stat == 1 ? Form("#LTd%s#GT", yz ? "Z" : "Y")
                                               : Form("#LT#sigma_{d%s}#GT", yz ? "Z" : "Y");
            const TString name = Form("L%d_HB%d_%s_%s", layer, hb, yz ? "dZ" : "dY", stat == 1 ? "mean" : "sigma");
            auto h = std::make_unique<TH2D>(name,
              Form("ITS layer %d, half-barrel %d: %s%s;Chip position along stave (increasing Z);Chip row (increasing #phi);%s [#mum]",
                   layer, hb, quantity.Data(), offset == 100000 ? " (IOB)" : "", quantity.Data()),
              nZ, -0.5, nZ - 0.5, nPhi, -0.5, nPhi - 0.5);
            h->SetDirectory(nullptr);
            for (auto* axis : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
              axis->SetTitleFont(42);
              axis->SetTitleSize(0.04);
              axis->SetLabelFont(42);
              axis->SetLabelSize(0.035);
            }
            h->GetYaxis()->SetLabelSize(0.03);
            h->GetYaxis()->SetTitleOffset(1.8);
            h->GetYaxis()->CenterTitle();
            h->GetZaxis()->SetTitleOffset(1.3);
            std::vector<double> fittedValues;
            for (int y = 0; y < nPhi; ++y) {
              int lr, stave, hs, mod, ch;
              mapping.expandChipInfoSW(rows[y][0].id, lr, stave, hs, mod, ch);
              // Label once per stave to keep the outer-layer pages readable;
              // the four chip rows remain separate bins, ordered in phi.
              h->GetYaxis()->SetBinLabel(y + 1, y % nRows == 0
                ? Form("S%d", stave - mapping.getFirstStavesOnLr(layer)) : " ");
              for (int x = 0; x < nZ; ++x) {
                const int bin = source->GetXaxis()->FindFixBin(rows[y][x].id);
                const double value = source->GetBinContent(bin) * 1.e4; // cm -> micrometres
                h->SetBinContent(x + 1, y + 1, value);
                h->SetBinError(x + 1, y + 1, source->GetBinError(bin) * 1.e4);
                if (nMads != 0 && std::isfinite(value) &&
                    (value != 0 || source->GetBinError(bin) != 0)) {
                  fittedValues.push_back(value);
                }
              }
            }
            // Missing/failed slice fits retain the source's zero content/error.
            if (!fittedValues.empty()) {
              auto median = [](std::vector<double> values) {
                std::sort(values.begin(), values.end());
                const auto n = values.size();
                return n % 2 ? values[n / 2] : 0.5 * (values[n / 2 - 1] + values[n / 2]);
              };
              const auto extrema = std::minmax_element(fittedValues.begin(), fittedValues.end());
              double low = *extrema.first, high = *extrema.second;
              if (nMads > 0) {
                const double center = median(fittedValues);
                for (auto& value : fittedValues) value = std::abs(value - center);
                const double mad = median(fittedValues);
                low = std::max(low, center - nMads * mad);
                high = std::min(high, center + nMads * mad);
              }
              if (low == high) {
                const double padding = std::max(1., std::abs(low)) * 1.e-6;
                low -= padding;
                high += padding;
              }
              h->SetMinimum(low);
              h->SetMaximum(high);
            }
            h->Draw("colz0");
            for (int y = nRows; y < nPhi; y += nRows) {
              auto line = std::make_unique<TLine>(-0.5, y - 0.5, nZ - 0.5, y - 0.5);
              line->SetLineStyle(3);
              line->Draw();
              lines.push_back(std::move(line));
            }
            out.cd();
            h->Write();
            plots.push_back(std::move(h));
            ++nMaps;
          }
        }
        if (!plots.empty()) {
          canvas.Print(pdfName);
          out.cd();
          canvas.Write(Form("canvas_L%d_HB%d", layer, hb));
        }
        canvas.Clear(); // remove pad references before deleting histograms/lines
      }
    }
    canvas.Print(pdfName + "]");
    out.Close();
    Info("drawChipResidualMaps", "Saved %d maps to %s and %s", nMaps, rootName.Data(), pdfName.Data());
  }
}
