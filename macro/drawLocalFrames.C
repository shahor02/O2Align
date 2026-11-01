// Lab reference followed by local rigid-body axes arranged parent -> child.
// Call as drawLocalFrames(<geomName>, std::vector<std::string> volNames,
// <azimuth = 35>, <elevation = 25>) to draw the local frames as seen in the lab
// for the whole hiearchy of the volume. E.g
// drawLocalFrames("o2sim_geometry-aligned.root",
//                     {"ITS/ITSULayer6/ITSUHalfBarrel0/ITSUStave6/ITSUHalfStave0/ITSUModule0/ITSUChip4",
//                      "ITS/ITSULayer6/ITSUHalfBarrel0/ITSUStave0/ITSUHalfStave0/ITSUModule0/ITSUChip4",
//                      "ITS/ITSULayer0/ITSUHalfBarrel0/ITSUStave3/ITSUChip4",
//                      "TRD/sm02/st3/pl0","TOF/sm03/strip35"})
// draws the doc/local_frames.pdf
//
// At zero angles: lab Z right, lab Y up, lab X into the screen.
// Azimuth rotates about lab Y from this view; elevation tilts the camera up.
// Both angles are in degrees.
// Display origins are separated for readability; physical translations are not
// drawn.

#include <DetectorsBase/GeometryManager.h>
#include <Framework/Logger.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

void drawLocalFrames(const std::string &volumeName, double azimuth = 35.,
                     double elevation = 25.) {
  auto *target = gPad;
  std::string volume{volumeName};
  while (!volume.empty() && volume.front() == '/')
    volume.erase(0, 1);
  while (!volume.empty() && volume.back() == '/')
    volume.pop_back();
  if (!gGeoManager->GetAlignableEntry(volume.c_str())) {
    std::cerr << "No alignable entry: " << volume << '\n';
    return;
  }
  struct Frame {
    std::string name;
    std::string label;
    TGeoHMatrix l2g;
  };
  std::vector<Frame> frames;
  std::size_t labelStart = 0;
  // Symbolic hierarchy, not the intervening mechanical nodes in the TGeo path.
  for (std::size_t end = volume.find('/');; end = volume.find('/', end + 1)) {
    const auto prefix = volume.substr(0, end);
    auto *entry = gGeoManager->GetAlignableEntry(prefix.c_str());
    if (entry) {
      auto *node = entry->GetPhysicalNode();
      if (!node)
        node = gGeoManager->MakePhysicalNode(entry->GetTitle());
      if (!node || !node->GetMatrix()) {
        std::cerr << "Cannot obtain local-to-lab matrix: " << prefix << '\n';
        return;
      }
      frames.push_back({prefix, prefix.substr(labelStart), *node->GetMatrix()});
      labelStart = prefix.size() + 1;
    } else {
      std::cerr << "Skipping symbolic prefix without an alignable entry: "
                << prefix << '\n';
    }
    if (end == std::string::npos)
      break;
  }

  // Identity local-to-lab transformation supplies the lab reference triplet.
  frames.insert(frames.begin(), Frame{"Lab frame", "Lab frame", TGeoHMatrix{}});

  const double az = azimuth * TMath::DegToRad();
  const double el = elevation * TMath::DegToRad();
  // Screen right/up; camera position direction is
  // (-cos(el)*cos(az), sin(el), cos(el)*sin(az)). No perspective distortion.
  const std::array<double, 3> right{{std::sin(az), 0., std::cos(az)}};
  const std::array<double, 3> up{{std::sin(el) * std::cos(az), std::cos(el),
                                  -std::sin(el) * std::sin(az)}};
  const int colors[3] = {kRed + 1, kGreen + 2, kBlue + 1};
  const char *labels[3] = {"X", "Y", "Z"};
  static unsigned serial = 0;
  ++serial;
  const double width = target->GetCanvas()->GetWw() * target->GetAbsWNDC();
  const double height = target->GetCanvas()->GetWh() * target->GetAbsHNDC();
  auto *heading = new TLatex;
  heading->SetNDC();
  heading->SetTextFont(42);
  heading->SetTextSize(0.035);
  heading->DrawLatex(
      0.015, 0.95,
      "Local rigid-body axes in lab orientation: parent #rightarrow child");
  heading->SetTextSize(0.026);
  heading->DrawLatex(0.015, 0.89,
                     Form("Orthographic lab view: azimuth %.1f#circ, elevation "
                          "%.1f#circ; origins separated",
                          azimuth, elevation));

  for (std::size_t i = 0; i < frames.size(); ++i) {
    target->cd();
    auto *pad = new TPad(Form("frame_%u_%zu", serial, i),
                         frames[i].name.c_str(), double(i) / frames.size(), 0.,
                         double(i + 1) / frames.size(), 0.83);
    pad->Draw();
    pad->cd();
    // Pad dimensions and range give equal horizontal/vertical pixel scales.
    const double aspect = width / (frames.size() * height * 0.83);
    const double halfX = 1.65 * (aspect > 1. ? aspect : 1.);
    pad->Range(-halfX, -halfX / aspect, halfX, halfX / aspect);
    auto *text = new TLatex;
    text->SetTextFont(42);
    text->SetTextAlign(22);
    text->SetTextSize(0.065);
    // Include non-alignable components since the preceding alignable ancestor.
    text->SetNDC();
    text->SetTextSize(0.055);
    text->DrawLatex(0.5, 0.88, frames[i].label.c_str());
    text->SetNDC(false);
    text->SetTextSize(0.065);
    std::cout << frames[i].name << "\n  local axes in lab (X,Y,Z):\n";
    for (int axis = 0; axis < 3; ++axis) {
      double local[3] = {0., 0., 0.}, lab[3];
      local[axis] = 1.;
      frames[i].l2g.LocalToMasterVect(local, lab);
      double x = 0., y = 0.;
      for (int k = 0; k < 3; ++k) {
        x += right[k] * lab[k];
        y += up[k] * lab[k];
      }
      std::cout << "  " << labels[axis] << ": (" << lab[0] << ", " << lab[1]
                << ", " << lab[2] << ")\n";
      text->SetTextColor(colors[axis]);
      if (std::hypot(x, y) < 1.e-8) {
        const double towardsCamera = -std::cos(el) * std::cos(az) * lab[0] +
                                     std::sin(el) * lab[1] +
                                     std::cos(el) * std::sin(az) * lab[2];
        text->DrawLatex(0., -0.22,
                        Form("%s %s", labels[axis],
                             towardsCamera > 0. ? "#odot" : "#otimes"));
        continue;
      }
      auto *arrow = new TArrow(0., 0., x, y, 0.006, "|>");
      arrow->SetLineColor(colors[axis]);
      arrow->SetFillColor(colors[axis]);
      arrow->SetLineWidth(3);
      arrow->Draw();
      text->DrawLatex(1.18 * x, 1.18 * y, labels[axis]);
    }
  }
  target->cd();
  target->Modified();
  target->Update();
}

TCanvas *cnv = nullptr;

void drawLocalFrames(const std::string &geomFile, std::vector<std::string> vols,
                     double azimuth = 35., double elevation = 25.) {
  if (!o2::base::GeometryManager::isGeometryLoaded()) {
    o2::base::GeometryManager::loadGeometry(geomFile);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGP(fatal, "failed to load geometry {}", geomFile);
    }
  }
  int nv = vols.size();
  if (!nv) {
    LOGP(warn, "vector of volumes to draw is empty");
    return;
  }
  cnv = new TCanvas("frames", "frames", 800, 1000);
  cnv->Divide(1, nv);
  int cnvCnt = 0;
  for (auto &vn : vols) {
    cnv->cd(++cnvCnt);
    drawLocalFrames(vn, azimuth, elevation);
  }
}
