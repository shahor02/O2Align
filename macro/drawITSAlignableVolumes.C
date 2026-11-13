// Draw ITS2 alignable frames. Lengths are in cm.
// root -l 'macro/drawITSAlignableVolumes.C+("o2sim_geometry-aligned.root",0,2,true)'
// Arguments: geometry file, first layer, last layer (-1 = last), ideal overlay.
// Plates use ideal chip envelopes in each alignable's own frame, flattened
// onto their mean local Y. Half-barrels use mean ideal sensor tracking X,
// exactly pi radians, and the ideal sensor Z envelope. Each frame is moved
// rigidly by its own alignable matrix (daughter misalignments do not resize it).
#include <DetectorsBase/GeometryManager.h>
#include <ITSBase/GeometryTGeo.h>
#include <TCanvas.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TLegend.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TView.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace itsAlignableFrames
{
using Point = std::array<double, 3>;
struct Volume {
  TGeoHMatrix ideal, aligned;
  int kind = -1, sensors = 0;
  Point lo{{1.e9, 1.e9, 1.e9}}, hi{{-1.e9, -1.e9, -1.e9}};
  double y = 0., radius = 0., cosPhi = 0., sinPhi = 0.;
};
constexpr int colors[] = {kRed, kCyan, kGreen, kBlue, kBlack};
int kindOf(const std::string& name)
{
  const auto leaf = name.substr(name.rfind('/') + 1);
  if (leaf.find("ITSUChip") == 0) return 0;
  if (leaf.find("ITSUModule") == 0) return 1;
  if (leaf.find("ITSUHalfStave") == 0) return 2;
  if (leaf.find("ITSUStave") == 0) return 3;
  if (leaf.find("ITSUHalfBarrel") == 0) return 4;
  return -1;
}
Volume makeVolume(const std::string& name, int kind)
{
  Volume v;
  v.kind = kind;
  auto* entry = gGeoManager->GetAlignableEntry(name.c_str());
  if (!entry || !gGeoManager->cd(entry->GetTitle()))
    throw std::runtime_error("Cannot navigate to " + name);
  v.aligned = *gGeoManager->GetCurrentMatrix();
  if (!o2::base::GeometryManager::getOriginalMatrix(name.c_str(), v.ideal))
    throw std::runtime_error("Cannot obtain original matrix for " + name);
  return v;
}
}

TCanvas* drawITSAlignableVolumes(const char* file = "o2sim_geometry-aligned.root", int firstLayer = 0, int lastLayer = -1, bool drawIdeal = false)
{
  using namespace itsAlignableFrames;
  // Loading another geometry after creating GeometryTGeo leaves its singleton
  // indexed against the old geometry. Require a fresh ROOT session instead.
  if (gGeoManager) {
    std::cerr << "Please run in a fresh ROOT session to load the requested geometry safely.\n";
    return nullptr;
  }
  o2::base::GeometryManager::loadGeometry(file);
  if (!gGeoManager) return nullptr;
  auto geom = o2::its::GeometryTGeo::Instance();
  if (lastLayer == -1) lastLayer = geom->getNumberOfLayers() - 1;
  if (firstLayer < 0 || lastLayer < firstLayer || lastLayer >= geom->getNumberOfLayers()) {
    std::cerr << "Invalid layer range: " << firstLayer << ".." << lastLayer << '\n';
    return nullptr;
  }
  std::map<std::string, Volume> volumes;
  try {
    for (int layer = firstLayer; layer <= lastLayer; ++layer) {
      for (int chip = geom->getFirstChipIndex(layer); chip <= geom->getLastChipIndex(layer); ++chip) {
        const std::string name = geom->getSymbolicName(chip);
        if (kindOf(name) != 0) throw std::runtime_error("Only planar ITS2 chips are supported: " + name);
        auto sensor = makeVolume(name, 0);
        auto* entry = gGeoManager->GetAlignableEntry(name.c_str());
        gGeoManager->cd(entry->GetTitle());
        auto* box = dynamic_cast<TGeoBBox*>(gGeoManager->GetCurrentVolume()->GetShape());
        if (!box) throw std::runtime_error("Chip shape is not a rectangular box: " + name);
        // Match GeometryTGeo::extractSensorXAlpha, using the ORIGINAL matrix.
        const auto* tr = sensor.ideal.GetTranslation();
        const double unitX[3] = {1., 0., 0.};
        double tangent[3];
        sensor.ideal.LocalToMasterVect(unitX, tangent);
        const double t = (tr[0] * tangent[0] + tr[1] * tangent[1]) /
                         (tangent[0] * tangent[0] + tangent[1] * tangent[1]);
        const double radius = std::hypot(tr[0] - t * tangent[0], tr[1] - t * tangent[1]);
        for (std::string parent = name; !parent.empty();) {
          const int kind = kindOf(parent);
          if (kind >= 0 && gGeoManager->GetAlignableEntry(parent.c_str())) {
            auto it = volumes.find(parent);
            if (it == volumes.end()) it = volumes.emplace(parent, makeVolume(parent, kind)).first;
            auto& v = it->second;
            Point center;
            v.ideal.MasterToLocal(tr, center.data());
            v.y += center[1];
            v.radius += radius;
            const double phi = std::atan2(center[1], center[0]);
            v.cosPhi += std::cos(phi);
            v.sinPhi += std::sin(phi);
            ++v.sensors;
            for (int ix : {-1, 1}) for (int iz : {-1, 1}) {
              const double local[3] = {box->GetOrigin()[0] + ix * box->GetDX(), box->GetOrigin()[1],
                                       box->GetOrigin()[2] + iz * box->GetDZ()};
              Point lab, p;
              sensor.ideal.LocalToMaster(local, lab.data());
              v.ideal.MasterToLocal(lab.data(), p.data());
              for (int k = 0; k < 3; ++k) {
                v.lo[k] = std::min(v.lo[k], p[k]);
                v.hi[k] = std::max(v.hi[k], p[k]);
              }
            }
          }
          const auto slash = parent.rfind('/');
          if (slash == std::string::npos) break;
          parent.resize(slash);
        }
      }
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return nullptr;
  }
  static unsigned serial = 0;
  auto* canvas = new TCanvas(Form("itsFrames%u", ++serial), "ITS alignable frames; coordinates in cm", 1200, 900);
  canvas->SetFillColor(kWhite);
  auto* view = TView::CreateView(1);
  Point low{{1.e9, 1.e9, 1.e9}}, high{{-1.e9, -1.e9, -1.e9}};
  auto line = [&](const std::vector<Point>& points, const Volume& v, bool ideal) {
    auto* poly = new TPolyLine3D(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      Point p;
      (ideal ? v.ideal : v.aligned).LocalToMaster(points[i].data(), p.data());
      poly->SetPoint(i, p[0], p[1], p[2]);
      for (int k = 0; k < 3; ++k) {
        low[k] = std::min(low[k], p[k]);
        high[k] = std::max(high[k], p[k]);
      }
    }
    poly->SetLineColor(colors[v.kind]);
    poly->SetLineStyle(ideal ? 2 : 1);
    poly->SetLineWidth(1);
    poly->SetBit(kCanDelete);
    poly->Draw();
  };
  std::array<int, 5> counts{};
  const double pi = std::acos(-1.);
  for (const auto& item : volumes) {
    const auto& v = item.second;
    ++counts[v.kind];
    for (int ideal = 0; ideal <= int(drawIdeal); ++ideal) {
      if (v.kind != 4) {
        const double y = v.y / v.sensors;
        line({{v.lo[0], y, v.lo[2]}, {v.hi[0], y, v.lo[2]}, {v.hi[0], y, v.hi[2]},
              {v.lo[0], y, v.hi[2]}, {v.lo[0], y, v.lo[2]}}, v, ideal);
      } else {
        const double r = v.radius / v.sensors;
        const double start = std::atan2(v.sinPhi, v.cosPhi) - pi / 2.;
        for (double z : {v.lo[2], v.hi[2]}) {
          std::vector<Point> arc;
          for (int j = 0; j <= 120; ++j) {
            const double phi = start + pi * j / 120.;
            arc.push_back({r * std::cos(phi), r * std::sin(phi), z});
          }
          line(arc, v, ideal);
        }
        for (double phi : {start, start + pi})
          line({{r * std::cos(phi), r * std::sin(phi), v.lo[2]},
                {r * std::cos(phi), r * std::sin(phi), v.hi[2]}}, v, ideal);
        if (!ideal) std::cout << item.first << ": mean ideal tracking X = " << r << " cm\n";
      }
    }
  }
  for (int k = 0; k < 3; ++k) {
    const double margin = std::max(1., 0.08 * (high[k] - low[k]));
    low[k] -= margin;
    high[k] += margin;
  }
  view->SetRange(low.data(), high.data());
  view->ShowAxis();
  int viewError = 0;
  view->SetView(40., 70., 0., viewError);
  auto* legend = new TLegend(0.01, 0.73, 0.25, 0.99);
  legend->SetHeader(Form("ITS layers %d-%d", firstLayer, lastLayer));
  legend->SetFillStyle(0);
  const char* labels[] = {"Sensors/chips", "Modules", "Half-staves", "Staves", "Half-barrels"};
  for (int kind = 0; kind < 5; ++kind) {
    auto* sample = new TPolyLine;
    sample->SetLineColor(colors[kind]);
    // Legend owns no sample: keep these tiny style objects alive with the canvas.
    sample->SetBit(kCanDelete);
    canvas->GetListOfPrimitives()->Add(sample);
    legend->AddEntry(sample, Form("%s (%d)", labels[kind], counts[kind]), "l");
  }
  if (drawIdeal) {
    auto* sample = new TPolyLine;
    sample->SetLineStyle(2);
    sample->SetBit(kCanDelete);
    canvas->GetListOfPrimitives()->Add(sample);
    legend->AddEntry(sample, "Ideal positions (same colors)", "l");
  }
  legend->SetBit(kCanDelete);
  legend->Draw();
  canvas->Modified();
  canvas->Update();
  std::cout << "Frames: " << volumes.size() << "; solid = aligned"
            << (drawIdeal ? ", dashed = ideal" : "") << '\n';
  return canvas;
}
