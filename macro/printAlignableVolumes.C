// root -l -b -q 'printAlignableVolumes.C("o2sim_geometry.root")'
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Group by detector prefix, including legacy O2 FT0/FV0/MFT names.
void printAlignableVolumes(const char* file = "o2sim_geometry.root")
{
  auto* geom = TGeoManager::Import(file);
  if (!geom) {
    std::cerr << "Cannot load geometry: " << file << '\n';
    return;
  }

  std::map<std::string, std::vector<std::string>> groups;
  for (int i = 0; i < geom->GetNAlignable(); ++i) {
    const auto* entry = geom->GetAlignableEntry(i);
    if (!entry) continue;
    const std::string name = entry->GetName();
    const auto first = name.find_first_not_of('/');
    const std::string trimmed =
      first == std::string::npos ? "" : name.substr(first);
    std::string detector = trimmed.substr(0, trimmed.find('/'));
    if (detector == "FT0A" || detector == "FT0C" ||
        detector.compare(0, 5, "0MOD_") == 0) {
      detector = "FT0";
    } else if (detector.compare(0, 3, "FV0") == 0) {
      detector = "FV0";
    } else if (detector == "MFT_0") {
      detector = "MFT";
    }
    groups[detector.empty() ? "(unnamed)" : detector].push_back(name);
  }

  std::cout << "Geometry: " << file << '\n'
            << "Alignable volumes: " << geom->GetNAlignable() << '\n';
  if (groups.empty()) {
    std::cout << "No alignable entries stored in this geometry.\n";
  }
  for (auto& group : groups) {
    std::sort(group.second.begin(), group.second.end());
    std::cout << "\n" << group.first
              << " (" << group.second.size() << ")\n";
    for (const auto& name : group.second) {
      std::cout << "  " << name << '\n';
    }
  }
}
