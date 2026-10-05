#include <iostream>
#include <fnmatch.h>

#include <TCollection.h>
#include <Framework/Logger.h>
#include "O2Align/ProcessingStats.h"

namespace o2::alignrs
{
namespace
{
bool matchesVolumePattern(const std::string& pattern, const std::string& symName)
{
  return pattern.empty() || fnmatch(pattern.c_str(), symName.c_str(), 0) == 0 ||
         fnmatch(("*" + pattern).c_str(), symName.c_str(), 0) == 0;
}
} // namespace

Long64_t ProcessingStats::Merge(TCollection* collection)
{
  if (collection == nullptr) {
    return 0;
  }

  Long64_t nMerged{0};
  TIter next(collection);
  while (auto* object = next()) {
    auto* other = dynamic_cast<ProcessingStats*>(object);
    if (other == nullptr || other == this) {
      continue;
    }
    nTF += other->nTF;
    nPV += other->nPV;
    nPVConstrAcc += other->nPVConstrAcc;
    nPVGBLAcc += other->nPVGBLAcc;
    nTrc += other->nTrc;
    nTrcGBLAcc += other->nTrcGBLAcc;
    nTrcSingleGBLAcc += other->nTrcSingleGBLAcc;
    nTrcMultiGBLAcc += other->nTrcMultiGBLAcc;
    for (const auto& [label, count] : other->nTrcByLabel) {
      nTrcByLabel[label] += count;
    }
    labelToSymName.insert(other->labelToSymName.begin(), other->labelToSymName.end());
    symNameToLabel.insert(other->symNameToLabel.begin(), other->symNameToLabel.end());
    ++nMerged;
  }
  return nMerged;
}

std::string ProcessingStats::asString() const
{
  return fmt::format("TF={}, PV={} PVConstrAcc={} PVGBLAcc={} trc={} trcGBLAcc={} trcSingleGBLAcc={} trcMultiGBLAcc={}",
                     nTF, nPV, nPVConstrAcc, nPVGBLAcc, nTrc, nTrcGBLAcc, nTrcSingleGBLAcc, nTrcMultiGBLAcc);
}

TH1F ProcessingStats::createStatHisto(const std::string& pattern) const
{
  constexpr int nCounters = 8;
  int nVolumeCounters{0};
  for (const auto& [label, count] : nTrcByLabel) {
    const auto nameIt = labelToSymName.find(label);
    if (nameIt != labelToSymName.end() && matchesVolumePattern(pattern, nameIt->second)) {
      ++nVolumeCounters;
    }
  }

  TH1F histogram("Stats", fmt::format("statistics: {}", pattern).c_str(), nCounters + nVolumeCounters, 0., nCounters + nVolumeCounters);
  histogram.SetDirectory(nullptr);
  const uint32_t counterValues[] = {nTF, nPV, nPVConstrAcc, nPVGBLAcc, nTrc, nTrcGBLAcc, nTrcSingleGBLAcc, nTrcMultiGBLAcc};
  const char* counterNames[] = {"nTF", "nPV", "nPVConstrAcc", "nPVGBLAcc", "nTrc", "nTrcGBLAcc", "nTrcSingleGBLAcc", "nTrcMultiGBLAcc"};
  int bin{1};
  for (size_t iCounter = 0; iCounter < nCounters; ++iCounter) {
    histogram.SetBinContent(bin, counterValues[iCounter]);
    histogram.GetXaxis()->SetBinLabel(bin, counterNames[iCounter]);
    ++bin;
  }
  for (const auto& [label, count] : nTrcByLabel) {
    const auto nameIt = labelToSymName.find(label);
    if (nameIt == labelToSymName.end() || !matchesVolumePattern(pattern, nameIt->second)) {
      continue;
    }
    histogram.SetBinContent(bin, count);
    histogram.GetXaxis()->SetBinLabel(bin, nameIt->second.c_str());
    ++bin;
  }
  return histogram;
}

void ProcessingStats::Print(Option_t* option) const
{
  LOGP(info, "Stats: {}", asString());
  const std::string pattern = option == nullptr ? std::string{} : std::string{option};
  if (!pattern.empty()) {
    LOGP(info, "Volume statistics for pattern '{}':", pattern);
    for (const auto& [label, count] : nTrcByLabel) {
      const auto nameIt = labelToSymName.find(label);
      if (nameIt != labelToSymName.end() && matchesVolumePattern(pattern, nameIt->second)) {
        LOGP(info, "\t{}: {}", nameIt->second, count);
      }
    }
  }
}

} // namespace o2::alignrs