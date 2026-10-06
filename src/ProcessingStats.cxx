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
    nTrcMultiPVAcc += other->nTrcMultiPVAcc;
    nTrcSinglePVAcc += other->nTrcSinglePVAcc;
    nTrcNoPVAcc += other->nTrcNoPVAcc;
    for (size_t iLayer = 0; iLayer < nITSLrPoints.size(); ++iLayer) {
      nITSLrPoints[iLayer] += other->nITSLrPoints[iLayer];
      nITSLrPointsOvl[iLayer] += other->nITSLrPointsOvl[iLayer];
    }
    for (size_t iLayer = 0; iLayer < nTRDLrPoints.size(); ++iLayer) {
      nTRDLrPoints[iLayer] += other->nTRDLrPoints[iLayer];
    }
    nTPCPoints += other->nTPCPoints;
    nTOFPoints += other->nTOFPoints;
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
  return fmt::format("TF={}, PV={} PVConstrAcc={} PVGBLAcc={} trc={} trcMultiPVAcc={} trcSinglePVAcc={} trcNoPVAcc={}",
                     nTF, nPV, nPVConstrAcc, nPVGBLAcc, nTrc, nTrcMultiPVAcc, nTrcSinglePVAcc, nTrcNoPVAcc);
}

TH1F ProcessingStats::createStatHisto(const std::string& pattern) const
{
  constexpr int nBaseCounters = 8;
  constexpr int nDetectorCounters = 7 + 7 + 6 + 2;
  int nVolumeCounters{0};
  for (const auto& [label, count] : nTrcByLabel) {
    const auto nameIt = labelToSymName.find(label);
    if (nameIt != labelToSymName.end() && matchesVolumePattern(pattern, nameIt->second)) {
      ++nVolumeCounters;
    }
  }

  constexpr int nCounters = nBaseCounters + nDetectorCounters;
  TH1F histogram("Stats", fmt::format("statistics: {}", pattern).c_str(), nCounters + nVolumeCounters, 0., nCounters + nVolumeCounters);
  histogram.SetDirectory(nullptr);
  const uint32_t counterValues[] = {nTF, nPV, nPVConstrAcc, nPVGBLAcc, nTrc, nTrcMultiPVAcc, nTrcSinglePVAcc, nTrcNoPVAcc};
  const char* counterNames[] = {"nTF", "nPV", "nPVConstrAcc", "nPVGBLAcc", "nTrc", "nTrcMultiPVAcc", "nTrcSinglePVAcc", "nTrcNoPVAcc"};
  int bin{1};
  for (size_t iCounter = 0; iCounter < nBaseCounters; ++iCounter) {
    histogram.SetBinContent(bin, counterValues[iCounter]);
    histogram.GetXaxis()->SetBinLabel(bin, counterNames[iCounter]);
    ++bin;
  }
  auto addArrayCounters = [&histogram, &bin](const char* name, const auto& counters) {
    for (size_t iElement = 0; iElement < counters.size(); ++iElement) {
      histogram.SetBinContent(bin, counters[iElement]);
      histogram.GetXaxis()->SetBinLabel(bin, fmt::format("{}[{}]", name, iElement).c_str());
      ++bin;
    }
  };
  addArrayCounters("nITSLrPoints", nITSLrPoints);
  addArrayCounters("nITSLrPointsOvl", nITSLrPointsOvl);
  addArrayCounters("nTRDLrPoints", nTRDLrPoints);
  histogram.SetBinContent(bin, nTPCPoints);
  histogram.GetXaxis()->SetBinLabel(bin++, "nTPCPoints");
  histogram.SetBinContent(bin, nTOFPoints);
  histogram.GetXaxis()->SetBinLabel(bin++, "nTOFPoints");
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
  LOGP(info, "Detector points: TPC={} TOF={}", nTPCPoints, nTOFPoints);
  for (size_t iLayer = 0; iLayer < nITSLrPoints.size(); ++iLayer) {
    LOGP(info, "\tITS layer {}: points={} overlap={}", iLayer, nITSLrPoints[iLayer], nITSLrPointsOvl[iLayer]);
  }
  for (size_t iLayer = 0; iLayer < nTRDLrPoints.size(); ++iLayer) {
    LOGP(info, "\tTRD layer {}: points={}", iLayer, nTRDLrPoints[iLayer]);
  }
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
