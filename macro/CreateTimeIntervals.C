#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>
#include <TROOT.h>
#include <TGrid.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "CommonDataFormat/TFIDInfo.h"
#include "Framework/Logger.h"
#endif

using TFIDInfo = o2::dataformats::TFIDInfo;

//_________________________________________________________________
// make sure the grid connection exists before an alien:// file is accessed
void checkGridConnection(const char* flName)
{
  if (!TString(flName).BeginsWith("alien://")) {
    return;
  }
  if (!(gGrid && gGrid->IsConnected()) && !TGrid::Connect("alien://")) {
    LOGP(fatal, "Failed to create a grid connection");
  }
}

//_________________________________________________________________
// load a chain of "tfidTree" trees from a single .root file or from a text file with 1 root file per line
TChain* loadTFIDChain(const char* inpData, const char* chName = "tfidTree")
{
  TChain* chain = new TChain(chName);
  if (!inpData) {
    return chain;
  }
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    checkGridConnection(inpData);
    LOGP(info, "Adding {}", inpData);
    chain->AddFile(inpData);
  } else {
    std::ifstream inpf(inpData);
    if (!inpf.good()) {
      LOGP(error, "Failed on input filename {}", inpData);
      delete chain;
      return nullptr;
    }
    TString flName;
    flName.ReadLine(inpf);
    while (!flName.IsNull()) {
      flName = flName.Strip(TString::kBoth, ' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {
        flName.ReadLine(inpf);
        continue;
      }
      flName = flName.Strip(TString::kBoth, ',');
      flName = flName.Strip(TString::kBoth, '"');
      checkGridConnection(flName.Data());
      LOGP(info, "Adding {}", flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  int n = chain->GetEntries();
  if (n < 1) {
    LOGP(warn, "Obtained chain is empty");
    delete chain;
    return nullptr;
  }
  LOGP(info, "Opened {} chain with {} entries", chName, n);
  return chain;
}

//_________________________________________________________________
void CreateTimeIntervals(const char* o2tfinfolist, const char* outName, int intervalSec, int startID = 0, int maxIntervalID = -1)
{
  if (intervalSec <= 0) {
    LOGP(error, "intervalSec must be positive, got {}", intervalSec);
    return;
  }
  TChain* chain = loadTFIDChain(o2tfinfolist);
  if (!chain) {
    LOGP(error, "Failed to load tfidTree chain from {}", o2tfinfolist);
    return;
  }
  TFIDInfo* tfid = nullptr;
  chain->SetBranchAddress("tfidinfo", &tfid);

  std::vector<TFIDInfo> tfids;
  long nent = chain->GetEntries();
  tfids.reserve(nent);
  for (long ient = 0; ient < nent; ient++) {
    chain->GetEntry(ient);
    if (!tfid || tfid->creation == -1UL || tfid->runNumber == -1U) { // ignore dummy/incomplete records
      continue;
    }
    tfids.push_back(*tfid);
  }
  delete chain;
  if (tfids.empty()) {
    LOGP(error, "No valid TFIDInfo records found in {}", o2tfinfolist);
    return;
  }
  // input timestamps come in random order
  std::sort(tfids.begin(), tfids.end(), [](const TFIDInfo& a, const TFIDInfo& b) {
    return a.creation == b.creation ? a.firstTForbit < b.firstTForbit : a.creation < b.creation;
  });

  const uint64_t intervalMS = uint64_t(intervalSec) * 1000;
  std::vector<std::pair<size_t, size_t>> intervals; // 1st / last index in tfids of each interval
  size_t iStart = 0;
  for (size_t i = 1; i < tfids.size(); i++) {
    // close the current interval if its duration would exceed intervalSec or if the run number changes
    if (tfids[i].creation - tfids[iStart].creation > intervalMS || tfids[i].runNumber != tfids[iStart].runNumber) {
      intervals.emplace_back(iStart, i - 1);
      iStart = i;
    }
  }
  intervals.emplace_back(iStart, tfids.size() - 1);

  int lastID = startID + int(intervals.size()) - 1;
  if (maxIntervalID >= 0 && lastID > maxIntervalID) {
    LOGP(fatal, "Reached interval ID {} exceeding maxIntervalID={}", lastID, maxIntervalID);
  }

  std::ofstream outf(outName);
  if (!outf.good()) {
    LOGP(fatal, "Failed to open output file {}", outName);
  }
  outf << "[\n";
  for (size_t iint = 0; iint < intervals.size(); iint++) {
    const auto& tfS = tfids[intervals[iint].first];
    const auto& tfE = tfids[intervals[iint].second];
    outf << "  {\n"
         << "    \"creationS\": " << tfS.creation << ",\n"
         << "    \"creationE\": " << tfE.creation << ",\n"
         << "    \"tfOrbitS\": " << tfS.firstTForbit << ",\n"
         << "    \"tfOrbitE\": " << tfE.firstTForbit << ",\n"
         << "    \"tfCountS\": " << tfS.tfCounter << ",\n"
         << "    \"tfCountE\": " << tfE.tfCounter << ",\n"
         << "    \"run\": " << tfS.runNumber << ",\n"
         << "    \"ID\": " << startID + int(iint) << "\n"
         << "  }" << (iint + 1 < intervals.size() ? "," : "") << "\n";
  }
  outf << "]\n";
  outf.close();
  LOGP(info, "Wrote {} intervals (IDs {} ... {}) of max {} s duration to {}", intervals.size(), startID, lastID, intervalSec, outName);
}
