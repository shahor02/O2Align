// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Filter TrackCosmicsExtended tree selecting only tracks with >= nITSCL ITS cluster candidates,
// Drop tfID branch as it is identical to the DPL TimingInfo 
// Run in an O2 environment:
//   root -l -b -q 'skimCosmExt.C()'
//   root -l -b -q 'skimCosmExt.C("_its2", 3, true, true)'
// Arguments: output suffix (empty replaces inputs), minimum ITS clusters,
// zero covariance matrices and both TPC dEdx objects, suppress empty TFs.
// Inputs in the current directory: o2_cosmics_full.root and optional o2_tfidinfo.root.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include "CommonDataFormat/TFIDInfo.h"
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#endif

void skimCosmExt(const char* outSuffix = "", unsigned int nITSCl = 3, bool noEmpty = true, bool zeroCovAndDEdx = false)
{
  using Cosmic = o2::dataformats::TrackCosmicsExtended;
  const char* inputName = "o2_cosmics_full.root";
  const char* tfidName = "o2_tfidinfo.root";
  const bool replaceInputs = !outSuffix || !*outSuffix;
  auto makeOutputName = [&](const char* name) {
    TString result(name);
    result.Insert(result.Length() - 5, replaceInputs ? TString::Format("_tmp_%d", gSystem->GetPid()).Data() : outSuffix);
    return result;
  };
  const auto outputName = makeOutputName(inputName);
  const auto tfidOutputName = makeOutputName(tfidName);
  std::unique_ptr<TFile> input(TFile::Open(inputName, "READ"));
  if (!input || input->IsZombie()) {
    throw std::runtime_error("Cannot open input file");
  }
  auto* tree = input->Get<TTree>("cosmicsFull");
  if (!tree || !tree->GetBranch("cosmics") || !tree->GetBranch("tfInfo")) {
    throw std::runtime_error("Expected cosmicsFull tree with cosmics and tfInfo branches");
  }
  std::unique_ptr<TFile> tfidInput, tfidOutput;
  TTree* tfidTree = nullptr;
  TTree* filteredTFID = nullptr;
  std::vector<o2::dataformats::TFIDInfo>* tfidInfo = nullptr;
  std::vector<o2::dataformats::TFIDInfo> filteredTFIDInfo;
  if (!gSystem->AccessPathName(tfidName)) {
    tfidInput.reset(TFile::Open(tfidName, "READ"));
    if (!tfidInput || tfidInput->IsZombie()) {
      throw std::runtime_error("Cannot open o2_tfidinfo.root");
    }
    tfidTree = tfidInput->Get<TTree>("tfidTree");
    tfidInput->GetObject("tfidinfo", tfidInfo);
    if (!tfidInfo && !tfidTree) {
      throw std::runtime_error("Expected tfidinfo vector and/or tfidTree in o2_tfidinfo.root");
    }
    if (tfidTree && !tfidTree->GetBranch("tfidinfo")) {
      throw std::runtime_error("Expected tfidinfo branch in tfidTree");
    }
    // Each available representation must have one record per input TF.
    if ((tfidTree && tfidTree->GetEntries() != tree->GetEntries()) ||
        (tfidInfo && tfidInfo->size() != static_cast<size_t>(tree->GetEntries()))) {
      throw std::runtime_error("TFID vector, TFID tree and cosmic tree must have matching entry counts");
    }
    if (tfidInfo) {
      filteredTFIDInfo.reserve(tfidInfo->size());
    }
    tfidOutput.reset(TFile::Open(tfidOutputName.Data(), "CREATE"));
    if (!tfidOutput || tfidOutput->IsZombie()) {
      throw std::runtime_error("Cannot create TFID output file (it may already exist)");
    }
    tfidOutput->cd();
    if (tfidTree) {
      filteredTFID = tfidTree->CloneTree(0);
      if (!filteredTFID) {
        throw std::runtime_error("Cannot clone tfidTree");
      }
    }
  }
  std::vector<Cosmic> buffer;
  auto* cosmics = &buffer;
  if (tree->SetBranchAddress("cosmics", &cosmics) < 0) {
    throw std::runtime_error("Cannot bind vector<TrackCosmicsExtended> branch");
  }
  // CREATE prevents accidental replacement of an existing filtered file.
  std::unique_ptr<TFile> output(TFile::Open(outputName.Data(), "CREATE"));
  if (!output || output->IsZombie()) {
    tree->ResetBranchAddresses();
    throw std::runtime_error("Cannot create output file (it may already exist)");
  }
  output->cd();
  // Build the retained schema explicitly. CloneTree can leave a split
  // tfID parent with dangling child pointers after pruning disabled leaves.
  o2::dataformats::CosmicsTFInfo tfInfoBuffer;
  auto* tfInfo = &tfInfoBuffer;
  if (tree->SetBranchAddress("tfInfo", &tfInfo) < 0) {
    throw std::runtime_error("Cannot bind CosmicsTFInfo branch");
  }
  auto* filtered = new TTree(tree->GetName(), tree->GetTitle());
  if (!filtered->Branch("cosmics", &cosmics) || !filtered->Branch("tfInfo", &tfInfo)) {
    throw std::runtime_error("Cannot create filtered cosmic branches");
  }
  Long64_t nInput = 0, nSelected = 0, nEmpty = 0;
  const auto nTF = tree->GetEntries();
  for (Long64_t entry = 0; entry < nTF; ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      throw std::runtime_error("Failed to read TF " + std::to_string(entry));
    }
    nInput += cosmics->size();
    std::vector<int> newIndex(cosmics->size(), -1);
    int count = 0;
    for (size_t i = 0; i < cosmics->size(); ++i) {
      if ((*cosmics)[i].clITS.size() >= nITSCl) {
        newIndex[i] = count++;
      }
    }
    // Compact in place, retaining the complete selected objects, especially
    // itsPatterns and clITS.pattEntry, raw TPC clusters, times and hit flags.
    for (size_t i = 0; i < newIndex.size(); ++i) {
      if (newIndex[i] < 0) {
        continue;
      }
      const auto dst = static_cast<size_t>(newIndex[i]);
      if (dst != i) {
        (*cosmics)[dst] = std::move((*cosmics)[i]);
      }
      auto& c = (*cosmics)[dst];
      const int original = c.duplicateOf;
      c.duplicateOf = original >= 0 && static_cast<size_t>(original) < newIndex.size() ? newIndex[original] : -1;
      if (zeroCovAndDEdx) {
        // setCov({}) zeros every element without touching track parameters.
        c.cosmic.setCov({});
        c.cosmic.getParamOut().setCov({});
        c.polished.setCov({});
        for (auto* tpc : {&c.tpcBottom, &c.tpcTop}) {
          tpc->setCov({});
          tpc->getParamOut().setCov({});
          tpc->setdEdx(o2::tpc::dEdxInfo{});
          tpc->setdEdxAlt(o2::tpc::dEdxInfo{});
        }
      }
    }
    cosmics->resize(count);
    nSelected += count;
    nEmpty += count == 0;
    if (noEmpty && count == 0) {
      continue;
    }
    if (tfidInfo) {
      filteredTFIDInfo.push_back((*tfidInfo)[entry]);
    }
    if (tfidTree) {
      if (tfidTree->GetEntry(entry) <= 0 || filteredTFID->Fill() < 0) {
        throw std::runtime_error("Failed to copy TFID entry " + std::to_string(entry));
      }
    }
    if (filtered->Fill() < 0) {
      throw std::runtime_error("Failed to write TF " + std::to_string(entry));
    }
  }
  output->cd();
  if (filtered->Write() <= 0) {
    throw std::runtime_error("Failed to write filtered tree");
  }
  filtered->ResetBranchAddresses();
  tree->ResetBranchAddresses();
  output->Close();
  if (tfidOutput) {
    tfidOutput->cd();
    if (tfidInfo && tfidOutput->WriteObjectAny(&filteredTFIDInfo, "std::vector<o2::dataformats::TFIDInfo>", "tfidinfo") <= 0) {
      throw std::runtime_error("Failed to write filtered TFID vector");
    }
    if (filteredTFID && filteredTFID->Write() <= 0) {
      throw std::runtime_error("Failed to write filtered TFID tree");
    }
    if (tfidTree) {
      filteredTFID->ResetBranchAddresses();
      tfidTree->ResetBranchAddresses();
    }
    tfidOutput->Close();
    tfidInput->Close();
  }
  input->Close();
  if (replaceInputs) {
    auto replaceFile = [](const TString& temporary, const char* original) {
      if (gSystem->CopyFile(temporary.Data(), original, true) != 0) {
        throw std::runtime_error("Failed to replace " + std::string(original) + "; output retained in " + temporary.Data());
      }
      gSystem->Unlink(temporary.Data());
    };
    replaceFile(outputName, inputName);
    if (tfidOutput) {
      replaceFile(tfidOutputName, tfidName);
    }
  }
  std::cout << " Folder: " << gSystem->pwd() << " Wrote " << (replaceInputs ? TString(inputName) : outputName) << ": " << (noEmpty ? nTF - nEmpty : nTF) << " TFs ("
            << nEmpty << (noEmpty ? " empty suppressed), " : " empty retained), ") << nSelected << '/' << nInput
            << " tracks selected; zeroCovAndDEdx=" << zeroCovAndDEdx;
  if (tfidOutput) {
    std::cout << " Wrote " << (replaceInputs ? TString(tfidName) : tfidOutputName);
  }
  std::cout << "\n";
}
