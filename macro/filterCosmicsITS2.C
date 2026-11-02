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
//   root -l -b -q 'filterCosmicsITS2.C()'
//   root -l -b -q 'filterCosmicsITS2.C("/path/input.root", 3, true)'
// The 2nd argument passes nCluster selection, 3d argument zeros covariance matrices and both TPC dEdx objects.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#endif

void filterCosmicsITS2(const char* inputName, int nITSCl = 2, bool zeroCovAndDEdx = false)
{
  using Cosmic = o2::dataformats::TrackCosmicsExtended;
  TString outputName(inputName);
  if (!outputName.EndsWith(".root")) {
    throw std::runtime_error("Input filename must end in .root");
  }
  outputName.Insert(outputName.Length() - 5, "_its2");
  std::unique_ptr<TFile> input(TFile::Open(inputName, "READ"));
  if (!input || input->IsZombie()) {
    throw std::runtime_error("Cannot open input file");
  }
  auto* tree = input->Get<TTree>("cosmicsFull");
  if (!tree || !tree->GetBranch("cosmics") || !tree->GetBranch("tfInfo")) {
    throw std::runtime_error("Expected cosmicsFull tree with cosmics and tfInfo branches");
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
  // Exclude tfID from the clone while retaining the other TF metadata.
  if (tree->GetBranch("tfID")) {
    tree->SetBranchStatus("tfID", 0);
  }
  auto* filtered = tree->CloneTree(0);
  if (!filtered) {
    throw std::runtime_error("Cannot clone cosmicsFull tree");
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
    // Fill even for an empty vector: one output entry per input TF.
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
  std::cout << "Wrote " << outputName << ": " << nTF << " TFs ("
            << nEmpty << " empty), " << nSelected << '/' << nInput
            << " tracks selected; zeroCovAndDEdx=" << zeroCovAndDEdx << '\n';
}
