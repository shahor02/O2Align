// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CallbacksPolicy.h"
#include "GlobalTrackingWorkflowHelpers/InputHelper.h"
#include "GlobalTrackingWorkflowHelpers/NoInpDummyOutSpec.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "O2Align/AlignmentSpec.h"
#include "TPCCalibration/CorrectionMapsOptions.h"
#include "TPCWorkflow/TPCScalerSpec.h"

using namespace o2::framework;
using GID = o2::dataformats::GlobalTrackID;
using DetID = o2::detectors::DetID;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  std::vector<o2::framework::ConfigParamSpec> options{
    {"disable-mc", o2::framework::VariantType::Bool, false, {"enable MC propagation"}},
    {"track-sources", VariantType::String, std::string{GID::ALL}, {"comma-separated list of track sources to use"}},
    {"detectors", VariantType::String, "ITS", {"comma-separated list of detectors in the alignment (clusters to load)"}},
    {"with-its3", VariantType::Bool, false, {"ITS3 alignment mode"}},
    {"output", VariantType::String, "", {"output steering"}},
    {"disable-root-input", VariantType::Bool, false, {"disable root-files input reader"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings ..."}}};
  o2::raw::HBFUtilsInitializer::addConfigOption(options);
  o2::tpc::CorrectionMapsOptions::addGlobalOptions(options);
  o2::itsmft::DPLAlpideParamInitializer::addITSConfigOption(options);
  std::swap(workflowOptions, options);
}
#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  o2::conf::ConfigurableParam::updateFromString(cfg.options().get<std::string>("configKeyValues"));
  const GID::mask_t allowedSourcesTrc = GID::getSourcesMask("ITS,TPC,ITS-TPC,ITS-TPC-TRD,ITS-TPC-TOF,ITS-TPC-TRD-TOF");
  const GID::mask_t allowedSourcesClus = GID::getSourcesMask("ITS,TPC,TRD,TOF");
  GID::mask_t srcTrc = allowedSourcesTrc & GID::getSourcesMask(cfg.options().get<std::string>("track-sources"));
  GID::mask_t dets = allowedSourcesClus & GID::getSourcesMask(cfg.options().get<std::string>("detectors"));
  const auto useMC = !cfg.options().get<bool>("disable-mc");
  const auto withITS3 = cfg.options().get<bool>("with-its3");
  const o2::alignrs::OutputEnum output(cfg.options().get<std::string>("output"));

  WorkflowSpec specs;

  bool requestCTPLumi = false;
  
  if (!output[o2::alignrs::OutputOpt::MilleRes]) {
    if (dets[GID::TPC]) { // the TPC cluster transformation needs the scalers
      auto sclOpt = o2::tpc::CorrectionMapsOptions::parseGlobalOptions(cfg.options());
      requestCTPLumi = sclOpt.requestCTPLumi;
      srcTrc = srcTrc | GID::getSourcesMask("CTP");
      specs.emplace_back(o2::tpc::getTPCScalerSpec(sclOpt));
    }    
    o2::globaltracking::InputHelper::addInputSpecs(cfg, specs, dets, srcTrc, srcTrc, useMC);
    o2::globaltracking::InputHelper::addInputSpecsPVertex(cfg, specs, useMC);
  } else {
    specs.emplace_back(o2::globaltracking::getNoInpDummyOutSpec(0));
  }
  
  specs.emplace_back(o2::alignrs::getAlignmentSpec(srcTrc, dets, useMC, withITS3, requestCTPLumi, output));

  o2::raw::HBFUtilsInitializer hbfIni(cfg, specs);
  return std::move(specs);
}
