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

#include <algorithm>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/StringUtils.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/Logger.h"
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
    {"track-sources", VariantType::String, std::string{GID::ALL}, {"comma-separated list of track sources to use, cosmic-ext for TrackCosmicExtended"}},
    {"detectors", VariantType::String, "ITS", {"comma-separated list of detectors or their hyphen-separated combinations, e.g. ITS,ITS-TRD,TOF-ITS: the union defines the detectors in the alignment (clusters to load), the tracks must be built by exactly the detectors of one of the listed combinations"}},
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
  const std::string PVDef = "pvertexer.useMeanVertexConstraint=false;pvertexer.meanVertexExtraErrSelection=0.2;pvertexer.iniScale2=100;pvertexer.acceptableScale2=10.;pvertexer.useTimeInChi2=false;";
  o2::conf::ConfigurableParam::updateFromString(PVDef); // first set proper defaults of PV refitting, the CL options will override them if provided
  o2::conf::ConfigurableParam::updateFromString(cfg.options().get<std::string>("configKeyValues"));
  const GID::mask_t allowedSourcesTrc = GID::getSourcesMask("ITS,TPC,ITS-TPC,ITS-TPC-TRD,ITS-TPC-TOF,ITS-TPC-TRD-TOF");
  const GID::mask_t allowedSourcesClus = GID::getSourcesMask("ITS,TPC,TRD,TOF");
  const auto srcTracks = cfg.options().get<std::string>("track-sources");
  bool isCosmicExt = srcTracks == "cosmic-ext";
  GID::mask_t dets;
  std::vector<DetID::mask_t> detCombinations;
  for (const auto& combo : o2::utils::Str::tokenize(cfg.options().get<std::string>("detectors"), ',')) {
    std::string comboList = combo;
    std::replace(comboList.begin(), comboList.end(), '-', ',');
    const GID::mask_t comboMask = allowedSourcesClus & GID::getSourcesMask(comboList);
    if (comboMask.none()) {
      continue;
    }
    dets |= comboMask;
    detCombinations.push_back(GID::getSourcesDetectorsMask(comboMask));
  }
  const auto useMC = !cfg.options().get<bool>("disable-mc");
  const auto withITS3 = cfg.options().get<bool>("with-its3");
  o2::alignrs::OutputEnum output(cfg.options().get<std::string>("output"));

  // rectify options
  if (output[o2::alignrs::OutputOpt::Regularize] && (output[o2::alignrs::OutputOpt::ProcessData] || output[o2::alignrs::OutputOpt::MilleData] || output[o2::alignrs::OutputOpt::VerboseGBL])) {
    LOGP(fatal, "Regularize cannot be combined with the options processing data (ProcessData, MilleData, VerboseGBL)");
  }
  if (output[o2::alignrs::OutputOpt::MilleRes] || output[o2::alignrs::OutputOpt::WeakModesDiag] || output[o2::alignrs::OutputOpt::WeakModesDiagIni] || output[o2::alignrs::OutputOpt::Regularize]) {
    output.reset(o2::alignrs::OutputOpt::ProcessData);
  }
  if (output[o2::alignrs::OutputOpt::MilleData] || output[o2::alignrs::OutputOpt::VerboseGBL]) {
    output.set(o2::alignrs::OutputOpt::ProcessData);
  }
  GID::mask_t srcTrc;
  srcTrc = allowedSourcesTrc & GID::getSourcesMask(srcTracks);

  WorkflowSpec specs;

  bool requestCTPLumi = false;
  
  if (output[o2::alignrs::OutputOpt::ProcessData]) {
    if (!isCosmicExt) {
      if (dets[GID::TPC]) { // the TPC cluster transformation needs the scalers        
        auto sclOpt = o2::tpc::CorrectionMapsOptions::parseGlobalOptions(cfg.options());
        requestCTPLumi = sclOpt.requestCTPLumi;        
        srcTrc = srcTrc | (GID::getSourcesMask("CTP") | GID::getSourcesMask("TPC"));
        specs.emplace_back(o2::tpc::getTPCScalerSpec(sclOpt));
      }    
      // the TPC refit needs the TPC tracks referred by the matches, even if standalone TPC is not processed
      const GID::mask_t srcTrcLoad = (dets[GID::TPC] && (srcTrc & GID::getSourcesMask("ITS-TPC,ITS-TPC-TRD,ITS-TPC-TOF,ITS-TPC-TRD-TOF")).any()) ? (srcTrc | GID::getSourcesMask("TPC")) : srcTrc;
      o2::globaltracking::InputHelper::addInputSpecs(cfg, specs, dets, srcTrc, srcTrcLoad, useMC);
      o2::globaltracking::InputHelper::addInputSpecsPVertex(cfg, specs, useMC);
    } else {
      o2::globaltracking::InputHelper::addInputSpecsCosmicsExtended(cfg, specs);
    }
  } else {
    specs.emplace_back(o2::globaltracking::getNoInpDummyOutSpec(0));
  }
  
  specs.emplace_back(o2::alignrs::getAlignmentSpec(srcTrc, dets, detCombinations, isCosmicExt, useMC, withITS3, requestCTPLumi, output));

  o2::raw::HBFUtilsInitializer hbfIni(cfg, specs);
  return std::move(specs);
}
