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

#ifndef O2_ALIGN_SPEC_H
#define O2_ALIGN_SPEC_H

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Framework/DataProcessorSpec.h"
#include "CommonUtils/EnumFlags.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include <vector>

namespace o2::alignrs
{

enum class OutputOpt : uint16_t {
  VerboseGBL = 0,
  ProcessData, // expect data input for processing
  MilleData,   // produce mille data (imposes ProcessData)
  MilleSteer,  // extract mille steering files (standalone or in the beginning of ProcessData)
  MilleRes,    // convert pede results to CCDB objects
  MisRes,
  WeakModesDiag, // with MilleRes: diagnose the global weak modes of the fitted corrections
  WeakModesDiagIni, // standalone: diagnose the global weak modes already in the initial (CCDB) alignment
  Regularize,       // write <DET>_<algParamsOutFile> with the common motion of the children transferred to their parents: of the MilleRes result, or else of the initial (CCDB) alignment
  Debug,
};
using OutputEnum = utils::EnumFlags<OutputOpt>;

o2::framework::DataProcessorSpec getAlignmentSpec(o2::dataformats::GlobalTrackID::mask_t srcTracks, o2::dataformats::GlobalTrackID::mask_t detectors, const std::vector<o2::detectors::DetID::mask_t>& detCombinations, bool cosmicsExt, bool useMC, bool withITS3, bool requestCTPLumi, o2::alignrs::OutputEnum out);
} // namespace o2::alignrs

#endif
