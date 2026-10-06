#ifndef O2_ALIGN_PROCESSING_STATS_H
#define O2_ALIGN_PROCESSING_STATS_H

#include <TObject.h>
#include <TH1F.h>
#include <array>
#include <map>
#include <string>

#include "O2Align/Label.h"

class TCollection;

namespace o2::alignrs
{
struct ProcessingStats : public TObject {
  uint32_t nTF{0};                  // processed TFs
  uint32_t nPV{0};                  // processed PVs, excluding the unassociated-track bucket
  uint32_t nPVConstrAcc{0};         // PVs accepted by the vertex-constraint refit
  uint32_t nPVGBLAcc{0};            // PVs with a successful composed (multi-track) GBL fit
  uint32_t nTrc{0};                 // selected tracks submitted for refit and GBL processing
  uint32_t nTrcMultiPVAcc{0};       // tracks accepted by a composed multi-track PV GBL fit
  uint32_t nTrcSinglePVAcc{0};      // tracks accepted by a single-track GBL fit with PV constraint
  uint32_t nTrcNoPVAcc{0};          // tracks accepted by a single-track GBL fit without PV constraint

  std::array<uint32_t, 7> nITSLrPoints{0}; // number of ITS points per layer contribuing
  std::array<uint32_t, 7> nITSLrPointsOvl{0}; // number of ITS overlap points per layer contribuing
  std::array<uint32_t, 6> nTRDLrPoints{0}; // number of TRD points per layer contribuing
  uint32_t nTPCPoints{0}; // number of TPC points contribuing
  uint32_t nTOFPoints{0}; // number of TOF points contribuing

  std::map<int, uint32_t> nTrcByLabel; // number of tracks accepted by GBL fit, indexed by the int label of the volume they belong to
  std::map<int, std::string> labelToSymName; // optional symbolic name of the volume corresponding to the int label
  std::map<std::string, int> symNameToLabel; // optional int label of the volume corresponding to the symbolic name
  
  Long64_t Merge(TCollection* collection);

  std::string asString() const;
  TH1F createStatHisto(const std::string& pattern) const;
  void Print(Option_t* option = "") const override;

  ClassDefOverride(ProcessingStats, 3);
};
} // namespace o2::alignrs

#endif
