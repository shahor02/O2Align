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

#ifndef O2_ALIGN_DETECTORTPC_H
#define O2_ALIGN_DETECTORTPC_H

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "DataFormatsTPC/VDriftCorrFact.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::gpu
{
class GPUParam;
class TPCFastTransformPOD;
} // namespace o2::gpu

namespace o2::alignrs
{

class DetectorTPC final : public Detector
{
 public:
  static constexpr int NStacks = 4; // IROC, OROC1, OROC2, OROC3

  DetectorTPC() : Detector(DetTPC) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;

  /// both must be provided by the caller before the 1st prepareTrack call
  void setCorrMaps(const o2::gpu::TPCFastTransformPOD* maps) { mCorrMaps = maps; }
  void setTPCParam(const o2::gpu::GPUParam* par) { mTPCParam = par; }
  const o2::gpu::TPCFastTransformPOD* getCorrMaps() const { return mCorrMaps; }
  const o2::gpu::GPUParam* getTPCParam() const { return mTPCParam; }

  /// The drift calibration already applied to the correction maps: the fitted VDRIFT correction is
  /// relative to it, hence it is needed to express the result wrt the reference drift velocity.
  void setVDrift(const o2::tpc::VDriftCorrFact& v) { mVDrift = v; }
  const o2::tpc::VDriftCorrFact& getVDrift() const { return mVDrift; }

  const std::string& getTimeSlotsJson() const final;

  /// Millepede label of the drift calibration DOFs of a given time slot. The drift is a property of
  /// the whole detector, hence the DOFs belong to the envelope and the slot ID plays the role of
  /// its volume ID: every slot is calibrated independently.
  static Label getVDriftLabel(int slotID = 0) { return Label(DetTPC, slotID, false, true); }
  Volume* getEnvelopeVolume() const { return mEnvelope; }

  /// Hand the TPC drift length to the drift calibration DOFs, which need it to convert the cluster Z
  /// into the drift length. Must be called once the correction maps are set and the DOF
  /// configuration is applied, i.e. after buildHierarchy.
  void finaliseCalib();

  /// index of the stack the padrow belongs to, -1 if the padrow is out of range
  static int getStack(int padrow)
  {
    for (int i = 0; i < NStacks; ++i) {
      if (padrow <= sStackMinMaxRow[i].second) {
        return i;
      }
    }
    return -1;
  }

  /// distance (in padrows) to the closest min/max padrow of the stack containing the padrow
  static int getDistanceToStackEdge(int padrow)
  {
    const int st = getStack(padrow);
    if (st < 0) {
      return -999;
    }
    return std::min(padrow - sStackMinMaxRow[st].first, sStackMinMaxRow[st].second - padrow);
  }

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;
  void onSlotChange(int slotID) final;

 private:
  static constexpr std::array<std::pair<int, int>, NStacks> sStackMinMaxRow = {
    std::pair<int, int>{0, 62}, std::pair<int, int>{63, 96}, std::pair<int, int>{97, 126}, std::pair<int, int>{127, 151}};

  const o2::gpu::TPCFastTransformPOD* mCorrMaps{nullptr};
  const o2::gpu::GPUParam* mTPCParam{nullptr};
  o2::tpc::VDriftCorrFact mVDrift{}; // drift calibration accounted by mCorrMaps
  Volume* mEnvelope{nullptr};    // fictitious top volume, owner of the drift calibration DOFs
  std::vector<Volume*> mSensors; // sector volumes, indexed by sector
};

} // namespace o2::alignrs

#endif
