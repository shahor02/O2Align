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

#include <cmath>
#include <algorithm>
#include <array>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "Framework/Logger.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITStracking/IOUtils.h"
#include "MathUtils/Utils.h"
#include "O2Align/DetectorITS.h"
#include "O2Align/Params.h"
#include "O2Align/SensorITS.h"
#include "O2Align/TrackFit.h"
#include "ITSBase/GeometryTGeo.h"

namespace o2::alignrs
{

void DetectorITS::prepareData(o2::globaltracking::RecoContainer* recoData)
{
  const auto clusITS = recoData->getITSClusters();
  const auto clusITSROF = recoData->getITSClustersROFRecords();
  const auto patterns = recoData->getITSClustersPatterns();
  auto pattIt = patterns.begin();
  const auto& params = Params::Instance();

  mITSPointsInfo.clear();
  mITSPointsInfo.reserve(clusITS.size());
  // The overlap search relies on the ITS2 chip mapping and topology dictionary. With ITS3 the IDs of
  // the outer barrel chips are shifted by the number of the IB tiles, hence it is not supported.
  const bool findOverlaps = params.ITSOverlapMargin > 0 && !mIsITS3;
  if (params.ITSOverlapMargin > 0 && mIsITS3) {
    static bool warned = false;
    if (!warned) {
      warned = true;
      LOGP(warn, "The search of the ITS overlaps is not supported with ITS3, ignoring ITSOverlapMargin={}", params.ITSOverlapMargin);
    }
  }
  if (findOverlaps && mOverlaps.empty()) { // static, compute once
    mOverlaps = o2::itsmft::ChipMappingITS{}.getOverlapsInfo();
  }
  if (findOverlaps) {
    mITSOvlClusRef.assign(clusITS.size(), -1);
    mITSOvlCandidateID.clear();
    mITSOvlCandidateID.reserve(clusITS.size());
  } else {
    mITSOvlClusRef.clear();
    mITSOvlCandidateID.clear();
  }

  auto geom = o2::its::GeometryTGeo::Instance();
  std::vector<int> edgeClusters;
  int rofCount = 0;
  struct ROFChipEntry {
    int rofCount = -1;
    int chipFirstEntry = -1;
  };
  std::array<ROFChipEntry, o2::itsmft::ChipMappingITS::getNChips()> chipROFStart{};

  for (const auto& rof : clusITSROF) {
    const int maxic = rof.getFirstEntry() + rof.getNEntries();
    edgeClusters.clear();
    for (int ic = rof.getFirstEntry(); ic < maxic; ++ic) {
      const auto& cls = clusITS[ic];
      const auto sensID = cls.getSensorID();
      const auto lay = geom->getLayer(sensID);
      float sigmaY2{0.f}, sigmaZ2{0.f};
      math_utils::Point3D<float> locXYZ;
      auto pattItCopy = pattIt;
      if (mIsITS3) {
        locXYZ = o2::its3::ioutils::extractClusterData(cls, pattIt, mIT3Dict, sigmaY2, sigmaZ2);
      } else {
        locXYZ = o2::its::ioutils::extractClusterData(cls, pattIt, mITSDict, sigmaY2, sigmaZ2);
      }
      sigmaY2 += params.extraClsErrYITS[lay] * params.extraClsErrYITS[lay];
      sigmaZ2 += params.extraClsErrZITS[lay] * params.extraClsErrZITS[lay];
      const auto gloXYZ = geom->getMatrixL2G(sensID) * locXYZ;
      auto trkXYZ = geom->getMatrixT2L(sensID) ^ locXYZ;
      double alpha = geom->getSensorRefAlpha(sensID);
      double x = geom->getSensorRefX(sensID);
      if (mIsITS3 && o2::its3::constants::detID::isDetITS3(sensID)) {
        trkXYZ.SetY(0.f);
        x = std::hypot(gloXYZ.x(), gloXYZ.y());
        trkXYZ.SetX(x);
        alpha = std::atan2(gloXYZ.y(), gloXYZ.x());
      }
      math_utils::bringToPMPid(alpha);
      o2::BaseCluster<float> clus(sensID, trkXYZ, sigmaY2, sigmaZ2, 0.f);
      auto& pointInfo = mITSPointsInfo.emplace_back();
      pointInfo.lr = lay;
      pointInfo.label = Label(mDetIdx, sensID, true);
      pointInfo.x = x;
      pointInfo.alpha = alpha;
      pointInfo.cluster = clus;
      if (findOverlaps) {
        int row = 0, col = 0;
        o2::itsmft::SegmentationAlpide::localToDetectorUnchecked(locXYZ.X(), locXYZ.Z(), row, col);
        int drow = row < o2::itsmft::SegmentationAlpide::NRows / 2 ? row : o2::itsmft::SegmentationAlpide::NRows - row - 1;
        if (drow * o2::itsmft::SegmentationAlpide::PitchRow < params.ITSOverlapMargin) {
          pointInfo.cluster.setBit(row < o2::itsmft::SegmentationAlpide::NRows / 2 ? DetectorITS::EdgeFlags::LowRow : DetectorITS::EdgeFlags::HighRow);
          if (params.ITSOverlapEdgeRows > 0) {
            auto pattID = cls.getPatternID();
            drow = cls.getRow();
            if (pattID != itsmft::CompCluster::InvalidPatternID) {
              if (!mITSDict->isGroup(pattID)) {
                const auto& patt = mITSDict->getPattern(pattID);
                if (row > o2::itsmft::SegmentationAlpide::NRows / 2) {
                  drow = o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() - 1);
                }
              } else {
                o2::itsmft::ClusterPattern patt(pattItCopy);
                drow = row < o2::itsmft::SegmentationAlpide::NRows / 2 ? drow - patt.getRowSpan() / 2 : o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() / 2 - 1);
              }
            } else {
              o2::itsmft::ClusterPattern patt(pattItCopy);
              if (row > o2::itsmft::SegmentationAlpide::NRows / 2) {
                drow = o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() - 1);
              }
            }
            if (drow < params.ITSOverlapEdgeRows) {
              pointInfo.cluster.setBit(DetectorITS::EdgeFlags::Biased);
            }
          }
          if (!pointInfo.cluster.isBitSet(DetectorITS::EdgeFlags::Biased)) {
            if (chipROFStart[sensID].rofCount != rofCount) {
              chipROFStart[sensID].rofCount = rofCount;
              chipROFStart[sensID].chipFirstEntry = edgeClusters.size();
            }
            edgeClusters.push_back(ic);
          }
        }
      }
    }
    for (auto ic : edgeClusters) {
      auto& cl = mITSPointsInfo[ic].cluster;
      const int sensID = cl.getSensorID();
      const auto ovl = mOverlaps[sensID];
      int ovlCount = 0;
      for (int ir = 0; ir < DetectorITS::OVL::NSides; ++ir) {
        if (ovl.rowSide[ir] == DetectorITS::OVL::NONE) {
          continue;
        }
        const int chipOvl = ovl.rowSide[ir];
        if (chipROFStart[chipOvl].rofCount == rofCount) {
          auto oClusID = edgeClusters[chipROFStart[chipOvl].chipFirstEntry];
          while (oClusID < static_cast<int>(mITSPointsInfo.size())) {
            const auto oClus = mITSPointsInfo[oClusID].cluster;
            if (oClus.getSensorID() != chipOvl) {
              break;
            }
            if (oClus.isBitSet(ovl.rowSideOverlap[ir]) && !oClus.isBitSet(DetectorITS::EdgeFlags::Biased) && std::abs(oClus.getZ() - cl.getZ()) < params.ITSOverlapMaxDZ) {
              if (!ovlCount) {
                mITSOvlClusRef[ic] = mITSOvlCandidateID.size();
              }
              mITSOvlCandidateID.push_back(oClusID);
              ++ovlCount;
            }
            ++oClusID;
          }
        }
      }
      cl.setCount(std::min(127, ovlCount));
    }
    ++rofCount;
  }
}

bool DetectorITS::prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack)
{
  const auto& params = Params::Instance();
  const bool allowOverlaps = params.ITSOverlapMaxChi2 > 0.f;
  const auto gidITS = ids[GTrackID::ITS];
  const auto gidITSAB = ids[GTrackID::ITSAB];
  const bool useITS = gidITS.isIndexSet();
  const bool useITSAB = !useITS && gidITSAB.isIndexSet();
  if (!useITS && !useITSAB) {
    return false;
  }
  const auto* itsTrack = useITS ? &recoData->getITSTrack(gidITS) : nullptr;
  auto trFitOut = itsTrack ? convertTrack<double>(itsTrack->getParamIn()) : convertTrack<double>(recoData->getTrackParam(resTrack.gid));
  auto trFitInw = itsTrack ? convertTrack<double>(itsTrack->getParamOut()) : convertTrack<double>(recoData->getTrackParam(resTrack.gid));
  auto prop = o2::base::PropagatorD::Instance();
  std::array<FrameInfoExt*, 7> frameArr{};
  std::array<FrameInfoExt*, 7> overlapArr{};
  std::array<int, 7> clusIDArr{};
  clusIDArr.fill(-1);
  std::array<TrackD, 7> trkOutAt{};
  std::array<bool, 7> hasTrkOutAt{};

  auto findBestOverlap = [&](const FrameInfoExt& frame, int clusID, const TrackD& tr) -> FrameInfoExt* {
    if (clusID < 0 || !frame.cluster.getCount() || mITSOvlClusRef.empty()) {
      return nullptr;
    }
    int clusIDtoCheck = mITSOvlClusRef[clusID];
    if (clusIDtoCheck < 0) {
      return nullptr;
    }
    int bestClusID = -1;
    float bestChi2 = params.ITSOverlapMaxChi2;
    for (int iov = 0; iov < frame.cluster.getCount(); ++iov) {
      if (clusIDtoCheck >= static_cast<int>(mITSOvlCandidateID.size())) {
        break;
      }
      const int clusOvlID = mITSOvlCandidateID[clusIDtoCheck++];
      auto trPropOvl = tr;
      const auto& frameOvl = mITSPointsInfo[clusOvlID];
      if (!prop->propagateToAlphaX(trPropOvl, nullptr, frameOvl.alpha, frameOvl.x, false, params.maxSnp, params.maxStep, 1, o2::base::PropagatorD::MatCorrType::USEMatCorrNONE)) {
        continue;
      }
      const auto ovlChi2 = static_cast<float>(trPropOvl.getPredictedChi2Quiet(frameOvl.cluster));
      if (ovlChi2 < bestChi2) {
        bestChi2 = ovlChi2;
        bestClusID = clusOvlID;
      }
    }
    return bestClusID >= 0 ? &mITSPointsInfo[bestClusID] : nullptr;
  };

  int nPoints = 0;
  auto registerCluster = [&](int clusID) {
    auto& curInfo = mITSPointsInfo[clusID];
    if (curInfo.cluster.getBits() && curInfo.cluster.isBitSet(DetectorITS::EdgeFlags::Biased)) {
      return;
    }
    frameArr[curInfo.lr] = &curInfo;
    clusIDArr[curInfo.lr] = clusID;
    ++nPoints;
  };

  if (useITS) {
    const auto itsClRefs = recoData->getITSTracksClusterRefs();
    const int nCl = itsTrack->getNClusters();
    for (int i = 0; i < nCl; i++) { // clusters are ordered from the outermost to the innermost
      registerCluster(itsClRefs[itsTrack->getClusterEntry(i)]);
    }
  } else {
    const auto& trkITSABref = recoData->getITSABRef(gidITSAB);
    const auto abTrackClusIdx = recoData->getITSABClusterRefs();
    const int nCl = trkITSABref.getNClusters();
    const int clEntry = trkITSABref.getFirstEntry();
    for (int icl = 0; icl < nCl; ++icl) { // ITSAB clusters are stored from inner to outer layers
      registerCluster(abTrackClusIdx[clEntry + icl]);
    }
  }
  if (nPoints < (useITS ? params.minITSCls : 2)) {
    return false;
  }

  if (allowOverlaps) {
    resetTrackCovariance(trFitOut);
    resetTrackCovariance(trFitInw);
    o2::track::TrackParD trkOutRef, *refLinOut = nullptr;
    if (params.useStableRef) {
      refLinOut = &(trkOutRef = trFitOut);
    }
    for (int i = 0; i < 7; ++i) {
      if (!frameArr[i]) {
        continue;
      }
      if (!prop->propagateToAlphaX(trFitOut, refLinOut, frameArr[i]->alpha, frameArr[i]->x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
        return false;
      }
      trkOutAt[i] = trFitOut;
      hasTrkOutAt[i] = true;
      if (!trFitOut.update(frameArr[i]->cluster)) {
        return false;
      }
      if (refLinOut) {
        refLinOut->setY(frameArr[i]->cluster.getY());
        refLinOut->setZ(frameArr[i]->cluster.getZ());
      }
    }
    o2::track::TrackParD trkInwRef, *refLinInw = nullptr;
    if (params.useStableRef) {
      refLinInw = &(trkInwRef = trFitInw);
    }
    for (int i = 7; i--;) {
      if (!frameArr[i]) {
        continue;
      }
      if (!prop->propagateToAlphaX(trFitInw, refLinInw, frameArr[i]->alpha, frameArr[i]->x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
        return false;
      }
      if (hasTrkOutAt[i] && frameArr[i]->cluster.getCount()) {
        auto smooth = interpolateTrackParCov<double>(trkOutAt[i], trFitInw);
        if (!smooth.isValid()) {
          return false;
        }
        if (!smooth.update(frameArr[i]->cluster)) {
          return false;
        }
        overlapArr[i] = findBestOverlap(*frameArr[i], clusIDArr[i], smooth);
      }
      if (!trFitInw.update(frameArr[i]->cluster)) {
        return false;
      }
      if (refLinInw) {
        refLinInw->setY(frameArr[i]->cluster.getY());
        refLinInw->setZ(frameArr[i]->cluster.getZ());
      }
    }
  }
  const size_t nFramesIni = resTrack.info.size(); // N frames before adding new ITS frames
  for (int i = 0; i < 7; ++i) {
    if (!frameArr[i]) {
      continue;
    }
    resTrack.info.push_back(*frameArr[i]);
    if (overlapArr[i]) {
      resTrack.info.push_back(*overlapArr[i]);
    }
  }
  // the overlap frames are interleaved with their host layer frames, order them along the track
  std::stable_sort(resTrack.info.begin() + nFramesIni, resTrack.info.end(), [](const auto& a, const auto& b) {
    return a.x < b.x;
  });
  // continue the fit of the seed prepared in advance outward, w/o resetting its cov.matrix
  return resTrack.continueFitOutward(static_cast<int>(nFramesIni));
}


// ITS2 hierarchy: ITS -> half-barrels -> staves -> half-staves -> modules -> chips. With ITS3 the
// half-barrels of the 3 inner layers are the sensors themselves, with the tiles as pseudo children
// carrying the measurements.
Volume::Ptr DetectorITS::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  uint32_t gLbl{0};
  const uint32_t det = mDetIdx;
  auto geom = o2::its::GeometryTGeo::Instance();
  auto isLayITS3 = [this](int lr) { return mIsITS3 && lr < 3; };
  std::unordered_map<std::string, Volume*> sym2vol;
  auto addVolume = [&](Volume* parent, const char* symName) {
    auto* vol = parent->addChild(symName, gLbl++, det, false);
    sym2vol[vol->getSymName()] = vol;
    return vol;
  };
  auto root = std::make_unique<Volume>(geom->composeSymNameITS(), gLbl++, det, false);
  sym2vol[root->getSymName()] = root.get();
  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {
    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      auto* volHB = addVolume(root.get(), geom->composeSymNameHalfBarrel(ilr, ihb, isLayITS3(ilr)));
      if (isLayITS3(ilr)) {
        volHB->setSensorId((2 * ilr) + ihb);
        continue;
      }
      for (int ist = 0; ist < geom->getNumberOfStaves(ilr) / 2; ist++) {
        auto* volSt = addVolume(volHB, geom->composeSymNameStave(ilr, ihb, ist));
        for (int ihst = 0; ihst < geom->getNumberOfHalfStaves(ilr); ihst++) {
          auto* volHSt = addVolume(volSt, geom->composeSymNameHalfStave(ilr, ihb, ist, ihst));
          for (int imd = 0; imd < geom->getNumberOfModules(ilr); imd++) {
            addVolume(volHSt, geom->composeSymNameModule(ilr, ihb, ist, ihst, imd));
          }
        }
      }
    }
  }
  int lay = 0, hba = 0, sta = 0, ssta = 0, modd = 0, chip = 0;
  for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
    geom->getChipId(ich, lay, hba, sta, ssta, modd, chip);
    const Label lbl(det, ich, true);
    const std::string parSym = isLayITS3(lay) ? geom->composeSymNameHalfBarrel(lay, hba, true) : (modd < 0 ? geom->composeSymNameStave(lay, hba, sta) : geom->composeSymNameModule(lay, hba, sta, ssta, modd));
    const auto parIt = sym2vol.find(parSym);
    if (parIt == sym2vol.end()) {
      LOGP(fatal, "did not find parent {} for chip {}", parSym, ich);
    }
    Volume* sensor{nullptr};
    if (isLayITS3(lay)) {
      sensor = parIt->second->addChild<SensorIT3>(geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip, true), lbl);
      sensor->setPseudo(true);
    } else {
      const int nch = modd < 0 ? geom->getNumberOfChipsPerStave(lay) : geom->getNumberOfChipsPerModule(lay);
      sensor = parIt->second->addChild<SensorITS>(geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip % nch), lbl);
    }
    sensor->setSensorId(ich);
    sensorMap[lbl] = sensor;
  }
  return root;
}

namespace
{
/// free-RB / free-calib / total counts of one hierarchy kind (half-barrel, stave, ...) at a layer
struct DOFCounts {
  int nRB{0}, nCal{0}, nTot{0};
  void add(const Volume* v)
  {
    if (v == nullptr) {
      return;
    }
    ++nTot;
    nRB += v->getRigidBody() && v->getRigidBody()->nFreeDOFs() > 0;
    nCal += v->getCalib() && v->getCalib()->nFreeDOFs() > 0;
  }
};

std::ostream& operator<<(std::ostream& os, const DOFCounts& c)
{
  return os << c.nRB << '/' << c.nCal << '/' << c.nTot;
}

/// comma-separated names of the free DOFs of a single volume, "none" if it has none
std::string freeDOFNames(const Volume* v)
{
  std::string s;
  auto append = [&s](const char* tag, const DOFSet* d) {
    if (d == nullptr) {
      return;
    }
    for (int i = 0; i < d->nDOFs(); ++i) {
      if (d->isFree(i)) {
        if (!s.empty()) {
          s += ",";
        }
        s += tag;
        s += d->dofName(i);
      }
    }
  };
  append("", v->getRigidBody());
  append("CAL:", v->getCalib());
  return s.empty() ? "none" : s;
}

} // namespace

std::string DetectorITS::reportDOFSummary() const
{
  if (mTopVolume == nullptr) {
    return {};
  }
  auto geom = o2::its::GeometryTGeo::Instance();
  auto isLayITS3 = [this](int lr) { return mIsITS3 && lr < 3; };
  std::unordered_map<std::string, Volume*> sym2vol;
  mTopVolume->traverse([&sym2vol](Volume* v) { sym2vol[v->getSymName()] = v; });
  auto find = [&sym2vol](const std::string& symName) -> Volume* {
    const auto it = sym2vol.find(symName);
    return it == sym2vol.end() ? nullptr : it->second;
  };

  std::ostringstream oss;
  oss << "ITSenvelope: " << freeDOFNames(mTopVolume) << '\n';

  int lay = 0, hba = 0, sta = 0, ssta = 0, modd = 0, chip = 0;
  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {
    const bool isIB3 = isLayITS3(ilr);
    DOFCounts halfBarrels, staves, halfStaves, modules, sensors;
    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      halfBarrels.add(find(geom->composeSymNameHalfBarrel(ilr, ihb, isIB3)));
      if (isIB3) {
        continue; // the half-barrel IS the sensor; its tiles are counted below as "sensors"
      }
      for (int ist = 0; ist < geom->getNumberOfStaves(ilr) / 2; ist++) {
        staves.add(find(geom->composeSymNameStave(ilr, ihb, ist)));
        for (int ihst = 0; ihst < geom->getNumberOfHalfStaves(ilr); ihst++) {
          halfStaves.add(find(geom->composeSymNameHalfStave(ilr, ihb, ist, ihst)));
          for (int imd = 0; imd < geom->getNumberOfModules(ilr); imd++) {
            modules.add(find(geom->composeSymNameModule(ilr, ihb, ist, ihst, imd)));
          }
        }
      }
    }
    for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
      geom->getChipId(ich, lay, hba, sta, ssta, modd, chip);
      if (lay != ilr) {
        continue;
      }
      const int nch = isIB3 ? -1 : (modd < 0 ? geom->getNumberOfChipsPerStave(lay) : geom->getNumberOfChipsPerModule(lay));
      const std::string symChip = isIB3 ? geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip, true) : geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip % nch);
      sensors.add(find(symChip));
    }
    oss << "Layer " << ilr << ": halfbarrels: " << halfBarrels;
    if (!isIB3) {
      oss << ", staves: " << staves << ", halfstaves: " << halfStaves << ", modules: " << modules;
    }
    oss << ", sensors: " << sensors << '\n';
  }
  return oss.str();
}

} // namespace o2::alignrs
