// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include <cmath>
#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Framework/Logger.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITStracking/IOUtils.h"
#include "O2Align/DetectorITS.h"
#include "O2Align/Params.h"
#include "O2Align/SensorITS.h"
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
  mOverlaps = o2::itsmft::ChipMappingITS{}.getOverlapsInfo();
  if (params.ITSOverlapMargin > 0) {
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
      pointInfo.x = x;
      pointInfo.alpha = alpha;
      pointInfo.cluster = clus;
      if (params.ITSOverlapMargin > 0 && (!mIsITS3 || lay > 2)) {
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
            if (oClus.getSensorID() != sensID) {
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

Volume::Ptr DetectorITS::buildHierarchyITS(Volume::SensorMapping& sensorMap)
{
  uint32_t gLbl{0}, det{0};
  auto geom = o2::its::GeometryTGeo::Instance();
  Volume *volHB{nullptr}, *volSt{nullptr}, *volHSt{nullptr}, *volMod{nullptr};
  std::unordered_map<std::string, Volume*> sym2vol;
  auto root = std::make_unique<Volume>(geom->composeSymNameITS(), gLbl++, det, false);
  sym2vol[root->getSymName()] = root.get();
  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {
    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      volHB = root->addChild(geom->composeSymNameHalfBarrel(ilr, ihb), gLbl++, det, false);
      sym2vol[volHB->getSymName()] = volHB;
      int nstavesHB = geom->getNumberOfStaves(ilr) / 2;
      for (int ist = 0; ist < nstavesHB; ist++) {
        volSt = volHB->addChild(geom->composeSymNameStave(ilr, ihb, ist), gLbl++, det, false);
        sym2vol[volSt->getSymName()] = volSt;
        for (int ihst = 0; ihst < geom->getNumberOfHalfStaves(ilr); ihst++) {
          volHSt = volSt->addChild(geom->composeSymNameHalfStave(ilr, ihb, ist, ihst), gLbl++, det, false);
          sym2vol[volHSt->getSymName()] = volHSt;
          for (int imd = 0; imd < geom->getNumberOfModules(ilr); imd++) {
            volMod = volHSt->addChild(geom->composeSymNameModule(ilr, ihb, ist, ihst, imd), gLbl++, det, false);
            sym2vol[volMod->getSymName()] = volMod;
          }
        }
      }
    }
  }
  int lay = 0, hba = 0, sta = 0, ssta = 0, modd = 0, chip = 0;
  for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
    geom->getChipId(ich, lay, hba, sta, ssta, modd, chip);
    Label lbl(det, ich, true);
    Volume* parVol = sym2vol[modd < 0 ? geom->composeSymNameStave(lay, hba, sta) : geom->composeSymNameModule(lay, hba, sta, ssta, modd)];
    if (!parVol) {
      LOGP(fatal, "did not find parent for chip {}", ich);
    }
    int nch = modd < 0 ? geom->getNumberOfChipsPerStave(lay) : geom->getNumberOfChipsPerModule(lay);
    auto* chipVol = parVol->addChild<SensorITS>(geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip % nch), lbl);
    chipVol->setSensorId(ich);
    sensorMap[lbl] = chipVol;
  }
  return root;
}

Volume::Ptr DetectorITS::buildHierarchyIT3(Volume::SensorMapping& sensorMap)
{
  uint32_t gLbl{0}, det{0};
  auto geom = o2::its::GeometryTGeo::Instance();
  Volume *volHB{nullptr}, *volSt{nullptr}, *volHSt{nullptr}, *volMod{nullptr};
  std::unordered_map<std::string, Volume*> sym2vol;
  auto root = std::make_unique<Volume>(geom->composeSymNameITS(), gLbl++, det, false);
  sym2vol[root->getSymName()] = root.get();
  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {
    const bool isLayITS3 = (ilr < 3);
    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      volHB = root->addChild(geom->composeSymNameHalfBarrel(ilr, ihb, isLayITS3), gLbl++, det, false);
      sym2vol[volHB->getSymName()] = volHB;
      if (isLayITS3) {
        volHB->setSensorId((2 * ilr) + ihb);
        continue;
      }
      int nstavesHB = geom->getNumberOfStaves(ilr) / 2;
      for (int ist = 0; ist < nstavesHB; ist++) {
        volSt = volHB->addChild(geom->composeSymNameStave(ilr, ihb, ist), gLbl++, det, false);
        sym2vol[volSt->getSymName()] = volSt;
        for (int ihst = 0; ihst < geom->getNumberOfHalfStaves(ilr); ihst++) {
          volHSt = volSt->addChild(geom->composeSymNameHalfStave(ilr, ihb, ist, ihst), gLbl++, det, false);
          sym2vol[volHSt->getSymName()] = volHSt;
          for (int imd = 0; imd < geom->getNumberOfModules(ilr); imd++) {
            volMod = volHSt->addChild(geom->composeSymNameModule(ilr, ihb, ist, ihst, imd), gLbl++, det, false);
            sym2vol[volMod->getSymName()] = volMod;
          }
        }
      }
    }
  }
  int lay = 0, hba = 0, sta = 0, ssta = 0, modd = 0, chip = 0;
  for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
    geom->getChipId(ich, lay, hba, sta, ssta, modd, chip);
    const bool isLayITS3 = (lay < 3);
    Label lbl(det, ich, true);
    if (isLayITS3) {
      Volume* parVol = sym2vol[geom->composeSymNameHalfBarrel(lay, hba, true)];
      if (!parVol) {
        LOGP(fatal, "did not find parent for chip {}", ich);
      }
      auto* tile = parVol->addChild<SensorIT3>(geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip, true), lbl);
      tile->setPseudo(true);
      tile->setSensorId(ich);
      sensorMap[lbl] = tile;
    } else {
      Volume* parVol = sym2vol[modd < 0 ? geom->composeSymNameStave(lay, hba, sta) : geom->composeSymNameModule(lay, hba, sta, ssta, modd)];
      if (!parVol) {
        LOGP(fatal, "did not find parent for chip {}", ich);
      }
      int nch = modd < 0 ? geom->getNumberOfChipsPerStave(lay) : geom->getNumberOfChipsPerModule(lay);
      auto* chipVol = parVol->addChild<SensorITS>(geom->composeSymNameChip(lay, hba, sta, ssta, modd, chip % nch), lbl);
      chipVol->setSensorId(ich);
      sensorMap[lbl] = chipVol;
    }
  }
  return root;
}

} // namespace o2::alignrs
