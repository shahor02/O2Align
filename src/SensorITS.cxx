// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

#include <nlohmann/json.hpp>

#include "Framework/Logger.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "O2Align/SensorITS.h"
#include "O2Align/DOFSet.h"
#include "ITSBase/GeometryTGeo.h"

namespace o2::alignrs
{

void SensorITS::defineMatrixL2G()
{
  // the chip volume is not the measurment plane, need to correct for the epitaxial layer
  const auto* chipL2G = mPN->GetMatrix();
  mL2G = *chipL2G;
  double delta = itsmft::SegmentationAlpide::SensorLayerThickness - itsmft::SegmentationAlpide::SensorLayerThicknessEff;
  TGeoTranslation tra(0., 0.5 * delta, 0.);
  mL2G *= tra;
}

void SensorITS::defineMatrixT2L()
{
  double locA[3] = {-100., 0., 0.}, locB[3] = {100., 0., 0.}, gloA[3], gloB[3];
  mL2G.LocalToMaster(locA, gloA);
  mL2G.LocalToMaster(locB, gloB);
  double dx = gloB[0] - gloA[0], dy = gloB[1] - gloA[1];
  double t = (gloB[0] * dx + gloB[1] * dy) / (dx * dx + dy * dy);
  double xp = gloB[0] - (dx * t), yp = gloB[1] - (dy * t);
  double alp = std::atan2(yp, xp);
  o2::math_utils::bringTo02Pid(alp);
  mT2L.RotateZ(alp * TMath::RadToDeg()); // mT2L before is identity and afterwards rotated
  const TGeoHMatrix l2gI = mL2G.Inverse();
  mT2L.MultiplyLeft(l2gI);
}

void SensorITS::computeJacobianL2T(const double* posLoc, Matrix66& jac) const
{
  jac.setZero();
  Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> rotT2L(mT2L.GetRotationMatrix());
  Eigen::Matrix3d skew, rotL2T = rotT2L.transpose();
  skew << 0, -posLoc[2], posLoc[1], posLoc[2], 0, -posLoc[0], -posLoc[1], posLoc[0], 0;
  jac.topLeftCorner<3, 3>() = rotL2T;
  jac.topRightCorner<3, 3>() = -rotL2T * skew;
  jac.bottomRightCorner<3, 3>() = rotL2T;
}

// The caller (Volume::writeMillepedeResults) only invokes this when getCalib() has free DOFs.
bool SensorITS::MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const
{
  using json = nlohmann::json;
  const auto* cal = getCalib();
  const int id = getSensorId();
  const auto calibLbl = getLabel().asCalib();

  if (cal->type() == DOFSet::Type::Legendre) {
    const auto* leg = dynamic_cast<const LegendreDOFSet*>(cal);
    const int order = leg->order();
    static const std::vector<std::vector<double>> emptyM;
    const auto* injMPtr = &emptyM;
    if (inj) {
      if (auto injIt = inj->matrix.find(id); injIt != inj->matrix.end()) {
        injMPtr = &injIt->second;
      }
    }
    const auto& injM = *injMPtr;
    json matrix = json::array();
    int idx = 0;
    for (int i = 0; i <= order; ++i) {
      json row = json::array();
      for (int j = 0; j <= i; ++j) {
        uint32_t raw = calibLbl.raw(idx);
        auto it = labelToValue.find(raw);
        double fitted = it != labelToValue.end() ? it->second : 0.0;
        double ref = (i < static_cast<int>(injM.size()) && j < static_cast<int>(injM[i].size())) ? injM[i][j] : 0.0;
        row.push_back(fitted - ref);
        ++idx;
      }
      matrix.push_back(row);
    }
    entry["matrix"] = matrix;
    return true;
  }

  if (cal->type() == DOFSet::Type::Inextensional) {
    const auto* inexSet = static_cast<const InextensionalDOFSet*>(cal);
    const int maxN = inexSet->maxOrder();
    static const InjectedMisalignment::Inextensional emptyI;
    const auto* injIPtr = &emptyI;
    if (inj) {
      if (auto injIt = inj->inextensional.find(id); injIt != inj->inextensional.end()) {
        injIPtr = &injIt->second;
      }
    }
    const auto& injI = *injIPtr;

    json inexEntry;
    json modesObj = json::object();
    for (int n = 2; n <= maxN; ++n) {
      int off = InextensionalDOFSet::modeOffset(n);
      std::array<double, 4> injCoeffs = {0., 0., 0., 0.};
      if (injI.modes.contains(n)) {
        injCoeffs = injI.modes.at(n);
      }
      json modeArr = json::array();
      for (int k = 0; k < 4; ++k) {
        uint32_t raw = calibLbl.raw(off + k);
        auto it = labelToValue.find(raw);
        double fitted = it != labelToValue.end() ? it->second : 0.0;
        modeArr.push_back(fitted - injCoeffs[k]);
      }
      modesObj[std::to_string(n)] = modeArr;
    }
    inexEntry["modes"] = modesObj;

    // alpha
    uint32_t rawAlpha = calibLbl.raw(inexSet->alphaIdx());
    auto itA = labelToValue.find(rawAlpha);
    inexEntry["alpha"] = (itA != labelToValue.end() ? itA->second : 0.0) - injI.alpha;

    // beta
    uint32_t rawBeta = calibLbl.raw(inexSet->betaIdx());
    auto itB = labelToValue.find(rawBeta);
    inexEntry["beta"] = (itB != labelToValue.end() ? itB->second : 0.0) - injI.beta;

    entry["inextensional"] = inexEntry;
    return true;
  }

  LOGP(warn, "MP2JSON_Calib: unhandled calib type {} for volume {}", static_cast<int>(cal->type()), getSymName());
  return false;
}

void SensorIT3::defineMatrixL2G()
{
  mL2G = *mPN->GetMatrix();
}

void SensorIT3::defineMatrixT2L()
{
  double locA[3] = {-100., 0., 0.}, locB[3] = {100., 0., 0.}, gloA[3], gloB[3];
  mL2G.LocalToMaster(locA, gloA);
  mL2G.LocalToMaster(locB, gloB);
  double dx = gloB[0] - gloA[0], dy = gloB[1] - gloA[1];
  double t = (gloB[0] * dx + gloB[1] * dy) / (dx * dx + dy * dy);
  double xp = gloB[0] - (dx * t), yp = gloB[1] - (dy * t);
  double alp = std::atan2(yp, xp);
  o2::math_utils::bringTo02Pid(alp);
  mT2L.RotateZ(alp * TMath::RadToDeg());
  const TGeoHMatrix l2gI = mL2G.Inverse();
  mT2L.MultiplyLeft(l2gI);
}

void SensorIT3::computeJacobianL2T(const double* posLoc, Matrix66& jac) const
{
  jac.setZero();
  Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> rotT2L(mT2L.GetRotationMatrix());
  Eigen::Matrix3d skew, rotL2T = rotT2L.transpose();
  skew << 0, -posLoc[2], posLoc[1], posLoc[2], 0, -posLoc[0], -posLoc[1], posLoc[0], 0;
  jac.topLeftCorner<3, 3>() = rotL2T;
  jac.topRightCorner<3, 3>() = -rotL2T * skew;
  jac.bottomRightCorner<3, 3>() = rotL2T;
}

} // namespace o2::alignrs
