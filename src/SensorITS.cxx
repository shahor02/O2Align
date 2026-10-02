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
#include "Framework/Logger.h"

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

namespace
{
/// azimuth of the foot of the normal from the beam axis on a flat chip, as in
/// GeometryTGeo::extractSensorXAlpha: the local X axis lies in the chip plane
double getFlatChipAlpha(const TGeoHMatrix& l2g)
{
  double locA[3] = {-100., 0., 0.}, locB[3] = {100., 0., 0.}, gloA[3], gloB[3];
  l2g.LocalToMaster(locA, gloA);
  l2g.LocalToMaster(locB, gloB);
  double dx = gloB[0] - gloA[0], dy = gloB[1] - gloA[1];
  double t = (gloB[0] * dx + gloB[1] * dy) / (dx * dx + dy * dy);
  double xp = gloB[0] - (dx * t), yp = gloB[1] - (dy * t);
  double alp = std::atan2(yp, xp);
  o2::math_utils::bringTo02Pid(alp);
  return alp;
}
} // namespace

void SensorITS::defineMatrixT2L()
{
  setT2LFromAlpha(getFlatChipAlpha(mL2G));
}

// The caller (Volume::writeMillepedeResults) only invokes this when getCalib() has free DOFs.
bool SensorITS::MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, const InjectedMisalignment* inj, nlohmann::json& entry) const
{
  using json = nlohmann::json;
  const auto* cal = getCalib();
  const int id = getSensorId();
  const auto calibLbl = getCalibLabels().front(); // the ITS3 deformations are not time-sliced

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
        const double fitted = getFittedValue(labelToValue, calibLbl.raw(idx));
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
        modeArr.push_back(getFittedValue(labelToValue, calibLbl.raw(off + k)) - injCoeffs[k]);
      }
      modesObj[std::to_string(n)] = modeArr;
    }
    inexEntry["modes"] = modesObj;

    inexEntry["alpha"] = getFittedValue(labelToValue, calibLbl.raw(inexSet->alphaIdx())) - injI.alpha;
    inexEntry["beta"] = getFittedValue(labelToValue, calibLbl.raw(inexSet->betaIdx())) - injI.beta;

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
  // The local origin of the tube segment is the centre of the circle, hence the flat-chip formula
  // does not apply: take the reference frame of the tile from the geometry, which treats ITS3
  // separately (GeometryTGeo::extractSensorXAlpha). It is used only as a reference, the derivatives
  // being computed in the tracking frame of every cluster, see hasFixedTrackingFrame.
  setT2LFromAlpha(o2::its::GeometryTGeo::Instance()->getSensorRefAlpha(getSensorId()));
}

} // namespace o2::alignrs
