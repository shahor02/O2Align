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
  //const auto* chipL2G = mPN->GetMatrix();
  //mL2G = *chipL2G;
  // RS: using mPN for ITS2 is not correct since this provides the chip matrix rather than the embedded sensor layer
  // Use the same method as in GeometryTGeo::fillCache
  mL2G = o2::its::GeometryTGeo::Instance()->extractMatrixSensor(getSensorId());
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
// The values are those of labelToValue (any reference is subtracted from it by the caller).
bool SensorITS::MP2JSON_Calib(const std::map<uint32_t, double>& labelToValue, nlohmann::json& entry) const
{
  using json = nlohmann::json;
  const auto* cal = getCalib();
  const auto calibLbl = getCalibLabels().front(); // the ITS3 deformations are not time-sliced

  if (cal->type() == DOFSet::Type::Legendre) {
    const auto* leg = dynamic_cast<const LegendreDOFSet*>(cal);
    const int order = leg->order();
    json matrix = json::array();
    int idx = 0;
    for (int i = 0; i <= order; ++i) {
      json row = json::array();
      for (int j = 0; j <= i; ++j) {
        row.push_back(getFittedValue(labelToValue, calibLbl.raw(idx)));
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
    json inexEntry;
    json modesObj = json::object();
    for (int n = 2; n <= maxN; ++n) {
      const int off = InextensionalDOFSet::modeOffset(n);
      json modeArr = json::array();
      for (int k = 0; k < 4; ++k) {
        modeArr.push_back(getFittedValue(labelToValue, calibLbl.raw(off + k)));
      }
      modesObj[std::to_string(n)] = modeArr;
    }
    inexEntry["modes"] = modesObj;
    inexEntry["alpha"] = getFittedValue(labelToValue, calibLbl.raw(inexSet->alphaIdx()));
    inexEntry["beta"] = getFittedValue(labelToValue, calibLbl.raw(inexSet->betaIdx()));
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
