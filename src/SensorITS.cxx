// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

#include "Framework/Logger.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "O2Align/SensorITS.h"
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
