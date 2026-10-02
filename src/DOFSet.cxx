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

#include "O2Align/DOFSet.h"
#include "O2Align/AlignmentTypes.h" // RSTODO this is just for legendrePols, get rid of it
#include "O2Align/Detector.h"       // for the detector index of the derivative context
#include "ITS3Base/SpecsV2.h" // RSTODO try to get rid of this dependence
#include "CommonConstants/MathConstants.h"
#include "DataFormatsTPC/Constants.h"
#include <format>

namespace o2::alignrs
{

double* DOFSet::gParVals = nullptr;
double* DOFSet::gParErrs = nullptr;

void DOFSet::validateDerivativeOutput(Eigen::Ref<Eigen::MatrixXd> out) const
{
  if (out.rows() != 2 || out.cols() != nDOFs()) { // rows: Y, Z of the measurement
    throw std::invalid_argument(std::format("Derivative buffer shape {}x{} does not match expected 2x{}", out.rows(), out.cols(), nDOFs()));
  }
  out.setZero();
}

std::string LegendreDOFSet::dofName(int idx) const
{
  int i = 0;
  while ((i + 1) * (i + 2) / 2 <= idx) {
    ++i;
  }
  int j = idx - (i * (i + 1) / 2);
  return std::format("L({},{})", i, j);
}

std::string InextensionalDOFSet::dofName(int idx) const
{
  if (idx == alphaIdx()) {
    return "alpha";
  }
  if (idx == betaIdx()) {
    return "beta";
  }
  int n = (idx / 4) + 2;
  int sub = idx % 4;
  static constexpr const char* subNames[] = {"a", "b", "c", "d"};
  return std::format("{}_{}", subNames[sub], n);
}

  
void RigidBodyDOFSet::fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const
{
  validateDerivativeOutput(out);

  // We need the slopes of the locally straight track wrt the normal to the sensor, y' = dy/dx and
  // z' = dz/dx, while the version below computes dy/ds and dz/ds (s = path length), underestimating
  // e.g. the z-row derivatives wrt TX and RY by a factor 1/sqrt(1+tgl^2), ~0.7 at tgl = 1:
  // const double csp = 1. / std::sqrt(1. + (ctx.tgl * ctx.tgl));
  // const double uP = ctx.snp * csp;
  // const double vP = ctx.tgl * csp;
  const double uP = ctx.dydx; // snp / sqrt(1 - snp^2)
  const double vP = ctx.dzdx; // tgl / sqrt(1 - snp^2)

  out(0, TX) = uP;
  out(0, TY) = -1.;
  out(0, RX) = ctx.trkZ;
  out(0, RY) = ctx.trkZ * uP;
  out(0, RZ) = -ctx.trkY * uP;

  out(1, TX) = vP;
  out(1, TZ) = -1.;
  out(1, RX) = -ctx.trkY;
  out(1, RY) = ctx.trkZ * vP;
  out(1, RZ) = -ctx.trkY * vP;
}

std::pair<double, double> LegendreDOFSet::getSensorPhiBorders(int sensorID, double radius)
{
  const bool isTop = sensorID % 2 == 0;
  const double gapPhi = std::asin(o2::its3::constants::equatorialGap / 2. / radius);
  return {o2::math_utils::to02Pid(((isTop ? 0. : 1.) * o2::constants::math::PI) + gapPhi),
          o2::math_utils::to02Pid(((isTop ? 1. : 2.) * o2::constants::math::PI) - gapPhi)};
}

double LegendreDOFSet::getSensorPhiWidth(int sensorID, double radius)
{
  const auto [phiBorder1, phiBorder2] = getSensorPhiBorders(sensorID, radius);
  const double width = phiBorder2 - phiBorder1;
  return (width < 0.) ? width + o2::constants::math::TwoPI : width;
}

std::pair<double, double> LegendreDOFSet::computeUV(double gloX, double gloY, double gloZ, int sensorID, double radius)
{
  const double phi = o2::math_utils::to02Pid(std::atan2(gloY, gloX));
  const auto [phiBorder1, phiBorder2] = getSensorPhiBorders(sensorID, radius);
  const double u = (((phi - phiBorder1) * 2.) / (phiBorder2 - phiBorder1)) - 1.;
  const double v = ((2. * gloZ + o2::its3::constants::segment::lengthSensitive) / o2::its3::constants::segment::lengthSensitive) - 1.;
  return {u, v};
}

  
void LegendreDOFSet::fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const
{
  validateDerivativeOutput(out);
  if (ctx.sensorID < 0 || ctx.layerID < 0) {
    throw std::invalid_argument("LegendreDOFSet requires an ITS3 measurement context");
  }

  const double gloX = ctx.measX * std::cos(ctx.measAlpha);
  const double gloY = ctx.measX * std::sin(ctx.measAlpha);
  const auto [u, v] = computeUV(gloX, gloY, ctx.measZ, ctx.sensorID, o2::its3::constants::radii[ctx.layerID]);
  const auto pu = o2::alignrs::legendrePols(mOrder, u);
  const auto pv = o2::alignrs::legendrePols(mOrder, v);
  const double phiWidth = getSensorPhiWidth(ctx.sensorID, o2::its3::constants::radii[ctx.layerID]);

  // same intergration as `evaluateLegendreShift' but now for each order separateley
  Eigen::VectorXd arcMismatch = Eigen::VectorXd::Zero(nDOFs());
  if (std::abs(u) > o2::constants::math::Almost0) {
    const auto& x = GaussNodes;
    const auto& w = GaussWeights;
    const double mid = 0.5 * u;
    const double half = 0.5 * u;
    for (int iq = 0; iq < 8; ++iq) {
      const double up = mid + (half * x[iq]);
      const auto puQ = o2::alignrs::legendrePols(mOrder, up);
      int idx = 0;
      for (int i = 0; i <= mOrder; ++i) {
        for (int j = 0; j <= i; ++j) {
          arcMismatch[idx] += w[iq] * puQ[j] * pv[i - j];
          ++idx;
        }
      }
    }
    arcMismatch *= 0.5 * phiWidth * half;
  }

  int idx = 0;
  for (int i = 0; i <= mOrder; ++i) {
    for (int j = 0; j <= i; ++j) {
      const double basis = pu[j] * pv[i - j];
      out(0, idx) = (ctx.dydx * basis) + arcMismatch[idx];
      out(1, idx) = ctx.dzdx * basis;
      ++idx;
    }
  }
}

void InextensionalDOFSet::fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const
{
  validateDerivativeOutput(out);
  if (ctx.layerID < 0) {
    throw std::invalid_argument("InextensionalDOFSet requires an ITS3 measurement context");
  }

  const double r = o2::its3::constants::radii[ctx.layerID];
  const double phi = std::atan2(r * std::sin(ctx.measAlpha), r * std::cos(ctx.measAlpha));
  const double z = ctx.measZ;

  for (int n = 2; n <= mMaxOrder; ++n) {
    const double sn = std::sin(n * phi);
    const double cn = std::cos(n * phi);
    const double n2 = static_cast<double>(n * n);
    const int off = modeOffset(n);

    out(0, off + 0) = -(z / r) * (n * sn + ctx.dydx * n2 * cn);
    out(1, off + 0) = -cn - ctx.dzdx * (z / r) * n2 * cn;

    out(0, off + 1) = (z / r) * (n * cn - ctx.dydx * n2 * sn);
    out(1, off + 1) = -sn * (1. + ctx.dzdx * (z / r) * n2);

    out(0, off + 2) = -cn + ctx.dydx * n * sn;
    out(1, off + 2) = ctx.dzdx * n * sn;

    out(0, off + 3) = -sn - ctx.dydx * n * cn;
    out(1, off + 3) = -ctx.dzdx * n * cn;
  }

  out(0, alphaIdx()) = z / r;
  out(1, alphaIdx()) = -phi;

  out(0, betaIdx()) = -phi - ctx.dydx;
  out(1, betaIdx()) = -ctx.dzdx;
}

void TPCVDriftDOFSet::fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const
{
  validateDerivativeOutput(out);
  if (ctx.detID != Detector::DetTPC || ctx.volID < 0) {
    throw std::invalid_argument(std::format("TPCVDriftDOFSet requires a TPC measurement context, got det {} vol {}", ctx.detID, ctx.volID));
  }
  if (mZLength <= 0.) {
    throw std::invalid_argument("TPCVDriftDOFSet was not given the TPC drift length, see setZLength");
  }
  // the side is defined by the sector, not by the sign of the cluster Z, which for a cluster
  // reconstructed with a wrong track time may end up on the wrong side of the central electrode
  const double side = (ctx.volID < o2::tpc::constants::MAXSECTOR / 2) ? 1. : -1.;
  const double driftLength = mZLength - side * ctx.measZ;
  // only the Z row: the drift affects neither the pad direction nor, to 1st order, the
  // space-charge correction of the cluster
  out(1, VDRIFT) = side * driftLength;
  out(1, DRIFTOFF) = side;
}

void MeanVertexDOFSet::fillDerivatives(const DerivativeContext& /*ctx*/, Eigen::Ref<Eigen::MatrixXd> /*out*/) const
{
  throw std::logic_error("MeanVertexDOFSet derivatives are set by the mean vertex constraint, not by the measurement chain");
}

}  // namespace o2::alignrs
