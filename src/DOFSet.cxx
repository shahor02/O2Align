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

#include "Align/DOFSet.h"

namespace o2::alignrs
{

void DOFSet::validateDerivativeOutput(Eigen::Ref<Eigen::MatrixXd> out)
{
  if (out.rows() != 3 || out.cols() != getNDOFs()) {
    throw std::invalid_argument(std::format("Derivative buffer shape {}x{} does not match expected 3x{}", out.rows(), out.cols(), dofSet.nDOFs()));
  }
  out.setZero();
}

void RigidBodyDOFSet::fillDerivatives(const DerivativeContext& ctx, Eigen::Ref<Eigen::MatrixXd> out) const
{
  validateDerivativeOutput(out);

  const double csp = 1. / std::sqrt(1. + (ctx.tgl * ctx.tgl));
  const double uP = ctx.snp * csp;
  const double vP = ctx.tgl * csp;

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
  
}  // namespace o2::alignrs
