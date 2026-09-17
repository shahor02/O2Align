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

#include <string>
#include <format>

#include "DetectorsBase/Propagator.h"
#include "MathUtils/Utils.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Params.h"

namespace o2::alignrs
{

std::string FrameInfoExt::asString() const
{
  return std::format("[{}] Sensor={} Layer={} X={} Alpha={} y={} z={}", label.asString(), cluster.getSensorID(), lr, x, alpha, cluster.getY(), cluster.getZ());
}

bool Track::fitTrack(int frameStart, int frameStop, bool cropOnFailure, bool reset)
{
  const int nFrames = static_cast<int>(info.size());
  if (frameStart < 0 || frameStart >= nFrames || frameStop < 0 || frameStop >= nFrames) {
    return false;
  }
  const auto& params = Params::Instance();
  auto prop = o2::base::PropagatorD::Instance();
  auto trFit = track; // fit a copy: the seed must be preserved if the fit fails
  if (reset) {
    trFit.resetCovariance();
    trFit.setCov(trFit.getQ2Pt() * trFit.getQ2Pt() * trFit.getCov()[14], 14);
  }
  const int step = frameStart > frameStop ? -1 : 1;
  const int frameEnd = frameStop + step; // 1st slot past the range
  const bool crop = cropOnFailure && step > 0; // cropping is meaningful for the outward fit only
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trFit);
  }
  float chi2 = 0.f;
  for (int i = frameStart; i != frameEnd; i += step) {
    const auto& frame = info[i];
    if (!frame.isValid()) { // invalid point
      continue;
    }
    if (!prop->propagateToAlphaX(trFit, refLin, frame.alpha, frame.x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      if (crop) {
        info.resize(frameStart);
      }
      return false;
    }
    const auto& cluster = frame.cluster;
    chi2 += static_cast<float>(trFit.getPredictedChi2Quiet(cluster));
    if (!trFit.update(cluster)) {
      if (crop) {
        info.resize(frameStart);
      }
      return false;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(cluster.getY());
      refLin->setZ(cluster.getZ());
    }
  }
  track = trFit;
  kfFit.chi2 = reset ? chi2 : kfFit.chi2 + chi2;
  // RSTODO fill the ndf / chi2Ndf info during the inward refit
  return true;
}

bool Track::continueFitOutward(int frameStart)
{
  if (frameStart >= static_cast<int>(info.size())) {
    return true; // no new frames were appended: nothing to fit
  }
  return fitTrack(frameStart, static_cast<int>(info.size()) - 1, true, false);
}

bool Track::updateWithVertex(const o2::dataformats::VertexBase& vtx)
{
  if (info.empty()) {
    return false;
  }
  auto& frame = info.front(); // the slot prebooked for the vertex point
  frame.lr = FrameInfoExt::Invalid; // flagged as valid only by a successful update
  auto trDCA = track; // find the DCA on a copy, the fit itself starts from the innermost measured point
  auto prop = o2::base::PropagatorD::Instance();
  if (!prop->propagateToDCA(vtx, trDCA, prop->getNominalBz())) {
    return false;
  }
  double ca{0}, sa{0};
  frame.alpha = static_cast<float>(trDCA.getAlpha());
  o2::math_utils::bringToPMPi(frame.alpha);
  o2::math_utils::sincosd(frame.alpha, sa, ca);
  frame.x = static_cast<float>(vtx.getX() * ca + vtx.getY() * sa); // the vertex position in the track frame
  frame.cluster = o2::BaseCluster<float>(-1, frame.x, static_cast<float>(-vtx.getX() * sa + vtx.getY() * ca), vtx.getZ(),
                                         0.5f * (vtx.getSigmaX2() + vtx.getSigmaY2()), // the vertex is round in the transverse plane
                                         vtx.getSigmaZ2(), 0.f);
  frame.lr = FrameInfoExt::Vertex;
  if (!fitTrack(0, 0, false, false)) {
    frame.lr = FrameInfoExt::Invalid;
    return false;
  }
  return true;
}

// RSTODO temporary here
std::vector<double> legendrePols(int order, double x)
{
  std::vector<double> p(order + 1);
  p[0] = 1.;
  if (order > 0) {
    p[1] = x;
  }
  for (int n = 1; n < order; ++n) {
    p[n + 1] = ((2 * n + 1) * x * p[n] - n * p[n - 1]) / (n + 1);
  }
  return p;
}

} // namespace o2::alignrs
