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

#include "DetectorsBase/Propagator.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Detector.h"
#include "O2Align/Params.h"
#include "O2Align/TrackFit.h"

namespace o2::alignrs
{

bool Detector::continueFitOutward(Track& resTrack, size_t nFramesIni)
{
  const auto& params = Params::Instance();
  auto prop = o2::base::PropagatorD::Instance();
  auto trFit = resTrack.track; // fit a copy: the seed must be preserved if the fit fails
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trFit);
  }
  float chi2 = 0.f;
  for (auto frameIt = resTrack.info.begin() + nFramesIni; frameIt != resTrack.info.end(); ++frameIt) {
    const auto& cluster = frameIt->cluster;
    if (!prop->propagateToAlphaX(trFit, refLin, frameIt->alpha, frameIt->x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      resTrack.info.resize(nFramesIni);
      return false;
    }
    chi2 += static_cast<float>(trFit.getPredictedChi2Quiet(cluster));
    if (!trFit.update(cluster)) {
      resTrack.info.resize(nFramesIni);
      return false;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(cluster.getY());
      refLin->setZ(cluster.getZ());
    }
  }
  resTrack.track = trFit;
  resTrack.kfFit.chi2 += chi2;
  // RSTODO fill the ndf / chi2Ndf info during the inward refit
  return true;
}

} // namespace o2::alignrs
