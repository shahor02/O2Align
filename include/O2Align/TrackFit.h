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

#ifndef O2_ITS3_ALIGN_TRACKFIT
#define O2_ITS3_ALIGN_TRACKFIT

#include <Eigen/Dense>

#include "ITSBase/GeometryTGeo.h"
#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"

namespace o2::alignrs
{
using Mat51 = Eigen::Matrix<double, 5, 1>;
using Mat55 = Eigen::Matrix<double, 5, 5>;
using TrackD = o2::track::TrackParCovD;

template <typename T, typename F>
track::TrackParametrizationWithError<T> convertTrack(const track::TrackParametrizationWithError<F>& trk)
{
  if constexpr (std::is_same_v<T, F>) {
    return trk;
  }
  track::TrackParametrizationWithError<T> dst;
  dst.setX(trk.getX());
  dst.setAlpha(trk.getAlpha());
  for (int iPar{0}; iPar < track::kNParams; ++iPar) {
    dst.setParam(trk.getParam(iPar), iPar);
  }
  dst.setAbsCharge(trk.getAbsCharge());
  dst.setPID(trk.getPID());
  dst.setUserField(trk.getUserField());
  for (int iCov{0}; iCov < track::kCovMatSize; ++iCov) {
    dst.setCov(trk.getCov()[iCov], iCov);
  }
  return dst;
}

// Both tracks must be at the same (alpha, x).
// Returns the interpolated track.
template <typename T>
o2::track::TrackParametrizationWithError<T> interpolateTrackParCov(
  const o2::track::TrackParametrizationWithError<T>& tA,
  const o2::track::TrackParametrizationWithError<T>& tB)
{
  auto res = tA;
  if (!tA.isValid() || !tB.isValid() || tA.getAlpha() != tB.getAlpha() || tA.getX() != tB.getX()) {
    res.invalidate();
    return res;
  }
  auto unpack = [](const std::array<T, track::kCovMatSize>& c) {
    Mat55 m;
    for (int i = 0, k = 0; i < 5; ++i) {
      for (int j = 0; j <= i; ++j, ++k) {
        m(i, j) = m(j, i) = (double)c[k];
      }
    }
    return m;
  };
  Mat55 cA = unpack(tA.getCov());
  Mat55 cB = unpack(tB.getCov());
  Eigen::LLT<Mat55> lltA(cA), lltB(cB);
  if (lltA.info() != Eigen::Success || lltB.info() != Eigen::Success) {
    res.invalidate();
    return res;
  }
  Mat55 wA = lltA.solve(Mat55::Identity());
  Mat55 wB = lltB.solve(Mat55::Identity());
  Mat55 wTot = wA + wB;
  Eigen::LLT<Mat55> lltTot(wTot);
  if (lltTot.info() != Eigen::Success) {
    res.invalidate();
    return res;
  }
  Mat55 cTot = lltTot.solve(Mat55::Identity());
  Mat51 pA, pB;
  for (int i = 0; i < 5; ++i) {
    pA(i) = tA.getParam(i);
    pB(i) = tB.getParam(i);
  }
  Mat51 pTot = cTot * (wA * pA + wB * pB);
  // build result - same alpha/x as inputs
  for (int i = 0; i < 5; ++i) {
    res.setParam(pTot(i), i);
  }
  for (int i = 0, k = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j, ++k) {
      res.setCov(static_cast<T>(cTot(i, j)), k);
    }
  }
  return res;
}

/// Reset the covariance of a seed before a refit: the default diagonal errors, with a 100% error on q/pt
inline void resetTrackCovariance(TrackD& trk)
{
  trk.resetCovariance();
  trk.setCov(trk.getQ2Pt() * trk.getQ2Pt() * trk.getCov()[14], 14);
}

} // namespace o2::alignrs

#endif
