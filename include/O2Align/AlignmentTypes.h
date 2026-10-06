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

#ifndef O2_ALIGN_TYPES_H
#define O2_ALIGN_TYPES_H

#include <string>
#include <vector>
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "DataFormatsITS/TrackITS.h"
#include "O2Align/Label.h"

namespace o2::alignrs
{

struct Measurement final {
  double dy = 0.f;
  double dz = 0.f;
  double sig2y = 0.f;
  double sig2z = 0.f;
  double phi = 0.f;
  double z = 0.f;
  ClassDefNV(Measurement, 1)
};

struct FrameInfoExt final {
  enum {
    Invalid = -2, // the point is not usable
    Vertex = -1   // the point is the primary vertex
  };
  int8_t lr = Invalid;            // detector-specific layer-like index, -1 = vtx, -2 = invalid point
  Label label{};                  // label of the sensitive volume this point belongs to
  float x{-999.f};                // X of the measurement in the tracking frame
  float alpha{-999.f};            // rotation angle of the tracking frame
  o2::BaseCluster<float> cluster; // cluster info
  bool zFromTrack{false};         // the Z is not measured but taken from the track (TRD): no Z derivatives wrt the alignment
  std::string asString() const;
  bool isValid() const { return lr > Invalid; }
  bool isVertex() const { return lr == Vertex; }
  ClassDefNV(FrameInfoExt, 3)
};

struct FitInfo final {
  float chi2Ndf{-1}; // Chi2/Ndf of track refit
  float chi2{-1};    // Chi2
  int ndf{-1};       // ndf
  ClassDefNV(FitInfo, 1)
};

struct Track {
  o2::dataformats::VtxTrackIndex gid; // global track ID
  o2::its::TrackITS its;           // original ITS track
  o2::track::TrackParCovD track;   // prepared track state
  FitInfo kfFit;                   // kf fit information
  FitInfo gblFit;                  // gbl fit information
  std::vector<Measurement> points; // measurment point
  std::vector<FrameInfoExt> info;  // frame info, owned by the track (detectors append to it)

  /// KF refit over the frames in the inclusive [frameStart, frameStop] slot range of `info`,
  /// starting from the state stored in `track` (which must be defined at the frameStart point).
  /// The range is traversed inward if frameStart > frameStop, invalid frames are skipped.
  /// The fit is done on a copy, hence `track` and `kfFit` are modified only on success:
  /// on success `track` becomes the state at the frameStop point and the chi2 of the fitted
  /// points is added to `kfFit.chi2`.
  /// reset: reset the seed covariance and zero `kfFit.chi2` before the fit, otherwise the fit
  ///        continues from the seed as is.
  /// cropOnFailure: on a failure of an outward fit drop the frames from the frameStart slot on.
  bool fitTrack(int frameStart, int frameStop, bool cropOnFailure = false, bool reset = false);

  /// Continue the KF fit outward, using the frames appended by the caller to `info` starting
  /// from the slot frameStart. The caller must append them in the outward direction. On failure
  /// the new frames are dropped and both `track` and `kfFit` stay intact.
  bool continueFitOutward(int frameStart);

  /// Update the track with the vertex point: the vertex is stored, in the frame of the track DCA
  /// to it, in the frame prebooked by the caller in the info[0] slot, then the track is fitted to
  /// it from the innermost measured point. On failure the info[0] slot is flagged as Invalid and
  /// both `track` and `kfFit` stay intact.
  /// covScale: scale of the vertex covariance used in the KF update only, the vertex point stored
  /// in info[0] keeps the true one. With covScale << 1 the updated state passes (almost) through
  /// the vertex, as required for the reference of a composed trajectory with the common vertex.
  /// The chi2 of the vertex point added to `kfFit.chi2` is then that of the scaled covariance.
  bool updateWithVertex(const o2::dataformats::VertexBase& vtx, double covScale = 1.);

  ClassDefNV(Track, 2)
};

struct TrackSlopes {
  double dydx{0.};
  double dzdx{0.};  
  static TrackSlopes computeTrackSlopes(double snp, double tgl)  {
    const double csci = 1. / std::sqrt(1. - (snp * snp));
    return {.dydx = snp * csci, .dzdx = tgl * csci};
  }
  ClassDefNV(TrackSlopes, 1)
};

std::vector<double> legendrePols(int order, double x); // RSTODO temporary here

  
} // namespace o2::alignrs

#endif
