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

#ifndef O2_ALIGN_PARAMS_H
#define O2_ALIGN_PARAMS_H

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "DetectorsBase/Propagator.h"

namespace o2::alignrs
{

struct Params : public o2::conf::ConfigurableParamHelper<Params> {
  // Track selection
  float minPt = 0.5f;        // minimum pt required  
  int minITSCls = 7;         // minimum number of ITS clusters
  float maxITSChi2Ndf = 1.2; // maximum ITS track chi2 // RSTODO not used, remove?

  int minGBLPoints = 4;      // minimum number of GBL points per track
  
  bool writeLocalAlignParams = true; // create output AlignParams in local convention.

  // ITS overlap handling
  float ITSOverlapMargin = 0.15;     // consider for overlaps only clusters within this margin from the chip edge (in cm)
  float ITSOverlapMaxChi2 = 16;      // max chi2 between track and overlapping cluster
  float ITSOverlapMaxDZ = 0.3;       // max difference in Z for clusters on overlapping ITS chips to consider as candidate for a double hit
  int ITSOverlapEdgeRows = 1;        // require clusters to not have pixels closer than this distance from the edge

  // TPC
  int minTPCClusters = 10;            // discard tracks with less clusters
  int minTPCPadRow = 6;               // min TPC pad-row to account
  int maxTPCPadRow = 146;             // max TPC pad-row to account
  int maxTPCRowsCombined = 1;         // allow combining clusters on so many rows to a single cluster
  int discardEdgePadrows = 3;         // discard padrow if its distance to stack edge padrow < this
  float discardEdgePadDepth = 2.5;    // discard clusters too close to the sector edge (in pad units)

  // TRD
  bool applyXORTRD = false;           // apply XOR in TRD tranfformer
  int minTRDTracklets = 3;         // min TRD tracklets to accept the track
  float TRDCorrDVT = 0.f;          // correction to Vdrift*t (equivalent to a shift in X at which Y is evaluated)
  float TRDNonRCCorrDzDtgl = 0.f;  // correction in Z proportional to tgl for non-row-crossing tracklets

  // TOF
  int minTOFClusters = 1; // min TOF clusters to accept track

  // propagation opt
  double maxSnp = 0.85;
  double maxStep = 2.0;
  // o2::base::PropagatorD::MatCorrType matCorrType = o2::base::PropagatorD::MatCorrType::USEMatCorrTGeo;
  o2::base::PropagatorD::MatCorrType corrType = o2::base::PropagatorD::MatCorrType::USEMatCorrLUT;

  // PV constraint
  int usePVConstraintMinTracks = 5; // minimum number of tracks to use PV constraint
  // steering of multi-track PV constraint:
  // 0: no multi-track PV constraint, 
  // 1: use multi-track PV constraint (at most 1st maxTracksPerMultiTrackPV) and allow single tracks with eventual PV constraint for the rest
  // 2: use multi-track PV constraint (at most 1st maxTracksPerMultiTrackPV) and use the rest w/o PV constraint at all
  // 3: multi-track PV constraint only (at most 1st maxTracksPerMultiTrackPV), no single-track at all
  int useMultiTrackPVConstraint = 1; 
  int maxTracksPerMultiTrackPV = 15;       // maximum number of tracks to use for multi-track PV constraint, 
  bool refitVtxFollowsMV = true;           // per-track PV constraint: the refitted vertex moves with the correction of the mean vertex position (same global derivatives as the prior), otherwise it anchors the frame of the starting geometry
  bool scaleMVPriorWithNTracks = true;     // per-track PV constraint: scale the luminous region covariance of the mean vertex prior by the number of tracks of the vertex, for the prior to count once per collision
  float vtxMultiTrackRefCovScale = 1e-2f;  // multi-track PV constraint: scale of the vertex covariance in the KF update defining the reference state at the vertex point, for the references of all tracks to pass (almost) through the common vertex. The stored vertex point keeps the true covariance
  std::string MVTimeSlotsJson = ""; // json file with mean vertex calibration intervals in ms

  // TPC drift calibration
  std::string VDTimeSlotsJson = ""; // json file with TPC drift calibration intervals in ms

  int verbose = 0; // verbosity level
  bool useStableRef = true; // use input tracks as linearization point
  float meanPtB0 = 0.3f;    // impose this pT on tracks for B=0 (ignore if <=0)
  float maxChi2Ndf = 10;    // maximum Chi2/Ndf allowed for GBL fit
  float minMS = 1e-6f;      //  minimum scattering to account for
  // A step whose material would otherwise be lumped into a single thin scatterer at its far end is
  // instead split into two (one extra scatterer-only point at the radius midpoint) when both the
  // crossed material and the radial gap are large, e.g. the ITS-TPC transition: lumping biases the
  // lever arm GBL uses between the surrounding measurements. Must have splitMSThreshold >= minMS.
  float splitMSThreshold = 5e-3f;   // minimum scattering angle (rad) of the full step to warrant splitting
  float splitMSStepMinDX = 5.f;     // minimum radial gap (cm) of the full step to warrant splitting

  // per chip extra error
  float extraClsErrYITS[7] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  float extraClsErrZITS[7] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  // extra systematic errors of the non-ITS measurements
  float extraClsErrYTPC = 0.f;
  float extraClsErrZTPC = 0.f;
  float extraClsErrYTRD = 0.f;
  float extraClsErrZTRD = 0.f;
  float extraClsErrYTOF = 0.f;
  float extraClsErrZTOF = 0.f;

  // misalignment simulation
  bool doMisalignmentLeg = false;  // simulate Legendre deformation on ITS3 layers
  bool doMisalignmentRB = false;   // simulate rigid body misalignment on ITS3 layers
  bool doMisalignmentInex = false; // simulate in-extensional deformation on ITS3 layers
  std::string misAlgJson;          // JSON file with deformation and/or rigid body params

  // DOF configuration (JSON file defining which volumes have which DOFs)
  std::string dofConfigJson; // if empty, no DOFs are configured

  // Ridder options
  int ridderMaxExtrap = 10;
  double ridderRelIniStep[5] = {0.01, 0.01, 0.02, 0.02, 0.02};
  double ridderMaxIniStep[5] = {0.1, 0.1, 0.05, 0.05, 0.05};
  double ridderShrinkFac = 2.0;
  double ridderEps = 1e-16;

  // MillePede output
  bool volumeStatistics = true; // include per-volume accepted-track counters in the statistics output
  std::string statFile = "procStat";
  std::string milleBinFile = "mp2data";
  std::string milleConFile = "mp2con.txt";
  std::string milleParamFile = "mp2param.txt";
  std::string milleTreeFile = "mp2tree.txt";
  std::string milleResFile = "millepede.res";
  std::string milleResOutJson = "result.json";

  // conversion of the fitted rigid-body corrections to o2::detectors::AlignParam (MilleRes stage)
  std::string algParamsOutFile = "alignment.root"; // if not empty, write the combined alignment of each detector to <DET>_<algParamsOutFile>, to be applied to the ideal geometry

  O2ParamDef(Params, "AlignParams");
};
} // namespace o2::alignrs

#endif
