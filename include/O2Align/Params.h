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
  // Track selection for alignment
  float minPt = 0.4f;        // minimum pt required  
  int minITSClsTrack = 7;    // minimum number of ITS clusters in the full track
  int minITSClsAB = 4;       // minimum number of ITS clusters in the AB track

  // Track selection for PV refitting
  float minPtPV = 0.3f;      // minimum pt required for tracks used in PV refitting
  int minITSClsPV = 4;       // minimum number of ITS clusters in the track used for PV refitting
  // PV constraint
  int usePVConstraintMinTracks = 5; // minimum number of tracks to use PV constraint
  // steering of multi-track PV constraint:
  // 0: no multi-track PV constraint, 
  // 1: use multi-track PV constraint (at most 1st maxTracksPerMultiTrackPV) and allow single tracks with eventual PV constraint for the rest
  // 2: use multi-track PV constraint (at most 1st maxTracksPerMultiTrackPV) and use the rest w/o PV constraint at all
  // 3: multi-track PV constraint only (at most 1st maxTracksPerMultiTrackPV), no single-track at all
  // 4: all contributors in multi-track PV constraints, split round-robin into groups of at most maxTracksPerMultiTrackPV,
  //    each a separate composed trajectory with its own free vertex; the non-contributors w/o PV constraint, as in 2
  int useMultiTrackPVConstraint = 1; 
  int maxTracksPerMultiTrackPV = 15;       // maximum number of tracks to use for multi-track PV constraint (mode 4: per group), 
  bool refitVtxFollowsMV = true;           // per-track PV constraint: the refitted vertex moves with the correction of the mean vertex position (same global derivatives as the prior), otherwise it anchors the frame of the starting geometry
  bool scaleMVPriorWithNTracks = true;     // per-track PV constraint: scale the luminous region covariance of the mean vertex prior by the number of tracks of the vertex, for the prior to count once per collision
  float vtxMultiTrackRefCovScale = 1e-2f;  // multi-track PV constraint: scale of the vertex covariance in the KF update defining the reference state at the vertex point, for the references of all tracks to pass (almost) through the common vertex. The stored vertex point keeps the true covariance
  std::string MVTimeSlotsJson = ""; // json file with mean vertex calibration intervals in ms


  int minGBLPoints = 4;      // minimum number of GBL points per track
  float cosmMaxChi2Cl = 30.f;      // cosmics: max chi2 wrt the reference track of an ITS cluster, TRD tracklet or TOF cluster to be accepted (ambiguity resolution)
  
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


  // TPC drift calibration
  std::string VDTimeSlotsJson = ""; // json file with TPC drift calibration intervals in ms

  int verbose = 0; // verbosity level
  bool useStableRef = true; // use input tracks as linearization point
  float meanPtB0 = 0.4f;    // impose this pT on tracks for B=0 (ignore if <=0)
  float meanPtB0Cosm = 2.f; // impose this pT on cosmic tracks for B=0 (ignore if <=0)
  float maxChi2Ndf = 15;    // maximum Chi2/Ndf allowed for GBL fit
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

  // misalignment simulation (closure test), see doc/MisalignmentTest.md
  std::string misAlgJson;                    // rule file naming the volumes, DOFs and values of the misalignment to inject, none if empty
  std::string misAlgRefMilleRes;             // millepede.res of a control run (same data and DOF config, no misalignment) subtracted at MilleRes stage in addition to the injected values: needed for the closure on real data
  std::string closureOutJson = "closure.json"; // MilleRes stage with misAlgJson: sensor-level closure summary

  // debug output: fraction of the tracks (randomly selected) written with their points, GBL residuals and kinks
  // to the "gbltrack" tree of debug_alg[_<lane>].root (DbgTrack, see AlignmentDbg.h). <=0: none, >=1: all
  float debugTracksFraction = -1.f;

  // DOF configuration (JSON file defining which volumes have which DOFs)
  std::string dofConfigJson; // if empty, no DOFs are configured

  // Ridder options
  int ridderMaxExtrap = 10;
  double ridderRelIniStep[5] = {0.01, 0.01, 0.02, 0.02, 0.02};
  double ridderMaxIniStep[5] = {0.1, 0.1, 0.05, 0.05, 0.05};
  double ridderShrinkFac = 2.0;
  double ridderEps = 1e-16;

  // MillePede output
  // weight the hierarchical mean-of-children constraint (Volume::writeChildrenMeanConstraints) by each
  // contributor's accumulated track count instead of an equal 1/N split, so a low-statistics edge
  // volume (e.g. an IB end chip) isn't pulled to match higher-statistics siblings in that equation.
  // A per-volume/per-DOF "weightedMean" rule in dofConfigJson overrides this default and, on its own,
  // requests the same statistics regardless of this switch. The steering files are always written from
  // endOfStream(), not init(), so the counts are available by then either way:
  // - with ProcessData, they are this run's own, accumulated as its TFs are processed;
  // - without it (steering files regenerated standalone), see statInFile below.
  bool useStatWeightedMeanConstraints = false;
  // when regenerating the steering files standalone (MilleSteer w/o ProcessData, e.g. to try a more
  // restricted dofConfigJson against mille data already recorded with a wider-open one), there is no
  // data pass of its own to accumulate statistics from: if useStatWeightedMeanConstraints is set, or
  // dofConfigJson has at least one "weightedMean" rule asking for true, give here the path of a
  // ProcessingStats ROOT file (statFile of a previous ProcessData run) to load each volume's track
  // count from instead. Left empty, or if the file can't be read or doesn't cover every contributor of
  // a weighted equation, falls back to an equal-split write for that equation, with a warning.
  std::string statInFile = "";
  // coefficients of the "globalModes"/"groupModes" equations (Volume ModeEquationOptions, see
  // doc/DOFConfig_rules.md#how-the-mode-coefficients-are-built)
  bool modeSensorProjection = true; // project each mode on the displacement a DOF gives to the sensor box corners of its subtree; false: legacy, the mode evaluated at the volume's origin, each volume weighted 1
  bool groupModesRelToRoots = true; // groupModes: constrain the mode relative to the free rigid motion of each group volume (needs modeSensorProjection)
  std::string statFile = "procStat";
  std::string milleBinFile = "mp2data";
  std::string milleConFile = "mp2con.txt";
  std::string milleParamFile = "mp2param.txt";
  std::string milleTreeFile = "mp2tree.txt";
  std::string milleResFile = "millepede.res";
  std::string milleResOutJson = "result.json";

  // conversion of the fitted rigid-body corrections to o2::detectors::AlignParam (MilleRes stage)
  std::string algParamsOutFile = "alignment.root"; // if not empty, write the combined alignment of each detector to <DET>_<algParamsOutFile>, to be applied to the ideal geometry
  // Regularize output option: transfer to every selected TGeo alignable parent the common motion of its
  // alignable children (robust L1 rigid-motion fit of the displacements of their bounding box corners),
  // deepest parents first, keeping the placement of every sensor
  std::string regularizeParents = ""; // regex (DOF-config style, implicit arbitrary prefix, '|' for alternatives) of the symbolic names of the parents to regularize, all if empty
  int regularizeMinChildren = 2;      // skip the parents with fewer alignable children

  // WeakModesDiag output (MilleRes stage): <weakModesOutBase>.json and <weakModesOutBase>.root
  std::string weakModesOutBase = "weakModes";

  O2ParamDef(Params, "AlignParams");
};
} // namespace o2::alignrs

#endif
