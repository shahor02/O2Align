# Cosmic tracks (`--track-sources cosmic-ext`)

The input is the `TrackCosmicsExtended` of the O2 `cosmics-cluster-collector` (`o2_cosmics_full.root`, option
`--cosmics-full-infile`): the raw clusters of both legs are inside, no cluster/track input of the collision data is needed.
There is no primary vertex, the mean vertex is not used. Needed options in addition: `--disable-mc` and `--hbfutils-config <root file>`
with the `tfidinfo` vector (`std::vector<TFIDInfo>`, e.g. filled from the `tfID` branch of the input), otherwise the CCDB
objects are queried for the current time.

## What is done per cosmic (`AlignmentSpec::processCosmicExt`)

1. Skipped: `duplicateOf >= 0`, cosmics with the points of only one leg.
2. Each requested detector provides the points of each leg (`Detector::prepareCosmicLeg`, leg 0 = bottom, 1 = top), the track-dependent
   quantities (cluster errors, TRD tilt correction, ambiguity resolution) being evaluated with the reference track of the cosmic:
   `polished` (refit at the TOF time) if available, otherwise the matcher's `cosmic`; both fly upward.
   - ITS: per leg and layer the cluster with the best chi2 wrt the reference (`AlignParams.cosmMaxChi2Cl`); in addition the best
     cluster looking like the other hit of the overlap with it (`DetectorITS::isOverlapPair`, `ITSOverlapMaxChi2`) is added.
     Clusters biased by the chip edge are discarded. ITS3 is not supported.
   - TPC: only the `Attached` clusters, transformed with the ideal map (no distortion corrections) and the drift velocity and t0
     of the TF stored in the input, at the time of the cosmic (TOF time if any) common for both legs; merged to the superclusters
     as for the collision tracks (the shared cluster state is not accounted).
   - TRD, TOF: per leg and layer (TOF: per leg) the point closest to the reference. TOF errors are those of a pad.
3. A detector contributes if the number of its measurements in both legs together reaches the minimum of the collision tracks
   (`minITSClsTrack`, `minTPCClusters`, `minTRDTracklets`). If an obligatory detector fails the cosmic is lost, otherwise the
   detector is dropped. The set of the contributing detectors must be accepted by `--detectors` (combinations, as for the collisions).
4. KF fit as in `MatchCosmics::refitWinners`: a track parametrisation cannot describe the motion opposite to the X axis of the tracking
   frame, therefore the reference track is fitted downward (inverted) through the bottom leg, outward and back inward, then inverted
   again (`TrackParCov::invert`) and continued upward through the top leg. The reference state is at the highest point
   (largest Y of the lab frame) of the top leg. The points are ordered in Y, top leg first; `Track::invertAt` is the first point
   of the bottom leg.
5. GBL: the trajectory goes from the highest point down along the top leg (still in the upward parametrisation), the track is
   inverted at the step to the first bottom point (`fillGBLPoints`, the transport jacobian of this step includes the inversion,
   the scatterer of this step is not split) and continues along the bottom leg.

Energy loss is propagated with the default sign (along the parametrisation), not along the flight as in `MatchCosmics`.
