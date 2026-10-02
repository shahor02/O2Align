---
name: gbl
description: Reference for the GBL (General Broken Lines) library and how O2Align wraps it — GblPoint/GblTrajectory/GblMeasurement semantics, jacobian and sign conventions, scatterers, composed trajectories with inner transformations, global/local derivatives, and how these map onto buildPointGlobals, toGBLOrder, Label packing, and the vertex-constraint trick in AlignmentSpec.cxx. Load whenever a question touches GBL points, trajectories, globals/labels, scatterers, composed/vertex fits, or jacobian conventions in this codebase.
---

# GBL (General Broken Lines) reference

This is reference material, not a workflow — read the relevant section(s) below, then answer or
edit code directly. It exists because GBL knowledge is split across the external library headers,
this project's wrapping code, and undocumented rationale; re-deriving all three from scratch each
time is the expensive part this skill avoids.

Library headers live at (version may differ, same layout):
`~/alice/sw/ubuntu2404_x86-64/GBL/V03-01-04-19/include/{GblPoint,GblTrajectory,GblMeasurement,GblData}.h`.
Go there for exact signatures; this file gives the semantics and the O2Align-specific mapping.

## Core model

A `GblTrajectory` is an ordered list of `GblPoint`s connected by propagation jacobians. Each point
may carry:
- **A jacobian from the previous point** (5×5, passed to the `GblPoint` constructor) — propagates
  the local curvilinear track parameters from the previous point to this one. The identity jacobian
  marks the first point of a (sub-)trajectory.
- **0+ measurements** (`addMeasurement`) — up to 5D, residual + precision (inverse covariance) in
  the *measurement* system; GBL diagonalizes the precision matrix internally and transforms
  residuals/derivatives accordingly, so callers don't need to pre-diagonalize.
- **A scatterer** (`addScatterer`, 2D kink angles) — thin scatterer at this point, or
  `addThickScatterer` for a 4D thick one. `scatDim` is 0/2/4 (none/thin/thick).
- **Global derivatives** (`addGlobals(labels, derivatives)`) — one column per external
  (Millepede) parameter, residual-vs-parameter. **Convention: derivatives are of the *prediction*,
  not the residual** — i.e. the negative of d(residual)/d(parameter) — because Millepede's
  linear model is `measurement = prediction(params) + noise`. Get this sign backwards and the fit
  still converges but to the wrong sign of correction.
- **Local derivatives** (`addLocals`) — per-trajectory-only parameters, not shared with Millepede
  (not used in this codebase as of this writing; globals are used throughout instead).

A **simple trajectory** is `GblTrajectory(pointList, flagCurv, flagU1dir, flagU2dir)`.

A **composed trajectory** —
`GblTrajectory(vector<pair<vector<GblPoint>, Eigen::MatrixXd>>&)` — stitches several independent
sub-trajectories (e.g. one per track) together through a shared set of *external* parameters (e.g.
a common vertex position). Each pair is `(points, innerTransformation)`: `innerTransformation` maps
the external parameters onto the local offset parameters **at the first point of that sub-list**.
The first point is therefore the anchor of the stitching and must exist in the point list even if
it carries no measurement of its own — see "Vertex-constrained composed trajectories" below. An
overload additionally takes `(extDerivatives, extMeasurements, extPrecisions)` to attach a genuine
external measurement (e.g. a prior) on the shared parameters themselves, on top of the stitching.

## O2Align's wrapping layer

- **`toGBLOrder`** ([AlignmentSpec.cxx](../../../src/AlignmentSpec.cxx)) reorders the ALICE track
  jacobian (its own parameter order) into GBL's curvilinear order before it's handed to the
  `GblPoint` constructor. Any time a jacobian looks transposed or permuted vs. what you'd expect
  from the ALICE track model, check this conversion first.
- **`buildPointGlobals()`** is the single implementation of the measurement-leaf-to-root walk that
  produces the `(labels, derivatives)` pair for `addGlobals`: base derivative in the tracking frame
  → sensor local frame via `computeJacobianL2T` → transported level-by-level with each child's
  `getJL2P()`, picking up columns at every volume with an active `RigidBodyDOFSet`. `fillGBLPoints`
  calls it once per non-vertex frame; the vertex frame bypasses it entirely (`frame.isVertex()`),
  since the mean vertex has no hierarchy position — see `DetectorPVT` in
  [CLAUDE.md](../../../CLAUDE.md).
- **Millepede labels** are packed via [Label.h](../../../include/O2Align/Label.h):
  `DOF(7)|CALIB(1)|ID(19)|SENS(1)|DET(3)`, 31 usable bits — **one bit short of the full 32** because
  GBL treats labels as a signed `int` and the MSB is reserved. `DET` uniqueness is the only thing
  that keeps per-detector volume-counter restarts from colliding.
- **Scattering precision placement**: `fillGBLPoints` lumps the material crossed *since the
  previous point* into a thin scatterer *at the current point* (`msErr > mParams->minMS && ip < np
  - 1`), not split between the two ends. `lt.clearFast()` resets the accumulated path only once it's
  actually been charged to a scatterer.
- **Measurement precision diagonalization is GBL's job, not ours** — `getMeasurementPrecision`
  builds the raw 2×2 precision from the cluster covariance; don't pre-diagonalize it before calling
  `addMeasurement`.

## Vertex-constrained composed trajectories (`buildGBLVertex`)

The common-primary-vertex fit (`AlignmentSpec::buildGBLVertex`) is the trickiest use of composed
trajectories in this codebase:

1. Each contributing track's `info[0]` slot is pre-booked for the vertex point by
   `Track::updateWithVertex` (projects the vertex onto the track's DCA frame, stores it as a
   pseudo-cluster with the propagated covariance).
2. `fillGBLPoints(track, 0, /*skipFirstMeas=*/true, points)` builds that track's points **with no
   measurement added at the vertex point** (`addMeas = !(skipFirstMeas && points.empty())`) — the
   point still exists, still gets the propagation jacobian and (if it crossed material) a
   scatterer, but no `addMeasurement` call. **This is deliberate, not a bug**: the vertex position
   is now a *shared external parameter* across all tracks via the composed trajectory, and adding
   the refitted-vertex measurement there as well would double-count information already folded
   into the per-track covariance by step 1.
3. `computeVertexTransformation(track)` builds that track's `innerTransformation` — the 2×3
   (local-offset-vs-vertex-XYZ) Jacobian at the DCA frame, same math as the `trans` matrix inside
   `Track::updateWithVertex`.
4. For the *first* used track only, `addMeanVertexPrior(...)` attaches the luminous-region prior
   `V ~ N(mu, Sigma_lumi)` as an actual measurement on `points.front()` — this is the one place a
   "vertex point" *does* get a measurement, but it's the prior on the shared mean-vertex parameter,
   not the per-track vertex residual skipped in step 2. If `mPVT->getPositionLabels()` is non-empty,
   this measurement also carries global derivatives `-trans` (sign flipped vs. the usual convention
   because here it's the *measurement*, not the prediction, that moves with the parameter — see the
   comment at that call site).
5. `GblTrajectory(pointsAndTrans)` stitches all tracks' point lists through their individual
   `innerTransformation`s; fitting it determines the shared vertex offset plus each track's own
   curvature/scattering parameters simultaneously, making the vertex constraint exact by
   parameterization rather than by an extra penalty term.

A consequence worth remembering: because the vertex point is skipped in `addMeas`, it's also absent
from `resTrack.points` (`fillGBLPoints` only `emplace_back`s there when `addMeas` is true) — so
`resTrack.points` and `resTrack.info` are **not index-aligned** for vertex-constrained tracks. Don't
assume `points[i]` corresponds to `info[i]` without checking which frames had `addMeas == true`.

## Quick sanity checks when something's wrong

- Fit converges but corrections have the wrong sign → check whether a global derivative was added
  as d(residual)/d(param) instead of d(prediction)/d(param) (should be negated), or vice versa for
  a measurement-side quantity like the mean-vertex prior.
- Composed/vertex fit silently drops a track → `fillGBLPoints` returning false, or
  `points.size() < 2` in `buildGBLVertex` (a track needs at least the vertex point plus one real
  measurement to be usable).
- Labels collide across detectors → check the `DET` field assignment, not `ID`/`SENS`; every
  detector restarts its own volume counter from 0 ([Label.h](../../../include/O2Align/Label.h)).
- A point with no measurement and no scatterer appearing in `VerboseGBL` output isn't necessarily
  wrong — it may be a composed-trajectory anchor point (see above) or a frame whose material was
  below `mParams->minMS`.
