# DOF Configuration Rules

`Volume::applyDOFConfig` traverses the hierarchy and applies matching rules to non-pseudo volumes. A `defaults` object, if present, is applied first to every volume; rules then run in order, and a later matching clause replaces the earlier DOF set of the same kind. Patterns are matched against symbolic volume names, with an implicit `*` prepended. Since `*` can match `/`, use exact-depth patterns when targeting a hierarchy level.

## Rigid-Body DOFs

The `rigidBody` clause frees an arbitrary subset of a volume's six DOFs (`TX,TY,TZ,RX,RY,RZ`, in
the volume's local frame), all the others of that volume staying fixed; a volume not matched by any
rule keeps none free. Free several at once by listing them together, e.g. the two in-plane shifts
plus the in-plane rotation of every ITS stave:

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?", "rigidBody": ["TX", "TZ", "RZ"] }
```

`"rigidBody": "all"` frees all six, `"rigidBody": "fixed"` (or `false`) fixes all six — the latter
is only needed to override a `defaults` clause or an earlier rule, since no DOF is free by default.
As with every other DOF-set clause, a later matching rule replaces the DOF set of an earlier one
rather than merging with it, so listing a subset is not cumulative across rules.

## Mean Vertex

The mean vertex volume is `PVT/meanVertex`. Its calibration DOFs are the position `X, Y, Z` and the
beam-line slopes `SlopeX, SlopeY` (dX/dZ, dY/dZ), all 5 defaulting to free. Each of the 5 is
independent of the others — any subset may be free while the rest stay fixed, with no all-or-nothing
grouping. A fixed DOF still enters the fit at its current prior value (from the `MeanVertexObject`
delivered by the CCDB for the active calibration slot), it is simply not adjusted.

Enable all 5 with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex" } }
```

Fix all 5 with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex", "fixed": true } }
```

Fix only the slopes, keeping X, Y, Z free:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex", "fix": ["SlopeX", "SlopeY"] } }
```

As with any calibration DOF set, an external measurement of an individual DOF (e.g. a value derived
from a prior calibration) can be added via the `measurement` clause, see
[Measurements of Rigid-Body and Calibration DOFs](#measurements-of-rigid-body-and-calibration-dofs)
below; it goes in a separate rule following the one that frees/fixes the DOFs.

## Automatic vs. Optional Constraints at a Parent/Children Boundary

`writeRigidBodyConstraints` decides, for every rigid-body DOF of every volume, whether a constraint
over its children is needed at all, and whether it is written automatically or only on request. The
two cases that matter in practice:

1. **Parent DOF free, and some descendant has *any* rigid-body DOF free** — not necessarily the same
   one (a chain of free levels, with possibly several non-rigid-body/fixed levels looked through in
   between — see `writeChildrenMeanConstraints`). `addChildrenMeanTerms` builds the parent's DOF `iDOF`
   constraint from *every* free DOF `jDOF` of each contributing child, weighted by the full jacobian
   element `j(iDOF, jDOF)`, not only `j(iDOF, iDOF)`: a child's `RX`, say, contributes to the parent's
   `TY` constraint too, with a coefficient set by the lever arm between them. This is not an
   approximation — `getJP2L()`/`getJL2P()` is the same fixed 6x6 matrix per child used by
   `buildPointGlobals` to chain every leaf measurement up through the hierarchy, so a parent's free DOF
   and a differently-named free DOF of a child really can be exactly indistinguishable once propagated
   through it: a genuine degeneracy of the fit. The constraint "weighted mean of the children's
   movement (transported to the parent's frame) vanishes" is written **automatically**, with no DOF
   config needed beyond freeing both levels:
   ```json
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?", "rigidBody": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?", "rigidBody": ["TY"] }
   ```
   Here the half-barrel's `TY` and the mean `TY` of its staves are tied together automatically; a stave
   free only in, say, `RZ` instead would still contribute to that same constraint via `j(TY, RZ)`.

2. **Parent DOF fixed, descendants' DOF free.** There is no parameter at the parent level to be
   degenerate with, so **no constraint is written by default**: the children (and the hierarchy below
   them) are free to drift together in that DOF, which may be a genuinely unconstrained or only
   weakly constrained mode of the fit (e.g. nothing anchors the overall `TY` of a detector whose
   envelope is fixed but whose staves are all free in `TY`). Use `pinChildrenMean` to request the same
   mean-zero constraint explicitly even though the parent DOF is fixed — this **does not** turn the
   parent DOF into a fit parameter, it only pins the common mode of what's below it to the parent's
   (fixed) position:
   ```json
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?", "rigidBody": "fixed", "pinChildrenMean": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?", "rigidBody": ["TY"] }
   ```
   Without the `pinChildrenMean` clause, this configuration leaves the mean `TY` of the staves
   unconstrained; with it, the mean is pinned to the (fixed) half-barrel position.

A third, unproblematic case is a free parent DOF with **no** free descendant DOF below it at all (see
[Aligning a Subtree as One Rigid Body](#aligning-a-subtree-as-one-rigid-body)): there is nothing to be
degenerate with, so no constraint is needed or written, automatically or otherwise, and the parent
DOF is still fully observable from the leaf measurements.

None of this applies to calibration DOFs (e.g. the mean vertex, TPC drift): they are never organized
in a parent/children hierarchy of their own, so there is no equivalent automatic constraint — only the
explicit `measurement` clause (previous section) adds a prior to a calibration DOF.

## Pinning the Common Mode of Children

`Volume::writeRigidBodyConstraints` forces, for every *free* rigid-body DOF of a parent, the
weighted mean of its active children's movement, transported to the parent's frame and projected
onto that DOF, to vanish — summed over *every* free DOF of each child (not only the same-named one),
weighted by the corresponding element of the child-to-parent jacobian; see
[Automatic vs. Optional Constraints](#automatic-vs-optional-constraints-at-a-parentchildren-boundary)
for why a child's differently-named free DOF is included on exactly the same footing. This breaks a
real degeneracy: moving the parent and moving all its children together in the same direction are
otherwise indistinguishable. A child without any free rigid-body DOF (no `rigidBody` rule, or a
`fixed` one) is looked through: its own children, with the jacobians chained, take its place, and so
on down to the nearest level having a free DOF. The walk stops at the first descendant with *any*
free DOF, so free DOFs of a grandchild along a direction its (partially free) parent keeps fixed
are not included in the mean — avoid such configurations. When the parent DOF is fixed instead, that degeneracy doesn't arise,
so no such constraint is written by default — the children are then free to collectively drift in
that DOF, which may be a genuinely unconstrained (or only weakly constrained) mode of the fit.

The `pinChildrenMean` clause asks for the same mean-zero constraint on a DOF of the parent even
though that DOF is fixed there, to remove such a mode explicitly. It takes a boolean (all DOFs), the
string `"all"`, or an array of DOF names (`TX,TY,TZ,RX,RY,RZ`), and — like `rigidBody`/`calib` — a
later matching rule replaces the flags of an earlier one rather than merging with them:

```json
{ "match": "TRD/sm[0-9][0-9]", "rigidBody": "fixed", "pinChildrenMean": ["TY", "RZ"] }
```

With the chambers free, this pins the mean `TY` shift and `RZ` rotation of the chambers of every
supermodule to the (fixed) supermodule position, without making the supermodule itself a fit
parameter. The flags don't depend on the volume's `rigidBody` clause: they also work on a volume with
no rigid-body DOF set at all. A pinned DOF that no free descendant DOF contributes to (nothing free below,
or none of the free DOFs feeds it through the chained jacobians) is skipped with a warning.
`pinChildrenMean` is ignored, with a warning, on a volume that is not rigid-body alignable (the root,
or an envelope without its own geometry such as `TRD_envelope`), so the common mode of the
supermodules of a detector can't be pinned this way.

## Measurements of Rigid-Body and Calibration DOFs

The `measurement` clause supplies an external measurement (e.g. from a survey, or a value carried
over from a prior calibration) of individual DOFs of the matching volumes — rigid-body DOFs,
calibration DOFs, or both at once (each name is resolved against whichever of the volume's two DOF
sets has it; the two never share a name, e.g. `TX..RZ` vs. `X,Y,Z,SlopeX,SlopeY` for the mean
vertex). It is written to `mp2con.txt` as a Millepede `Measurement` record (`Measurement value
sigma`, followed by the label of the DOF with coefficient 1), i.e. a soft constraint `p = value ±
sigma` which, unlike a `Constraint`, does not fix the parameter. It is an object `DOF name ->
[value, sigma]` (or `{"value": v, "sigma": s}`), in the units of the DOF (cm for `TX,TY,TZ` and the
mean-vertex `X,Y,Z`; rad for `RX,RY,RZ`; dimensionless for `SlopeX,SlopeY`), relative to the geometry
used by the fit:

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?", "rigidBody": ["TX", "TY"] },
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel?", "measurement": { "TX": [0.01, 0.005], "TY": {"value": -0.02, "sigma": 0.01} } }
```

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex" } },
{ "match": "PVT/meanVertex", "measurement": { "SlopeY": [0.0002, 0.0001] } }
```

A `rigidBody`/`calib` clause *defines* a DOF set (listed DOFs free, all the others fixed), while
`measurement` only adds a prior to DOFs that are already free. A rule containing `measurement`
together with `rigidBody` or `calib` is therefore rejected (fatal): put the measurement in its own
rule, after the one defining the DOFs.

As for `pinChildrenMean`, no DOF set is created and a later matching rule replaces all the
measurements of both DOF sets of the volume (even if the new rule only touches one of them).
`sigma` must be positive. A measurement of a DOF which is not free at the time of writing, as well as
the clause on a volume having neither a rigid-body nor a calibration DOF set at all, is ignored with
a warning (under an exact, non-wildcard pattern) or silently skipped (under a wildcard pattern, e.g.
a `measurement` rule meant for a different kind of volume).

## Aligning a Subtree as One Rigid Body

The opposite situation to `pinChildrenMean` is freeing a parent's rigid-body DOF while giving none
of its descendants a `rigidBody` rule at all. This is not auto-disabled and needs no mean-of-children
constraint: the parent-to-child jacobian is applied unconditionally down to every leaf measurement
regardless of whether a descendant carries a DOFSet, so the parent's free DOF is fully observable on
its own. A constraint is only needed to resolve the degeneracy of a DOF being free at *both* a parent
and a descendant level simultaneously (the case of the previous section); with nothing free below
the parent, at any depth, there is no such degeneracy. Use this to align, say, a half-barrel as a single rigid body
while keeping its staves/modules/chips fixed to it, without needing to give them any DOF config.

## TPC Drift

The TPC drift calibration belongs to `TPC_envelope` and is not enabled by default. Enable both `VDRIFT` and `DRIFTOFF` with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift" } }
```

Fix both with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift", "fixed": true } }
```

## ITS Alignable Match Strings

The patterns below are based on the output from the `macro/printAlignableVolumes.C` running over Run3 TGeometry root file. 
The ITS envelope is the `ITS` alignable. The layer name is a path component, not a separate alignable entry. In the match patterns, layer IDs use `[0-2]` for the inner barrel and `[3-6]` for the outer barrel. Half-barrel, half-stave, module, and chip IDs in these paths are single digits, matched with `?`. Inner-barrel stave IDs are single digits; outer-barrel stave IDs can be one or two digits, so both exact-width forms are listed there. Avoid a trailing `*` when matching a specific level: `fnmatch` allows it to match `/` and descendants too.

### Inner Barrel (Layers 0–2)

The inner-barrel hierarchy has half-barrels, staves, and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?` |
| Stave | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?` |
| Chip | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?/ITSUChip?` |

### Outer Barrel (Layers 3–6)

The outer-barrel hierarchy adds half-staves and modules between staves and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?` |
| Stave, one-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?` |
| Stave, two-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]` |
| Half-stave, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?` |
| Half-stave, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?` |
| Module, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?/ITSUModule?` |
| Module, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?/ITSUModule?` |
| Chip, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?/ITSUModule?/ITSUChip?` |
| Chip, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?/ITSUModule?/ITSUChip?` |

For example, replace a `match` value with one of the strings above in a rule such as `{ "match": "...", "rigidBody": ["TX", "TY", "TZ"] }`. ITS3 uses different names, such as `ITS3CarbonForm`, and is not selected by these patterns.

## TRD Alignable Match Strings

TRD has a dummy `TRD_envelope` (not a geometry alignable), followed by alignable supermodules and chambers. The geometry log uses two digits for supermodule IDs and one digit each for stack and plane IDs:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TRD_envelope` |
| Supermodule | `TRD/sm[0-9][0-9]` |
| Chamber | `TRD/sm[0-9][0-9]/st?/pl?` |

## TOF Alignable Match Strings

TOF has a dummy `TOF_envelope` (not a geometry alignable), followed by alignable supermodules and strips. Both IDs in the geometry names are formatted as two digits:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TOF_envelope` |
| Supermodule | `TOF/sm[0-9][0-9]` |
| Strip | `TOF/sm[0-9][0-9]/strip[0-9][0-9]` |

The envelopes are not rigid-body alignable; rigid-body clauses on them are ignored. Calibration clauses can be configured on these synthetic volumes.