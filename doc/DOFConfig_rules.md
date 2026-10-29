# DOF Configuration Rules

The name of the JSON file for alignment configuration is passed to the workflow by AlignParams.dofConfigJson configurable. Below the rules to build such file are described. The misalignment to inject for a closure test (`AlignParams.misAlgJson`) is given by a file of rules of the same style, see [MisalignmentTest.md](MisalignmentTest.md). Note that any JSON input 
`"#<something>" : "your comment"` or `"_comment<somethig>" : "your comment"` are ignored as user comments (avoid repeting `<something>` to not repeat the same key twice).

`Volume::applyDOFConfig` traverses the hierarchy and applies matching rules to non-pseudo volumes. A `defaults` object, if present, is applied first to every volume; rules then run in order, and a later matching clause replaces the earlier DOF set of the same kind. Patterns are matched against symbolic volume names, as regular expressions, also retried with an implicit leading `.*`. Since `.*` can match `/`, use exact-depth patterns when targeting a hierarchy level.

## Rigid-Body DOFs

The `rigidBody` clause frees an arbitrary subset of a volume's six DOFs (`TX,TY,TZ,RX,RY,RZ`, in
the volume's local frame), all the others of that volume staying fixed; a volume not matched by any
rule keeps none free. Free several at once by listing them together, e.g. the two in-plane shifts
plus the in-plane rotation of every ITS stave:

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TX", "TZ", "RZ"] }
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
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TY"] }
   ```
   Here the half-barrel's `TY` and the mean `TY` of its staves are tied together automatically; a stave
   free only in, say, `RZ` instead would still contribute to that same constraint via `j(TY, RZ)`.

   The equation is written only if some contributor really realises the parent DOF, i.e. the largest
   `|j(iDOF, jDOF)|` over the contributors' free DOFs is at least `Volume::MinChildrenCoupling` (1e-2).
   Otherwise the parent DOF is not degenerate with the children (case 3 below) and the equation would
   consist only of the leakage of small frame misrotations. Example: OB chips free only in `TX` under a
   module free in `TX,TY,TZ`. The module `TY`/`TZ` equations would be its `TX` equation scaled by
   10⁻⁷–10⁻⁵, and pede stops on the singular constraint matrix. The same cut applies to `pinChildrenMean`.

2. **Parent DOF fixed, descendants' DOF free.** There is no parameter at the parent level to be
   degenerate with, so **no constraint is written by default**: the children (and the hierarchy below
   them) are free to drift together in that DOF, which may be a genuinely unconstrained or only
   weakly constrained mode of the fit (e.g. nothing anchors the overall `TY` of a detector whose
   envelope is fixed but whose staves are all free in `TY`). Use `pinChildrenMean` to request the same
   mean-zero constraint explicitly even though the parent DOF is fixed — this **does not** turn the
   parent DOF into a fit parameter, it only pins the common mode of what's below it to the parent's
   (fixed) position:
   ```json
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": "fixed", "pinChildrenMean": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TY"] }
   ```
   Without the `pinChildrenMean` clause, this configuration leaves the mean `TY` of the staves
   unconstrained; with it, the mean is pinned to the (fixed) half-barrel position.

A third, unproblematic case is a free parent DOF with **no** free descendant DOF realising it (see
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

### `pinChildrenMeanNoCross`: the same pin without the jacobian cross-terms

The transported form above sums *every* free DOF of each child weighted by a jacobian element, and
those elements carry the lever arm from the child to the parent. That is the correct operator when
the parent's DOF is a fit parameter — the degeneracy genuinely is "parent pose vs transported
children pose". It is the **wrong** operator for pinning a collective drift, because the fit can
satisfy it with whichever contributing DOF is cheapest in χ², and a weakly-observable rotation with a
large lever arm is by far the cheapest. Measured in the ITS (see
[AlignmentWeakModes_handover.md](AlignmentWeakModes_handover.md)): in `DOF TY for …/ITSUStave0` the two
half-staves moved −119 and −75 µm, and the equation was satisfied to 5·10⁻⁴ µm because each took an
RZ of ±7.5 mrad — rotation about the sensor normal, which displaces the measured coordinates by
almost nothing but enters the equation with a coefficient of 0.645 against 0.5 for TY. The pin did not
prevent the 97 µm drift; it manufactured a 0.43° rotation to pay for it. At the OB half-barrel the
coefficient ratio is 1.724 vs 0.042 — a factor of 41, the stave radius in cm.

`pinChildrenMeanNoCross` writes the same equation with the cross-terms dropped: the plain mean of the
children's *own* same-named DOF, which no other DOF of the child can pay for.

```json
{ "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]", "rigidBody": "fixed",
  "pinChildrenMeanNoCross": ["TX", "TY", "TZ", "RX", "RY"] }
```

It takes the same values as `pinChildrenMean` (a bool, `"all"`, or a list of DOF names), has the same
"replace, don't merge" semantics, and its flags are independent — so both clauses may appear on one
volume and each writes its own equation. **In practice, use one or the other per volume**: over the
same set of contributors the two span largely overlapping subspaces, and requesting both at a level
with 5-6 free DOFs per child makes the constraint set rank-deficient (pede then reports a rank below
the number of constraint equations). Keep the transported form where the parent's DOF is free, and
use the cross-term-free form to pin a fixed parent's children.

Two things it deliberately cannot do, because they are not collective *drifts* of any single DOF:

- a global rigid **shift** or **tilt** of the whole branch is a cos φ / sin φ pattern in each child's
  local DOFs, so its plain mean vanishes and this clause is blind to it (the transported form does see
  it — that is what the cross-terms are for);
- any deformation mode — use [`globalModes`](#constraining-global-deformation-modes).

A pinned DOF that no descendant has free is skipped with a warning, which is the expected message
when, say, the stave RZ below an OB half-barrel is fixed.

## Weighting the Mean by Statistics

Both forms above split the mean equally (`1/N`) over the contributors, regardless of how much data
actually constrains each one. That is usually harmless, but not when the contributors' own statistics
differ by an order of magnitude or more — e.g. an ITS IB end chip (chip0/chip8 of a stave) sees far
fewer hits than the chips in the middle, since only the high-|tgλ| tail of tracks reaches that far
along the stave. The equal-weight mean still ties that chip's DOF to its well-determined siblings on
exactly the same footing as if it had just as much data, which can leave its own, genuinely different
correction under-applied: the DOF's formal pede error looks deceptively small (set mostly by its
high-statistics neighbours through the shared constraint) while its best-fit value only partially
reflects what its own, scarce hits call for.

`AlignParams.useStatWeightedMeanConstraints` (default `false`) switches every mean-of-children
equation — automatic or `pinChildrenMean`/`pinChildrenMeanNoCross`, transported or cross-term-free —
from the equal `1/N` split to one weighted by each contributor's own accumulated statistics
(`Volume::getDataCounter()`, the number of tracks that actually contributed to that volume), normalised
to sum to 1 within that equation. It requires per-volume track counts from an actual data pass. The
steering files (`mp2tree.txt`/`mp2con.txt`/`mp2param.txt`) are always written from `endOfStream`
(never from `init`), precisely so these counts — whichever of the two ways below they come from — are
available by the time `mp2con.txt` is written:
- With `ProcessData` (i.e. `MilleData,MilleSteer`), the counts are this run's own, accumulated as its
  TFs are processed and final once `endOfStream` is reached.
- With `MilleSteer` requested standalone (no data at all — the usual way to try a more restricted
  `dofConfigJson` against mille data already recorded with a wider-open one, reusing the same
  `mp2data.bin`), there is no data pass of its own. Point `AlignParams.statInFile` at the
  `ProcessingStats` ROOT file (`AlignParams.statFile`) a *previous* `ProcessData` run wrote (needs
  `AlignParams.volumeStatistics=true`, the default) to load its counts back instead. `statInFile` is
  consulted whenever it is needed — either `AlignParams.useStatWeightedMeanConstraints` is `true`, *or*
  `dofConfigJson` has at least one `weightedMean` rule asking for `true` on some DOF, since that
  per-rule override takes effect regardless of the global switch and would otherwise silently fall back
  to an equal split for lack of counts. Without a usable `statInFile` in either case, the write falls
  back to an equal split, with a warning.

Pipelined (`--pipeline`) runs use lane 0's own counts only, not summed with the other lanes': every
weight is normalised within its own equation, so a common scale across all its contributors cancels
out, and lane 0 alone is already equivalent to the full dataset for this purpose, given the TFs are
evenly distributed across lanes.

The `weightedMean` clause overrides this global default per volume and per DOF, the same "replace,
don't merge" way as the other clauses — a bool (all DOFs), `"all"`, or an array of DOF names (every
other DOF of that volume is then explicitly *not* weighted, like `pinChildrenMean`'s array form):

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "weightedMean": ["RY"] }
```

A DOF never touched by any matching `weightedMean` rule keeps falling back to
`AlignParams.useStatWeightedMeanConstraints`.

**Prefer the per-rule `weightedMean` clause over the global switch.** The global switch changes the
default for *every* mean-of-children equation in the hierarchy at once, including ones whose
contributors are not comparable in the way the feature was designed for. A mean-of-children equation
sums contributors that are siblings under the same parent, and the statistics-weighted split is only
meaningful when those siblings' track counts are proportional to how well each one's own correction is
determined — true for same-size, same-role volumes like the chips of one stave (the motivating case:
an IB end chip's far fewer hits versus its middle-chip siblings). It stops being meaningful once a
single equation mixes contributors of different granularity or role, which the global switch
reaches just as readily as the intended ones:
- `ITS`'s own `pinChildrenMeanNoCross` sums the two IB half-barrels *and* the six OB half-barrels in
  one equation. Their track counts differ by detector-geometry, not by alignment-quality, reasons
  (far more OB sensors/hits than IB) — weighting by them does not express "trust the better-determined
  side more," just "trust whichever side happens to have more silicon."
- A `pinChildrenMeanNoCross` whose matched contributors span deliberately asymmetric roles (e.g. a
  config giving one branch much coarser segmentation than its sibling) has the same problem: counts
  there differ by design, not by how well-determined the correction is.
- Even a structurally uniform equation can be swayed by something unrelated to alignment precision —
  e.g. instantaneous luminosity or trigger-selection differences slicing an otherwise symmetric
  detector's track yield unevenly run-by-run — silently shifting which contributor the mean constraint
  trusts, in a way that has nothing to do with this run's alignment quality.

None of this is diagnosed or validated end-to-end; the feature was designed and tested against the
specific IB-chip-in-a-stave case. Turning on `AlignParams.useStatWeightedMeanConstraints=true` applies
the same reasoning everywhere at once, including to equations where it hasn't been checked and may do
the opposite of what's intended. Scope it instead with a `weightedMean` rule matched only to the
volume(s) actually diagnosed (e.g. the IB stave's mean-of-chips equation), leaving every other
equation's default untouched; `statInFile` is consulted whenever any rule anywhere needs it (see
above), so a scoped rule gets the same track-count loading as the global switch without its reach.
Since it changes only the *weights* of an existing equation, not which DOFs are free or how many
equations are written, it can't by itself turn a
rank-deficient or missing constraint into a valid one, or vice versa.

## Constraining Global Deformation Modes

`pinChildrenMean` and the automatic mean constraint are both *rigid-body* constraints: 6 scalar
equations per parent, removing the net translation and rotation of that parent's set of free
descendants and nothing else. A coherent deformation of the barrel as a whole — a uniform radial
scaling, an elliptical squashing of each layer, a twist — is not a net rigid motion of any single
parent, carries no net rigid component, and therefore passes straight through it, however
`pinChildrenMean` is set. Such modes are also (very nearly) null directions of the track χ², so
Millepede takes large amplitudes in them for a negligible χ² gain. This is a different operator, not
a different setting of the same one.

The `globalModes` clause flags the volumes whose free rigid-body DOFs take part in one global
equation per mode, written by `Volume::writeGlobalModeConstraints` (called once on the root from
`writeRigidBodyConstraints`), forcing the total amplitude of that mode over the whole flagged set to
vanish:

```json
{ "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+",
  "globalModes": ["radialScaling", "elliptical", "twist"] }
```

| mode | displacement field | equations |
|---|---|---|
| `radialScaling` | Δr = ε·r (uniform transverse dilation, absorbed by the 1/p<sub>T</sub> scale) | 1 |
| `zScaling`      | Δz = ε·z | 1 |
| `curl`          | Δφ = ε (rigid rotation about the beam) | 1 |
| `twist`         | Δφ = ε·z (torsion along the beam) | 1 |
| `elliptical`    | Δr = ε·cos2φ and Δr = ε·sin2φ (r-independent amplitude) | **2** |
| `bowing`        | Δr = ε·z² (radial sagitta along the beam) | 1 |
| `radialHarmonic` | Δr = ε·r·cos(ℓφ) and ε·r·sin(ℓφ), ℓ = 1…`GlobalHarmonicLMax` (4) | **8** |
| `tangentialHarmonic` | rΔφ = ε·r·cos(ℓφ) and ε·r·sin(ℓφ), same ℓ range | **8** |
| `globalTilt`    | rigid rotation of the branch about the global X and about the global Y axis | **2** |
| `radialPolar`   | Δr = ε·r^p·t^k, t = z/r, (p,k) ∈ {(0,0),(0,1),(0,2),(1,1),(1,2)} | **5** |
| `tangentialPolar` | rΔφ = ε·r^p·t^k, same (p,k) | **5** |
| `longitudinalPolar` | Δz = ε·r^p·t^k, (p,k) ∈ {(0,1),(0,2),(1,0),(1,2)} | **4** |
| `globalShift`   | rigid translation along the global X, Y, Z axis (`globalShift_x/_y/_z`) | **3** |

The two harmonic families are the generic form of the azimuthal modes, with the amplitude explicitly
proportional to r because that is the measured shape: the fitted ε = amplitude/r of the surviving ITS
modes is coherent to 2-6 % across the IB and the OB, two independently mounted structures spanning
r = 2.3 → 39 cm. At ℓ = 0 they reduce exactly to `radialScaling` and `curl`, which is why their orders
start at ℓ = 1, and `tangentialHarmonic2` is the companion that `elliptical` — a purely radial field —
does not cover. `globalTilt` supplies the two transverse partners of `curl`; it is needed as an
explicit constraint because a tilt shows up as a child TZ pattern ∝ −r·cos φ whose mean vanishes both
over a half-shell and over the full circle, so no `pinChildrenMean*` of the TZ DOF can see it.

The three polar families are the modes that depend on the polar angle rather than on φ. A track from
the vertex is at z ≈ r·t at radius r, t = tgλ, so a field r·f(t) is absorbed by its φ (rφ component),
its tgλ (z component) or, with B ≠ 0, its curvature (radial component): these p = 1 members are null
directions whatever the vertex does. A field f(t) (p = 0) shifts the track's d0 or z0 instead, which
only the vertex constraint can see. The members which are already modes of their own are left out:
radial (1,0) is `radialScaling`, tangential (1,0) is `curl`, longitudinal (1,1) = z is `zScaling`, and
longitudinal (0,0) is a rigid z shift. Longitudinal (1,0), Δz = ε·r, is the telescope mode, and radial
(0,1), Δr = ε·z/r, is the cone that showed up as the IB chip-TY dipole. Like `zScaling` these families
have little or no handle on stave-level DOFs (all at z ≈ 0, see "How the mode coefficients are
built"): they matter as soon as DOFs distributed in z (modules, chips) are free, which is exactly when
they must be requested. Note that `bowing`
(Δr = ε·z², one ε for all layers) is *not* a null direction; its null counterpart is
`radialPolar_r0t2`, Δr = ε·z²/r². The equations are named `<family>_r<p>t<k>`, e.g.
`radialPolar_r0t1`. They are only the axisymmetric, low-order members of a family that is
really a function of (φ, t): with free DOFs at module/chip granularity a fit can still take large
amplitudes in, e.g., Δz = r·g(φ, t) with g not covered here (measured in `ITS_RB_v3`, see
AlignmentWeakModes_handover.md section 6). The `rayModes` clause below covers the whole family.

The clause takes a boolean (all equations), the string `"all"`, or an array of names, each of which may
be a mode name (all its equations), a single order of a harmonic family (`radialHarmonic3`), or a
single equation (`elliptical_cos2phi`, `tangentialHarmonic2_sin2phi`, `globalTilt_y`,
`radialPolar_r0t1`). Like
`rigidBody`/`calib`/`pinChildrenMean` a later matching rule replaces the flags of an earlier one.
It creates no DOF set and frees nothing: it only adds terms for DOFs that are already free, so it
must follow the rule that frees them. Flagging several nested levels at once is legitimate and is not
double counting: the equation then constrains the total amplitude summed over all of them, which is
the degenerate quantity, regardless of which level realises it.

### How the mode coefficients are built

The coefficients come from the volumes' own `getL2G()`, not from nominal radii. Two definitions exist,
selected by `AlignParams.modeSensorProjection` (`Volume::ModeEquationOptions`, set in
`AlignmentSpec::writeMilleSteerFiles`):

- **Sensor projection (default, `true`).** The equation is the projection of the displacement of the
  sensors on the mode. Each sensor is sampled at the 8 corners **p** of its box (the box of its geometry
  volume, placed in its L2G frame; its origin alone if it has none). A unit value of the free DOF *k* of
  a volume with origin **O** and local axes **e**<sub>k</sub> moves **p** by **d**(**p**) =
  **e**<sub>k</sub> (`TX,TY,TZ`) or **e**<sub>k</sub> × (**p** − **O**) (`RX,RY,RZ`), and its coefficient is
  Σ<sub>p</sub> **u**(**p**)·**d**(**p**) over the sample points of the sensors *of its subtree*, i.e.
  **U**·**e**<sub>k</sub>, resp. (**M** − **O** × **U**)·**e**<sub>k</sub>, with **U** = Σ **u**(**p**) and
  **M** = Σ **p** × **u**(**p**), sums which add up over the subtrees (`Volume::collectProjectedModeTerms`).
  Since the sensor displacement is the sum of the contributions of all levels, the equation is the same
  whichever level realises the motion: a half-barrel shift gives the same total at stave, half-stave or
  module level (`globalShift_y` over Lr3 HB0 in `ITS_RB_v17`: coefficient sums −7.69, −7.69, −7.71 per
  level), and a volume counts in proportion to its sensors. The rotations get their lever arm from the
  sensors' extent: the rotation of a stave is weighted by its length, that of a chip about its own centre
  by its few cm. This is the metric of the B amplitude of `WeakModesDiag` (a fit of the sensor
  displacements), up to the sampling of the sensor rotations.
- **Origin evaluation (legacy, `false`).** For every free DOF, the projection of the mode's field
  **u**(x) (for `TX,TY,TZ`) or of the local rotation **ω**(x) = ½∇×**u** it induces (for `RX,RY,RZ`),
  evaluated at the volume's *origin*, onto the corresponding local axis; every volume counts with
  weight 1, whatever its size. The rotation part is what makes `twist` constrain the OB staves' `RY`
  rather than their translations. With several nested levels flagged, the equation then weights the
  levels by their number of volumes rather than by the displacement they cause (the same Lr3 HB0 shift:
  −7.69, −15.4, −61.6 at stave, half-stave and module level), so its content depends on which level
  carries the motion. It reproduces the `mp2con.txt` written before the option existed exactly.
  Metric A of `WeakModesDiag` still uses these coefficients.

Each equation is homogeneous, so it is normalised to a largest coefficient of 1, and terms below 1e-9
of that are dropped as rounding residues of the L2G matrices.

A mode with no handle on the flagged level is skipped with a warning rather than silently. With the
origin evaluation, `zScaling`, `bowing` and the polar families project on a volume's DOFs with
coefficients ∝ z (z²) of the volume's origin, which vanish for an ITS stave — a stave is centred at
z = 0 — so these need free DOFs that are *distributed* in z (module or chip level) to be observable at
all. With the sensor projection the modes even in z (`bowing`, `radialPolar_r0t2`, ...) do reach the
stave DOFs, through Σ z² of the stave's sensors; the odd ones (`zScaling`, `radialPolar_r0t1`, ...)
still cancel over a stave symmetric in z. `globalModes` is ignored, with a warning, on a volume that
is not rigid-body alignable.

`globalTilt` is nearly the same equation as a `pinChildrenMean` (transported form) of the parent's
`RX`/`RY`, so requesting both makes the constraint set rank-deficient. Prefer `globalTilt`: the pin's
contributor walk stops at the first descendant with any free DOF, which in the ITS excludes the 48 IB
staves, while `globalTilt` flagged at the stave level covers all 192.

Adding a mode to `Volume::GlobalMode` propagates automatically to the `WeakModesDiag` output and to
`macro/drawWeakModes.C`, which resolve the mode columns by name; a `static_assert` after the class
guards `Volume::NGlobalEq` against going out of step.

`globalShift` (Δ**x** = ε·**e**, **e** = x̂, ŷ or ẑ) is the translation partner of `globalTilt`. Over a whole
detector the top-level pins remove it already; it exists for the `groupModes` clause below. In
`WeakModesDiag` its B amplitude duplicates the rigid shift of the fit (like `curl` and `globalTilt` the
rigid rotation) and carries no information.

### Per-group mode equations: `groupModes`

A `globalModes` equation is one equation over the *whole* flagged set, whatever the rules: flagging
each layer with its own rule still gives a single equation over all of them. A deformation which
cancels between the parts of the flagged set therefore escapes it. In `ITS_RB_v15`/`v16` the OB bowed
radially by +12/+13 µm on Lr3/Lr4 and −15…−17 µm on Lr6 at the stave ends, which meets the single
`bowing` equation by cancellation between layers. The `groupModes` clause writes the equations per
group instead:

```json
{ "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]", "groupModes": ["globalShift_x", "globalShift_y"] },
{ "match": "ITS/ITSULayer3/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+",
  "groupModes": { "group": "OB-Lr3", "modes": ["bowing"] } }
```

- **Group.** With a plain mode list (the same syntax as `globalModes`) every matched volume is a group
  of its own: the first rule gives 2 equations for each of the 8 OB half barrels. With
  `{"group": "<name>", "modes": <list>}` all volumes given the same name form one group: the second rule
  gives one bowing equation over all the staves of Lr3. A group name can collect volumes from several
  rules, which is how a group can span a level that has no volume of its own (the ITS has no layer
  volume: the half barrels sit directly under `ITS`). All members of a group must request the same
  modes.
- **Group names** are free labels with no convention: they are not parsed nor matched as patterns, only
  compared as strings, and serve to collect the volumes and to name the equations in `mp2con.txt`
  (`OB-Lr3` could be any other string). Membership is decided only by the `match` of the rules. An
  unnamed group is keyed by the volume's symName, so do not choose a name equal to the symName of a
  volume (e.g. `ITS/ITSULayer3/ITSUHalfBarrel0`), which would merge it with that volume's own group.
- **Terms.** The equation sums the free rigid-body DOFs of the *whole subtree* of every volume of the
  group, at all levels, with the same coefficients as `globalModes` (see "How the mode coefficients are
  built"). Nothing has to be
  flagged below, and the matched volume itself needs no free DOF (a fixed half barrel is fine).
  Summing every level is what makes it impossible to escape: a `pinChildrenMean*` sees only the nearest
  free descendants, so a motion it pins can reappear one level down. The OB half-barrel "scissor" did
  exactly that under the stave-level pin, through module TX + RY; `globalShift` per half barrel removes
  the net translation whichever level carries it.
- **Relative to the free motion of the group volumes** (`AlignParams.groupModesRelToRoots`, default
  `true`, needs `modeSensorProjection`). Before the projection, the mode is replaced, for each volume of
  the group, by its residual after a least-squares fit, over the sample points of that volume's sensors,
  of the rigid motions of the volume's own *free* DOFs (`Volume::removeFreeRigidMotion`). The free DOFs
  of the group volumes then drop out of the equation exactly, and it constrains only the deformation of
  their content. For the per-layer stave bow of `ITS_RB_v17` (staves free in TX, TY, TZ, RY) the field
  z²·r̂ acts as (z² − ⟨z²⟩)·r̂ per stave: the radial offset of each stave, which the data measure,
  stays free, while the module TY coefficients sum to ≈ 0 per half stave. Without it, the projected
  bow would load the stave TY with Σ z² of its sensors and pin the stave radii together with the bow;
  the origin evaluation escaped this only because a stave origin sits at z = 0, and was centred at module
  level only through the equal-weight child-mean constraints. A group volume without free DOFs (the
  fixed OB half barrels of the `globalShift` groups) is left unchanged. A mode which is entirely a free
  rigid motion of the group volumes (e.g. `globalShift` on free half barrels) leaves nothing to
  constrain and is skipped with a warning; fix the DOF or set the switch to `false` to constrain it.
- **Semantics.** A volume carries one `groupModes` clause: a later matching rule replaces the modes and
  the group name of an earlier one. To put a volume in two kinds of group, use different levels, as
  above (half barrels for the shift, staves for the per-layer bow). The equations are written by
  `Volume::writeGroupModeConstraints` after the global ones, named `Group mode <eq> over <group>`, with
  the same normalisation and pede grouping anchor.
- **Rank.** Group equations must not tile a global one: if the groups of a mode together cover exactly
  the volumes of a `globalModes` set requesting the same mode, the global equation is the sum of the
  group ones, an exact rank defect. Drop the mode from the global list (`ITS_RB_v17` removes `bowing`
  there). Likewise `globalShift` on a volume is almost the transported `pinChildrenMean` of its TX/TY
  at the first level, so use one or the other. Per-group members of the `rayModes` span (`zScaling`, the
  p = 1 polar members) inherit the degeneracy warning of the global ones.
- **Measured directions.** A group equation on a direction the data do measure costs χ². Keep the
  groups to the weak combinations (net shifts of a part, its bow), not to shapes such as a per-layer
  constant radius (`radialPolar_r0t0`), which dz vs tgλ measures.
- **Not diagnosed yet.** `WeakModesDiag` still reports only the modes over each detector branch, not per
  group.

## Constraining Ray Modes

A track from the vertex crosses every layer at the same azimuth φ and the same t = z/r (its tgλ). A
sensor displacement r·F(φ, t)·**e**, with **e** = φ̂, ẑ or r̂, is therefore absorbed by the track's φ,
tgλ or (B ≠ 0) curvature *for any function F*: a function space of null directions, not a short list
of global modes. Rigid-body DOFs at stave level, all centred at z ≈ 0, cannot realise it, but free DOFs
distributed in z (modules, chips) can.

The `rayModes` clause flags a branch, typically the detector envelope. The sensors under it are
binned in (φ, t) cells, and per cell and component one equation is written:

  Σ over the sensors s of the cell of (**e**·**d**ₛ)/rₛ = 0

Here **d**ₛ is the linearised displacement of the sensor centre by every free rigid-body DOF on its
path to the root. The equation states that the mean angular displacement of the cell, seen from the
vertex, vanishes. A common-ray motion of the cell is removed. So is one other direction per cell,
the mean angular motion of its sensors, which for a single layer would be a measurable motion; hence
cells whose sensors span less than `rSpanMin` in radius are skipped. The per-DOF displacement is
built by finite differences from the same `AlignParam` matrix the corrections are applied with, so
it follows the fit's rotation convention by construction.

```json
{ "match": "ITS", "rigidBody": "fixed",
  "rayModes": { "nPhi": 24, "nT": 16, "tMax": 4.0, "components": ["rdphi", "dz", "dr"] } }
```

| key | default | meaning |
|---|---|---|
| `nPhi`, `phi0` | 24, 0 | azimuthal cells over [φ0, φ0 + 2π) |
| `nT`, `tMax` | 16, 4.0 | cells in t = z/r over [−tMax, tMax]; sensors beyond go to the edge cells |
| `rSpanMin` | 1.0 cm | minimal radial span of a cell's sensors |
| `rEdges` | [] | optional radial boundaries; each range gets its own cells (this forbids, e.g., a relative IB/OB ray motion, which tracks do measure) |
| `rMin`, `rMax` | 0, ∞ cm | only sensors within this radial range are binned |

The cells use the nominal t = z/r, i.e. a vertex at z = 0. A track from z_v crosses the layer at
radius r at t + z_v/r. With σ(z_v) ≈ 5 cm this shifts t by about ±2 at the ITS IB but only by
±0.13-0.25 in the OB. So at the IB a cell is not a set of rays, and its equation is not a pure
null direction. Use `rMin` (e.g. 10 cm for the ITS) to bin only the OB; see
AlignmentWeakModes_handover.md section 6 for the measured effect.
| `components` | all | any of `rdphi`, `dz`, `dr` |

`true` takes the defaults and `false` removes the flag. A flagged ancestor takes precedence over a
flagged descendant.

**Do not combine with the global modes that lie in the ray span.** That excludes `radialScaling`,
`curl`, `zScaling`, `globalTilt`, the two harmonic families, and the p = 1 members of the polar
families: each is, up to the discretisation, a weighted sum of the cell equations, so requesting it
as well makes the constraint set nearly degenerate. Still meaningful alongside `rayModes` are `twist`
(∝ r², not a ray field), `elliptical`, and the p = 0 polar members (∝ F(t), broken only by the vertex).

Like the global modes, every equation carries the pede grouping anchor (see
`Volume::addPedeGroupingAnchor`).

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
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": ["TX", "TY"] },
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "measurement": { "TX": [0.01, 0.005], "TY": {"value": -0.02, "sigma": 0.01} } }
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
a warning (under a literal pattern, without regex metacharacters) or silently skipped (under any other pattern, e.g.
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
The ITS envelope is the `ITS` alignable. The layer name is a path component, not a separate alignable entry. In the match patterns, layer IDs use `[0-2]` for the inner barrel and `[3-6]` for the outer barrel. Half-barrel, half-stave, module, and chip IDs in these paths are single digits, matched with `[0-9]`. Inner-barrel stave IDs are single digits; outer-barrel stave IDs can be one or two digits, so both exact-width forms are listed there. Avoid a trailing `.*` when matching a specific level: it also matches `/` and descendants.

### Inner Barrel (Layers 0–2)

The inner-barrel hierarchy has half-barrels, staves, and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]` |
| Stave | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]` |
| Chip | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]/ITSUChip[0-9]` |

### Outer Barrel (Layers 3–6)

The outer-barrel hierarchy adds half-staves and modules between staves and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]` |
| Stave, one or two digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]+` |
| Half-stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]+/ITSUHalfStave[0-9]` |
| Module | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]+/ITSUHalfStave[0-9]/ITSUModule[0-9]` |
| Chip | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]+]/ITSUHalfStave[0-9]/ITSUModule[0-9]/ITSUChip[0-9]+` |

For example, replace a `match` value with one of the strings above in a rule such as `{ "match": "...", "rigidBody": ["TX", "TY", "TZ"] }`. ITS3 uses different names, such as `ITS3CarbonForm`, and is not selected by these patterns.

## TRD Alignable Match Strings

TRD has a dummy `TRD_envelope` (not a geometry alignable), followed by alignable supermodules and chambers. The geometry log uses two digits for supermodule IDs and one digit each for stack and plane IDs:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TRD_envelope` |
| Supermodule | `TRD/sm[0-9]+` |
| Chamber | `TRD/sm[0-9]+/st[0-9]/pl[0-9]` |

## TOF Alignable Match Strings

TOF has a dummy `TOF_envelope` (not a geometry alignable), followed by alignable supermodules and strips. Both IDs in the geometry names are formatted as two digits:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TOF_envelope` |
| Supermodule | `TOF/sm[0-9]+` |
| Strip | `TOF/sm[0-9]+/strip[0-9]+` |

The envelopes are not rigid-body alignable; rigid-body clauses on them are ignored. Calibration clauses can be configured on these synthetic volumes.
