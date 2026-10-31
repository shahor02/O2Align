# Regularization of an alignment

## What it does

A hierarchical alignment has a gauge freedom: a common rigid motion of all children of a parent
can equally be carried by the parent, with the same sensor placement in the lab. The fit does not
always choose the split sensibly. Children often end up with similar, same-sign and relatively large
corrections, while the parent carries little or something compensating. Examples are the chips of
an IB stave all shifted by ~70 µm, or IB staves at −3.8 mm compensated by their half-barrel. This
happens because the hierarchical mean constraint fixes only the net motion of the *increments* of
one pede iteration, not the accumulated alignment.

The `Regularize` output option rewrites an alignment a posteriori so that:

- every selected parent absorbs the **common rigid motion of its children**, estimated robustly
  so that outliers do not bias it;
- the children keep only what is not a rigid motion of the parent;
- the placement of every sensor in the lab is **unchanged**, to numerical precision.

Tracking is therefore not affected; only the bookkeeping is. The rewritten objects are easier to
read and interpret, and they are a better starting point when later iterations fix or constrain
some levels.

## Running it

`Regularize` can be combined with any `--output` option that does not process data (`MilleRes`,
`MilleSteer`, `WeakModesDiag`, `WeakModesDiagIni`). Combining it with `ProcessData`, `MilleData` or
`VerboseGBL` is fatal. It writes the usual `<DET>_<AlignParams.algParamsOutFile>`
(`ITS_alignment.root` by default):

- `--output MilleRes,Regularize` regularizes the combined result, i.e. the initial alignment plus
  the pede corrections. Only the regularized file is written.
- `--output Regularize` without `MilleRes` regularizes the initial (CCDB) alignment itself. No
  `millepede.res` is needed.

```bash
o2-dev-alignment-workflow -b --run --output Regularize --detectors ITS --disable-mc \
  --hbfutils-config <some>/o2_tfidinfo.root \
  --condition-remap "file://<dir>=ITS/Calib/Align,GLO/Config/GeometryAligned" \
  --configKeyValues "AlignParams.regularizeParents=ITSUStave[0-9]+|ITSUModule[0-9]+"
```

**The aligned geometry loaded (`GLO/Config/GeometryAligned`) must be the one obtained by applying
the input alignment (`ITS/Calib/Align`) to the ideal geometry.** The ideal matrices are derived from
the two as G<sup>id</sup> = (P<sup>p</sup>)<sup>−1</sup>·G<sup>p</sup>. With a mismatched pair they
are wrong, and the closure check described below would not notice it, because it is
self-consistent by construction.

### Parameters (`AlignParams.*`)

| parameter | default | meaning |
|---|---|---|
| `regularizeParents` | `""` | regex on the symbolic names of the parents to process, DOF-config style (also matched with an implicit arbitrary prefix), alternatives separated by `\|`. Empty means every alignable having children |
| `regularizeMinChildren` | 2 | skip the parents with fewer alignable children |
| `writeLocalAlignParams` | true | output convention (local or global deltas), independent of the convention of the input |

Notes on the selection:

- Parents are processed **deepest first**, for example chips→module, then modules→half-stave,
  half-staves→stave, staves→half-barrel and so on. The common part of one level is thus passed
  further up, and it stops at the shallowest selected parent.
- With everything selected, a part of the corrections ends up in the detector top volume (`ITS`).
  This does not matter for the sensors, but it moves the detector as a whole relative to the others
  in the bookkeeping. Restrict `regularizeParents` if you do not want that.

## Parents and children

Parents and children are taken from the **TGeo alignable hierarchy**, not from the `Volume` tree of
the fit. The child of an alignable is its nearest alignable TGeo descendant. The alignables
considered are those at or below the non-virtual volumes of the detector branch. The tool therefore
works on any alignment object, whatever DOF configuration produced it. Consequences:

- A child **without** an alignment object takes part with zero correction relative to its parent:
  it follows the parent rigidly.
- To keep the sensors in place, every child of a processed parent gets an explicit object, and so
  does the parent itself if it had none. The output may therefore contain more objects than the
  input; new ones have alignable ID −1, which `GeometryManager::applyAlignment` does not use.
- Non-alignable TGeo volumes inside a parent, such as supports, move with the parent. They have no
  influence on tracking.

## Formalism

The notation is that of the derivation in [Detector.cxx](../src/Detector.cxx):

- G<sup>id</sup><sub>j</sub> is the ideal L2G matrix of volume j.
- P<sub>j</sub> is its cumulative global delta, so that G<sub>j</sub> = P<sub>j</sub>·G<sup>id</sup><sub>j</sub>.
- D<sub>j</sub> = P<sub>j</sub>·P<sub>a(j)</sub><sup>−1</sup> is its global delta, a(j) being its
  nearest ancestor carrying an object.

The operation acts on the cumulative deltas, which is what makes it independent of the
local/global convention of both input and output. The steps are:

1. The input objects, of either convention and combined with the pede corrections if `MilleRes` is
   given, are turned into the cumulative deltas P<sup>c</sup> of the volumes having an object
   (`Detector::MP2AlignParams`, steps 1–3).
2. For each selected parent p, the motion of child i relative to p, expressed in the frame of p, is
   E<sub>i</sub> = G<sub>p</sub><sup>−1</sup>·D<sub>i</sub>·G<sub>p</sub>, with
   D<sub>i</sub> = P<sub>i</sub>·P<sub>p</sub><sup>−1</sup>. A point x of the child, in parent
   coordinates, is displaced by u = E<sub>i</sub>(x) − x.
3. The parent is moved by a rigid motion e in its own frame, G<sub>p</sub> → G<sub>p</sub>·e, that is
   P<sub>p</sub> → G<sub>p</sub>·e·G<sub>p</sub><sup>−1</sup>·P<sub>p</sub>. The P<sub>i</sub> of the
   children are **kept**, so their placement, and with it every grandchild and every sensor, is
   unchanged. Only the objects of p and of its direct children change:
   - in the global convention, D<sub>p</sub>' = P<sub>p</sub>'·P<sub>a(p)</sub><sup>−1</sup> and
     D<sub>i</sub>' = P<sub>i</sub>·P<sub>p</sub>'<sup>−1</sup>;
   - in the local convention, d<sub>i</sub>' = L<sub>i</sub><sup>−1</sup>·e<sup>−1</sup>·L<sub>i</sub>·d<sub>i</sub>,
     with L<sub>i</sub> the ideal placement of the child in the parent frame.
4. The modified P<sup>c</sup> are converted back to objects in the requested convention, ordered
   parents first.

### Estimating the common motion

The motion e is written to first order as x → x + T + ω × x. It is estimated from a point cloud:
the **8 corners of the TGeo bounding box of every child**, at their ideal place in the parent
frame. For an assembly, such as an ITS stave, the box is computed explicitly, because TGeo leaves it
at 0 until it is requested. (T, ω) minimizes

  Σ<sub>i,k</sub> ‖ u<sub>ik</sub> − T − ω × x<sub>ik</sub> ‖

the sum of the Euclidean norms (not of the squares) of the residual displacements of the corners.
This is the L1 / geometric-median counterpart of a least-squares rigid fit. Children moving
differently from the majority affect it like outliers affect a median, not like they affect a
mean.

How it is solved (`fitRobustRigidMotion` in [AlignGeomUtils.cxx](../src/AlignGeomUtils.cxx)):

- The minimum is found by iteratively reweighted least squares. The first iteration is the plain
  least-squares fit; afterwards each point gets the weight 1/max(‖r‖, 10<sup>−7</sup> cm).
- The 6×6 normal system is solved by complete orthogonal decomposition. A combination not
  constrained by the points, such as a rotation about a line that all points lie on, is left at 0
  and stays in the children.
- The fit is linearized. e is therefore built as an exact rotation (Rodrigues) and composed, and
  the step is repeated, at most 5 times, until it vanishes (10<sup>−14</sup> rad, 10<sup>−12</sup> cm).

### Why not the componentwise median of the child parameters

The most obvious definition is to make the median of every child DOF (TX…RZ) vanish. That is
achievable: take ω as the median of the children's rotations, then T as the median of their origin
displacements after removing the lever arm ω × r<sub>i</sub>. It goes wrong, however, when the
children share a rotation **about their own axes** that no rigid motion of the parent produces.

Example: all the staves of one half-barrel may carry RZ ≈ −22 mrad (a "turbine" pattern).
Transferring it to the half-barrel as a rotation about the beam axis forces lever-arm compensations
of ω·r ≈ 22 mrad × 19.6 cm onto the staves. They end up with ±4.7 mm translations instead of ±1.8 mm.
The sensors are still in place, but the children's corrections get larger, not smaller.

The bounding-box field fit sees what such a pattern actually does to the volumes: a turn of a stave
about its own long axis displaces its corners by only ω × half-width. The −21 mrad RZ therefore
stays in the staves, and their translations drop to 0.1–0.4 mm. The price is that the medians of
the individual child parameters are close to zero rather than exactly zero.

Each child carries the same 8 points, whatever its size, so it counts once.

## Checks and logs

- **Closure.** After the conversion to objects, these are reapplied, parents first, to the ideal
  geometry. The placement of every sensor, meaning every alignable without alignable children, is
  compared with the placement before the regularization. The result is logged at info level, or at
  error level if the maximum deviation exceeds 10<sup>−9</sup> cm in translation or
  10<sup>−11</sup> in a rotation-matrix element. A typical value is 2·10<sup>−10</sup> cm.
- **Summary (info).** The log gives the number of parents processed, the largest transferred shift
  and rotation, and the mean corner displacement of the children relative to their parents before
  and after.
- **Per parent (debug).** For each parent the log gives the transferred T and ω and the mean child
  displacement before and after.

To inspect the result, compare the input and output vectors object by object. Remember that the
`AlignParam` angles (psi, theta, phi) are in **radians**.

## Code map

| piece | where |
|---|---|
| option `OutputOpt::Regularize`, consistency check | [AlignmentSpec.h](../include/O2Align/AlignmentSpec.h), [alignment-workflow.cxx](../src/alignment-workflow.cxx) |
| invocation, standalone or with `MilleRes` | `AlignmentSpec::writeAlignParams`, `AlignmentSpec::endOfStream` ([AlignmentSpec.cxx](../src/AlignmentSpec.cxx)) |
| combined cumulative deltas, extra output entries, closure check | `Detector::MP2AlignParams` and helpers ([Detector.cxx](../src/Detector.cxx)) |
| alignable hierarchy, robust fit, transfer | `geomutils::collectAlignablesBelow`, `collectAlignableChildren`, `transferChildrenCommonMotion` ([AlignGeomUtils.cxx](../src/AlignGeomUtils.cxx)) |
| parameters | `regularizeParents`, `regularizeMinChildren` ([Params.h](../include/O2Align/Params.h)) |
