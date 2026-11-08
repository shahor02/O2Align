# To do

## ITS OB chips without data (entries)

Some OB chips have no data in the collision tracks, mostly the stave-end chips (Lr5 module 0 and module 6, at
|tgλ| ≈ 2). Their alignment parameters cannot be determined by the fit:
- pede fixes them through `entries 55`: 851 parameters in C1v19, 847 of them chip TX;
- 2.7k OB chips have fewer than 100 tracks in the checkResid samples.

These chips keep the v10n chip TY/TZ offsets, which are the large ones: their net z deviates from the module median by
210 µm (MAD). See `history.md` (2026-11-07, v10n OB chip TY/TZ). Freeing chip TY/TZ in the fit does not change them,
so a dedicated treatment is needed.

- Count the entries per parameter: add `printcounts 2` to the pede steering. The count from the binary files, which
  `entries` is compared to, then becomes the last column of `millepede.res`.
- Decide, from the counts, whether a lower `entries` cut (about 20, with chip priors) is safe for the chips
  constrained only by cosmics.
- For chips with no (or too few) entries, put them back at their ideal position within their parent (module) after
  the fit, rather than keeping the inherited offsets. Options for implementing this in `Regularize`
  ([Regularization.md](Regularization.md)):
  - A minimum-statistics threshold (e.g. `AlignParams.regularizeMinSensorStat`, in tracks from `ProcessingStats` /
    `statInFile`, or in pede entries from `millepede.res` with `printcounts`). A sensor below it is placed at its
    ideal position in the parent frame, so its local delta becomes the identity. The parent is moved with the rigid
    motion fitted to its other sensors. Unlike the present `Regularize`, which keeps every sensor in place, this
    moves the reset sensors in the lab. It must therefore be a separate switch, reported per sensor.
  - Choice of the DOFs to reset (e.g. only TY, TZ), keeping the others (e.g. TX), since the fit may determine some
    DOFs and not others.
  - Exclude the sensors below the threshold from the robust rigid-motion fit of the parent's children, so that their
    offsets do not pull the transferred common motion.
  - Alternative without a post-fit step: per-chip `measurement` priors centred on the value that cancels the
    inherited deviation from the ideal placement. This only acts on parameters pede does not fix, and needs the
    module-mean constraints weighted by statistics (`weightedMean`), so the priors do not fight the constraint.

## Regularize: common child rotations grow when the children translations are large

With `Regularize` applied to the v10n alignment (`AlignParams.regularizeParents` including `ITSUModule[0-9]+`),
the common rotation of the chips of a module gets larger rather than smaller. For example, on
`ITS/ITSULayer6/ITSUHalfBarrel1/ITSUStave23/ITSUHalfStave1`:

| module | chip rotation, v10n | chip rotation, regularized | mean chip \|t\|, v10n → regularized |
|---|---|---|---|
| Module0 | (−0.016, +0.039, +0.000) mrad | (−0.20, −0.25, +1.45) mrad | 210 → 159 µm |
| Module1 | (+0.032, +0.051, +0.000) mrad | (−1.32, −0.27, +2.16) mrad | 167 → 94 µm |

The angles are the local pitch, roll and yaw of the chips. In each module, all chips carry the same rotation, before
and after. The sensor placement is unchanged (closure verified).

**Cause.** The motion `e` transferred to the module is fitted to the corner displacements of its chips'
bounding boxes ([Regularization.md](Regularization.md#estimating-the-common-motion)). In v10n these displacements come
almost entirely from the chip translations, which are 100–300 µm and differ from chip to chip and between the two chip
rows. The pattern looks like a mrad-level rotation of the module. The common chip rotation (16–51 µrad) is a rotation
about each chip's own centre. Over a chip's half-size it moves the corners by less than 1 µm, and it is not a rigid
motion of the module anyway (the "turbine" case of [Regularization.md](Regularization.md)). The sensors are kept in
place, so every chip of the module gets the same counter-rotation `e⁻¹`. The common chip rotation therefore becomes
about `−ω(e)`.

**Check.** A least-squares fit of the module motion from the v10n chip deltas and the real chip boxes was made with
the translations only and with the rotations only. This fit is the first IRLS iteration. The values are in the module
frame:

| module | ω, fit with translations only | ω, fit with rotations only | observed change of the chip rotation (sign flipped) |
|---|---|---|---|
| Module0 | (−0.59, +0.21, +2.34) mrad | (+0.001, +0.003, 0) mrad | (−0.19, +0.29, +1.45) mrad |
| Module1 | (−1.52, +0.31, +2.16) mrad | (−0.002, +0.004, 0) mrad | (−1.35, +0.32, +2.16) mrad |

The translations reproduce the transferred rotation. For Module0, the L1 iterations then move away from the
least-squares starting point. The fit minimizes the summed corner displacement in cm, which counts a mrad chip
rotation as almost free.

**Options, if the children's rotations should stay small:**

- Leave the level out: exclude `ITSUModule[0-9]+` from `regularizeParents`, so the chip→module level is not touched.
  No code change is needed.
- A translation-only transfer for selected parents: fit with ω = 0 and transfer only T. This needs a new parameter,
  e.g. a regex `AlignParams.regularizeShiftOnlyParents`.
- A penalty on the children's rotations: add a term λ²·Σᵢ‖Δωᵢ‖² to the robust fit, Δωᵢ being the change of child i's
  own rotation implied by `e`, with λ a length (e.g. the child half-size or a few cm). A rotation then costs as much as
  the corner shift it would produce at the distance λ, instead of almost nothing.
