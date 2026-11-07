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
