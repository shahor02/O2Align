# Handover: pede fails to converge on B=0 multi-track vertex records

## Context

- **Data:** pp run 562862 (LHC25ab), B=0. GBL records are written by `o2-dev-alignment-workflow` and fitted with pede (Millepede-II).
- **What works:** single-track records, where each track is refitted with the vertex as one measurement (`testNVSingle`), converge.
- **What fails:** multi-track records, where all tracks of a vertex form one composed GBL trajectory with a common vertex (`buildGBLVertex`, `AlignParams.useMultiTrackPVConstraint > 0`), stop with `Too many rejects (>33.3%) - stop`.

### Test outputs

All are in `/data/devAlg/`.

| Directory | Test |
|---|---|
| `testNVSingle` | single-track reference (converges) |
| `testNVMulti` | multi-track, original code |
| `testNVMultiSingle` | mixed: multi-track plus single tracks for the rest |
| `testNVMulti_sign1` | multi-track with the inner transformation sign flipped (generation log `rr1.log`) |
| `testNVMulti_nooutlier` | test 4.1: `outlierdownweighting` and `dwfractioncut` disabled |
| `testNVMulti_ITSEnvFix` | test 4.2: ITS envelope DOFs fixed, down-weighting on |

### Configuration of these runs

- **Parameters** (`mp2param.txt`): 1242 ITS parameters plus the mean vertex X/Y/Z (PVT), all free.
- **DOF config:** `doc/ITS_RB_from_staves.json`, with the ITS, HalfBarrel and Stave levels free.
- **Mean-vertex prior** (from CCDB): position (−0.0542, −0.0047, −0.085) cm; σx ≈ 51 µm, σy ≈ 46 µm, σz ≈ 5 cm.

## Findings

### F1. The vertex sign conventions and the zeroed curvature at B=0 are correct

- `computeVertexTransformation` returns the derivative of the predicted (Y,Z) at the vertex point with respect to the vertex position. For a straight track through vertex V: `Y = V_Y + y'·(frame.x − V_X)`.
- `addMeanVertexPrior` puts the mean vertex on the measurement side, so its global derivatives are `−trans`. This is correct.
- `makeZeroFieldInnerTransformation` builds a 5 × (3 + 2n) matrix:
  - the q/pt row is null;
  - the snp and tgl entries are unity;
  - `trans` fills rows 3–4.

  GBL builds `innerTransDer = matLocalToFit · innerTrans`, so the curvature is fixed to 0 at every point.
- Flipping the sign of every inner transformation (commit be81ebe, `innerTrans = -innerTrans` in `buildGBLVertex`) only replaces V by −V. The record chi² values stay identical and only the sign of the fitted vertex becomes wrong. The `testNVMulti_sign1` histograms match `testNVMulti`.
- The flip is no longer in the working tree, which is correct. If be81ebe survives on any branch, it should be reverted.

### F2. The divergence comes from a free gauge mode combined with an inconsistency in pede's outlier weights

**The gauge mode.** A global ITS translation or rotation, together with a shift of the mean vertex, is absorbed by the common vertex fitted in each record. Single-track records pin this mode, because there the refitted vertex is a measurement with no global derivatives.

**The pede inconsistency** (`pede.f90`, local fit near lines 4700–5060 and global update near lines 5110–5220):

- The last local-fit iteration uses Cauchy weights (iterations > 3).
- The global update ("fourth loop") uses Huber weights and skips the b correction, commented "after local fit it is zero!".
- Cauchy weights are smaller than Huber weights, so the subtracted Schur term is over-estimated. This produces negative curvature and a spurious b component, both along directions the local parameters can absorb.
- Combined with the gauge mode, the result is an indefinite, inconsistent linear system.

**Evidence:**

- **Original and `sign1` runs:**
  - The first solve, made before any down-weighting, is clean: 19 MINRES iterations, rnorm 649 → 12.
  - Every later solve stops after 2 iterations with `istop=1`, rnorm ≈ ‖b‖ ≈ 1.5E+04 and ynorm ≈ 1.5E+04.
  - The expected decrease `dfcn_exp` (≈ 0.27E+09) is larger than fcn itself, and fcn rises: 0.11339E+08 → 0.15682E+08 → 0.19452E+08.
- **Test 4.1 (down-weighting off):** the second solve is consistent (‖b‖ = 11.8, ynorm 16).
- **Test 4.2 (ITS envelope fixed, down-weighting on):**
  - fcn is flat: 0.11340E+08 → 0.11344E+08 → 0.11340E+08.
  - In the later solves rnorm is 5e-3 to 8e-4 and ynorm is 12 to 4, with MINRES `istop=5` (iteration limit). This is the same benign behaviour as `testNVSingle`.

Fixing the gauge cures the divergence. The `IEEE_DENORMAL` floating-point warning from pede is a symptom of this problem, not a separate issue.

### F3. The remaining abort is pede's chi² cut, not a fitting problem

- After alignment, multi-track records have chi²/Ndf = **3.14 ± 0.99** at Ndf ≈ 240, the same with or without down-weighting. The 3σ cut at that Ndf is chi²/Ndf ≈ 1.3.
- Single-track records go from 4.1 at iteration 0 to 1.33 after alignment, against a 3σ cut of ≈ 2.4 at Ndf ≈ 14.
- The pede cut schedule (`pede.f90` near lines 3545–3557) cannot avoid the abort:
  - the factor is reduced by `sqrt` at each iteration;
  - it is set to 1.0 once it falls below 1.5;
  - it is also set to 1.0 after any iteration with 0 rejects;
  - `chirem` is ≥ 1.

  The cut therefore always ends at the bare 3σ quantile, and tuning `chisqcut` cannot help.
- The excess is about 430 chi² per vertex, i.e. ≈ 23 per track for its 2 pointing degrees of freedom. These are pulls of ≈ 3.4σ to the common vertex.
- All DOFs are free, so residual misalignment cannot explain this. **The pointing errors at the vertex are underestimated by about 2–3×.**

### F4. The likely cause of F3 is multiple scattering underestimated at B=0

- The scattering angle goes as 1/p, so its contribution to chi² goes as 1/p². The error scale is set by RMS(1/p) = sqrt(⟨1/p²⟩), not by 1/⟨pT⟩.
- Example estimate:
  - spectrum dN/dpT ∝ pT·exp(−pT/T) with T = 0.3 GeV (⟨pT⟩ = 0.6 GeV);
  - low cutoff ≈ 0.05 GeV, since at B=0 there is no curvature cut;
  - then ⟨1/pT²⟩ = E1(pmin/T)/T² ≈ 1.37/0.09 ≈ 15 GeV⁻², giving an effective momentum of **≈ 0.26 GeV**.
- Assuming pT = 0.6 GeV therefore underestimates the scattering sigma by ≈ 2.3× on average (≈ 5× in chi²), and by 3–5× or more for soft tracks. Using p = pT·cosh η instead of pT reduces this only slightly within the ITS acceptance.
- The vertex pointing error is dominated by scattering in the beam pipe and L0. The common vertex forces about 19 tracks through one point, so multi-track records feel underestimated scattering much more than single-track records. In single-track records the vertex is a single weak measurement, and the per-layer kinks absorb part of the excess.
- The wide spread of true momenta also explains the broad, flat multi-track chi²/Ndf distribution, from 1 to 10 with heavy tails.
- In the current tree `AlignParams.meanPtB0 = 0.3` ([Params.h](../include/O2Align/Params.h); applied in `process()` of [AlignmentSpec.cxx](../src/AlignmentSpec.cxx)). This goes in the right direction but has not been tested yet.

### F5. Per-track reference states did not pass through the common vertex (implemented, uncommitted)

**The problem.** The Kalman-updated reference state of each track does not pass exactly through the common vertex. The mismatch hides part of the tension between tracks.

**The change in the working tree:**

- `Track::updateWithVertex(vtx, covScale)` ([AlignmentTypes.cxx](../src/AlignmentTypes.cxx)) updates the track with the vertex covariance scaled by `covScale`. The vertex point stored in `info[0]` keeps the true covariance.
- `constrainWithVertex(vtx, ivref, useCommonVertex, resTracks)` ([AlignmentSpec.cxx](../src/AlignmentSpec.cxx)) applies the new parameter `AlignParams.vtxMultiTrackRefCovScale` (default 1e-2) only when a common vertex is used. Otherwise the scale is 1.

**Expected effect.** Some hidden freedom is removed, so the multi-track chi² should rise a little.

**Side effect.** `kfFit.chi2` then contains the vertex chi² computed with the scaled covariance. It is used for diagnostics only.

## Recommendations

1. **Fix the gauge** whenever the multi-track vertex constraint is used with a free mean vertex. Either fix the ITS envelope translations and rotations, or fix the mean vertex. Put this in the DOF config. Ideally also add a warning in the code when both are free and `useMultiTrackPVConstraint > 0`.
2. **Next test:**
   - Regenerate the records with the ITS envelope fixed, `meanPtB0 = 0.3` (try 0.25 as well) and the new `vtxMultiTrackRefCovScale`.
   - Check that the multi-track chi²/Ndf ends near 1 and that pede converges with down-weighting on.
   - A conservative (low) momentum is preferable. High-momentum tracks then come out with chi²/Ndf < 1, which pede tolerates; too high a momentum gets the soft tracks rejected.
3. **If chi²/Ndf is still well above 1:**
   - Add a separate scale factor for the B=0 scattering in `Params`.
   - Revisit the ITS cluster extra errors (`extraClsErrY/ZITS`, currently 10 µm).
   - Compare the vertex pulls per track with the per-layer residuals, to separate scattering from intrinsic resolution.
4. **Housekeeping:**
   - Uncommitted changes are in `AlignmentTypes.h/.cxx`, `Params.h`, `AlignmentSpec.cxx` and `DetectorITS.cxx`. Stray files `gdb.txt` and `macro/applyAlignmentRed.C` are also in the tree. Commit selectively, not with `git add -A`.
   - `chi2Ndf` and `ndf` are still not filled in `Track::fitTrack` (an existing `RSTODO`).
