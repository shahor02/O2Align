# Half-barrel-level change between pp (May 2025) and PbPb (Nov 2025) B+ data

Date: 2026-11-21

## Inputs
checkResid (PV1) of two B+ runs, each with two alignments:

| period | run | best/IBEdge | V10n |
|---|---|---|---|
| early (pp) | LHC25ab 562837 | `/data0/devAlg/tests/best/IBEdge/checkResid_LHC25ab_562837_hman_AlgRS_PV1.root` (tree in `reducedFit/IBEdge/`) | `/data0/devAlg/tests/v10n/checkResid_LHC25ab_562837_hman_AlgV10n_PV1.root` |
| late (PbPb) | LHC25as 568664 | `/data1/PbPb/newAlg/IBEdgeRS4/checkResid_LHC25as_568664_hman_IBEdgeRS4_PV1.root` | `/data1/PbPb/newAlg/v10n/checkResid_LHC25as_568664_hman_v10n_PV1.root` |

Control: LHC25ah 564644 with best/IBEdge (`best/IBEdge/checkResid_LHC25ah_564644_hman_AlgRS_PV1.root`).

PDFs: `/data0/devAlg/tests/best/IBEdge/checkResid_PV1_hman.pdf` (early), `/data1/PbPb/newAlg/checkResid_PV1_hman.pdf` (late).

Observed in the late PDF:
1. dY vs φ gets worse at 1.8<φ<3 and 5<φ<6, correlated between Vtx, Lr0 and Lr2.
2. Lr1 dY vs φ has a slope at 0<φ<π and a wave at φ>π.
3. Lr4 dY vs φ has a slope at 0<φ<π.
4. The IB-OB sub-track differences vs φ (p.9) get worse, with structure at the half-barrel level.

## Method
- **Means.** Every slice of every TH2 is re-fitted with a Gaussian (`FitSlicesY`). The stored `_1` histograms cannot be
  used because they are hadd sums. The quantity analysed is the late − early difference of the fitted means, per
  alignment.
- **Per-chip histograms.** `dy_LrN_{0,1}` / `dz_LrN_{0,1}` are per HB, with x = chip ID. Chip φ and z are the mean
  `points.phi` / `points.z` per `points.sens`, taken from the res tree.
- **Chip ranges.** Chips per stave are 9 (IB), 112 (Lr3/4) and 196 (Lr5/6). The first chip of each layer is
  0, 108, 252, 432, 3120, 6480 and 14712.
- **Rigid fit per layer and HB:** dy = −tx·sinφ + ty·cosφ + r·ωz − z(ωx cosφ + ωy sinφ), dz = tz + r(ωx sinφ − ωy cosφ).
  The values are *apparent* motions, i.e. leave-one-out responses, not true displacements.
- **Vertex grouped by IB stave.** From the trees: the vertex dy of 7-cluster tracks with pT > 0.7 and |dy| < 60 µm is
  grouped by the Lr0/Lr1/Lr2 sensor that the same track crosses (`points[0]` = vertex, `points[1..3]` = Lr0..Lr2). In
  the hman histograms the vertex φ is the track direction, so this grouping is the only way to bin the vertex residual
  by IB stave.

## Results

### The change does not depend on the alignment
- In every bin, (late − early)<sub>IBEdge</sub> − (late − early)<sub>V10n</sub> is below 0.5 µm. This holds for dy vs φ
  of Vtx and Lr0–6, and for the IB-OB parameters.
- So the change is in the data or the detector, and the two alignments differ from each other in the same way in both
  periods.

### It looks geometric and comes after 564644
- **Flat in pT.** The late − early dy is flat in pT, e.g. Lr1 +0.7 µm and Lr4 −3.7 µm in every bin up to 9 GeV/c.
  Material or multiple scattering would not give that.
- **Not the PbPb vertex.** The vertex-point error is nearly the same in pp and PbPb: median 24 vs 23 µm at pT > 1, and
  12.2 vs 10.7 µm at pT > 4. A more precise PbPb vertex therefore cannot be what reveals an old misalignment.
- **Step after 564644.** 564644 − 562837 agrees within about 1 µm in the IB and within 8 µm in IB-OB par0, against
  ±50 µm for 568664. The only exception is Lr4 HB0, with ty = +4.7 µm. The change therefore happened after 564644.

### Per layer and HB (late − early, IBEdge; V10n agrees)
| | HB0 (0<φ<π) | HB1 (π<φ<2π) |
|---|---|---|
| Vtx (φ profile) | ty −3.3 | tx −2.7, ty +3.6 |
| Lr0 | tx −1.4, ty +1.4 | tx +2.9 |
| Lr1 | **ty +2.2** (rms 1.76→0.35) | tx −1.7, **ty −1.5** |
| Lr2 | tx −2.0, ty +1.0 | tx +4.3 |
| Lr3 | tx +1.4, ty +2.6 | tx +6.3, ty −2.0 |
| **Lr4** | **ty −12.9** (rms 13.0→8.8) | tx +3.7, ty +5.3 |
| Lr5 | tx +3.9 | tx +7.7, ty +4.8 |
| Lr6 | ≲2 | ≲3.5 |

Units are µm. The fitted tilts ωx and ωy are at most 15 µrad (under 0.6 µm at the stave ends). The dz changes are a few
µm: Lr3/Lr5 HB0 have tz ≈ +4 against Lr4. The vertex dz vs tgl gained a slope of about +5 µm per unit tgl, which is
radial-like and not yet understood.

### Lr4 HB0: one rigid shift along Y (item 3)
- Per stave, the change runs linearly from −17 µm at φ = 0.1 to +11 µm at φ = 3.0. Lr3 HB0 shows the opposite-sign
  leakage (+2.6).
- Lr3 has an extra −8…−10 µm on its last staves (φ 5.9–6.2). Lr5 HB1 has −6…−13 µm at φ 3.2–4.3.
- Lr5 and Lr6 also have a stave-to-stave scatter of ±5 µm.

### IB: the two halves move relative to the vertex, Lr1 against Lr0/Lr2 (items 1, 2)
**Expected leave-one-out response.** If an IB half moves rigidly by d relative to the vertex and the OB:
- the vertex sees about −d;
- Lr0 sees a fraction of +d (it lies between the vertex and Lr1);
- Lr2 sees a fraction of +d (it lies between Lr1 and Lr3);
- Lr1 sees about 0.

At 1.8<φ<3 the change is Vtx positive, Lr0 and Lr2 negative, and Lr1 small (−0.6…−1.6 µm), which is this pattern. The
vertex ty is −3.3 in HB0 and +3.6 in HB1: a vertical displacement of the top half relative to the bottom one of about
7 µm, seen through the vertex that both halves share. Lr1 moving against Lr0/Lr2 (ty +2.2 in HB0, −1.5 in HB1) adds a
per-HB IB sagitta of 1–2 µm.

**Vertex dy grouped by the Lr2 stave crossed** (late − early, µm, stat err ≈ 0.12):

| | per Lr2 stave, in increasing φ | rms → after rigid fit (tx, ty, const) |
|---|---|---|
| HB0 | −2.4 −1.4 −0.8 −0.9 +0.3 −0.1 +1.1 +0.7 +2.0 +4.1 | 1.76 → 0.68 |
| HB1 | −1.2 −2.4 −0.5 −1.7 −1.6 −0.7 +2.1 +3.7 +4.1 +3.0 | 2.35 → 0.92 |

A rigid per-HB motion explains about 85 % of the variance. The ±1 µm left over is stave-level:
- **HB0 edge stave at φ ≈ π.** It goes beyond the rigid fit in every layer: Vtx +4.1, Lr0 −2.7, Lr1 −1.4, Lr2 −4.3.
  This may be a non-rigid opening of the half-shell at its edge.
- **HB1 at 5.1–6.1.** The change there looks more like a step than a cosine.

Per-layer stave rms after the per-HB fit:

| | HB0 | HB1 |
|---|---|---|
| Lr0 | 0.3 µm | 0.7 µm |
| Lr1 | 0.3 µm | 0.7 µm |
| Lr2 | 0.8 µm | 1.6 µm |

The largest stave-level structure is Lr2 HB1 at 3.3–4.9 (−3.0, −0.6, +1.6, +2.2 µm), outside the windows of item 1.

### IB-OB (item 4) is the IB sagitta magnified
- **It is a curvature change.** dyIOB of the OB layers (the IB+vertex sub-track extrapolated to Lr3…Lr6) has ty =
  ±140, ±209, ±398, ±512 µm with the sign of the HB. Divided by r² this is a constant 0.35 µm/cm², i.e. a pure curvature
  change of the IB sub-track.
- **Size.** It corresponds to Δ(q/pT) ≈ ±0.05 GeV⁻¹ (IB-OB par4), which gives par0 ≈ ±55 µm at 12 cm (ty −53 in HB0,
  +44 in HB1).
- **Cause.** A sagitta of about 1.3 µm over the 3.8 cm of the IB is enough. No large IB-vs-OB displacement is needed:
  item 4 is the same IB change as items 1 and 2.
- **φ correlation.** The φ pattern of par4 follows the Lr1 dy change bin by bin.

### Structural limitation of the current configs
In the v21-type configs the OB HBs are fixed, with `pinChildrenMeanNoCross` and `groupModes globalShift_x/y` over their
staves. A stave-level fit therefore cannot absorb an OB HB move such as the one of Lr4 HB0.

## Tests
1. **When did it happen.** Run checkResid with best/IBEdge on runs between 564644 and 568664: the last pp/OO before
   PbPb, the first and a late PbPb run, and a B− run if there is one. Track the per-HB amplitudes of the tables above
   run by run.
   - A step at a given date points to a mechanical or thermal event (field ramp or polarity flip, an access, a cooling
     or power change).
   - A step exactly at the switch to PbPb points to PbPb conditions.
2. **Rule out PbPb-specific effects.** Split 568664 by occupancy or vertex multiplicity. A geometric change does not
   depend on it.
3. **Forward test.** Add the fitted HB motions to the IBEdge geometry and run checkResid on 562837. It must reproduce
   the late − early pattern, in particular the r² growth of IB-OB.

## Realignment
**Stage 1: HB only** ([ITS_RB_HBonly.json](ITS_RB_HBonly.json)).
- The 14 per-layer HBs get TX, TY, TZ and RZ, with priors of 30 µm and 1 mrad.
- The mean vertex is free.
- `ITS` is fixed with `pinChildrenMeanNoCross` [TX, TY, TZ, RZ] over the 14 HBs, which removes the global gauge.
- That is about 60 parameters.

How to run it:
1. Record `mp2data` from 568664 (and neighbouring PbPb runs) with the geometry remapped to best/IBEdge, so the start
   values are 0. Use the usual mvbias and `useMultiTrackPVConstraint=2` settings.
2. Run pede (inversion is fine at this size), then `MilleRes`, then checkResid on 568664.
3. Expect the Lr4 HB0 slope, the Lr1 sawtooth, the top/bottom flip of the vertex, the HB-sign pattern of IB-OB
   par0/2/4 and most of item 1 to disappear.

**Stage 2: IB staves**, only if the ~1 µm stave-level remainder at φ ≈ π and 5–6 survives stage 1. Add:
```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "rigidBody": ["TX"] },
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "measurement": { "TX": [0.0, 0.0005] } }
```
- This assumes stave TX is the rφ direction, as it is for the chips; check that against the jacobian before running.
- With the HBs free, the automatic stave-mean constraints keep the split between HB and stave defined.
- Running stage 1 first tells whether the edge effect is a single stave or the half-shell deforming.
