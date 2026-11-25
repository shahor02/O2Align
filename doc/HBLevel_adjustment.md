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
- `ITS` is fixed with `pinChildrenMeanNoCross` [TX, TY, TZ, RZ] over the 14 HBs. Because HB1's frame is rotated by π
  about z, this removes the global TZ and RZ gauge but not TX and TY (see
  [DOFConfig_rules.md](DOFConfig_rules.md#pinchildrenmeannocross-the-same-pin-without-the-jacobian-cross-terms)).
  [ITS_RB_HBonly_v3c.json](ITS_RB_HBonly_v3c.json) uses `pinChildrenMean`, which pins the common lab motion.
  [ITS_RB_HBonly_v3d.json](ITS_RB_HBonly_v3d.json) adds `pinChildrenMeanNoCross` [TX, TY, RX, RY], which also pins the
  layer-averaged top−bottom relative motion, not constrained by collision data.
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

## Result of stage 1 (`solPbHBA0`) and the weak-mode constraints
Fit: `/data0/devAlg/tests/solPbHBA0`, [ITS_RB_HBonly.json](ITS_RB_HBonly.json) on the PbPb mp2data of
`Mille_Pb_MP2/IBEdgeS4R` (568664, `useMultiTrackPVConstraint=2`, `maxTracksPerMultiTrackPV=6`), checkResid on 568664.

### The solution is one radial family
- The HB TX/TY are a few µm. The RZ carry the solution: r·RZ is about +10 µm on Lr0–Lr2 (both HBs), then −29, −51,
  −113 and −153 µm on Lr3–Lr6.
- Layer-averaged, rφ(r) = d + a·r + c·r² with residuals below 0.5 µm on all 7 layers:
  - d = +10.9 µm, a constant tangential offset: every track's d0 changes by d, which only the vertex sees;
  - a = 12.8 µrad, a global rotation (gauge);
  - c = −0.109 µm/cm², a sagitta: Δ(q/pT) ≈ −0.015 GeV⁻¹ for both charges. In the units of the 2026-11-22 history
    entry this is k = −10.9, against the B0-measured IBEdge sagitta of +1.38.
- `drawHBShift.C` (`3dview.png`) shows exactly this: the OB HB0 move by +16…+91 µm in x and the HB1 by the opposite,
  growing with r. It is a rotation pattern, not an IB vs OB translation.
- `weakModes.pdf` misattributes it. Its `curl` (17.8 % of the power) duplicates the rigid rotation. With no r² term in
  the basis, the sagitta appears as `tangentialPolar_r0t0` = +94 µm and `tangentialPolar_r0t2` = −139 µm. The `sagitta`
  mode (rΔφ = ε·r²) has been added to `Volume::GlobalMode` for this. Rerunning `WeakModesDiag` on the solution with it
  gives `sagitta` = −174 µm at 39.7 cm (c = −0.110 µm/cm²), `tangentialPolar_r0t0` = +10.5 µm and a curl of 3.8 µm. The
  rms after the modes drops from 8.1 to 1.5 µm.

### Effect on the residuals (new − IBEdge, both on 568664)
| Vtx | Lr0 | Lr1 | Lr2 | Lr3 | Lr4 | Lr5 | Lr6 | IB-OB par0 | IB-OB par4 |
|---|---|---|---|---|---|---|---|---|---|
| −9.0 | +0.4 | −0.7 | +2.3 | −8.5 | +11.6 | +6.0 | +2.2 | −52…−64 | +0.11 GeV⁻¹ |

(µm; constant term of the dy vs φ fit.)
- The refitted PV cannot follow a tangential offset common to all tracks. The vertex residual moves by −d, flat in pT
  and tgl, and pulls the IB sub-track. Hence the IB-OB jump and the layer-dependent leakage: Lr3 is −11 µm at
  pT < 1 and about 0 at pT > 2.
- The small structures did improve. Against the early data, the dy-vs-φ rms is 0.9/0.7/1.5 µm on Lr0/1/2 (before the
  refit 1.5/1.7/2.2), and the Lr4 HB0 cos term is −0.4 µm (before −3.8).

### Why pede took it: a change of dy with incidence
- Mille replay of the solution on the same data, with the profiled Hessian of the 14 HB RZ projected on d, a and c:
  - Removing d, a and c costs Δχ² = +27.8k on the records without a vertex (σ_d = 0.10 µm, σ_c = 0.0015 µm/cm²).
  - The vertex records gain 1.1k from it and alone would want d ≈ 4 µm. They are 5.4k of 1.9M records.
  - So the family is not a null direction of these data.
- The reason is that a layer RZ also rotates every stave plane, which changes dy with the crossing point on the stave,
  i.e. with incidence. The res trees show such a change, fitted as dy = c0 + k·tanα with tanα from `points.snp`:

  | | Lr0 | Lr1 | Lr2 | Lr3 | Lr4 | Lr5 | Lr6 |
  |---|---|---|---|---|---|---|---|
  | late − early k (µm) | −2.9 | −4.1 | −4.3 | −4…−6 | +5±2 | −4±4 | +20…+30 |
  | solPbHBA0 − early k (µm) | −1.4 | −1.6 | −0.4 | +0.7 | | | |

  - It is the same for both charges and for cluster sizes 1–2, 3–4 and 5+ (cluster-size fractions are unchanged too),
    so it looks geometric rather than a cluster-position effect.
  - In the IB it is what a radial (normal) offset of about 3–4 µm of the stave measurement planes would produce. It is
    not linear: Lr1 rises from −0.4 to +1.6 µm at snp ≈ 0.11, then falls back to 0 at 0.39.
  - The OB tanα range is only ±0.13.
- HB DOFs can absorb such a change only through RZ, and any rφ translation outside the helix-null family is expensive.
  The solution therefore lands exactly in d + a·r + c·r².

### Constraints
Two variants, both checked by regenerating `mp2con.txt` (`--output MilleSteer`) on the IBEdge geometry.

**[ITS_RB_HBonly_v2.json](ITS_RB_HBonly_v2.json): the rotation of every layer is pinned.** Run it with
`AlignParams.groupModesRelToRoots=false`; otherwise the curl, being a free rigid motion of the group volumes, is skipped.
- The ITS has no layer volume, so the pin is a named `groupModes` curl group per layer (`Lr0` … `Lr6`) over its two
  HBs. Each equation is the layer's rotation about the beam axis: the HB RZ plus a small TX term, since the HB origin
  is off the axis.
- This removes d, a, c and the other 4 layer rotations, keeping only the HB0−HB1 RZ difference per layer.
  - The late − early layer rotations are only 2–4 µm on Lr3–Lr5.
  - With B+ data only, the sagitta must stay at its IBEdge value, since only B0 data constrain it.
  - RZ leaves the `ITS` pin, which the 7 layer equations already contain.
- `globalModes: ["radialHarmonic1", "tangentialHarmonic1"]` on the HBs removes the transverse shift ∝ r. It is near-null
  for tracks from the vertex, and was 3–4.5 µm at 40 cm in solPbHBA0 (IB TX about +1.6, OB about −1 µm).
- The RZ prior is tightened from 1 to 0.1 mrad. The old prior allowed 150 µm at Lr6.

**[ITS_RB_HBonly_v2b.json](ITS_RB_HBonly_v2b.json): only the helix family is constrained.** It keeps the `ITS` RZ pin
(a) and adds `globalModes: ["tangentialPolar_r0t0", "sagitta"]` (d and c) to the harmonic-1 modes. The 4 other layer
rotations stay free. A rotation pattern outside the helix family is measured by the rφ translation it implies, so it
can absorb little of the incidence change.

With these constraints the IB incidence slope stays in the residuals. To absorb it, free the IB staves along their
normal (stave local Y for ITS volumes; check against the jacobian), leaving the mean over each layer free, since a
uniform radial offset is seen through incidence. Before that, check whether the slope follows individual staves or is
common to all of them.
