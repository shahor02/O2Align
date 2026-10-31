# Handover: residual patterns after the cluster-x fix (`solA2v5F`)

Written for: the person continuing the ITS alignment tuning (knows O2Align and checkResid).
Date: 2026-10-21. State: 2026-10-21. Everything below is read from `/data0/devAlg/tests/solA2v5F` (config `doc/ITS_RB_v5.json`),
compared with `solA2v5` (before the fix) and with V10n (`/data0/devAlg/tests/v10n`). Nothing was re-run for this note
except reading existing outputs.

## 1. What was fixed, and why it explained the IB per-chip pattern

`DetectorITS::prepareData` stored the **nominal sensor reference x** (`getSensorRefX`) as the x of every cluster. The
fit propagated the track to that x and used it as the derivative pivot. checkResid and `refitITSTrack` use the cluster's own
tracking x (`cls.getX()`). The two coincide only if the chip plane is at constant x. In V10n the IB chip planes are
pitched about the local x axis: dx_trk/dz_loc = +2.3 / +0.6 / −1.3 mrad (mean) for Lr0/1/2, rms 1.7 / 0.8 / 1.1 mrad
(`macro/chipPitch.C` on `Pede_v10n/GLO/Config/GeometryAligned`). Over a 3 cm chip the plane moves by up to ~35 µm in x,
so the prediction at the wrong x leaves dy = tan(φ_inc)·θ·z_loc and dz = tan(Φ_z)·θ·z_loc. The fit absorbed that into a
common chip RY (+480 / +104 / −234 µrad) and the resulting sawtooth.

Evidence (kept as reusable tools): the Mille rows alone reproduce the fitted RY as a leave-one-out slope of dy vs z_loc,
and one global term `tx·z_loc` per layer replaces all nine per-chip RY (`macro/milleChipSlope.py`). The fix is one line:
`double x = trkXYZ.X();` in `prepareData` (ITS3 tiles keep their own branch).

Result of the fix (`solA2v5F`): mean chip RY per layer −37 / −58 / +26 µrad (was +482 / +104 / −234); IB dy residuals are
flat to <1 µm in every projection; the Lr0 sawtooth is gone.

## 2. How to read the numbers here (one trap)

`macro/checkResidMedians.C` dumps per-bin **medians** of the 2D residual histograms. They agree with the PDF. The
fitted-mean histograms `*_1` stored in the v5F `_hman_` roots are ~7× larger than the projection medians (e.g. Vtx dz at
tgλ=−1.2: `_1` = 31 µm, median 4.0 µm); the V10n roots do not show this. Cause not investigated, probably the post-process
of the run with `nBinsPhi/Z/Tgl = 100/100/60`. Use the PDF or the medians, not `_1`, from these files.
pp B0 (562862) and OO B+ (564356) are independent samples; they agree on every feature below, so these are geometry
or reconstruction features, not statistics.

## 3. Remaining patterns

### p.2 — dz(Vtx) vs φ and IB dz vs φ
Size: ≲ 3 µm in the vertex and IB, to be compared with ±10 µm (V10n) before. The OB layers no longer have the ±40 µm
per-stave spikes of V10n (Lr3, Lr4, Lr6 are now within ±3–5 µm; Lr5 still scatters by up to ±10–15 µm from stave to stave).

| feature | where | amplitude |
|---|---|---|
| Vtx dips | φ = 0.97, 1.92, 4.12 (also 5.06) | −3 µm (OO and pp equal) |
| Vtx ramp across the wrap | φ = 0.03 → 0.45: +2.4 → −1.1; φ = 5.9–6.25: −1.2 … −1.8 | 3–4 µm |
| IB ramp inside one stave | φ = 0.03 → 0.28: Lr0 −1.3 → +1.5, Lr2 −1.4 → +1.9, **Lr1 opposite sign** (+1.3 → −1.0) | ±1.5 µm |
| same alternating pattern | φ ≈ 4.4 (Lr0 −1.0, Lr1 +0.8, Lr2 −1.6), φ ≈ 5.6–5.8 | ±1–2 µm |

Reading:
- The vertex dips sit at the same φ where V10n had its OB per-stave z errors (φ ≈ 1, 2, 4, 5.7, see the checkResid skill,
  section 6). V10n shows −10 µm there, v5F −3 µm: the OB staves at those φ are corrected but not completely. The vertex
  residual is an extrapolation of the 7 layers, so an OB stave TZ error of 2–3 µm shows up there with a lever arm ~1.
  Not verified which staves: that needs the `res` tree (per-sensor, `points.sens`), which is not in the v5F directory.
- The IB ramp inside one stave (3 µm over Δφ=0.25, i.e. ~0.6 cm of x_loc) means dz vs x_loc = −θ with θ ≈ 5·10⁻⁴: an
  in-plane rotation of that stave of ~500 µrad, opposite in the middle layer. In v5 the **stave RY is not free** (only the
  chips have it), and the mean chip RY is not pinned: a stave-level RY remains unconstrained at the 100–500 µrad level.
  The alternating Lr0/Lr2 vs Lr1 sign is the leave-one-out leakage of a z-shift of the neighbours.

### p.4 — dz vs Z in the IB
- |Z| < 9 cm: rms 0.9 (Lr0), 0.5 (Lr1), 0.7 µm (Lr2); a 3 cm sinusoid of 0.8–0.95 µm in Lr0 only (chip pitch in z).
- |Z| > 9.5 cm (end chips c0, c8), identical in OO and pp:
  - Lr0: an oscillation with the 3 cm chip period, amplitude ±3–6 µm (e.g. +1.9, −5.6 at Z=10.0, 10.65; +1.3, −6.1 at 12.15, 13.35).
  - Lr1, Lr2: monotonic, odd in Z: dz ≈ −3 … −4.5 µm (Lr1) and −5 … −8 µm (Lr2) at Z = +13, and +2 … +3 µm at Z = −13.
  - V10n has the same shape (−10 µm for Lr2 at +13): the fit did not remove it.
- The end chips are also where the fitted chip RY stays large: c8 = −561 / −491 / −119 µrad (Lr0/1/2) and c1 = −402 /
  −231 / −125, the other chips are within ±230 (`result.json`, mean over staves), against statistical errors of 20–40 µrad.
- Hits at |Z| > 10 in Lr1/2 come from tracks with tgλ ≳ 1.5 and a vertex displaced by ≳ 5 cm, so incidence-angle
  effects are mixed in. The dz vs tgλ of the same layers is small (see below, −1 … −1.8 µm at tgλ = 1.9), smaller than the
  −5 … −8 µm at the Z edge, so the Z-edge effect is mostly positional.

### p.7 — dy vs tgλ
- IB (Vtx, Lr0–Lr2): < 1 µm everywhere, except Lr1 at tgλ = −1.9 / +1.9 (+1.9 µm).
- OB, **even in tgλ**, the same in V10n and v5F (so not fixed by the fit):

  | tgλ | −1.9 | −1.7 | −1.5 | −1.3 | 0 | +1.3 | +1.5 | +1.7 | +1.9 |
  |---|---|---|---|---|---|---|---|---|---|
  | Lr3 (OO) | −4.8 | −10.5 | −4.4 | +0.4 | −0.9 | +2.1 | −2.9 | −7.3 | −3.3 |
  | Lr4 | −6.0 | −6.1 | −4.3 | −1.7 | −0.3 | −1.0 | −3.5 | −6.2 | −6.0 |
  | Lr5 | −16.3 | −14.1 | −5.3 | −1.0 | +0.7 | +3.0 | −2.8 | −10.7 | −20.9 |
  | Lr6 | −4.3 | +7.9 | +3.1 | −1.5 | +1.1 | −1.3 | +3.7 | +6.6 | −6.8 |

  Inside |tgλ| < 1.3 the OB curves are flat to ~1–2 µm.
  The structure starts at |tgλ| ≈ 1.4, where Z = r·tgλ reaches the stave ends (Lr3: 27 cm, Lr5: 48 cm).

### p.8 — dz vs tgλ
- Vtx: odd linear slope of **−3.4 µm per unit tgλ** (OO) / −2.4 (pp): +4.1 at −1.3, −3.6 at +1.3, flattening beyond |tgλ| 1.3,
  then −6 … −12 µm at +1.7 … +1.9.
- IB: Lr1 −0.7, Lr2 −1.2 µm per unit; Lr0 ≈ 0 (V10n had −2.5 / +3 / −5 µm at tgλ = 2). So the IB is now radially consistent to ~1–3 µm.
- OB **odd and rising like |tgλ|³**: Lr3 −8 (−1.9) … +30 µm (+1.9); Lr5 −30 … +57 µm (OO) / +63 µm (pp); Lr4 +6 … +17;
  Lr6 −4 (both ends, OO), −11 … −16 µm (pp). Inside |tgλ| < 1.4 they are within ±5 µm with a wiggle of ~1.2 period in tgλ in
  Lr3/Lr5 (e.g. Lr5: +4.2, −4.9, +2.6, −3.9, +6.1 at −1.3 … +0.9, i.e. a period of Z ≈ 40 cm, i.e. about two OB modules of 21 cm).
  The same as V10n.

## 4. Which deformations could produce them

Rules used: dz = tgλ·ΔR (radial shift of a layer) ⇒ linear in tgλ; with ΔR(z) ∝ z² and z = r·tgλ ⇒ dz ∝ tgλ³;
dy = Δ(rφ) even in tgλ ⇒ rφ displacement even in z; layer-to-layer differences only, because the leave-one-out residual
sees the shape of a layer relative to the others and is blind to a coherent shift of all (checkResid skill, section 3).

1. **dz vs tgλ ∝ tgλ³ in Lr3/Lr5 (odd, p.8): radial bowing, Δr(layer) ∝ z² with a layer-dependent amplitude.** The weak-mode
   diagnostic of v5F still shows `bowing` = +83 µm (B) at the barrel end, `radialPolar_r1t2` = −78 µm, `zScaling` +17 µm,
   `radialPolar_r0t0` +15 µm (before the fix: 71, −67, 16, −): all large and none is in the `globalModes` list of v5
   (A ≈ 2 µm, unconstrained). Candidate: add `bowing` (not in the span of `rayModes`, Δr ∝ r²t²) and `zScaling`. Do not
   add the p = 1 polar terms while `rayModes` is on (nearly degenerate, `doc/WeakModes.md`).
2. **dy vs tgλ, even (p.7): an rφ displacement even in z with layer-dependent amplitude**, i.e. a `tangentialPolar` r·t²
   type (rΔφ = ε·r·t², p = 1) or a stave bending in the transverse plane (a "banana": Δ(rφ) ∝ z² of the whole stave or
   module chain). `tangentialPolar_r0t2` is constrained, `r1t2` is only constrained through `rayModes` (common to all layers),
   so the layer-dependent part is unconstrained. The same sign for Lr3, Lr4, Lr5 (−6 … −20 µm at |tgλ| = 1.9) with an opposite
   Lr6 (+8, +7 then −7) fits a layer-alternating bending (the leave-one-out leakage ratio ~ −1/2).
   The alternative is a non-geometric cause: cluster position bias vs incidence angle (pattern dictionary at large cluster
   size), which would be the same for V10n and the fit, as observed. The two cannot be separated from these plots.
3. **dz vs tgλ, wiggle of ~40 cm period in Lr3/Lr5 and the Lr5/6 shapes at |tgλ|<1.4**: module-scale effects: each OB module
   is a rigid body in the fit (TX, TY, TZ, RY with TZ prior 30 µm) but its length or chip pitch in z is not a DOF (no
   scale). A z-scale error per module gives a sawtooth in dz vs Z of the module length (21 cm) and, through z = r·tgλ, the period
   in tgλ of 21/r ≈ 0.6–0.9 (r = 24–34 cm); the observed 1.2 is about twice that, so a module-pair effect cannot be excluded.
4. **Vtx dz vs tgλ linear −3.4 µm/unit:** an IB-vs-OB radial offset of ~3 µm (all IB layers see −0.7 … −1.2). In the
   v5 `globalModes` the corresponding radial modes (`radialPolar_r0t0`, `r0t1`) are present, but `radialScaling` (−20 µm) and
   `radialPolar_r0t0` (+15 µm) remain. Expected to be small; low priority.
5. **IB end chips (p.4): IB end-of-stave effect**: c8/c0 TZ (z of the end chips) and RY. Hypotheses: (a) statistics at the
   acceptance edge (the formal errors at c8 are the largest, 40 µrad, but the value is 14σ); (b) the incidence angle of
   the tracks reaching the end chips (tgλ ≳ 1.5 and displaced vertices); (c) end-of-stave material not in the
   propagation (the IB end has the flex tail). Not decided.

## 5. Suggested next steps (cheapest first)

1. **Constraints:** `ITS_RB_v7.json` = v5 + `bowing` and `zScaling` in `globalModes` (check that the weak-mode A goes to 0 and
   that p.8 Lr3/Lr5 cubic rise moves, not just the diagnostic bars).
2. **Stave RY in the IB** (the v6 variant, `ITS_RB_v6.json`) now makes sense again: v6 was tested before the fix, where it could
   only redistribute the artefact. Rerun it on the fixed code and compare the ramp in p.2 (φ = 0.03 … 0.28).
3. **Chip RY need:** the original reason to free the IB chip RY was the V10n in-chip slope seen in checkResid. After the fix,
   rerun with chip RY fixed (v5 without RY on chips) and compare `dz vs φ/Z` and the IB `dy vs Z`: if nothing changes, RY on the
   chips is not needed (chip RY rms is still 280–350 µrad against 20–40 µrad errors, which may be real).
4. **Separate geometry from incidence-angle bias in p.7/p.8:** cut on |tgλ| < 1.3 and rerun the OB alignment, or split the
   residuals by (tgλ, Z) with the `res` tree: a bias that depends only on tgλ is the same at all Z, a geometry effect follows Z = r·tgλ + z_v.
5. **Identify the OB staves behind the vertex dips (p.2)** with the tree (`points.sens` at φ = 0.97, 1.92, 4.12) and look at their TZ.
6. Re-run the `closure_chipRY.sh` closure on the new build as a control (should still give slope −1).

## 6. Files

- `src/DetectorITS.cxx`: the fix (`x = trkXYZ.X()`).
- `macro/chipPitch.C`: pitch of the chip planes of any GeometryAligned snapshot.
- `macro/milleChipSlope.py`: chip RY and incidence terms from the Mille rows (`python3 -I macro/milleChipSlope.py mp2data_0.bin mp2tree.txt 40000`; pure Python, ~15 s for 40000 records; skips the large multi-track records).
- `macro/checkResidMedians.C`: robust per-bin medians of a checkResid `_hman` file.
- `macro/injectChipRY.C`, `macro/chipRYSummary.py`, `doc/closure_chipRY.sh`: the chip-RY closure test (derivative chain verified: slope −0.9 to −0.83).
- `doc/ITS_RB_v6.json`: v5 + free IB stave RY.
