# IB end chips and chip bowing: dy-vs-Z and dz-vs-Z patterns

Written for: whoever continues the ITS alignment tuning (knows O2Align, Millepede-II and checkResid).
Date: 2026-10-25 (last updated). State: 2026-10-25. Supersedes the versions of 2026-10-22/23 ("rphi-only, non-rigid term", which was an
artefact, section 5) and of 2026-10-24 ("the end-chip priors block the fit", which was wrong, section 6).
Sections 9–10 (restart test from IBv0, IBv1–IBv3, v16, other checkResid features) added on 2026-10-25.

## 1. Summary

- **End chips (dy vs Z, Lr0 |Z| > 10 cm).** The pattern is a rigid displacement of the IB end chips:
  - a normal shift plus a tilt about the chip's rφ axis (RX);
  - strongest on Lr0 chip0, the z>0 end; on chip8 the shift varies from stave to stave without a coherent
    pattern.

  The IB Mille data are consistent with rigid chips.
- **Why earlier fits missed it.**
  - v8ns, v10 and v11 kept chip RX fixed.
  - Freeing it (v13) was not enough: pede's `globalModes` constraints, written once "over ALICE" with the
    IB layers included, forbid the corresponding IB deformation.
  - The chip priors and pede's outlier treatment were tested and are not the cause.
- **Fix, first step.** `doc/ITS_RB_IBv0.json` (= v13b with `globalModes` restricted to the OB) moves the end
  chips, and the dy-vs-Z residuals improve (`/data/devAlg/solA2IBv0`, section 7). The changes are confined
  to the end chips, apart from a common ~9 µm translation of the IB plus mean vertex relative to the OB.
- **Still constrained in IBv0.** A pede restart from the IBv0 solution with only the 96 end chips free gains
  Δχ² −14k (Lr0 c0: TY ≈ −25 µm, RX −30…−40 µm at the edge), and checkResid confirms the move (section 9).
  Three things blocked it:
  - the IB half-barrel `pinChildrenMeanNoCross: ["RX"]`, which forbids a cone (removing it in IBv1 gave
    −5.9k);
  - the σ 5 µm chip TY priors acting on gauge-dependent chip-level values (relaxing them in IBv2 gave −3.8k
    more);
  - the stave `pinChildrenMeanNoCross: ["RZ"]` (removing it in IBv3 gave −7.4k more).

  `doc/ITS_RB_IBv3.json` reaches the restart-test level on the end chips, with the beam slopes free and the
  meaningless mean-vertex prior removed. It introduces stave rolls of ±0.3–0.9 mrad that checkResid does not
  see. `doc/ITS_RB_v16.json` merges it with the v15 OB.
- **Chip bowing (dz vs Z "V" shapes).** The remaining V shape inside each chip comes from chips bent out of
  plane along z:
  - normal displacement ∝ w², sagitta about −1…−3 µm;
  - chip-to-chip spread 2–6 µm.

  No rigid DOF can absorb it; it needs a new per-chip deformation DOF (section 8).
- **The handout.** The pooled RY "endChipTilt" of the (since deleted) end-chip handout of 2026-10-22 rotates about the wrong
  axis: in the local frame RY is the in-plane rotation.

## 2. checkResid with cluster size and snp (v8ns trees)

Trees: `/data/devAlg/solA2v8ns/checkResid_LHC25ab_562862_AlgRS_PV1.root` (pp, B0) and
`checkResid_LHC25ae_AlgRS_PV1.root` (OO, B+), with `points.clSize` and `points.snp` added.

- **No cluster-size dependence.** On Lr0 chip0 the dy-vs-local-z profile is the same for every `clSize`
  from 1 to more than 7. The two datasets agree.
- **snp is the hit position.** Within a chip, snp and the rφ hit position across the chip (v) are
  degenerate for primary tracks. In B0, snp minus the position-based angle y/x is a constant per sensor
  (±0.08, sign set by the half barrel) with only 0.005 rms. That offset is the beam position, at
  Y ≈ −1.8 mm. So a dependence on snp is a dependence on where the track crosses the chip.
- **Shape on Lr0 chip0.** The residual behaves like dy ≈ A(u) + B(u)·v, where u is the local z pointing
  towards the stave end. B grows towards the outer end: dy runs from −2 to +25 µm across v at the outer
  end, but only from −10 to −2 µm at the inner end. Read as a normal displacement δX(u)·tanφ_inc, this
  gives δX ≈ 13 → 28 → 62 µm across chip0 (B+): a shift plus a tilt.
- **Fine binning (0.3 cm).** On the z>0 side, Lr0 shows a per-chip sawtooth: dy rises within each chip and
  drops at every chip boundary. Chip0 goes −6.5 → +13 µm, chip1 −2 → +2 and chip2 −1 → +1 µm. The z<0
  side has a weaker, mirrored version.

## 3. Hit displacement or track bias?

The Lr0 residual is an extrapolation of the Lr6→Lr1 fit, since the outward prediction holds only the
vertex, which gives no direction.

1. **pT and charge (B+).** On the chip0 outer end, dy0 is +11.6, +11.5, +10.1 and +11.6 µm for the four
   combinations of charge and pT (< 0.8 / > 0.8 GeV/c).
   - This excludes a bias from curvature or from the field model, and an offset that differs between Lr1
     and Lr2: the weight of such an offset at Lr0 changes about 2× between the pT bins.
   - It does not exclude a common shift of Lr1 and Lr2 relative to Lr0. That reaches Lr0 with a weight
     of about −1 at any momentum (−1.03 at 0.6 GeV/c, −1.07 at 3 GeV/c in the toy).
2. **The OB has almost no weight.** A toy of the checkResid leave-one-out Kalman combination (straight
   lines, or a helix with a loose curvature prior for B≠0, with multiple scattering and approximate
   material) gives the OB almost no weight in the IB residuals at these momenta.
3. **The vertex anchors the attribution.** The mean residuals of the same tracks (dyv, dy0, dy1, dy2) were
   solved for the hit offsets (δ0, δ1, δ2) with the toy response matrix.
   - Only offsets relative to the vertex can be determined.
   - With a fixed vertex, the Lr0 end-chip hits come out displaced by +15…+40 µm at the outer end of
     chip0, and δ1, δ2 are smaller.
   - Forcing δ0 = 0 needs a "vertex offset" that depends on where the track crosses the Lr0 chip, even
     after averaging over all 12 staves. No PV translation can do that.

The vertex "vortex" for |zv| > 11 cm tracks (a φ-independent dyv of −20…−97 µm) is the same effect seen
from the vertex. The Mille fit (section 4), which uses the vertex, attributes it to the end chips, and the
IBv0 result confirms that (section 7).

## 4. Mille fit (mean vertex free)

Data: `/data/devAlg/mille{B0,BPlus,BMinus}/mp2data_0.bin`, recorded with all six chip DOFs; the same files
as the v13…IBv0 pede runs.
- **Parameters:** the IB chips with all six DOFs (y and z rows tied), the IB staves and half barrels, and
  the mean vertex of every time slot. The mean vertex comes out at (+0.32, −1.82) mm, as in pede.
- **Other globals:** all other globals stay at their recorded reference, or at the v13b values where
  stated.
- **Locals:** the local parameters are profiled out per record.

Δχ² (B0 / B+ / B−). The gain expected by chance equals the number of parameters.

| over the rigid model (6 DOFs per chip) | parameters | B0 | B+ | B− |
|---|---|---|---|---|
| full y/z split | 5184 | 4.6k | 10k | 9.6k |
| y-only normal β (∝ tanφ_inc, y rows only) | 432 | 2.3k | 6.3k | 6.3k |
| out-of-plane bow ∝ w² (y and z rows) | 432 | 10.2k | 18.7k | 18.2k |
| y-only normal × w | 432 | 0.6k | 1.0k | 0.9k |

| over the v8ns set (TX, TY, TZ, RY per chip) | parameters | B0 | B+ | B− |
|---|---|---|---|---|
| RX, all chips | 432 | 31.3k | 68.2k | 70.0k |
| RX, end chips only | 96 | 3.6k | 7.7k | 7.4k |
| RZ, all chips | 432 | 19.3k | 36.1k | 36.6k |
| RX + RZ | 864 | 48.2k | 101k | 103k |

- **Rigid chips are enough,** apart from the bow (section 8).
- **The RZ pattern.** RZ shows a coherent twist along Lr0, from +1.5 mrad on chip0 to −1.3 mrad on chip8.
  Pede finds the same in every run.
- **Staves and half barrels are exactly redundant with the chips.** Adding all 324 of their parameters to
  the chip-only fit gains Δχ² = 0, as it should.

Lr0 chip0, stave mean (rms), all six chip DOFs free, B+:

| TX | TY | TZ | RX (normal shift at w = 1.5 cm) | RY (at 1.5 cm) | RZ |
|---|---|---|---|---|---|
| +11 (13) µm | −37 (24) µm | +1 (5) µm | −28 (15) µm | −6.5 (5) µm | +1.5 (1.1) mrad |

The pattern is the same in all three fields. Chip8 has a large stave-to-stave scatter instead of a
coherent shift (TY rms 34 µm, TZ rms 24 µm).

## 5. Artefact: the y-only normal term with the mean vertex fixed

A first version of the accumulator kept every non-IB global fixed, including the mean vertex. The vertex
records measure the mean vertex relative to a reference at (0, 0), so each vertex was pulled by about
1.8 mm (~40σ), purely in the transverse plane. That bends the tracks in rφ only. The fit could absorb it
only with a y-only chip term, and produced:
- β: Δχ² 80–150k;
- full y/z split: 88–101k;
- a per-stave first harmonic in the stave azimuth, ∝ R·cos α, twice as large in B0 as in B±.

All of this disappears with the mean vertex free (section 4). `macro/milleSplit.cxx` now frees the mean
vertex, and the rule is noted in `CLAUDE.md`.

## 6. Why v13 and v13b did not move the end chips

**v13** (`doc/ITS_RB_v13.json`, `/data/devAlg/solA2v13`) is v8/v8ns with these changes:

| level | change |
|---|---|
| IB chips | free RX and RZ (all 6 DOFs), priors 3 mrad; the v8 priors TX 10, TY 5, TZ 10 µm kept |
| IB staves | free RX (prior 0.5 mrad); `pinChildrenMeanNoCross: ["RZ"]` |
| IB half barrels | `pinChildrenMeanNoCross: ["RX"]` |
| removed | the stave `weightedMean` rule (v8ns used equal-weight constraints) |

- **Result.** dz vs Z of Lr1/Lr2 and the vertex dy at Z > 6 cm improved. The Lr0 dy vs Z did not: per-bin
  medians at the z>0 end went +13 → +10 µm, at the z<0 end +6 → +10 µm, and the sawtooth remained.
- **What pede did.** On Lr0 chips 1–7, pede agrees with the free Mille fit in every DOF to within a few µm.
  On chip0 it gives TY +0.7 µm and RX +12.8 µm, against −37 / −28 µm in the free fit.

**Hypotheses tested and rejected:**

| test | outcome |
|---|---|
| v13b: end-chip priors relaxed to TX 30, TY 50, TZ 30 µm | values unchanged within 0.5 µm; χ² −3k out of 29.2M |
| `solA2v13bNOCHI`: no record cut, no down-weighting (0 rejected records) | values unchanged |
| my fit with all records kept (no χ²/ndf cut) | still wants the end-chip shift |
| my fit with the OB and all other globals at the v13b values | still wants the end-chip shift |

For the last row the χ² cost of pede's end chips stays at Δχ² 8–16k per field.

**Direct check of the objective.** A record-level evaluator applies a parameter set, refits the locals of
each Mille record and sums the χ²:
- **Same objective as pede.** At pede's v13b solution it gives 3.9751×10⁷ for the three files, against
  pede's reported 3.9753×10⁷.
- **Pede is not at the data minimum.** With my chip solution, the total is 54k lower (B0 −9.8k,
  B+ −24.6k, B− −20.0k).

So pede's v13b solution is not the χ² minimum of the data. Pede solves by matrix inversion, so the
difference must come from what it adds to the data: priors or constraints.

**The constraints are the cause.**
- **They reproduce pede.** Imposing pede's constraints from `mp2con.txt` exactly (via a null-space
  projection) reproduces pede's end chips: chip0 TY +3 µm, RX +11 µm.
- **Without the priors** they still hold chip0 at TY −3 µm, RX +10 µm, at a real data cost of 2.9k in B0.
  So the priors are not the cause.
- **Each hierarchical constraint type alone is a pure gauge fix.** This holds for Σ chip RZ = 0 per stave
  and for the transported stave TX/TY/TZ and RX/RY means. Each can be satisfied inside the data null space.
- **The combination is not.** Combined, the constraints have rank 328, of which only 324 lie in the data
  null space. So 4 equations restrict directions that the data measure.
- **Where the 4 come from.** The excess appears when the 11 `globalModes` equations "over ALICE" (twist,
  bowing, zScaling, radial/tangential/longitudinal polar) are combined with the stave, half-barrel and
  ITS-level pins. It is a collective dependence, not one bad equation: removing any single global mode,
  any of the ITS-level NoCross pins, or one of many IB half-barrel equations lowers the excess by 1.
- **Ray modes reach the IB too.** The `rayModes` equations touching the IB are "r-range 0" cells (18–127
  sensors each). They are part of the same constraint set.

Each `globalModes` equation forces one pattern to zero over the whole flagged set. With `ITS/.+` and the
per-layer IB rules, the IB layers sit in the same equations as the OB. An IB deformation that the vertex
and the IB tracks measure is then forbidden, unless the OB pays for it with the opposite one.

All other constraints kept and no priors, B0, cost relative to the free optimum:

| global modes dropped | χ² cost | Lr0 chip0 TY / RX |
|---|---|---|
| none (as pede) | 2886 | −3 / +10 µm |
| any single mode | 1.6k–2.8k | mostly little change |
| `radialPolar_r0t2` + `tangentialPolar_r0t2` | 875 | −39 / −6 µm |
| the four r0t1/r0t2 radial and tangential modes | 538 | −29 / −13 µm |
| all | 127 | −51 / −23 µm |

Two caveats on this test. My first KKT solve with `lstsq` silently dropped the constraints, because the
normal-equation entries (~10¹⁰) swamp the constraint coefficients (~1). Only the null-space projection is
reliable. And dropping every mode also frees directions that really are unconstrained: chip0 drifts beyond
the free-fit value. Hence the IBv0 choice: keep the modes for the OB only.

## 7. IBv0: global modes for the OB only (`doc/ITS_RB_IBv0.json`, `/data/devAlg/solA2IBv0`)

**Config.** IBv0 = v13b with two changes:
- `ITS/.+` replaced by `ITS/ITSULayer[3-6]/.+`, with the same 10-mode list;
- the three per-layer IB `globalModes` rules deleted.

**Pede.** It used the standard cuts and down-weighting. There is no rank defect, and the final χ² is
29.203M against 29.216M for v13b (−13k). The dy-vs-Z residuals of the IB end chips improve.

**Chip positions, IBv0 − v13b.** These come from `algdump.root` (`pos` is the full aligned geometry), after
removing one common IB translation:

| stave mean | normal shift at centre / outer end | rφ shift |
|---|---|---|
| Lr0 chip0 | +12.0 / +19.5 µm | −8.2 µm |
| Lr0 chip8 | +7.5 / +14.6 µm | −5.6 µm |
| Lr1 chip0, Lr2 chip0 | +4.7 / +4.5 µm | −6.1 / −4.3 µm |
| Lr1 chip8, Lr2 chip8 | below 1 µm | −4.5 / −4.3 µm |
| Lr0–Lr2 chips 1–7 | ≤ 1.6 µm (≤ 0.7 µm on chips 2–6) | ≤ 2 µm |

- **Common translation.** The whole IB moved by (−7.3, −4.9, +0.6) µm, the same in all layers and half
  barrels. The mean vertex moved by (−8.2, −5.5) µm in X, Y in every slot. The IB and the vertex thus moved
  together by about 9 µm relative to the OB, whose sensors moved 2 µm on average. Relaxing the IB global
  modes allows this. Watch it in later iterations.
- **Tilts, z scale, twist** per layer and half barrel: changes ≲ 1.3 µm at the stave ends.
- **Chip-level parameters** look very different: all Lr0 chip TY shifted by −12…−30 µm, the staves by
  +13 µm. This is bookkeeping between hierarchy levels: the transported stave constraint ties the mean chip
  TY to the chips' RX·z lever, and it cancels in the positions. Compare chip positions, not
  `millepede.res` values.

## 8. Chip bowing: the dz-vs-Z "V" shapes

In the IBv0 PDF, page 4, the IB dz vs Z has a 3 cm periodic structure (one chip) of 1–3 µm, with maxima
near the chip boundaries. It is the same in V10n and MP2 and in both datasets. The checks below use the
v8ns trees, because `solA2IBv0` has only the PDFs.

**It depends on the track angle.** On the middle chips (2–6), dz inside the chip:
- is flat for |tanλ| < 0.5;
- has a U or Λ shape for |tanλ| > 0.5 that roughly flips sign with tanλ. On Lr1, for example, one sign
  gives +1.4 … −0.5 … +2.5 µm across the chip, and the other the inverted shape, +0.7 … −2.0 µm.

That is dz = δX(w)·tanλ with a normal displacement δX ∝ w²: a chip bent out of plane along its length. An
in-plane z effect (pitch, segmentation, cluster position) would not depend on tanλ, and an incidence
effect alone would not repeat with the chip period.

**Fit per chip.** The fit uses dz = a + b·w + tanλ·(c0 + c1·w + c2·w²) and the same form in dy with tanφ.
The table gives the sagitta, the normal shift at the chip ends minus that at the centre (c2·1.5²), for
chips 1–7:

| | from dz, B0 | from dz, B+ | corr(dz-, dy-based) over chips |
|---|---|---|---|
| Lr0 | −2.3 ± 0.5 µm | −1.8 ± 0.5 µm | +0.45 / +0.64 |
| Lr1 | −2.6 ± 0.6 µm | −2.5 ± 0.6 µm | +0.85 / +0.81 |
| Lr2 | −1.3 ± 0.5 µm | −1.2 ± 0.4 µm | +0.44 / +0.50 |

- **Spread.** The chip-to-chip rms is 2–6 µm.
- **dy partner.** dy sees the same per-chip pattern, with a weaker lever (⟨tanφ⟩ ≈ 0.25).
- **Mille.** The Mille fit finds the bow as the only significant non-rigid shape: Δχ² 10k / 19k / 18k for
  432 parameters.
- **Sign.** Which way the chips bend has not been checked against the sign conventions.

**What it would take to fit.** No rigid DOF can absorb a w² shape. It needs a per-chip deformation
parameter: a normal displacement ∝ P₂(w) along the chip, entering both the y and z derivatives. The
deformation DOF sets (`LegendreDOFSet`) exist only for ITS3. Calibration-type derivatives are taken only
from the direct parent of the measurement leaf, so an IB chip bow needs a code change, not a config rule.

## 9. IBv0 still holds the end chips: restart test from the IBv0 solution

**Starting point.** In the IBv0 checkResid trees the IB dz shows single-chip outliers of −3…−5 µm, all on
end chips: Lr0 63 (st7 c0), 72 (st8 c0), 107 (st11 c8); Lr2 297, 351, 360, 378 (c0), 404, 413, 431 (c8).
On top of that most end chips share a pattern: c0 dz ≈ −1…−2 µm, c8 ≈ +1…+2 µm, dz/dtgλ −2…−4 µm (a normal
offset rather than TZ).

**Method.** A pede fit of the same `mp2data_0.bin` (B0, B−, B+) with every label fixed at its IBv0
`millepede.res` value (presigma −1), except the mean vertex and the 6 DOFs of selected chips, which start
from IBv0 and are free. There is no constraint file: with everything else fixed there is no gauge freedom
left. The steering repeats IBv0's (`chisqcut 20 10`, `outlierdownweighting 4`, `dwfractioncut 0.15`,
`method inversion 5 0.1`). This needs no new Mille data, because the data were recorded against the
reference geometry IBv0 solved for. Each run takes minutes.

Gradients of the data χ² vanish at an unconstrained optimum, and gauge constraints never cost χ². Any χ²
gained by freeing parameters on top of the fixed IBv0 solution is therefore paid by IBv0's constraints
acting on measured directions.

**Results** (files in `/data/devAlg/chipTest/<run>/`; shifts run − IBv0; rotations as the shift at the chip
edge, RX/RY at |w| = 1.5 cm, RZ at |v| = 0.75 cm):

| run | free | χ² | Δχ² vs `mvonly` |
|---|---|---|---|
| `mvonly` | mean vertex X, Y, Z (18) | 29 183 092 | — (reproduces IBv0: 27 067 rejects) |
| `null` | + c4 of the 10 staves of the targets | | −20 (60 params) |
| `targets` | + the 10 outlier chips above | | −850 (60) |
| `ends` | + all 96 IB c0/c8 chips | 29 169 098 | **−13 994** (576) |
| `mvslope` | `mvonly` + SlopeX/Y | | −264 (12) |
| `endsslope` | `ends` + SlopeX/Y | 29 168 833 | −14 259 (588) |

- **null:** all moves ≤ 0.5 µm.
- **targets:** the single chips move mostly in the normal direction (63: TY −6, RX −11 µm; 72: TY −12, RX
  −23 µm; 107: RZ +13 µm), TZ only ~0.3 µm. Freeing outliers one by one captures little.
- **ends:** a coherent move of all end chips, largest at c0 (the +z end), decreasing with layer radius:

  | | TY (normal) | RX at chip edge | TZ |
  |---|---|---|---|
  | Lr0 c0 | −8…−28 µm | −13…−44 µm | +1…+4 µm |
  | Lr1 c0 | −4…−14 µm | −3…−23 µm | +1…+2 µm |
  | Lr2 c0 | −1…−8 µm | −5…−19 µm | +1…+3 µm |
  | c8 | −5…+15 µm | mixed | |

  TZ of 63, 72, 107 comes out +3.6, +3.6, +4.8 µm, matching their checkResid dz. The outliers are the
  visible part of the collective move.
- **Beam slope.** The mean-vertex slopes (fixed at 0 in IBv0) are significant: SlopeX −0.03…−0.06 mrad in
  all slots, SlopeY +0.02…+0.03 mrad (slots 0–3) and −0.06…−0.07 mrad (slots 4–5), errors ~0.007 mrad. They
  are worth only −264 and are independent of the end chips: the `endsslope` gain is the sum of the two, and
  the chip shifts are unchanged to 0.1 µm. The beam tilt does not explain the end chips.

**checkResid validation** (`endsslope` geometry, pp B0 tree only; per chip, RMS over chips, IBv0 →
endsslope):

| | dz | dz/dtgλ | dy | dy/dtgλ |
|---|---|---|---|---|
| Lr0 c0 | 1.6 → 1.1 µm | 2.3 → 1.1 | 1.2 → 0.4 | 1.8 → 0.6 |
| Lr0 c8 | 1.7 → 1.0 | 1.9 → 1.4 | 1.1 → 0.7 | 1.5 → 0.6 |
| Lr1 c0 | 1.4 → 1.3 | 1.3 → 1.1 | 0.6 → 0.4 | 1.0 → 0.9 |
| Lr2 c0 | 1.5 → 1.1 | 1.5 → 1.3 | 0.6 → 0.5 | 1.5 → 0.8 |
| Lr2 c8 | 1.3 → 0.8 | 1.5 → 1.2 | 0.6 → 0.5 | 1.2 → 0.6 |
| chips 1–7 | unchanged | unchanged | unchanged | unchanged |

- Former outliers: 63 dz −3.3 → −1.7, 72 −3.7 → −1.3, 107 −4.0 → −0.9, 431 −2.4 → −0.3; Lr1 c0
  (153, 162, 171) barely change.
- dz vs Z, outermost bins: Lr0 −1.1/+1.6 → −0.2/−0.1, Lr1 −2.3 → −0.5, Lr2 −1.6 → −0.7 µm. Inner bins, the
  vertex residuals and the core widths are unchanged.

The move is real. The RX tilt compensates the large TY over the chip, so checkResid does not see the
−25 µm normal shift as a residual. The end chips were not absorbing a Mille-specific artefact.

**Suspected blocker: the IB half-barrel `pinChildrenMeanNoCross: ["RX"]`.** NoCross pins the plain mean of
the staves' own RX per half barrel, i.e. a common stave tilt in their own frames: a cone, normal shift ∝ z,
with opposite signs at the two stave ends. The other constraints are excluded or unlikely:
- `rayModes` has `rMin` 10 cm and no longer touches the IB;
- the priors are too weak (end chips σ 30–50 µm; the move costs ≲ 1 per chip);
- the stave RZ pin is a roll, unrelated to this shape.

**Test config `doc/ITS_RB_IBv1.json`** = IBv0 without the half-barrel RX pin and with the mean-vertex
SlopeX/SlopeY free. Removing the pin opens no gauge freedom: the half-barrel RX stays fixed, so no automatic
hierarchical constraint is written for it, and a whole half-barrel tilt is still measured against the OB.
The risk is a weakly measured IB-wide cone (`radialPolar_r0t1`); the stave RX prior (0.5 mrad, 67 µm at the
stave end) limits it only loosely. What to check:
- χ² against IBv0's 29 202 684: about −264 from the slopes, and at least ~14k more if the pin was the
  blocker;
- constraint rank = number of equations, no near-zero eigenvalues;
- Lr0 c0 as in `endsslope` (TY ≈ −25 µm, RX −30…−40 µm at the edge), small changes elsewhere;
- the IB `radialPolar_r0t1` amplitude in the `WeakModesDiag` output. If it runs away, add an IB-only
  `globalModes: ["radialPolar_r0t1"]` rule, separate from the OB one. If the χ² gain is much smaller than
  14k, the next suspects are the ITS-level NoCross over the IB and OB half barrels, then the stave RZ pin.

**IBv1 and IBv2 results** (`/data/devAlg/solA2IBv1`, `solA2IBv2`). Pede's total χ² contains the
`measurement` priors of `mp2con.txt` (7374 equations in IBv0/IBv1, exactly the ndf excess over the
restart runs). Only the data χ² = total − prior χ² (the prior part re-evaluated from `millepede.res`) is
comparable with the restart test:

| | total χ² | prior χ² | data χ² | Δ data vs IBv0 |
|---|---|---|---|---|
| IBv0 | 29 202 684 | 19 588 | 29 183 096 | — |
| IBv1 (no HB RX pin, slopes free) | 29 197 757 | 20 517 | 29 177 240 | −5 856 |
| IBv2 (IBv1 + relaxed chip priors, no MV prior) | 29 177 081 | 3 619 | 29 173 462 | −9 634 |
| `endsslope` (restart test) | | | 29 168 833 | −14 263 |

- **The prior χ² of IBv0 was dominated by two meaningless terms.**
  - The chip TY priors, σ 5 µm on chips c1–c7, cost 8.5k (9.4k in IBv1). They act on chip-level values,
    which are hierarchy bookkeeping relative to the stave. In IBv1 the stave-mean chip TY runs −17 → −32 µm
    along Lr0. The stave's transported constraint is gauge for the data, but not once the chip-level
    values carry priors: the end-chip move then pays prior χ².
  - The mean-vertex Y prior, 0 ± 50 µm against the beam at −1.82 mm, costs 1.33k per slot (8k in total).
    Its pull on Y is ~0.1 µm, so it changes nothing.
- **`doc/ITS_RB_IBv2.json`** = IBv1 with all IB chip priors set to the former end-chip values (TX 30, TY 50,
  TZ 30 µm, RX and RZ 3 mrad) and the mean-vertex `measurement` prior removed. Its constraint rank is full
  (2800 of 2800).
- **checkResid, per-chip means from the `_hman` 2D histograms**, RMS over chips, in µm:

  | | IBv0 | IBv1 | IBv2 | `endsslope` |
  |---|---|---|---|---|
  | Lr0 c0 dz, pp | 1.7 | 1.2 | 1.1 | 1.1 |
  | Lr0 c8 dz, pp | 1.6 | 1.3 | 1.2 | 1.2 |
  | Lr2 c0 dz, pp | 1.5 | 1.2 | 1.4 | 1.1 |
  | Lr2 c8 dz, OO | 1.4 | 1.3 | 1.1 | 0.8 |
  | Lr0 c8 dz, OO | 1.4 | 1.1 | 1.0 | 0.8 |
  | middle chips | unchanged | unchanged | unchanged | unchanged |

  On Lr0 and Lr1 IBv2 is close to `endsslope`. The Lr2 end chips lag: pp dz of 297, 360 and 413 is still
  −2.5 µm, against −1.3…−2.6 µm in `endsslope`.
- **The remaining 4.6k.** IBv2 recovers two thirds of the restart-test gain. The IB priors now cost only
  1.05k, and would cost 0.73k at the `endsslope` values, so they are not the remaining gap. `endsslope` need
  not satisfy IBv2's non-gauge pins (the ITS-level NoCross, the stave RZ NoCross), so part of the gap is
  their cost. Pede also stops at its 5-iteration limit with an expected ΔF of ~380 in the last iteration,
  so it is not fully converged; that is worth well under 1k.

**IBv3: stave RZ pin removed (`doc/ITS_RB_IBv3.json`, `/data/devAlg/solA2IBv3`).** IBv3 = IBv2 without
the IB stave `pinChildrenMeanNoCross: ["RZ"]`. The stave RZ stays fixed, so the common RZ of a stave's
chips, a roll of the stave about its long axis, becomes a fit parameter constrained only by the data and
the chip RZ priors. The constraint rank is full (2752 of 2752).

- **χ².** Total 29 170 192, prior 4 088, data **29 166 104**: −7.4k against IBv2 and −17.0k against IBv0.
  That is 2.7k below the restart-test bound, which is possible because IBv3 also moves DOFs the restart
  test kept fixed. The prior χ² rises mainly through the Lr0 HB0 stave TX prior (σ 5 µm), from 145 to 505.
- **checkResid per chip** (from the `_hman` 2D histograms) is now at the `endsslope` level, and Lr2
  improves where IBv2 stalled:

  | RMS over chips (µm) | IBv2 | IBv3 | `endsslope` |
  |---|---|---|---|
  | Lr2 c8 dz, pp | 1.2 | 1.0 | 0.9 |
  | Lr2 c0 dy, OO | 1.3 | 1.0 | 1.0 |
  | Lr2 c8 dy, OO | 1.1 | 0.8 | 0.7 |
  | Lr0 c8 dy, OO | 1.1 | 0.9 | 0.9 |

  Chips 351 and 413 are fixed (pp dz −2.4/−2.5 → −1.1/−1.3 µm). Chips 297, 360 and 378 keep −2.7…−2.9 µm
  in pp while they are normal in OO (−1.0…−1.5 µm). One geometry for both datasets cannot remove such a
  pp-only offset, so it is probably not alignment.
- **IBv3 introduces stave rolls.**
  - In `millepede.res` the stave-mean chip RZ goes from exactly 0 (IBv2, pinned) to 1.6 mrad rms on Lr0
    (HB0 −0.9…−2.3, HB1 +0.3…+1.9 mrad), 0.6 on Lr1 and 0.45 on Lr2. The errors are 0.08–0.12 mrad.
  - In the aligned geometry (`algdump.root`, IBv3 − IBv0; IBv1 has none) the real per-stave rolls are
    smaller: ±0.3–0.9 mrad, alternating in sign from stave to stave. The difference is hierarchy
    bookkeeping, as for the chip TY.
  - checkResid does not see them. dy vs the position across the stave changes by ≤ 0.4 µm per stave, even
    on the strongest rolls (Lr0 st2/st9/st11 +0.8, Lr1 st3 +0.95 mrad), in pp and OO, split by charge and
    pT as well. A 0.8 mrad roll alone should give an edge-to-centre change of ~1.5 µm on Lr0. The pooled
    profiles, the vertex jigsaw and the stave-edge excesses are unchanged.
  - So the rolls neither improve nor spoil checkResid. Mille fits them from data checkResid does not
    contain (the B− sample, tracks that are not vertex contributors, low pT), or they are weakly determined.
    Keep them under observation.
- **IB–OB translation.** In `algdump`, IBv3 − IBv0 moves the IB sensor centres by (0, −3.0, −4.5…−5.1) µm.
  The OB moves by y −0.8…+0.5 µm and z 0.0 µm, and the mean-vertex Y follows the IB (−3.7 µm). This adds to
  IBv0's ~9 µm. The checkResid IB–OB comparisons hardly see it: IBv3 − IBv0 dPar1 (Z at 12 cm) changes by
  −0.9 µm (pp) / −1.2 µm (OO), dPar0 has no first harmonic (≤ 1 µm), and `dyIOB`/`dzIOB` change by
  < 0.5 µm. The IB+PV sub-track carries the event vertex, which moves with the IB. The relative IB–OB(–mean
  vertex) position is thus close to a weak mode for these observables, which is consistent with its drift
  from config to config. It does not affect the tracking residuals, but it matters for the absolute geometry
  (e.g. TPC matching).

**Conclusion and next step.** IBv3 is the IB baseline. `doc/ITS_RB_v16.json` = v15 with the IBv3 IB block
and mean-vertex treatment, global modes for the OB only (`ITS/ITSULayer[3-6]/.+`); the OB rules are those
of v15 unchanged. Run it on the existing Mille data (all 45 files):
- the `ITS_RB_all.json` recording already holds every derivative v16 uses (IB chip RX/RZ, mean-vertex
  slopes);
- the linearisation is not a limit at these sizes (tens of µm, ≲ 1 mrad);
- a recording against a new reference would turn every prior into a prior on the increment.

Record once more at the end, against the final geometry and again with the wide-open `ITS_RB_all.json`,
as a convergence check (near-zero corrections expected). That run also tests whether the recording-time
selection made with the v10 geometry (GBL χ² cut, vertex refit) biased the tracks through the end chips.
For the v16 run:
- **Solver.** With ~32k free OB parameters pede needs MINRES, which hit its 2000-iteration limit in v11.
  Check `istop`, and that the IB end chips match IBv3 (chip positions, per-chip checkResid).
- **OB module priors** (30 µm) also act on gauge-dependent values below free staves and half-staves.
  Check the prior χ² per group; relax them as in IBv2 if they cost thousands.

**`MilleRes` bug found on the way (fixed in `f9ae6ee`).** `Volume::readMillepedeResults` kept only labels
with presigma ≥ 0, so parameters fixed at a non-zero preset value were dropped. For the restart test this
gave the reference geometry plus the freed chips instead of IBv0 plus the freed chips. All labels are kept
now; consumers treat missing labels as 0, so ordinary runs are unaffected. Two limits remain: corrections are
applied only to volumes with at least one free rigid-body DOF in the `MilleRes` config (use the config of the
fit that produced the values), and the calibration JSON writes only calibration DOFs free in that config.

## 10. Other IBv0 checkResid features (not end chips)

- **Vertex jigsaw in dY_vtx vs φ.** It has the Lr0 stave pitch: dY_vtx rises linearly across each Lr0 stave
  by ~2 µm (pp) to ~3 µm (OO), the same in both half barrels, while Lr0's own dy is flat across the stave.
  Lr2 dy shows the opposite ramp (+0.8 → −0.9 µm pp). A ramp in the position across the stave is a ramp in
  incidence angle, the signature of a normal offset of the measurement plane. It is unchanged by the
  `endsslope` geometry. Not understood.
- **Outlier bins at φ ≈ π (and weaker at 2π) in dy vs φ, Lr0 and Lr2.** They come from the outermost
  ~1.5 mm in rφ of the last HB0 stave of each layer at the half-barrel gap: edge minus centre +8…+13 µm on
  Lr0 st5, −2…−6 µm on Lr1 st7, +6…+15 µm on Lr2 st9, uniform along z, the same in pp and OO and independent
  of cluster size. The alternating sign is the leave-one-out signature of one layer. A real deformation of
  that size on the edge staves only is considered improbable; the cause is open. (An earlier attribution to
  tracks crossing between half barrels was wrong: those are only 1–10 % of the bin.)

## 11. Open points

- **v16.** Run `doc/ITS_RB_v16.json` (v15 OB + IBv3 IB) on the 45-file set; check the MINRES
  convergence, the IB end chips against IBv3 and the OB module prior χ² (section 9).
- **IB baseline.** IBv3 (section 9): no IB half-barrel RX pin, no stave RZ pin, relaxed chip priors, no
  mean-vertex prior, beam slopes free.
- **Stave rolls.** IBv3's rolls (±0.3–0.9 mrad in the geometry) are invisible in checkResid; watch whether
  they stay stable in v16 and in the convergence recording.
- **pp-only end-chip dz** on Lr2 chips 297, 360, 378 (−2.7…−2.9 µm in pp, normal in OO).
- **Priors on chip-level values.** Chip priors act on gauge-dependent values. Keep them loose wherever the
  parent has free DOFs with a transported constraint.
- **Beam slopes.** Free SlopeX/SlopeY in all further configs, and drop the mean-vertex `measurement` prior.
- **IB–OB translation.** Watch the IB + mean-vertex translation relative to the OB in `algdump` (~9 µm in xy
  at IBv0, plus (0, −3, −5) µm at IBv3). It is nearly invisible to the checkResid IB–OB comparisons. If it
  keeps drifting in the combined fit, constrain it, e.g. with a pin of the IB half barrels relative to the
  OB.
- **Convergence recording** with the final geometry and `ITS_RB_all.json`.
- **Chip bow DOF.** Implement it (section 8); the expected gain is Δχ² ≈ 10–19k per field.
- **Vertex jigsaw and the HB-gap edge strip** (section 10).

## 12. Tools

- **In the repository:** `macro/milleSplit.cxx`, `macro/milleSplitModels.py` and `macro/sensDump.C`,
  documented in `doc/milleSplit.md`. The mean vertex is free since 2026-10-24.
- **Session scratchpad only:**
  - the accumulators: `milleSplit3.cxx` (mean vertex free), `milleSplit4.cxx` (+ IB staves and half
    barrels), `milleSplit5/6.cxx` (non-target globals fixed at a pede solution, optional χ² cut);
  - `milleEval.cxx`: the record-level χ² evaluator for a given parameter set;
  - the fits: `fit6.py` – `fit16.py`, including `fit12.py` (pede's constraints imposed by null-space
    projection, with or without priors, selectable by constraint type or global mode) and `fit15/16.py`
    (constraint ranks inside the data null space);
  - the checkResid scripts: `flatIB2.C`, `d1.py`/`d2.py` (chip bow), `toy.py`/`toyB.py` (Kalman response
    toys), `c1.py` – `c4.py` (vertex-anchored inversions).

  Copy them into `macro/` if they are needed again.
- **Restart test (`/data/devAlg/chipTest/`):** `mkChipTest.py` writes, per run directory, `param.txt` (all
  labels at their IBv0 values, the chosen chips and the mean vertex free), `steer.txt` and `chips.txt`;
  `cmpChipTest.py <runs>` prints χ², Δχ² vs `mvonly` and the per-chip DOF changes. Session scratchpad:
  `o1.py` – `o13.py` (per-chip checkResid tables, tgλ fits, vertex profiles, edge checks).
