# Handover: the rφ ∝ r² sagitta of the IBEdge ITS and the attempts to remove it

Written for: whoever continues the ITS alignment (a person or Claude). Assumes O2Align, Millepede-II, checkResid and the
`rayModes` constraints are familiar.
Date: 2026-11-23. State: IBEdge (`/data0/devAlg/tests/best/IBEdge`) is the reference. Corrections v0–v4 live in
`/data0/devAlg/tests/sag/`. The best candidate so far is v3R, but it leaves inner barrel (IB) jumps the user does not
accept. The open decision is the gauge (see §7).
Chronological details and all numbers are in `doc/history.md`, entry "2026-11-22: rφ ∝ r² sagitta…". This note is the
digest.

## 1. Goal and constraints
- **Goal:** an ITS that is as bias-free as possible at all φ and tgl, with an extrapolation to TRD/TOF that has the right
  curvature scale and the least bias.
  - The absolute TRD/TOF alignment does not matter: they will be realigned to the ITS.
- **Designated data:**
  - Collisions: `/data0/devAlg/tests/Mille_IBEdgeS4R/coll_LutPt08ITSSEED_ITS`, ITS points only, with B0, B0Bis, B+,
    B+Bis and B−. It was recorded with
    `--condition-remap file:///data0/devAlg/tests/best/IBEdge/=GLO/Config/GeometryAligned`, so every solution is
    relative to IBEdge.
  - Cosmics: `…/cosmHP4LUT_ITS`. The same cosmic data processed with v10n are in `/data0/devAlg/tests/Mille_V10N`.
  - The `Mille_TRDTOF_C1V21` tree is unreliable: ignore it.
  - B− is for validation only.
  - Do not use TRD/TOF points to fit the ITS sagitta.
- **Recording DOFs:** `ITSnochipTRDTOF_RB.json` has no chip DOFs. The lowest levels are the IB stave and the OB module.
  The IB stave TY/TZ/RX/RY/RZ labels exist, but they are fixed.

## 2. The effect
- **Model:** each B0 track is fitted with ITS points only as v = a + b·u + k·u·r: u along the track, v transverse, r the
  point radius. k = 1 (in units of 10⁻⁶ cm⁻¹) is a 1.7 µm sagitta over the ITS, or about 1.3·10⁻³ (GeV/c)⁻¹ in q/pT at
  0.5 T. Tool: `sag/tools/sagITS2.C`.
- **IBEdge B0 at 100 % (28.3M tracks, |tgl| < 1):**
  - mean k = 1.38 ± 0.02;
  - n=2 harmonic 1.9 at −11° (maxima near 170° and 350°), n=1 0.8;
  - maximum 7.3 at 350–360°;
  - k rises with |tgl| and has an odd part, larger at tgl > 0. Near 150–170° and 330–350° it reaches 13–19 at
    |t| ≈ 1.9.
  - Map: 36 φ × 16 t over |t| < 2, in `sag/v4/kmapX.txt`.
- **It is a smooth field, not a displaced stave.** In leave-one-layer-out residuals with k free, every stave and z bin is
  within 1–2 µm. So it is a coherent r² field over all layers.
- **B± hide it:** they absorb k into q/pT. B0 and cosmics see it.
- **Its effect on TRD/TOF:** it drives the B0 − B± gap of the TRD/TOF `pre` residuals in `solB1A4`. Correlation 0.92,
  about 0.9 mm per unit k.
- **Origin** (geometry difference decomposition, `sag/tools/geo/geoK.C`):
  - All of k was created in the step v10n → C1v21w3: Δk 1.06, n=2 1.41.
  - That step moved the OB rφ relative to the IB, with n=2 of 10–14 µm.
  - The same step changed the ray rotation b by 30–70 µrad.
  - v10n itself is predicted at k ≈ 0.4.

## 3. Gauge freedom
- **What tracks cannot see:** for each ray (φ, t = z/r), an rφ field a + b·r.
  - a is fixed by the vertex, up to a global translation.
  - b(φ, t) is an azimuthal remap.
  - B± see b only through a q/pT scale 1 + ∂b/∂φ.
- **What cosmics through the axis fix:** a(φ) + a(φ+π) and b(φ) − b(φ+π). An even-harmonic b, such as n=2, stays free.
  Only cosmics with large d0, or external detectors, are sensitive to b.
- **Consequence:** the correction δrφ = −k(φ,t)·r·(r − β) removes k for any β. β only chooses which layers carry it,
  i.e. it fixes b.
  - The "ob" gauge, β = ⟨r⟩_OB = 33.25 cm, keeps the mean OB ray rotation, so the B0 line extrapolation is unchanged.
  - The IB gets a large lever: r − β ≈ −30 cm.
- **Pede cannot hold k:** in the Mille model k is a quasi-weak mode, with Hp/Hu ≈ 10⁻⁴, as small as for b. The B0 MS
  kinks at `meanPtB0` = 0.45 GeV/c absorb it. Replays (`sag/tools/rep`) give the k each dataset requests:

  | dataset | k requested |
  |---|---|
  | B0 | 0.8 |
  | B+ | 2.35 |
  | B− | 0.0 |
  | cosmics | 0.3 |

  So k has to be imposed from the debug-tree map and pinned, using the `rayModes` `rdphi` cells with right-hand sides
  re-centred on the start values by `pinToStart.py`.

## 4. The top/bottom half-barrel step seen by cosmics
- **Fit:** each cosmic is fitted with line + k·u·r + a step between its upper and lower halves.
- **The step depends on d0:**
  - |d0| < 2 cm: about 0 (−0.3 ± 0.9 µm overall, at most +12 ± 5 µm at 140–160°);
  - |d0| = 5–20 cm: −27 to +42 µm, varying with φ (about −20 µm at 20–60°, +30…+40 µm at 140–160°);
  - the same in 2025 and 2026.
- **The collision side agrees with zero.** For a cosmic through the axis the step equals a(φ) + a(φ+π). The ITS-only
  collision line extrapolated to the fitted beam measures the same combination: (d0(φ) + d0(φ+π))/2. That is ≤ 1 µm in
  B0 and B+, and B− is ≤ 2 µm off. So the step is not a compromise between B0/B+/B−/cosmics in which B− is in tension.
- **A z-dependent field is ruled out.** A field c·z² would make the step grow with z0², and the data reject it.
- **It reacts to geometry that collisions cannot see.** The step at the same |d0| differs between IBEdge and v10n by
  15–26 µm in several φ bins. Large-d0 cosmics are the only ITS data sensitive to b.
- **Not explained yet.** It could be a b-type geometry (an even-b component in the half-barrel structure), or a cluster
  bias at large incidence. The Mille records carry no d0, so pede cannot cut on it.
- **For this reason the cosmics are commented out in the v3R steering.** Their large-d0 part could pull the free module
  TX toward the step. Also, there are only about 130k cosmic records against 40M collision records, and under v0 they
  had the largest χ² cost (+1.1 %).
- **Related:** a separate HB-level change exists in PbPb, at IB top/bottom ~7 µm (memory note `pbpb-late-hb-change`).

## 5. Corrections tried
- **Directory:** each version is in `/data0/devAlg/tests/sag/<v>`, with a README. Every `.res` is relative to the IBEdge
  recording.
- **How a geometry is made from a `.res`:** run MilleRes with the IBEdge geometry as the initial one. The user then makes
  the geometry, runs checkResid against v10n, v2 and IBEdge, and produces B0 debug trees.

| version | construction | result |
|---|---|---|
| v0 (`sag/v0`, `millepede_sagCorr.res`) | ob gauge; IB stave TX gets chip-averaged k; OB module TX | B0: k 1.36 → 0.07, harmonics ≤ 0.25. B± worse in checkResid: IB Lr2, Lr4 stave 12, OB \|t\| > 1. B+ ≡ B− (corr 0.93–1.0), so this is visible geometry. Causes: chip-t-averaged k on rigid IB staves; high-t noise × r(r − β). |
| v1 (`sag/v1`, `_comb.res`) | gauge β = 18.3 cm | corrections up to −115 µm; dropped |
| v0R (`sag/v0R`) | reduced fit on v0 | not to be used (inherits the v0 defects) |
| v2 (`sag/v2`, `mkSagCorr2.py 4`) | k smoothed with a degree-4 polynomial in t per φ bin, 1-2-1 in φ; IB k track-weighted over t | halves the v0 cost; OB dips remain at the +z module boundaries (Lr3 z ≈ 20, Lr5 z ≈ 32, Lr4 z > 25) |
| v3 (`sag/v3R/base_millepede_sagCorr3.res`, `mkSagCorr3.py`) | v2 plus OB module TX + RY fitted to the chip-level field (Δx_loc = TX + RY·z_loc) | within-module misfit Lr3 4.4 → 1.3 µm |
| v3R (`sag/v3R`) | reduced pede on v3 (details below the table) | IB rms TX change 0.3–0.5 µm |
| v4 (`sag/v4`, `mkSagCorr4.py`, `model.py`) | hybrid gauge (details below the table) | see §7 |

- **v3R details:**
  - Free: OB module TX + RY, IB stave TX, and the mean vertex.
  - 212 rdphi pins: OB 24 φ × 8 t, IB 24 φ × 1 t, `rSpanMin` 0.05.
  - Priors: TX 10 µm and RY 1 mrad, both centred on v3.
  - Fitted with B0, B0Bis, B+ and B+Bis. Cosmics and B− are commented out.
- **v4 details:** the n ≤ 1 φ harmonics of k keep β = 33.25 cm; n ≥ 2 use β = Σr²/Σr over the IB chips = 3.45 cm. So
  the IB carries only n ≤ 1, and the OB carries the n ≥ 2 rotation.

Replay noPV ⟨χ²/ndf⟩ (`o2-dev-mille-replay -res`):

| | IBEdge | v0 | v2 | v3 |
|---|---|---|---|---|
| B0 | 0.91954 | 0.92184 | 0.92066 | 0.92055 |
| B+ | 1.00792 | 1.01343 | 1.01174 | 1.01014 |
| B− | 1.00879 | 1.01548 | 1.01350 | 1.01179 |

- v3R and v4 have not been replayed.
- The cost is quadratic in the correction scale (λ-scan minimum at −0.07), so it comes from the visible part
  (granularity), not from the weak k itself.

## 6. checkResid of v3R: what is left
- **The OB module-boundary dips of v2 are gone.**
- **IB jumps, as in v2.** checkResid p.9 shows IB−OB par0 at 12 cm dipping by up to −50 µm at φ ≈ 2.1, 2.7, 5.3 and 5.85,
  in B+ and B− only (B0 is flat). par2 and the field-odd par4 (±0.06) dip at the same φ, and Vtx/Lr2 dy jump on p.1.
  - Per-chip change in Lr2 (pp.29/31): a constant offset per stave.

    | stave | offset |
    |---|---|
    | HB1 stave 6 | −1.55 µm |
    | HB0 stave 1 | +0.97 µm |
    | HB1 stave 8 | +0.9 µm |
    | HB0 stave 6 | −0.7 µm |

    The slope along z is ≤ 0.3 µm per 10 cm and dz is unchanged. The structure along z inside the staves is already in
    IBEdge.
- **Cause: the IB stave corrections are not collinear.**
  - For B≠0 the IB sub-track is a parabola through the vertex and Lr0–Lr2. Extrapolating it to 12 cm magnifies a sub-µm
    non-collinearity by about 30.
  - A rigid IB stave takes the field at its own φ (30°, 22.5° and 18° pitch). Where k(φ) varies, the three layers of one
    track are not on one parabola.
  - The non-collinearity predicted from the v3R stave TX alone gives a Δpar0 that correlates at 0.73 with the measured
    B± change (slope 0.65).
  - It is the stave-boundary sawtooth of the n=2 itself, not only the steep peaks. IB with n ≤ 2 of k still gives
    13–24 µm rms; n ≤ 1 gives about 6 µm.
- **Why pede leaves it:**
  - The 15° IB rdphi cells pin the IB rotation pattern. For example, one cell holds only Lr0 stave 3 and Lr2 stave 6.
  - In full-track fits the effect is a 0.5 µm misfit, so pede moved the IB staves by ≤ 0.5 µm.
- **Things that do not work** (model scripts in `sag/v4`, checks in `sag/tools/cr`):
  - **Tuning the IB stave TX to cancel Δpar0:** not robust to the vertex resolution. An optimum at σ_PV = 15 µm gives
    84 µm rms at 60 µm.
  - **IB stave radial shift:** it carries the field gradient across the stave and gives about 6 µm IB−OB, but it needs
    ΔR of 3.4/4.8/6.1 µm rms (Lr0/1/2, max 15 µm). That adds a dz-vs-tgl slope of the same size per stave. IBEdge has
    0.4–0.7 µm per unit tgl, uncorrelated, so dz would get 5–10× worse.
  - **IB stave RZ:** the 17° tilt gives dz = RZ·x·tgl, odd across the stave. It needs 0.40–0.46 mrad rms (max 1.2) and
    halves IB−OB (16 → 8 µm rms). IBEdge has 0.10–0.17 mrad-equivalent, uncorrelated, so dz would get about 3× worse.
  - **IB chip DOFs:** not in the recording. The jumps are at stave level anyway.

## 7. The open choice: gauge v3 (ob) or v4 (hybrid)
Straight-track model (`sag/v4/model.py`), changes relative to IBEdge:

| | v3 gauge | v4 hybrid |
|---|---|---|
| B± IB−OB par0 mismatch | 17–27 µm rms | 6 µm, any σ_PV |
| vertex DCA change rms (max) | 1.7 (6.5) µm | 1.3 (4.1) µm |
| Lr2 leave-one-out rms | 0.9 µm | 0.5 µm |
| OB Lr5/6 drφ rms (max) | 2–10 (37) µm | 24–42 (150) µm |
| OB ray rotation ⟨drφ/r⟩_OB, max over \|t\| ≤ 1.6 | 57 µrad (33 at \|t\| ≤ 0.8) | 256 µrad (136) |
| q/pT scale error ∂b/∂φ, max | 2·10⁻⁴ | 1.1·10⁻³ (6·10⁻⁴ at \|t\| ≤ 0.8) |
| k left | 0.08 ± 0.26 | 0.06 ± 0.43 |
| full-track visible misfit | 0.5 µm | 0.65 µm |

- **The trade-off is unavoidable.** Shrinking the IB share of the correction (β → IB) removes the IB jumps but
  proportionally rotates the OB.
- **A t-averaged variant does not help enough.** Putting only the t-averaged n ≥ 2 part in the IB gauge still leaves
  33 µm rms at Lr6.
- **v4's OB moves grow with |t|.** v4 − v3 is about 80 µrad at t = 0, 160 at 0.8 and 300 at 1.6. The high-|t| k map
  noise is multiplied by r(r − 3.45).

**TRD/TOF test of v4 against v3** (`sag/tools/trd`; the `solB1A4` debug trees, where the ITS is IBEdge and TRD/TOF were
aligned with B0 + cosmics):
- The GBL seed is an ITS-only refit (`CollGBLSeedITSOnly`), so `pre` at TRD/TOF is the ITS-only extrapolation.
- `trdB.C` computes per track the seed response to the per-chip correction tables: a line for B0, a parabola for B±.
  `ana*.py` compare the maps (18 sectors × 4 tgl).
- χ²/ndf with per-sector constants removed (IBEdge / v3 / v4):

  | field | TRD | TOF |
  |---|---|---|
  | B0 | 20.0 / 20.3 / 22.2 | 17.0 / 17.1 / 17.9 |
  | B+ | 22.9 / 9.3 / 7.8 | 8.3 / 3.6 / 2.9 |
  | B− | 43.5 / 18.2 / 16.3 | 18.2 / 10.5 / 9.6 |

- **Both corrections close most of the B0 − B± gap:** 2.2–2.5 mm rms → 1.0 mm (v3) or 0.95 mm (v4).
- **B± favour v4 by 10–20 %.** The OB part of v4 − v3 has a fitted coefficient of 0.9–1.7 (±0.1–0.25). For a seed with
  free curvature this OB-vs-IB rotation acts as curvature.
- **B0 disfavours v4 by 5–10 %.** TRD gives −1.5 ± 0.3, TOF 0.4–1.2 ± 0.4. B0 `pre` keeps a 0.9 mm n=2 pattern in the
  phase of k that neither gauge removes, and v4 adds about 0.2 mm close to opposite to it.
- **The B0 verdict is partly circular:** TRD/TOF were aligned with B0 against IBEdge.
- **Not conclusive.** The decisive test, proposed but not run: refit only the TRD/TOF DOFs on the solB1A4 recording with
  the ITS labels fixed at v3, and then at v4. Compare the per-field χ² and the consistency of B0/B+/B− `pre`.

## 8. Suggested next steps / questions for review
1. Run the TRD/TOF refit with the ITS fixed at v3 and at v4 (§7). Then do the same with v3R.
2. Find a gauge between v3 and v4 that also respects TRD/TOF. For example, choose β(φ) per harmonic so that the B0
   `pre` n=2 is reduced, if that pattern turns out to be ITS b and not TRD/TOF.
3. Understand the large-d0 cosmic step (§4). It is the only ITS handle on b, and it might decide between the gauges.
   This needs a d0-aware cosmic selection. Options: a debug-tree based fit, or a recording that splits cosmics by d0.
4. Reduce the high-|t| noise of the k map (smoothing or a prior in t) before any OB-heavy gauge.
5. Validate v3R or a successor:
   - replay per dataset;
   - B0 debug: k flat, LOO z rms back to the IBEdge level;
   - checkResid;
   - a variant with cosmics in the steering.
6. Open question: v10n has its own predicted k ≈ 0.4 (n=1 1.5 at −32°). It has not been checked with v10n B0 data.

## 9. Files
- **Corrections, README per version:** `/data0/devAlg/tests/sag/{v0,v1,v0R,v2,v2R,v3R,v4}`.
- **v3R reduced-fit outputs:** `millepede.res`, `millepede_labeled.res` and `pede.log`, in `/data0/devAlg/tests/sag/v3R`.
  - Its checkResid PDFs: `checkResid_PV1_hman.pdf` (against v10n), `checkResidRefV2_PV1_hman.pdf` and
    `checkResidRefBest_PV1_hman.pdf`.
- **Analysis tools:** `/data0/devAlg/tests/sag/tools` (copied from the session scratch; relative paths kept).

| directory | contents |
|---|---|
| `sagITS2.C` | per-track k and the LOO TH3 from `gbltrack` |
| `f100/` | 100 % debug runs and k map (`kmapX.C`, `half.C`, `d0s.C`, `even.C`, `z0.C`, …) |
| `geo/` | geometry-difference k/b/a decomposition (`geoK.C`, `bray.C`) |
| `rep/` | replay projections of the sagitta/rotation modes |
| `corr/`, `corr2/` | `mkSagCorr*.py`, volume tables `vols.txt` / `volsOB.txt` |
| `cr/` | checkResid per-sensor analyses: `sensMean.C`, `ib.py`, `iobpred.py`, `chipz.py`, `dzTgl.C`, `dzX2.C`, `dzcmp.py` |
| `trd/` | TRD/TOF seed-response test (`trdB.C`, `ana.py`, `ana2.py`, `ana3.py`) |

- **Model of the visible effects for any gauge:** `/data0/devAlg/tests/sag/v4/model.py`, with `radial.py`, `rz.py` and
  `rz2.py`.
- **Per-chip correction tables:** `/data0/devAlg/tests/sag/v4/chipcorr.txt` holds v3 and v4 drφ in cm, by sensor ID.
