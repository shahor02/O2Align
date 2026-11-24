# Plan: sagitta correction with the inner barrel as the reference (v5)

Date: 2026-11-23. Follows `ITSsagittaFix.md` (same directory); read that first. Written as instructions for the next
session. Nothing below has been run except the `model.py` numbers in §2.

## 1. Why the gauge is not free: the argument to verify

The handover treats the ray rotation b(φ, t) as a free gauge and the v3/v4 difference as a trade-off. The claim of this
plan is that for n ≥ 1 in φ it is not free, because the IB staves are rigid.

- A φ-dependent b·r in the IB is a stretch along the stave width. No rigid stave move makes it. The nearest rigid
  moves are a tangential sawtooth (seen by B± as the IB−OB jumps) or a radial shift / RZ (seen as dz vs tgl).
- IBEdge has neither: its B± IB−OB par0 has no stave-pitch jumps, and its dz-vs-tgl stave slopes are 0.4–0.7 µm per
  unit tgl against the 3.4–6.1 µm the radial variant would need.
- So the error of IBEdge in the IB is collinear: the IB of IBEdge is right up to a global rigid motion. Any correction
  that moves the IB staves by b(φ)·r (v0–v3R) puts in an error that was not there. The jumps are that error, measured.
- With the IB and the vertex as the reference, the OB correction is fixed by the B0 data:
  δrφ_OB = −k(φ,t)·r·(r − β_IB). Its "OB ray rotation" (50 µrad rms, 256 max) is then the size of the error IBEdge has,
  not an induced one. It is the same size as the v10n ↔ IBEdge difference (30–70 µrad), which came from the OB moving
  against the IB by 10–14 µm (n=2).
- Independent support already in hand: the B± TRD/TOF `pre` maps want the OB part of v4 − v3 with a coefficient of
  0.9–1.7. The B0 verdict against it is circular (TRD/TOF aligned to IBEdge with B0).

Consequence: there is no rφ field that removes k, keeps the IB jump-free and leaves the OB where it is. What can be
done is to make the OB move no larger than the data require (§3) and to test it on data not used to build it (§4).

## 2. Model numbers for the IB share (`python3 model.py` in `/data0/devAlg/tests/sag/v4`)

| gauge | B± IB−OB par0 rms (σ_PV 15/30/60 µm) | OB rotation change rms |
|---|---|---|
| NLO=99 (v3) | 16.9 / 18.0 / 27.0 µm | 7.5 µrad |
| NLO=1 (v4) | 5.7 / 5.8 / 6.5 | 44 |
| NLO=0 | 4.9 / 5.0 / 5.3 | 47 |
| NLO=−1 (IB reference) | 4.9 / 4.9 / 4.9 | 50 |

n ≤ 1 in the OB gauge buys 6 µrad and costs IB collinearity, so use NLO = −1.

## 3. Task A (do first, about an hour): do the large-d0 cosmics want the OB move?

This is the test inside the ITS, on 2025 cosmics, and it may also explain the top/bottom step of handover §4.

- **Idea:** a cosmic with |d0| = 5–20 cm misses the IB and crosses the OB layers at different φ. The field
  −k(φ)·r² then does not look like `line + k·u·r`, and what is left goes into the fitted step. Order of magnitude:
  k(n=2) = 1.9, r·u ≈ 1500 cm² gives 28 µm times a geometry factor below 1. The observed step is −27…+42 µm. Nobody
  has computed what the known k field predicts for the step. v3 and v5 predict different steps, because they differ by
  an OB rotation that is even in φ and so invisible only at d0 = 0.
- **Do:**
  1. Per-chip tables for v3, v4 and v5 (NLO = −1): `chipcorr.txt` has v3 and v4; `sag/tools/trd/run1.sh` shows how it
     was made. Add v5 from `model.correction(-1)`.
  2. In the cosmic fitter (`sag/tools/sagITS2.C`, the `step`, `kStep`, `d0` branches used by `f100/d0s.C`), add the
     linear response of the same fit (line + k·u·r + step) to each table, as `trd/trdB.C` does for the seed. Inputs:
     `Mille_IBEdgeS4R/cosmHP4LUT_ITS/*2025*/debug_alg.root`, then 2026 as a cross-check.
  3. Tabulate, in the `d0s.C` binning (9 φ × 5 |d0|): observed step, step after each correction, and the regression
     coefficient of the observed step on the predicted removal.
- **Read:**
  - v5 coefficient near 1 and the 5–20 cm step rms down by a factor ≥ 2, with |d0| < 2 cm staying at 0: the step was
    the sagitta, the OB move is confirmed by cosmics, and cosmics can go back into the steering.
  - v3 does this and v5 does not: the argument of §1 fails somewhere; stop and report.
  - Neither: the step is something else. v5 still stands on §1, but cosmics stay out of the fit and the report must
    say that cosmic compatibility is not shown.
- Also check the same response for the v10n → IBEdge geometry difference against the 15–26 µm step difference between
  the two geometries. It validates the response code.

## 4. Task B: build v5 and its reduced fit

- **v5 correction** (`sag/v5`, README, relative to the IBEdge recording):
  - IB staves and the vertex untouched. Compare in `model.py` "IB exactly zero" with the NLO = −1 stave TX
    (≤ 0.3 µm) and keep whichever gives the smaller B± IB−OB mismatch.
  - OB module TX + RY fitted to the chip-level field −k̃·r·(r − 3.45), as `mkSagCorr3.py` does.
  - k̃: the map noise is multiplied by r(r − 3.45) ≈ 1400 cm² at Lr6, so smooth with the map errors as weights. Keep
    the peaks at the half-barrel junctions (the geometry decomposition reproduces 7.1 of 7.1 at 340–360°). Where a
    high-|t| cell is within 2σ of the smooth surface, use the surface. Report the max OB drφ and the rotation
    before and after; the aim is to bring the 256 µrad max down to what is significant.
- **v5R reduced fit** (as v3R, with these changes):
  - IB stave TX fixed. No IB rdphi pins.
  - Free: OB module TX + RY, mean vertex. OB rdphi pins re-centred on v5; priors as in v3R.
  - Data: B0, B0Bis, B+, B+Bis; add the 2025 cosmics if Task A confirms.
- **Validation, in this order:**
  1. `model.py` metrics and the replay per dataset (also for v3R and v4, which were never replayed).
  2. B0 debug: k flat, leave-one-out z rms at the IBEdge level.
  3. checkResid against IBEdge (the user makes the geometry and runs it): no IB jumps in B±, OB dips not back.
  4. `trd/trdB.C` + `ana*.py` with the v5 table, then the decisive test of the handover: TRD/TOF-only refit on the
     `solB1A4` recording with the ITS fixed at v3 and at v5, comparing B0/B+ `pre` consistency. B− is validation only.

## 5. What to report

- Whether §1 holds (Task A, and whether IBEdge B± IB−OB is jump-free at the 5 µm level).
- The v5 trade-off table in the format of handover §7, with a v5 column.
- Observations go to `doc/history.md`. No commit unless asked.

## 6. Review after the v4 checkResid (2026-11-23, Opus)
See the `doc/history.md` entry for the numbers.
- **§1 holds for n ≥ 2:** v4 already removes the IB jumps with n ≤ 1 left on the IB. The IB data alone do not decide
  between n ≤ 1 and IB-untouched (NLO = 1 vs −1); TRD/TOF and Task A must decide that.
- **A gap §4 does not cover:** at tgl ≳ 1, mostly on the +z side, B0 wants the k correction and B± reject it with the same
  size (Lr4/Lr5: −4…−8 µm in both fields).
  - This is the same in every gauge (v0/v2 show it too), so v5 will show it as well.
  - Smoothing k̃ will not fix it: it is a coherent B0/B± disagreement, not map noise.
  - A reduced fit trades it back. That is what v3R did, keeping the B0 high-tgl sagitta.
- **Add before Task B:**
  1. A variant with k(t) frozen at its |t| = 0.9 value beyond that (or tapered to 0). Run checkResid and see whether the
     B± OB +z dip disappears while the IB stays clean.
  2. Add multiple scattering to `model.py`: Kalman weights with x/X0/sinθ per layer and a realistic pT spectrum. Without
     it the "invisible to B±" statements fail at high tgl.
  3. Find out what the high +tgl B0 k is. It is odd in t and is also seen with OB points only. Candidates: a
     B0-specific effect at large incidence, or geometry that the B± refit absorbs differently. Cosmics at high tgl, few
     as they are, and the TRD/TOF B0 − B± gap versus tgl are the external handles.
- **Task A result (2026-11-23, Opus):** negative for v5.
  - The large-d0 cosmic step is not explained by the sagitta in any gauge.
  - Its pattern anti-correlates with the OB-rotation response: coefficient −1.4 to −3.7, χ² worse with v4 or v5 in both
    2025 and 2026. v3 is about neutral.
  - The method was validated on the v10n − IBEdge step difference (c = 1.27 ± 0.11).
  - Numbers are in `doc/history.md`, tools in `sag/tools/cosm`.
