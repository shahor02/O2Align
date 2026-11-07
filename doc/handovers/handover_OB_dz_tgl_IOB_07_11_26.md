# Handover: OB dz-vs-tgλ structures in checkResid (IOB vs standard)

Date: 2026-11-07

Data: `/data/devAlg/solA2v19/checkResid_PV1_hman.pdf` (standard) and `checkResid_PV1_iob_hman.pdf`
(OB clusters vs IB+PV sub-track extrapolation, IB clusters vs OB sub-track extrapolation), p.4 = dz vs Z,
p.8 = dz vs tgλ. Only p.4 and p.8 were looked at, by eye. Tree (`res`) not yet used.
Background: memory `ob-dz-tgl-incidence`, skill `checkResid`.

## Question from the user
The OB dz-vs-tgλ structures look to the user like rigid-body movements, with the outer two OB layers (Lr5, Lr6)
as culprits (IOB p.4 and p.8).

## Established in this session
- Standard residual of Lr5 = Lr5 vs (vertex+Lr0..Lr4); the inward prediction (Lr6 alone after a covariance reset)
  has negligible weight. So "Lr5 vs IB+L3+L4+PV" needs no separate test; Lr6 has no outer anchor.
- Standard p.8: Lr5 is within a few µm for |tgλ|<1.4, with a ±40 µm sigmoid only at |tgλ|>1.5, identical for V10n
  and MP2, all fields. Lr3 has the same-sign ±40 µm edge sigmoid, Lr4 (between) has none. Lr6 only ±10–20 µm.
  → Lr5/Lr6 are not displaced relative to Lr3/Lr4 at |tgλ|<1.4.
- IOB p.8: Lr3 small, Lr4 ~40–100 µm, Lr5 ~100, Lr6 150–250 µm; Lr6 dips near tgλ≈1.2 by −100 µm in B± vs −50 µm in B0
  (same geometry); V10n and MP2 differ at +z. IB panels (Lr0–2) flat to ±4 µm. p.4 and p.8 are the same tracks
  (tgλ ≈ (Z−z_PV)/R), not independent evidence.
- Hypothesis (back-of-envelope, unverified): the IOB signal is an OB-vs-(IB+PV) mismatch growing with R, not an
  Lr5/Lr6 fault. A common IB radial error Δ gives dz ≈ tgλ·Δ·R_OB/R_IB (×13 at Lr6, ×6.5 at Lr3); the linear part of the
  Lr6 curve (~ −55 µm per unit tgλ) needs Δ ≈ 4–5 µm. IB-internal residuals are blind to it (helix stays straight);
  the IB panels of the IOB plot would see it shrunk ×13, consistent with ±4 µm. Does not explain the tgλ≈1.2 dip
  or the B0/B± difference. A PV z error is an alternative.
- The nonlinear sigmoid at |tgλ|>1.5 is not a rigid shift (dR·tgλ is linear); a stave/chip tilt or bow (dR(Z)·tgλ)
  could produce it. Identical across alignments and fields points to reconstruction/incidence.

## Proposed next steps (not done)
1. From the `res` tree: IOB dz-vs-tgλ slope per OB layer vs R; if ∝ (R−R_IB) one IB or PV radial term explains it,
   if faster there is a genuine OB term.
2. `trIBOut` − `trOBInw` par3 (tgλ) and par1 (Z) at R=12 cm vs tgλ: a direct IB/OB slope mismatch.
3. Closure: inject ±5 µm on IB radii with `misAlgJson`; check it reproduces the IOB slope and leaves the standard
   residuals untouched.
4. Optionally free OB chip TY/RX/RY (or a Z-dependent mode) and see if the sigmoid shrinks.

## Tree check (2026-11-07, later): the hypothesis above is superseded
Trees are `/data/devAlg/solA2v19/checkResid_{LHC25ab_562837 (B+), LHC25ab_562862 (B0), LHC25ah_564644 (B−)}_AlgRS_PV1.root`,
MP2 alignment only. Means come from iterative ±2σ Gaussian core fits per tgλ bin. Slopes are linear fits for |tgλ|<1, in µm per unit tgλ.
The scratch macros (`iobSlopes.C`, `anaSlopes.C`, `dCtest.C`) are in the session scratchpad and are easy to redo.

| slot | std B+ | std B0 | std B− | IOB B+ | IOB B0 | IOB B− |
|---|---|---|---|---|---|---|
| vtx | −2.7 | −1.1 | −0.1 | −0.8 | −9.3 | −5.3 |
| Lr0/1/2 | −0.7/−0.6/−0.8 | −0.4/−0.4/−0.4 | −0.4/−0.1/−0.3 | +3.0/+3.3/+3.3 | −3.8/−3.2/−2.2 | −4.5/−3.3/−3.1 |
| Lr3 | −4.9 | −1.7 | −2.5 | −16 | −9.0 | −8.0 |
| Lr4 | +1.1 | −1.3 | −1.5 | −25 | −12.7 | −16 |
| Lr5 | −7.2 | −2.7 | −3.4 | −57 | −25.5 | −45 |
| Lr6 | +1.6 | −1.1 | −0.3 | −77 | −35.3 | −66 |
| IB−OB par1 (Z at R=12) | | | | +8.8 | +3.4 | +1.4 |
| IB−OB par3 (tgλ), 10⁻⁶ | | | | +61 | +71 | +52 |

- **In B≠0 the IOB OB curves are dominated by an estimator bias, not geometry.** The IB(+PV) sub-track has almost no
  curvature information: median σ(q/pT)/|q/pT| = 0.37 (pT>1.4), 0.42, 0.60 (pT<0.67), against 0.10 for the OB sub-track. The
  extrapolated z = z0 + tgλ·s(R) with s ≈ R + C²R³/24 is quadratic in C, so E[z_pred] is biased by +tgλ·σ_C²·R³/24 and
  dz = −tgλ·σ_C²R³/24.
  - Lr6 IOB slope by pT: B+ −127/−71/−54, B− −121/−59/−36 (pT <0.67 / 0.67–1.4 / >1.4).
  - At fixed pT, split by |C_IB/C_OB−1|: <0.1 → −56, 0.1–0.3 → −59, >0.3 → −231 (B+, pT<0.67); B− gives −48/−52/−240.
  - Geometry cannot depend on how well the IB sub-track measured the curvature, so the 100–250 µm and the deepening of the
    tgλ≈1.2 dip in B± vs B0 are this artefact.
  - par3 (tgλ is constant along the helix) has no such bias, and it is field-independent (52–71·10⁻⁶).
  - The IOB residual needs the IB sub-track's q/pT fixed to the full-track or OB value before it is usable in B≠0. This is a
    CheckResidSpec change, not done.
- **The B0 IOB is clean, and its geometric part is much smaller than the handover assumed.**
  - Lr3/Lr4 (−9.0/−12.7) agree with −(par1 + par3·(R−12)) = −8.8/−12.2. This is an IB+PV vs OB rz-slope mismatch of
    Δtgλ ≈ +7·10⁻⁵·tgλ (IB steeper). As a common IB radial shift through a vertex pivot that is Δ ≈ 2 µm, not 5. As an
    IB/OB relative radial scale it is 7·10⁻⁵.
  - Neither simple variant fits the vertex IOB value (−9.3, OB line extrapolated to r=0), so the model is not settled.
  - par3 vs tgλ is non-linear: it rises to +2·10⁻⁴ at tgλ≈1.25 and falls beyond 1.5, and V10n has a different shape (PDF
    p.12). It is alignment-dependent, so it is geometry, probably Z/stave-end related.
- **Lr5/Lr6 (the user's suspicion).** In B0, Lr5/Lr6 IOB exceed the Lr3–Lr4 line by about −5 and −11 µm/unit tgλ. That would
  mean a radial-type offset of order 5–10 µm of the outer layers relative to Lr3/Lr4. The standard Lr5/Lr6 residuals
  (−2.7/−1.1 in B0) do not confirm it at that size; with linear estimators they should show at least −5. A few-µm effect of
  Lr5/Lr6 is therefore possible, but it is at the level where core-fit means of very different-width distributions are not
  additive. It does not account for the large IOB structures.
- The PV is common to all tracks of an event, so it cannot make a per-track odd-in-tgλ term beyond the 1/Nc self-weight.
  It is dropped as a candidate.

Revised next steps:
- (a) In CheckResidSpec, give the IB sub-track the full-track (or OB) q/pT with a small error before the IOB extrapolation,
  or cut |C_IB/C_OB−1|<0.1, then redo the IOB PDF.
- (b) Study par3 vs (z_PV, tgλ) or vs the OB-hit Z, and compare V10n with MP2, to localise the 7·10⁻⁵ mismatch.
- (c) Use a 5 µm IB-radius closure only if (b) points to the IB; ~2 µm is the relevant size.
