# History of the alignment development

Chronological record of the observations that led to each constraint, DOF-config clause and code change, and of
problems found and solved on the way. The functional documentation (`DOFConfig_rules.md`, `WeakModes.md`, code
comments, …) describes only what the code does and how to use it; **why it is there** is recorded here.

How to maintain it: when a development is motivated by an observation, or a bug/degeneracy is found and fixed, add
an entry here (date, what was seen, what was changed) instead of narrating it in a code comment or a feature doc.
Entries are in chronological order; dates are those of the work (see `git log`). Current investigation state, to be
picked up by the next session, lives in `doc/handovers/`; superseded investigation notes in `doc/obsoleteNotes/`.

Run directories quoted below are `/data/devAlg/...` or `/data0/devAlg/tests/...`; configs are `doc/ITS_RB_*.json`.

---

## 2026-10-04 … 10-10: first ITS rigid-body fits, `ITSBarFix`, `StavesOnly`, `savStavesOnlyPinHB`

Context: incremental pede rounds on ITS (runs LHC25ae B=0.5 T and LHC25ab/562862 B=0) gave a factor ~5 in pede
χ²/ndf while `checkResid` residuals developed new patterns or got worse. All tests derive from
`ITS_RB_from_hstaves_IB.json` (the config that produced the Mille data), changing which DOFs are free at the
stave / half-stave / chip level and whether `pinChildrenMean` is used.

- **`ITSBarFix` (IB stave free in all 6 DOFs, IB chip free in TX,TY,TZ):** a strong pattern vs Z, correlated with
  the chip, on IB staves only. RX/RY of a stave transport to a child at offset z as a translation ∝ z, and a chip
  with translation but no rotation is the only thing that can absorb that lever arm, so a stave tilt and a
  linear-in-chip-index ramp of the nine chip TY are numerically almost indistinguishable. The automatic
  children-mean constraint fixes only the 0th moment, not this dipole. Confirmed in `millepede_labeled.res`: chip TY
  ramps almost linearly with chip index on every IB stave (Stave0/Lr0: +114 → −68 µm, chip0 → chip8), with a slope of
  the order of the stave's own RX. ⇒ the rule of thumb *one representation per physical deformation*: do not free
  a stave tilt together with the translation of its chips.
- **`StavesOnly` (only the stave level free, everything above fixed, no `pinChildrenMean`):** OB residuals worse than
  the initial geometry. A fixed parent with free descendants and no `pinChildrenMean` writes no constraint at all
  (`mp2con.txt` was empty), leaving the common rigid mode of all staves of a half-barrel unconstrained: mean stave TX
  and its RMS grew with the layer radius, −21 µm (Lr0) → −537 µm (Lr6), RMS up to 2.5 mm. ⇒ every fixed parent whose
  descendants are free needs `pinChildrenMean` (see the DOF-config documentation).
- **`savStavesOnlyPinHB` (as above + `pinChildrenMean:"all"` at the half-barrel):** the pin worked (Lr6 TX mean −537 →
  −81 µm, RMS 2.5 mm → 64 µm) but checkResid barely improved: a rigid stave captures only the stave-average part of
  the real misalignment, which has sub-stave structure (basis mismatch, not conditioning). pede's χ²/ndf is an
  in-sample number (42 M equations vs 1182 free parameters), so adding free parameters lowers it mechanically; the
  leave-one-out checkResid is the honest estimator and shows the OB geometry got worse because noisy stave
  corrections (poor OB observability) were applied on top of an already good geometry.
- **Rule matching moved from `fnmatch` to `std::regex`** (whole-name `regex_match`, implicit leading `.*`,
  patterns compiled once per rule). Numeric ID fields are written `[0-9]+`, which cannot match `/` and so cannot
  swallow a child level, and which removes the one-/two-digit OB-stave duplication; the old two-digit-only half-stave
  rule had left OB staves 0–9 without any half-stave freedom. All `ITS_RB_*.json` were converted (`?` → `[0-9]`).
  Inapplicable-rule warnings are suppressed unless the pattern has no regex metacharacters (`isLiteralPattern`).
- **Revised config `ITS_RB_from_hstaves_IB.json`:** ITS envelope fixed with its own `pinChildrenMean`
  (TX,TY,TZ,RZ); IB half-barrel freed in TX,TY,TZ,RZ with a 30 µm `measurement` prior (the two IB half-shells meet at
  a mechanical gap, so a hard zero-mean pin would have forbidden exactly the mode a gap-closure error produces; the
  prior regularises but lets the data move it); OB half-barrel fixed + `pinChildrenMean:"all"`; IB stave without
  RX,RY (the degenerate channel of `ITSBarFix`); IB chip TX,TY,TZ; OB stave / half-stave `all`.
- **Mille data production:** generate `MilleData` once with every DOF free (`[{"match": ".*", "rigidBody": "all"}]`,
  or one rule per level for ITS only; `ITS_RB_all.json`), and explore free/fixed/constraint variants downstream at the
  `MilleSteer`/pede stage. Trap found: `Volume::writeParameters` writes parameter lines only for volumes having a
  rigid-body DOF set, so labels of volumes matched by no rule were missing from `mp2param.txt` while their
  derivatives were in the binary, and pede treated them as variable (`solA0`: NVGB = 146016, void). Cure: a catch-all
  `{"match": ".*", "rigidBody": "fixed"}` as the *first* rule (later rules replace the DOF set), so every unmentioned
  volume is written as an `F` line (presigma −1).

## 2026-10-06: pede does not converge on B=0 multi-track vertex records

Multi-track records (all tracks of a vertex in one composed GBL trajectory with a common vertex,
`useMultiTrackPVConstraint>0`) stopped with "Too many rejects" on pp B=0 data (run 562862), unlike single-track
records. Findings (notes kept in `doc/obsoleteNotes/multiTrackGBLVertex_handover.md`): vertex sign conventions and the
zeroed curvature at B=0 are right; the divergence is a free gauge mode (mean vertex vs ITS envelope translations,
fix one of them) combined with an inconsistency in pede's outlier weights; the remaining abort is pede's χ² cut,
most likely from multiple scattering underestimated at B=0 (hence `meanPtB0`, a low imposed pT is safer); and the
per-track reference states did not pass through the common vertex, which `vtxMultiTrackRefCovScale` (KF update with
a scaled vertex covariance) now enforces.

## 2026-10-11: `solA1` (`ITS_RB_v0.json`): the first valid run, radial scaling

`solA1`: NVGB = 4134, 1084 constraints, full rank, all 157872 labels present in `mp2param.txt`; the regex migration
worked in production. checkResid still degraded at the OB, so the cause was not bookkeeping. All 1084 constraints
were satisfied to 1e-8. Decomposing the stave corrections into a coherent r-term plus per-half-barrel harmonics
(ℓ ≤ 2) left almost nothing:

- the stave TY (sensor normal ≈ radial) shows a **uniform radial scaling ε = −5.4·10⁻⁴**, linear in r over all seven
  layers (−13 µm at Lr0, −214 µm at Lr6); the genuinely local stave-by-stave part is 1–2 µm (IB) and 20–55 µm (OB).
  (Local axes: X = rφ, Y = sensor normal, Z = beam, from `SensorITS::defineMatrixL2G` / `getFlatChipAlpha`.)
- Why such modes survive: the hierarchical constraint is a *rigid-body* constraint (6 equations per parent: net
  translation and rotation of the child set); radial scaling, ℓ=2 squashing and twist carry no net rigid component and
  pass through, and no `pinChildrenMean` setting can fix that: it is the wrong operator. A uniform dilation of the
  transverse plane is an almost exact null direction of the track χ² (absorbed by 1/pT and the impact parameter);
  only a fixed external scale breaks it. It hurts checkResid at large radius through the φ-harmonic debris the pin
  forces in, plus a B=0 vs B=+0.5 T split (~15 µm) from the momentum-scale degeneracy.
- The mean vertex had escaped the catch-all (it carries a *calibration* DOF set) and was free with no prior: fitted
  X = +334 µm, Y = −1820 µm, Z = −1808 µm. ⇒ `ITS_RB_v1.json`: explicit `calib` clause freeing only X, Y, Z with priors
  (50 µm in X/Y, 1 mm in Z; the `measurement` must be a rule of its own after the `calib` rule). Also IB chip reduced to
  TX,TZ (the in-plane measured directions): the chip TY, 15–29 µm RMS over silicon known to a few µm, only absorbed noise
  and fed the `ITSBarFix` ramp.
- `mpSteer.txt` used `sparseMINRES`, which gives no parameter errors (error column 0.0000) and no spectrum; also the "×5
  χ²" was iteration 0 → 1 only, with the cut factor changing 20 → 10 across that step. ⇒ use `method inversion`, and
  `diagonalization` once to see the weak modes as the smallest eigenvalues.
- **`globalModes` clause.** First as a stop-gap, `macro/add_global_mode_constraints.py` (appends radial, curl,
  elliptical equations to an existing `mp2con.txt`, using nominal layer radii and φ from the stave index), then as
  `Volume::writeGlobalModeConstraints`, a sibling of `writeRigidBodyConstraints` that projects the analytic mode
  field **u**(x) (translations) and **ω** = ½∇×**u** (rotations) onto each free DOF using the volume's own `getL2G()`.
  Verified on the `v0` free set: 5 equations, rank 1089 (independent of the 1084 mean constraints, i.e. these modes lived
  in their null space); `radialScaling` coefficients ∝ r on stave TY plus TX terms on IB staves (the real IB stave tilt);
  `twist` lands only on OB stave RY; evaluated on `solA1`, Σcp/Σ|c| = +174 µm for `radialScaling`, −1.5 µm `curl`,
  +50/−20 µm the elliptical phases, −0.18 µm `twist` (never a problem). `zScaling` and `bowing` are skipped (with a
  warning) at stave level, since a stave is centred at z = 0 and their coefficient ∝ z.
- Priors cannot replace global-mode constraints: the honest mechanical tolerance (O(100 µm) for an OB stave radius) is
  far too loose to cap a 214 µm mode. Priors belong where external knowledge exists (IB half-barrel gap, beam spot,
  survey).

## 2026-10-12: `inversion` runs, the pin is paid off by lever arms; `ITS_RB_v2.json`

Four runs (`resOp/solA1Iv0`, `Iv1`, and `…GC` with the python constraints), `method inversion`:

- The constraints worked: radial scaling ε −5.4·10⁻⁴ → −5.6·10⁻⁸, at +0.29 % pede χ² (the proof that it was a null
  direction); `v1` ≈ `v0` (+1.0 % χ² for 444 fewer parameters).
- **The hierarchical pin is paid off by a lever-arm cancellation.** Printing one pin equation term by term (OB stave
  TY): half-staves moved −119 and −75 µm (mean −97 µm) and the equation was satisfied to 5·10⁻⁴ µm because each took
  an RZ of ±7.5 mrad, equal and opposite to 6 digits. The pin manufactured a 0.43° rotation to pay for the drift. The
  cancellation partner is RZ (rotation about the sensor normal, ~25 µm of unmeasured normal displacement over the
  stave half-width, weighted by the ~41 cm lever arm) at the envelope, OB half-barrel and OB stave; only IB stave →
  chip had no free rotation and was genuinely constrained, which is why the chip mode appeared as a ramp.
- **What survived** (ε = amplitude/r per layer, coherent across the independently mounted IB and OB, which a real
  deformation cannot be): tangential elliptical u_φ ℓ=2 (ε 1.9–2.5·10⁻⁴, 105 µm at Lr6), a **global rigid tilt about Y
  of −155 ± 20 µrad on all 14 half-barrels**, u_r ℓ=3, ℓ=4, u_φ ℓ=1; plus a localised OB clamshell (relative
  half-barrel motion ΔT_y = +36…+62 µm at Lr3–6, ΔR_x = −200 µrad) and the IB chip TY dipole (span 76/48/26 µm for
  Lr0/1/2, i.e. ∝ 1/r). The 0th moment is constrained, the 1st is not; with IB stave RX removed the mode moved into the
  measured directions as a z² bow of the chips (chip TZ +23 → −10.6 → +12.5 µm).
- **Implemented:**
  1. `pinChildrenMeanNoCross`: the mean equation with the jacobian cross-terms dropped, so only the children's own
     same-named DOF enters and no rotation can pay for it. Do not request both pin forms on one volume (14 of 849
     equations were rank-deficient when both were set at the envelope and the OB half-barrels). It cannot see a global
     rigid shift or tilt (cos φ / sin φ pattern, plain mean 0).
  2. Generic harmonic families `radialHarmonic`/`tangentialHarmonic` (u = ε·r·cos/sin ℓφ, ℓ = 1…4; the factor r is the
     measured shape: ε coherent across layers) and `globalTilt_x/y`; equations indexed by a flat index
     (`Volume::NGlobalEq`), flaggable one by one; propagated to `WeakModes.cxx` and `drawWeakModes.C`.
  3. `ITS_RB_v2.json`: IB chips fixed, IB stave RX/RY restored (the stave tilt is the physical parameter), stave and
     half-stave RZ fixed (buys ~25 µm of unmeasured displacement at a 41 cm lever arm in every constraint),
     `pinChildrenMeanNoCross` at the envelope and the OB half-barrels, the global-mode set at stave level. NVGB 2427,
     803 constraints (21 global modes), full rank.
  Validation of the new equations on `solA1Iv0GC`: Σcp/Σ|c| ≈ 0 where already constrained; −44/+49 µm
  `tangentialHarmonic2`, −41 µm `globalTilt_y`, +34 µm `radialHarmonic3`, −19/+19 µm `radialHarmonic4`.

## 2026-10-12: `solA1Iv2GC` (v2) in checkResid; what needs DOFs distributed in z

- The apparent growth of the OB dy-vs-φ scatter in the first v2 PDF was a binning mismatch (30 φ bins in the V10n
  hman files vs 100 in MP2). At equal binning v2 is better than V10n: dy Lr3–6 11.8/11.1/28.1/7.8 → 8.5/8.5/19.7/6.5 µm,
  dz 11.7/9.7/23.7/13.9 → 3.8/2.9/8.9/2.8 µm; the same in pp B0 and OO B+.
- Inside an OB half-stave, the per-chip dy RMS is still 12–24 µm, dominated by an offset between the two chip rows of a
  module (rows 7–13 vs 0–6: mean −10…−22 µm), reproduced in pp B0 and OO B+ (per-chip correlation 0.73–0.84); no v2 DOF
  can touch it (→ OB chip TX, see 2026-10-22).
- The large-scale Z/tgλ structures (checkResid p.3/4/7/8) are identical for V10n and v2 and in both datasets: the fit
  neither made nor tried to remove them, since below the stave everything is fixed and a rigid stave only makes
  structure constant or linear in z. Fitting every sensor's mean dy, dz on {1, z_loc, tgλ, z_loc·tgλ}: IB: the 3 cm
  sawtooth is an in-plane rotation of each chip, the same on all staves (Lr0 dy/dz_loc −3.5…+2.5 µm/cm for chips 0 → 8);
  end chips carry the stave-end deviations (dy −2.8/−4.9 µm on chip 8 of Lr1/Lr2, dz +5.6/+2.4 µm on chips 0/8 of Lr0);
  a normal displacement Δn gives dz ≈ Δn·tgλ with tgλ up to ±5 at the Lr0 end chips, so normal positions are measured
  at the stave ends and unmeasured at z ≈ 0. OB: features at the module scale (21 cm) on the outer modules, a module
  rotation about the normal (Lr3 dy +25 µm over 14 cm) and a bow-like odd dz vs tgλ (up to +45 µm at tgλ = 2, Lr3);
  only ~200 hits per chip, so per-chip fits are ill-conditioned. Lr5 dz is anomalously broad (RMS 350 µm vs 220–245 µm).
- ⇒ **v3 idea:** free IB chip `TX,TZ,RY` + `TY` with a ~5 µm prior, OB module `TX,TZ,RY` + `TY` with a 30 µm prior, drop
  stave/half-stave RX (a tilt would again be a TY ramp of the children).
- **New weak-mode family.** For a straight track from the vertex z = r·tgλ, so fields r·f(φ, tgλ) (rφ, z, r) are
  absorbed by the track's φ, tgλ, curvature (invisible at B=0), and f(tgλ) by z0, broken only by the vertex. The
  existing modes cover only the t-independent members, which sufficed while every free DOF sat at z = 0. ⇒
  `radialPolar`, `tangentialPolar`, `longitudinalPolar` (u = r^p·(z/r)^k, p ∈ {0,1}, k ≤ 2, duplicates of existing
  modes excluded; 14 equations, `GlobalPolarKMax`; rotation ω checked by finite differences), plus `zScaling`.
  `bowing` (Δr = εz², one ε for all layers) is not a null direction; its null counterpart is `radialPolar_r0t2`. The IB chip
  TY dipole is `radialPolar_r0t1` (Δr ∝ z/r).

## 2026-10-12: v3, two bugs (pede grouping crash, void `WeakModesDiag` fit B)

- **pede crash with long constraints** (`free(): invalid next size`, v3 with up to ~9000 terms per equation). A
  bounds-checked build of Millepede-II 4fd638f67f locates it in `grpcon` (`pede.f90:1902`, `globalParCons(0)`):
  constraints are grouped in the order of their lowest parameter index; a constraint joining several earlier groups
  relabels only its direct neighbours, which can leave an earlier group empty, and the redundancy check then indexes out
  of bounds. A replica of the algorithm gave v2: 436 groups, one empty (survived by luck); v3: 3 groups, two empty.
  Workaround: `Volume::writeGlobalModeConstraints` gives every mode equation the lowest free rigid-body label as an
  anchor with coefficient 1e-12 (not 0: pede's reader silently drops zero factors), so pede processes these equations
  first and forms a single group (`addPedeGroupingAnchor`, `minFreeRigidBodyLabel`). Any other long constraint (ray
  modes) needs the same anchor.
- **`WeakModesDiag` fit B was void** (also in v2): `Eigen::BDCSVD` on a design matrix exactly rank deficient by
  construction (`curl` and `globalTilt` duplicate the rigid rotation columns) returned NaN singular values (v2 rank 12,
  v3 rank 0). `JacobiSVD` gives rank 42/45, the three missing being exactly those duplicates. Metric B values in earlier
  `weakModes.json` are void; metric A was unaffected.
- **v3 fit** (`solA1Iv3`, 10454 parameters, 2010 constraints at rank): χ²/ndf 1.025 vs 1.12, objective −2.8 % vs v2, but
  the corrections grew a lot (sensor displacement RMS rφ IB 6–7 → 12–14 µm, rφ OB 20–45 → 33–85 µm, z OB 7–14 → 30–66 µm;
  Lr5 module TX up to 864 µm) while per-chip checkResid residuals are 2–4 µm (IB). Binning the sensors in (φ, t = z/r)
  cells and dividing by r: 87 % of the variance of dz/r is common to all seven layers (70 % already in v2 at a third of
  the amplitude), i.e. Δz = r·g(φ, t), and rφ/r is common within the IB (0.99) and within the OB (0.97): the **ray
  modes**, a whole function of (φ, t) per component, which short lists of polar/harmonic modes cannot cover. checkResid
  is blind to them by construction. Options considered: per-cell ray constraints, priors at the placement tolerance,
  data breaking the degeneracy (cosmics, tracks not from the vertex).

## 2026-10-12: v4, `rayModes` clause and priors; IB stave RY dropped

- **`rayModes` clause** on the detector envelope: sensors binned in (φ, t = z/r) cells (24 × 16 over |t| ≤ 4 gives 720
  equations); per cell and component (rφ, z, r) one equation says the cell's mean angular displacement seen from the
  vertex vanishes, with every free DOF on the sensor's path to the root entering through the sensor's linearised
  displacement (finite differences of the same `AlignParam` matrix used to apply corrections, so conventions match by
  construction). Cells whose sensors span less than 1 cm in radius (a single layer) are skipped, moving one layer is
  measurable. Must not be combined with global modes lying in the ray span (`curl`, `radialScaling`, harmonics,
  `globalTilt`, p = 1 polar terms): near-degenerate constraint set; v4 keeps only `twist` and the p = 0 polar terms.
- **Priors** (assumed tolerances, not survey numbers): IB chips TX,TZ 10 µm, TY 5 µm; OB modules 30 µm.
- **IB stave RY dropped:** with it free, the fit expressed a common in-plane rotation of a stave's chips as stave RY plus a
  compensating chip-TX ramp (chip TX RMS 32 µm against its 10 µm prior, correlation with RY·z 0.88), the same double
  representation as RX before; without it chip TX RMS 15 µm, χ² slightly lower with 48 fewer parameters.
- Result (`solA1Iv4`): 2655 constraints at rank, one pede group; χ²/ndf 1.030; all-ITS sensor displacement RMS 104 → 42
  µm (v2: 59), OB rφ 13–20 µm, OB module TX RMS 22 µm (max 91); cost 0.6 % of the χ² that v3 gained. Left: an
  anti-correlated Lr3-vs-Lr6 pattern in rφ/r and dr (−0.97, −0.80), zero mean per cell, sagitta-like (∝ r²), absorbed
  by curvature when B ≠ 0 (to be broken by B=0).
- **checkResid on v4:** the OB improves, but the IB and vertex panels get worse than v2 (Lr0 dy vs Z 1.7 → 4.6 µm,
  Vtx dz vs φ 1.8 → 5.4 µm, Lr6 dz vs Z 8.7 → 15.8 µm). The Lr0 dy-vs-Z sawtooth is the new IB chip RY: a common rotation
  per layer of +4.97/+1.15/−2.26·10⁻⁴ rad at Lr0/1/2, whereas the V10n in-chip dy slope that checkResid measures is
  +0.09/+0.02/−0.53 µm/cm. Ruled out: rotation sign/pivot convention (`AlignParam` is right-handed and active, as the
  derivatives), wrong application (`algdump.root` shows the chips rotated exactly as fitted, nothing compensating),
  the IB part of the ray cells (not a null direction since σ(z_vertex) = 4.9/5.4 cm; v5 = v4 + `rMin` 10 cm for `rayModes`
  gave the same χ² and chip parameters; `rMin`/`rMax` keys stay as options). The fit (GBL on Mille data) and the
  reconstruction (KF in checkResid) disagreed at IB chip level — **this was the cluster-x bug of 2026-10-20/21**.

## 2026-10-20 … 10-21: cluster-x fix

`DetectorITS::prepareData` stored the nominal sensor reference x (`getSensorRefX`) as the x of every cluster, the fit
propagating to that x; checkResid/`refitITSTrack` use the cluster's own tracking x. They coincide only if the chip plane
is at constant x, but the V10n IB chips are pitched about the local x (dx_trk/dz_loc = +2.3/+0.6/−1.3 mrad for Lr0/1/2,
up to ~35 µm over a chip), leaving dy = tan φ_inc·θ·z_loc and dz = tan Φ_z·θ·z_loc, absorbed by a common chip RY
(+480/+104/−234 µrad) and the sawtooth. Fix: `double x = trkXYZ.X();` in `prepareData` (ITS3 tiles keep their own branch);
chip RY per layer −37/−58/+26 µrad afterwards, IB dy residuals flat to < 1 µm. Details:
`doc/handovers/handover_patterns_afterFixX.md`.

## 2026-10-21: `milleSplit` first result invalid (mean vertex fixed)

First result of `macro/milleSplit*` on `/data/devAlg/mille{B0,BPlus,BMinus}` (251k/217k/235k records): split − rigid
Δχ² = 88 019 / 101 015 / 95 372 (2592 extra parameters); a single y-only incidence term β (∝ tan φ_inc) per sensor carried
80–95 % of it (β ≈ +8…+18 µm, per-sensor scatter 35–115 µm rms). **Invalid:** the tool kept the mean vertex fixed at (0,0)
while the 2025 beam is at about (+0.3, −1.8) mm, which pulls every vertex by ~40σ transversely and fakes an rφ-only chip
term. Fixed by freeing the mean-vertex calibration globals in every model (see also the memory note on offline Mille
fits). The checkResid view (Lr0 end-chip dy-vs-Z excess, vertex "vortex" at |zv|>12, the missing dz partner of
dy ∝ tan φ_inc) contradicted the premise of the 2026-10-22 end-chip handout (pooled RY "endChipTilt", deleted). A rerun
with the corrected tool, also without chip RZ in the rigid set, is pending.

## 2026-10-21 … 10-24: OB chip-row offset, HB scissor, IB end chips

Summaries only; the full story is in `doc/handovers/handover_OB_22_10_26.md` and
`doc/handovers/ib-endchip-rphi-incidence.md`.

- **v5…v10 (10-21/22):** ITS_RB_v5bow…v10. v10 (OB stave RZ free + per-layer `globalModes`) changed nothing: the extra
  `globalModes` lay inside the `rayModes` span. OB dy-vs-φ patterns are a rφ offset between the two chip rows of an OB
  module. Optional weighting of the automatic constraints by statistics (`useStatWeightedMeanConstraints`,
  `ProcessingStats`, 10-22).
- **v11 (10-22):** OB chip TX freed removes most of the chip-row offset (Lr3–6 rows). Open: chip TZ, TX/module-RY
  degeneracy, B0-only n = 2.
- **OB half-barrel Y scissor (10-23/24):** v12 (transported HB pin replacing the NoCross pins, module RY fixed, chip TZ
  free) was worse: the transported pin must complement the NoCross TX/TY pins, never replace them; OB modules move along
  the ray, so module TY acquires a z pattern with TZ ∝ TY·z/r. v14 pins the translations at the half-barrel with
  `pinChildrenMean` on top of the NoCross pins; the scissor escaped through the module RY cross-term (mean module TX paid
  by a fan of module RY ∝ z; half-stave pin's contributor walk stops at the staves, `collectMeanContributors`). v15 fixes
  module RY (kills the scissor at +641k χ², builds a per-layer OB bow the data want).
- **IB end chips (10-24/25):** a rigid shift + RX tilt of the end chips (Lr0 chip 0, z>0 end). v13/v13b (chip RX freed) did
  not move them: pede's `globalModes` equations, written over all of ITS including the IB, forbade the corresponding IB
  deformation. `ITS_RB_IBv0` restricts `globalModes` to the OB; IBv1…IBv3, v16 (= v15 OB + IBv3) follow. IB chip bowing
  shows up as dz V-shapes vs Z.

## 2026-10-25 … 10-26: group modes, sensor-projected coefficients, regularisation, v17/v18

- `groupModes` clause (per-group mode equations over whole subtrees at any level, named groups) and the `globalShift`
  mode: a bow that cancels between the parts of a flagged set escapes a single global equation (in v15/v16 the OB bowed
  radially by +12/+13 µm on Lr3/Lr4 and −15…−17 µm on Lr6 at the stave ends).
  Details: the OB scissor escaped a stave-level `pinChildrenMean` through module TX + RY, `globalShift` per half-barrel
  closes that route; v17 drops `bowing` from the global list (per-layer group equations would tile it) and uses
  `groupModes` bowing on the staves of each layer (staves free in TX,TY,TZ,RY), with `groupModesRelToRoots` keeping the
  data-measured stave radial offsets free.
- Mode-equation coefficients are by default projected on the displacement a DOF gives to the **sensor box corners of its
  subtree** and are root-relative (`AlignParams.modeSensorProjection`); the older "mode evaluated at the volume's origin"
  weights levels by their number of volumes rather than by the displacement they cause, and is blind to a volume whose
  origin sits at z = 0 (an ITS stave) while its sensors are not. Check: `globalShift_y` over Lr3 HB0 gives coefficient
  sums −7.69, −7.69, −7.71 at stave, half-stave and module level with the sensor projection (level-independent), against
  −7.69, −15.4, −61.6 with the origin evaluation; the latter reproduces `mp2con.txt` files written before the option existed.
- `Regularize` output option (10-25): transfer to each alignable parent the rigid motion robustly fitted to the
  displacements of its children's corners (see `Regularization.md`), motivated by the Layer-3 HB0 staves, all with
  RZ ≈ −22 mrad (a "turbine" pattern) in solA2IBv3.
  In solA2IBv3 transferring the median RZ to the half-barrel as a rotation about the beam axis forced lever-arm
  compensations ω·r ≈ 22 mrad × 19.6 cm onto the staves (±4.7 mm instead of ±1.8 mm); the bounding-box corner fit
  leaves the −21 mrad RZ in the staves and their translations drop to 0.1–0.4 mm (mean child-corner displacement
  179 → 68 µm with all levels selected).
- v17 (module RY free + `groupModes globalShift` per HB + per-layer bowing), v17N/v18 (projected coefficients): the
  scissor is killed; v18 preferred.

## 2026-10-26 … 10-31: multi-track PV modes, cosmics, new GBL, IB/OB kink

- Mode 4 for multi-track PV constraints (10-26); cosmic tracks read and added (10-27/29); misalignment tests
  (`misAlgJson`, 10-29); debug tracks (`debugTracksFraction`, 10-29); a workaround for a GBL `scatDataIndex` bug (10-29)
  was removed again once the fixed GBL was adopted (10-30/31).
- The ITS3 closure test formerly driven by `doMisalignmentLeg/RB` is replaced by `misAlgJson` with `--output MisRes` /
  `MilleRes` (same linearisation-level approximation for the Legendre mode), 10-29.
- v18 OB rdphi cells built an IB/OB n = 2 kink; v19 (IB `tangentialHarmonic` gauge) removes it at the price of pede χ² +51 M
  on the same data; v19s/v19rmc test the DOF set and combined rdphi cells (`rayModes` as an array of equation sets,
  `combineRanges`): OB-only dz,dr cells (rMin 10) plus combined rdphi cells over IB + OB (`rEdges [10]`) with IB
  `tangentialHarmonic` modes pin IB and OB separately against the relative-rotation (n = 2 kink) mode.
- v19 pp B0 checkResid shows an IB/OB kink (dPar2 n = 1, 2 up to ~230 µrad at |tgλ| > 1.2, sub-tracks crossing at
  R ≈ 12 cm). A linear response model of dPar2 to the hit rφ shifts, fitted on V10n/v18/v19/v19rmc, gives the OB
  ray rotation about 12 cm that removes it; applied to v19 as OB chip TX (`cure012`, 10-31) it brings pp dPar2 to
  ≤ 34 µrad but worsens OO B+ by the same amount, so the pp/OO kink difference is not geometric. Tools:
  `macro/kinkMap.C`, `macro/kinkResponse.C`, `macro/makeKinkCure.C`, `doc/resToCCDB.sh`
  (see `doc/handovers/kink_cure.md`). pede with the labels fixed at the cure costs Δχ² = +0.73 M (0.05 per record, equal
  in B0, B+, B−): the Mille data mildly prefer v19 independently of the field, the opposite of the pp checkResid. Removing the PV point
  from the checkResid IB sub-track (`checkresid.addPVAsCluster=false`, 11-01) changes the kink by ≤ 3 µrad for V10n,
  v19 and cure012: the kink is in the clusters, not set by the vertex. On pp B+ 562837, next to B0 562862 (11-03), with
  dPar2 constrained to the common IB/OB q/pT (`kinkMap.C` `curvCorr`, replacing the dPar4 regression), v19 shows no
  B0-like kink and the cure makes it worse: the B0 kink is B = 0 specific, not geometric. The pede Δχ² of the cure is
  ~150× the χ² the kink carries, so it does not test the kink. With more B+ statistics and pp B− 564644 (11-04), v19 is flat in B+ and
  nearly in B−, and cure − v19 is the same shift in all fields: the kink is a field-off (even in B) effect, the same for
  V10n and v19. The v18 kink that led to the v19 IB gauge was judged on B0, so v18 is to be rechecked on B±. v18 on B0/B+/B− (11-05): large B0
  kink, B± as flat as v19. The dPar2 response to the layer shifts differs by field (B0: OB direction from Lr3 − Lr4 only,
  pT fixed; B±: Lr5/Lr6 enter and the common q/pT absorbs curvature-like OB patterns), and folding the v18 − v19 shifts
  with each response reproduces the measured B0 and B± changes; dPar4 at B± flips sign with the field. So the B0 kink can
  be a geometric OB pattern that is not linear in r; the 11-04 "field-off only, non-geometric" reading is withdrawn.
- Cosmics: the ITS road of `TrackCosmicsExtended` is a candidate list (without a TOF time the window is the TPC one, so hits
  from many ROFs; `leg` is the half of the cosmic, not the half-barrel; ~95 % of IB and 5–20 % of OB road hits are off the
  track). `AlignParams.cosmiGoldOnly` (default true) keeps only the hits on the straight line of a golden cosmic
  (`CosmicsITSClean.h`, from the 569489 B=0 study: ~18 gold cosmics per CTF file with TOF required), 11-03.
- Unbiased (leave-one-out) residuals `dyUnb`, `dzUnb` with errors, and the cluster `sigYZ`, added to the `gbltrack` debug points
  (`DbgPoint` v2), 11-05. checkResid needs collision tracks with a PV; B = 0 cosmics crossing both half-barrels at non-radial
  angles test the ray-wise OB patterns (the IB/OB kink candidates) that radial collision tracks are blind to, and with few IB
  hits only one-out residuals are usable for them.
- pede on the B = 0 cosmics alone with the labels fixed at each solution (`/data0/devAlg/tests/cosmEval`, 11-06): all Mille
  solutions gain ~130k (4 %) over v10n, so the cosmics confirm them independently. v16 gains least. cure012 beats v19 by 14k (cut) / 32k (no cut),
  ~0.05–0.12 per record, against +0.05 per record for collisions: the cosmics favour the ray-wise OB rotation that collisions
  disfavour (`doc/handovers/kink_cure.md` 6f).
- Combined pede fit B0 + B+ + B− + cosmics in the v18 configuration (`solC1v18`, 11-06). The cosmics gain 127k over v18,
  mostly via stave-end OB chips outside the collision acceptance and a 4 µm IB + mean-vertex shift. Collision checkResid is
  unchanged to 0.01 µm, and the B0 kink stays, yet its pede χ² was ≈ 52 M above v18. Fixed-label evaluations with one pede
  (05.01.07, one file per field) gave v18 = v19 within 0.004 per record. Cause: the outlier down-weighting fix had been
  reverted on epn224, so v19, v19rmc, the pedeEval and C1v18 ran unfixed, while v18 (epn000) had it (χ² ratio 0.894).
  Refits with 05.01.07 (11-07): v19 is 155k *below* v18, so the IB gauge costs nothing. C1v18 is again ≈ v18 and the B0 kink
  stays; its cosmics gain 123k over v18 at a collision cost of 34k (−89.5k combined), so collisions and cosmics are
  compatible at weight 1 (`kink_cure.md` 6g, 6h). The combined fit in the v19 configuration (C1v19, 11-07) is 158k below
  C1v18 on the same records: its cosmics gain 151k over v19 at a collision cost of 54k (−248k combined vs v18). It keeps the v19 OB z-odd c1·z pattern and the end-of-acceptance B0 kink, while the cosmics
  halve the v19 OB s1·z term (`kink_cure.md` 6i). With the cosmics weighted ×3 (C1v19w3, 11-07) the s1·z term vanishes and
  c1·z drops by 15 %, but the B0 kink is unchanged (within 1.3σ), so the s1·z term is not its cause; the OB chip-level
  scatter grows (Lr5/Lr6 dz rms vs v18 23 µm instead of 13). Relative to C1v19 it gains 89k on cosmics but loses 148k
  on collisions (+59k combined at weight 1): the samples pull apart beyond weight 1, and C1v19 stays preferred.
- OB dz-vs-tgλ structures (11-07): the standard checkResid Lr5/Lr6 are flat (|tgλ|<1.4) while the IOB (IB+PV extrapolation) ones
  reach 100–250 µm at Lr5/Lr6, growing with R; Lr0–2 flat to ±4 µm. First guess was an IB radial error of ~5 µm. A check on the
  `res` tree (same day) found that in B≠0 most of the IOB signal comes from the estimator, not from geometry. The IB(+PV) sub-track
  has no curvature information (median σ(q/pT)/|q/pT| = 37–60%), and z extrapolated over 35 cm is quadratic in q/pT. Its mean is
  therefore biased by −tgλ·σ_C²·R³/24. At fixed pT, the Lr6 slope is −230 µm/unit tgλ for |C_IB/C_OB−1|>0.3 and −50 for <0.3.
  B0 (straight lines) is free of this bias: Lr3..Lr6 = −9/−13/−25/−35 µm/unit tgλ. Lr3/Lr4 follow from IB−OB Δtgλ ≈ +6–7·10⁻⁵·tgλ,
  which is field-independent and so geometric. Details are in the handover's "Tree check" section
  (`doc/obsoleteNotes/handover_OB_dz_tgl_IOB_07_11_26.md`, merged into `kink_cure.md` 6j on 11-09).
  - 11-09: the OB chip DOFs do not remove this IB−OB slope mismatch, they increase it (B0 par3 87 → 108·10⁻⁶ for A2v19 → A2v20,
    121 with cosmics ×3, 113 in C1v21w3; same ordering in B±). The B0 Lr5/Lr6 IOB excess over the Lr3–Lr4 line shrinks by
    ~5 µm/unit tgλ only in the fits with cosmics, i.e. through the projective component (`kink_cure.md` 6j).
  - 11-09, per-chip residual maps (`macro/drawChipResidualMaps.C`, C1v21w3): the chip-to-chip structure at large |Z| is
    mostly statistical (run-to-run correlation 0.1–0.35, χ²/n ≈ 2 at all |Z|). The significant part is stave-common and is
    the same with and without cosmics (A2v20 = C1v20w3). Candidate: the v17+ sensor-projected gauge (per-layer `bowing`,
    OB `radialPolar_r0t0/r0t2`, `longitudinalPolar_r0t1`) pins exactly the stave-common odd-in-tgλ dz shapes
    (`kink_cure.md` 6k).
  - 11-09, reduced refits on top of C1v21w3 (`/data0/devAlg/tests/reducedFit`, all fixed except a chosen subset) do not
    confirm the gauge candidate. Releasing those gauge equations with the OB modules free (S1, S2) changes the sigmoid
    little. Also freeing the end-module chips (S3) moves them by tens of µm and trades the odd part of the sigmoid for an even
    one, while the GBL kinks are only redistributed. The mean kink pulls are ≤ 0.06, so the sigmoid costs pede only a few k χ²
    per 2000 TFs: the fit has almost no lever on it (`kink_cure.md` 6k).
  - 11-10, the same refits on new collision data, with B0 tracks at an assumed 1.5 GeV and B± at pT > 2, with and without
    the B0 collisions (S*R, S*RR), give the same sigmoid and the same end-chip moves (correlation 0.85–0.94). Reducing
    the GBL scattering freedom does not give the sigmoid a χ² weight (`kink_cure.md` 6k).
  - 11-10, `AlignParams.tgLCut` was added to write edge-only Mille data (|tgλ| > 1.2, ≥ 4 ITS clusters, B± pT > 1.5).
    Boosting it with the bulk free (S3R) distorts the bulk. With only the OB stave ends free (S4R) it halves to removes the
    +tgλ edge sigmoid and most of the IOB edge kink. The central IOB dtgλ wave, the B+ vertex dz-vs-tgλ slope and the
    Lr0 |z| > 12 cm dY rise stay (`kink_cure.md` 6k).
  - 11-11, the IB z-stretch on C1v21w3 (ε = −1.2·10⁻⁴) removes the central IOB dtgλ wave but builds a dz(R = 12 cm)
    slope vs Z that the unstretched geometry does not have. The wave pivots near R ≈ 12 cm, so it is not an IB z-scale
    about z = 0 (`kink_cure.md` 6k).
  - 11-12, freeing only the IB end chips (chip0/8 TX, RY, RZ) on S4RBis with IB edge data halves the Lr0 dY rise at |z| > 12 cm
    without degrading other residuals. The fitted rφ shift ramps over the staves of each half barrel, instead of being
    uniform (`kink_cure.md` 6k).
  - 11-12, freeing chips 1/7 too (IBEdge2) gains nothing. What remains at the Lr0 stave ends scales with the rφ incidence:
    it is an out-of-plane displacement of the outer ~1 cm of the Lr0 end chips (≈ 7–14 µm), which TX/RY/RZ cannot
    absorb (`kink_cure.md` 6k).
  - 11-12, fixing the Lr5/6 Module1/5 (S5R) does not remove the Lr5 dz structure that S4R/S4RBis build at z ≈ +41…+47 cm,
    and costs +464k. The structure sits at tgλ ≈ 1.3–1.4, where tracks cross the inner part of the Lr3/4 end modules, so
    the edge correction itself, through the tilt of those modules, is the likely cause (`kink_cure.md` 6k).
  - 11-12, restoring the 30 µm module translation priors (S6R) leaves the 0.1–0.3 mm slides of the OB end modules along
    the IP ray almost unchanged (−10%) and all residuals within 0.5 µm. The slides are held by the data (≈ 10 µm per DOF),
    not by prior softness. Pede pays +34k, almost all of it in the prior term (`kink_cure.md` 6k).
  - 11-12, fixing the end-module translations at C1v21w3 (S7R) costs +1.17M, almost all in the track data, and loses
    the last part of the edge cure (Lr5 dz at tgλ = 1.9: 4 → 16 µm in B0). The 1.3–1.4 band stays, so the slides are
    needed and do not cause it. The cosmics agree: +152k against S4RBis (`kink_cure.md` 6k).
- v10n OB chip out-of-plane (TY) and longitudinal (TZ) corrections (11-07). They are carried unchanged into all our
  solutions, since OB chips are free only in TX, and many look unphysical (|TY|, |TZ| > 100 µm). All chips of a module share
  the same rotation, so only the chip translations scatter. Net chip displacements from the ideal geometry confirm that the
  scatter is real: it is not one level compensating another.
  - Chips without data (2.7k OB chips with < 100 tracks in the C1v19 pp B0 + B+ + B− checkResid samples, mostly at the stave
    ends) carry most of the outliers: their net z deviates from the module median by 210 µm (MAD).
  - For chips with data the within-module scatter is 33–36 µm rms in modules without such chips. It is 75 µm (TY) and
    126 µm (TZ) in modules with some.
  - Per-chip residuals do not correlate with it. The per-chip ⟨dz⟩ has an rms of 11–16 µm (3.5 µm error), and the
    regression on the v10n TZ deviation is 0.00 ± 0.02, also for the 1600 live chips with |TZ| > 100 µm. For TY (from the
    dz-vs-tgλ slope) the regression is ≤ 0.1, up to 0.3 in Lr3 HS0.
  - So the data support the v10n chip offsets of live chips at the 10–15 µm level, except possibly for coherent patterns at
    the acceptance edge, which the interpolation does not see. The values of chips without data cannot be corrected by any
    fit.
- v20 (11-08, `doc/ITS_RB_v20.json`: v19 plus OB chip TY, TZ, with priors of 30 µm and a stat-weighted module mean) was fitted
  with cosmics ×3 (C1v20w3). It is 4.84 M below C1v19w3 and rejects 32k fewer records. The OB dz-vs-tgλ sigmoid shrinks
  and sets in at larger |tgλ|; the Lr6 one, which was leakage of Lr5, disappears.
  - The chips move mainly along the direction projective from the IP, ΔZ = tgλ·ΔR: 30–45 µm rms for live chips and up to
    120 µm for chips without data. Collision tracks barely constrain this direction.
  - The combination the collision tracks measure, ΔZ − tgλ·ΔR, is mostly stave-specific (stave-to-stave rms 25–40 µm at
    |tgλ| > 1.4). That part of the sigmoid is therefore geometry, not a universal incidence bias.
  - The part common to all staves is odd in tgλ, ≈ 20–35 µm at |tgλ| ≈ 1.8, with the same sign at Lr5 and Lr6. It could
    still be a cluster-z incidence bias, since chip TY + TZ absorb any function of tgλ locally.
  - The B0 n = 2 kink is unchanged. The OB s1·z term reappears with the sign opposite to v19 (+26 µm at Lr5/Lr6 vs v18).
  - cosmEval: the cosmics gain 1.04 M (35 %, nocut −2.41 M) and the collisions 1.73 M vs C1v19w3 (−2.77 M combined at weight
    1). Both samples now agree, unlike in C1v19w3, so the chip moves (including the projective component) are supported by
    tracks not pointing to the IP. An independent check (v20 without cosmics) is pending.
  - (`kink_cure.md` 6i.)
  - The v20 steering files (C1v20w3 and A2v20) were regenerated without `AlignParams.statInFile`, so the `weightedMean`
    rule fell back to the equal split. A2v20 (v20, no cosmics): −1.78 M vs A2v19. It reproduces C1v20w3 on the combination
    the collisions measure, c = ΔZ − t0·ΔR (corr ≈ 0.9). Its projective component is half as large: that component in
    C1v20w3 is driven by the cosmics.
  - cosmEval of A2v20: cosmics −345k vs A2v19 without being in the fit, an independent confirmation of c. C1v20w3 gains a
    further 932k on the cosmics for +255k on the collisions (v19: −241k for +202k), so the cosmics fix the projective
    component at a small collision cost.
- v21 = v20 + OB chip RX/RY/RZ (`solC1v21w3`, cosmics ×3): −3.28 M vs C1v20w3, of which −2.92 M collisions and −0.12 M
  cosmics at weight 1. The gain shows in the checkResid widths and φ-uniformity of the OB dz (core σ −1…−3 µm, spread of
  the dz means vs φ halved), not in the layer-averaged means vs tgλ: the OB sigmoid is unchanged (its edges at Lr3/Lr6
  slightly worse). (`kink_cure.md` 6i.)

## 2026-11-11 … 11-13: material before the last point, energy-loss sign of cosmics, jacobian step, charge-less clusters

Context: PbPb cosmics with the field on (run 568041), TPC with TOF points at both ends, unbiased TOF residuals from the
`gbltrack` debug tree. The TOF y pulls were large (⟨pull²⟩ 5.6), twice as large on the bottom leg as on the top one;
cosmics with ITS points had ⟨pull²⟩ ≈ 1.4 with the same TOF hits.

- **Material before the last point dropped.** The last point of a trajectory carries no scatterer (`hasScatterer`
  requires `ip < np - 1`), and the step to it was split only above `splitMSThreshold` (5 mrad, i.e. low momenta), so
  the material on the way to it was lost: for a TOF point the TPC outer wall, the TRD and the TOF. 97 % of the cosmics
  had no scatterer between the last TPC point and the bottom TOF (the track error there 0.13 cm instead of 0.16 cm).
  The step to the last valid point is now always split when it crosses material, the scattering of the whole step on
  the mid point (in the pre-step frame, or in the frame of the last point if the mid point cannot be reached in the
  first); if neither works the step stays unsplit and is counted (`last-step MS lost` in the GBL summary). This also
  acts on collision tracks ending on TOF or TRD.
- **Energy-loss sign of cosmics.** The Propagator takes the sign of the energy-loss correction from the direction of x
  (outward = loss). The top leg of a cosmic is run in the upward parametrisation, so in the KF of the legs it lost
  energy on the way up and in the GBL reference it regained it on the way down. `Track::fitTrack` takes the sign
  (`eLossSign`); the KF of the legs loses energy outward along the bottom leg and regains it back inward and along the
  top leg; every GBL step of a cosmic (time order, top first) is a loss. Collision tracks are unchanged. Bottom-leg TOF
  ⟨pull²⟩ y at 2–5 GeV 5.7 → 5.0, GBL χ²/ndf 1.080 → 1.070.
- **Transport jacobian.** The numerical derivative starts with a perturbation of ±0.05 in snp; at |snp| ≈ 0.85 the
  perturbed track leaves maxSnp after the rotation to the next sector and the whole track was lost. The first step is
  now shrunk (down to `minStep`) before giving up, and a failure at a later extrapolation keeps the best estimate so
  far. Cosmics: propagation failures of the GBL stage 145 → 12 per 60 k TFs.
- **TPC superclusters of clusters without the charge.** Clusters read with qTot = 0 (ClusterNative written before the
  qTot → qTotPacked rename and read by a newer O2) gave a NaN supercluster position and no fit at all; they now get
  equal weights.
- Together, on one cosmics piece (TPC-TOF, ideal map, no selection): GBL fits 10 734 → 11 754, propagation failures
  2106 → 992 (the reference with the charge fix only).

## 2026-11-13: TRD pad tilt, cosmics across the central electrode, fixed DOFs in MilleRes, single-point bottom leg

Found in a code review of `ea32487`; each one checked on the laptop.

- **TRD pad tilt.** `makeTrackletPoint` set the tilt correction to 0 when the track is within one pad length of the
  pad-row centre, i.e. exactly where O2 (`GPUTRDTracker`, `TrackInterpolation`) applies it, and applied it only farther
  away. The condition is now O2's: within one pad length of the row centre and with a track z error σz² < l²/12 (the
  latter excludes 1 of 11 145 tracklets on a cosmics piece with the ideal map). The tilt-corrected y also depends on the z of the pad row, while the z row of the
  derivatives is zeroed for these points (`zFromTrack`): the y row now gets `-tilt * (z row)` where the correction is
  applied (`FrameInfoExt::tiltDz`). Checked by finite differences: the raw tracklet moved by a ±10 µm / ±10 µrad
  correction of each DOF of the chamber and of the supermodule (as `MilleRes` applies it), the measurement redone; all
  12 columns agree with the analytic ones at the track point (100 % of 1500 points within 1 %, TOF strip and
  supermodule as well). Cosmics with TRD on a calibrated TPC (4 pieces, 162 cosmics with TRD): TRD y ⟨pull²⟩ per layer
  1.35 / 1.22 / 0.88 / 1.35 / 1.62 / 1.13 → 1.23 / 1.15 / 0.82 / 1.29 / 1.60 / 1.13.
- **Cosmics crossing the central electrode.** With the ideal TPC map, the clusters of one pad row on both sides of the
  central electrode (sectors s and s + 18) have the same alpha and X. The zero step between them has an identity
  jacobian, rejected as "Near jacobian idendity": the whole cosmic was dropped (≈ 85 per cosmics piece, exactly the
  cosmics that link A and C). Such a frame is now a 2nd measurement on the previous GBL point (`DbgPoint::SamePoint`).
  When such frames end the track, the last GBL point is the first of them: it gets the split of the last step and no
  scatterer (GBL ignores a thin scatterer on the last point). One piece (TPC-TRD-TOF, ideal map): identity errors 83 → 0,
  fits 11 523 → 11 617; 105 same-point measurements, none on the last frame.
- **Fixed DOFs in MilleRes.** pede reports a fixed parameter at its preset value (`startFromMilleRes`), and
  `readMillepedeResults` keeps it, but the writers required a free DOF: a restart with the parents fixed exported the
  children only. Fixed DOFs with a non-zero value now go to `result.json`, the CCDB objects and the weak-mode
  diagnostics; nothing changes without such values. Test (TRD supermodules fixed at TX = 0.05 cm, chambers free):
  `result.json` 521 → 539 entries, the supermodules moved by 0.05 cm in the CCDB object, chambers unchanged.
- **Single-point bottom leg of a cosmic.** The step between the legs is not split, and with one bottom-leg point it is
  also the last step, whose point carries no scatterer: the interior material was lost (and not counted). It now goes
  on the last top-leg point, combined with the scatterer there. Those cosmics (0.6 %): χ²/ndf 5.46 → 4.33, ⟨pull²⟩ of
  the lone point 116 → 60.
- All four together on that piece (w01_0, TPC-TRD-TOF, ideal map): 11 523 → 11 611 fits, identity errors 83 → 0.

## 2026-11-14: TRD tracklet z is the pad-row centre

- `solB0A1` (B0 + B± data, TRD supermodules free, no constraints): `millepede_labeled.res` gave TY of the TRD supermodules
  of −1070…+1088 cm (e.g. sm00 −438, sm17 −839), while TOF came out small. TY of a supermodule is a shift along the beam.
  Cause: the TRD point took its z from the track (`zFromTrack`, since `0ee2cc9`), so the z residual was identically
  zero, its derivatives were zeroed, and a beam-axis displacement was seen only through the pad tilt (dy/dz = tan 2–3°).
  O2 (`TRDGlobalTrackingSpec.cxx`, Kalman update) uses the calibrated tracklet z, i.e. the pad-row centre, with σz² = l²/12.
  The point now has z = row centre (`traXYZ[2]`, which carries the chamber alignment) + `TRDNonRCCorrDzDtgl`·tgl, the z
  derivatives are kept, the tilt term on y is unchanged. `zFromTrack` removed; the KF refit also no longer shrinks σz with a
  no-information z update, and the ndf counts 2 coordinates per TRD point. Not yet rerun on `solB0A1`.
- Cross-check of the tracklet point with the O2 Kalman update (`TRDGlobalTrackingSpec.cxx`): the z of non-row-crossing
  tracklets is corrected with `getZCorrCoeffNRC()` (1.4 per unit tgl) in O2, our default was 0 — harmless with z from the
  track, a z bias once z is the row centre; `AlignParams.TRDNonRCCorrDzDtgl` now defaults to the O2 value (<-990). O2 (default
  `useAngularPull=3`) also increases σy with the angular pull and shifts y by `-pull*getCorrYDy(snp)`; both are now
  applied (`AlignParams.TRDUseAngularPull`, 0 = off). Still different from O2: chamber occupancy, pile-up y correction and error.
- The TRD tracklet error parameters (`trkltResRPhiIdeal`, `trkltResVsTanPhiMisalign`, pile-up ranges) come from the
  `GPU_rec_trd.*` configurable params, not from `AlignParams`: to reproduce the errors of the reconstruction, pass the same
  `GPU_rec_trd.*` keys to the workflow.
- `weakModes.pdf` of `solB0A1` showed TOF amplitudes of up to 4e7 µm (`bowing`, `radialPolar_r*t2`, `twist`), while the TOF
  corrections are a few mm and `rmsDisplacement` is 3.3 mm. TOF is a thin cylinder (r = 370–383 cm), so `bowing`,
  `radialPolar_r0t2`, `radialPolar_r1t2` (and `twist`/`tangentialPolar_r0t1`, `radialPolar_r0t1`/`radialPolar_r1t1`) are
  the same field (correlation ±1.0000, conditioning 6e-17): the fit with a 1e-8 SVD cut returned huge amplitudes cancelling
  each other. The SVD cut is now 1e-3 and the diagonal of the resolution matrix is reported per mode (`resolution`,
  `determined`); modes with resolution < 0.5 are flagged in the log and drawn in grey by `drawWeakModes.C`.


## 2026-11-19: TRD/TOF residuals of `solB0A1`, chamber DOFs for the next iteration

- `solB0A1` (sectors only, ITS fixed, B0+B±+cosmics) vs its start `best/IBEdge`, from the unbiased residuals of `trdTofResid.C`:
  TRD sector z offsets of 1–2 cm (sm03–06 +13…+17 mm, sm13–14 −11…−16 mm, all layers and fields) and the common −0.9 mm dy
  offset are gone; the chamber-in-sector dy-offset spread went 0.4–0.6 → 0.23 mm (stat. error 0.02–0.04 mm, field-consistent,
  corr 0.6–0.8), the per-sector dy-vs-tanφ slopes of sm12–16 (+1…+6 mm) went to ≤ 1 mm. Solution: TRD TX +6.6 mm and RY
  +2.4 mrad in every sector (a 2 mrad rotation about the beam), TY −9.1 mm·sinφ (a tilt), TZ −8.9 mm·sinφ − 4.6 mm·cos2φ with
  no matching rφ term (bottom sectors sm12–15 out by 9–15 mm).
- Left: L4/L5 dy-vs-tanφ slope −0.5…−1.5 mm at B± (radial-like, ≈ −0.5 mm at B0); chamber dz-centre spread 2.5–3.4 mm with
  poor field consistency; chamber dy offsets correlate with ±tan 2°·dz (alternating per layer, corr 0.56 at the start), i.e.
  chamber z enters y through the pad tilt. The apparent in-chamber dz-vs-tgl slope (−1…−2 cm/tgl, also in v10n) is an
  edge-bin effect: the inner tgl bins are flat. TOF: dy vs snp is a step at snp = 0 of ±2.5–3 mm at B± (negative slope at B0),
  end-slice dz offsets of ±1 cm at B± only — non-geometric; TOF sector RY came out +4.3 mrad on average with TX only +1 mm.
- pede rejected 28.5 % of the records (2.6 % in the IBEdge ITS fit, up to 45 % in B0 files), jumping from 28k to 1.46M when
  the cut factor reached 1. Unbiased pulls of accepted tracks: ITS 7–12 % beyond 3σ, TOF y 9–10 % at B± (core 1.1–1.3), TRD y
  core 0.4–0.6. The mean vertex was fixed at the reference (0,0) while the beam is at (+0.34, −1.82) mm.
- Cause of the rejections: the composed (multi-track PV) records. In the `Mille_IBEdgeS4R` debug trees (mode 2) 98 % (B0) / 89 %
  (B+) of the composed records exceed the 3σ χ² quantile (χ²/ndf ≈ 7, rising with ndf), the single tracks 38 % / 12 %; the
  accepted pede records have ⟨ndf⟩ 26–32, i.e. mostly single tracks. Their excess χ² sits in the kinks at ITS layer 0. The
  recorded prior is the CCDB mean vertex, (−0.54, −0.05) mm with σ ≈ 50 µm, while in the IBEdge ITS frame the beam (and the
  refitted vertices) is at ≈ (−0.20, −1.87) mm, which is the IBEdge pede correction (+0.34, −1.82) mm. With the mean-vertex
  labels fixed in `solB0A1` the 37σ pull remains in every composed record. The tail fraction per TRD chamber is uniform and
  uncorrelated with the chamber residuals. The ProcessData check runs (mode 0, per-track prior, no `mvbias`) carry the same pull.
- `solB1A1` (`ITSEnvTRDTOFch_RB.json`, mean vertex free): MV fitted to (+0.333, −1.813) mm per slot as in IBEdge; rejections
  22.3 % (28.5 %): B± 18–20 → 10.5–11.3 %, B0 42–44 → 35–37 %, cosmics 3.7 % in both. The MV pull was the charge-independent
  ±2σ kink pattern at ITS layer 0 of the composed tracks (sign change at φ ≈ −75°/+105°, i.e. a vertex pulled along φ ≈ 100°).
  The rest are single tracks (B0 36 %, B+ 12.5 % with χ²/ndf > 2.5 in the Mille debug trees) with kink pulls² of 10–19 at
  every scatterer (1.4–2.2 for the good ones): multiple scattering underestimated, at B0 that of the tracks softer than
  `meanPtB0` = 0.7 GeV; a composed B0 record with one such track is rejected as a whole.
- `doc/ITSEnvTRDTOFch_RB.json`: mean vertex free; TRD chambers free in TX, TY, TZ, RZ with priors 0 ± 1 mm / 1 mrad; chamber
  RX/RY kept fixed, since with a free sector the transported mean constraint gives them lever arms of ±240 cm (RX → sector TZ)
  and ±31 cm (RY → sector TX), through which weakly observed tilts could satisfy the constraint; TOF sector RX/RY with
  1 mrad priors, TOF strips fixed.

## 2026-11-19: `milleReplay`, and what `solB1A1` did with the TRD chambers

- Regenerating the `gbltrack` debug trees to validate a solution takes hours. `o2-dev-mille-replay` (`src/milleReplay.cxx`) replays the Mille
  binaries with the solution applied, with pede's local fit (Huber/Cauchy down-weighting) and rejections, and accumulates
  per parameter the profiled gradient and curvature, per input file. All the TRD/TOF binaries take about 3 min; the
  rejection fractions match pede's within 0.5 % (B0 file 36.1 vs 36.5 %, cosmics 3.7 vs 3.7 %).
- `ITSEnvTRDTOFch_RB.json` was written assuming the TRD chamber DOFs are in the TGeo chamber frame (`TX` rφ, `TZ` radial).
  They are in the sensor frame of `SensorTRD` (`TX` radial, `TY` rφ, `TZ` z): the freed `RZ` was the tilt about the beam,
  which enters the supermodule TX constraint with the layer lever arm, while the in-plane `RX` stayed fixed.
- `solB1A1` chambers: radial `TX` RMS 1.9 mm (up to 7.6 mm), rφ `TY` 1.4 mm, z `TZ` 2.2 mm (up to 9.8 mm), `RZ` 4.5 mrad
  (up to 19 mrad), with layer-linear means (`TX` −0.7 → +1.4 mm, `TY` −0.6 → +1.7 mm, `RZ` −3.1 → +4.4 mrad from L0 to L5);
  the 1 mm / 1 mrad priors do not hold well-measured chambers. Replayed per dataset, B0 and cosmics ask to undo the layer
  patterns (L0 radial +1.4 / +2.5 mm, L0 `RZ` +3.5 / +3.3 mrad) and B± for more (−0.2…−1.2 mm, −2.6…−2.9 mrad); per chamber,
  B0/cosmics would undo 20–32 % of the radial and `RZ` values, B± add 8–15 %; the rφ `TY` values are accepted by all four
  datasets (regression slope 0.00 ± 0.02). The chamber z `TZ` is rejected by B− only. The fixed in-plane `RX` would gain
  1.4–3k χ² per dataset, with shifts correlated between B0, B+ and cosmics (+0.58, +0.35).
- Common modes in tension between field-free and B≠0 data at the `solB1A1` optimum: TRD sector rφ (B0 +0.6 mm, B± −0.5/−0.9,
  Δχ² 6k/11k), TOF sector rφ (B0 −1.1, B± +1.0/+1.5, Δχ² 11k/15k), TRD sector radial (B− −2.5 mm, cosmics +4.5 mm, Δχ² 6.7k/7.2k),
  TOF strip normal (fixed; B0 −12.8 mm, B± +4…+6). B− chamber-level shifts do not correlate with any other dataset.
- `doc/ITSEnvTRDTOFch2_RB.json`: chambers free only in rφ `TY` and in-plane `RX` (priors 1 mm / 1 mrad), radial and z fixed.
- `macro/drawLocalFrames.C` drew the TGeo physical-node frames, so `doc/local_frames_TRD.pdf` showed the TRD chambers with
  the supermodule axes (X rφ, Z radial), not the sensor frame their DOFs are defined in. It now draws the DOF frame (from
  the `Sensor*` classes) solid and the TGeo frame dashed; the three `doc/local_frames_*.pdf` are regenerated.

## 2026-11-20: the B on/off rφ tension is TOF relative to TRD, and is already in IBEdge

- At `solB1A2` (`ITSEnvTRDTOFch2_RB.json`), single-group replays showed the common TRD and TOF sector rφ shifts requested
  with opposite signs by the B = 0 data (B0 collisions, cosmics) and the B ≠ 0 data, and B+ and B− agreeing (polarity-even).
  The TRD and TOF global rotations are correlated by 0.5–0.8, so a single-group shift mixes them. `o2-dev-mille-replay`
  now writes the profiled curvature matrix between groups (`_groupH.txt`) for joint fits.
- Joint fit of the TRD and TOF global rotations (sector TX = R·δφ, RY = δφ), 4 files per collision dataset, all cosmics:
  - at the IBEdge reference (recording with `applyTRDshift=false`, fitted mean vertex), every dataset asks for the TRD
    pad-shift rotation: +4.9 / +5.1 / +5.3 mm at R = 329.6 cm for B0 / B+ / B−, +6.2 mm for cosmics;
  - TOF relative to that: B0 −0.7 mm, B+ +1.9, B− +1.6, cosmics +0.5 (at R = 385 cm). At `solB1A2`: TRD +0.3 / +0.0 / +0.2 /
    +0.2, TOF −0.8 / +1.4 / +1.1 / −0.3. The tension is a 2.0–2.6 mm polarity-even TOF rφ offset between B on and off,
    present with the same size before the TRD/TOF fit.
- TRD sector TX and RY are well separated (correlation 0.2–0.3, layer lever arm ±32 cm); TOF RY is weakly measured
  (σ 0.3–0.6 mrad per 4 files) but correlated with TX by only 0.0–0.4.
- An ITS rφ mode ∝ r² (a charge-independent track curvature at B ≠ 0) does not explain it. Imposing it shifts the B on
  requests in the ratio (r_TRD/r_TOF)² ≈ 0.74 and the B0 ones by much less. Removing the B on/off difference would need
  +3…+5 µm at ITS L6 for TRD and −28…−40 µm for TOF. The data themselves ask for −11 µm (B0) and −16 µm (cosmics).
- `gbltrack` debug trees of `Mille_IBEdgeS4R` (10 % of the tracks, IBEdge reference). Each collision dataset is a single run
  (B0 562862, B+ 562837, B− 564644), and the tracks have no TPC points. 57 % of the TOF tracks (B+) are ITS+TOF only: at B ≠ 0
  their curvature is fixed by the TOF point itself, so a TOF rφ shift is nearly degenerate with it, and their unbiased TOF
  residual binned in the KF q/pT is self-correlated (a ±6 mm wiggle around q/pT = 0). At B = 0 they measure TOF against ITS
  directly: B0 +0.2…+0.4 mm, cosmics −0.7 mm.
- The TOF measurement is the centre of the main contributing pad (σ 7.5 mm), so the TOF residual of ITS+TRD+TOF tracks is a
  25 mm box smeared by the 3.7 mm (B0) or 4.6 mm (B+) prediction error. At B0 the box is flat (−4.9…+21 mm, centre +8.2).
  At B+ and B− the right edge is the same within ≈ 0.5 mm, but the box carries an extra narrow peak at 0…+5 mm. The medians
  then differ by 3 mm, the 5–95 % centres by 0.3–0.8 mm, and the skew (84−50)−(50−16) is +3 mm. The skew is the same for
  both polarities and at every incidence angle (+2…+4 mm already at |snp| < 0.05, where B0 has +0.1), and it is absent in
  cosmics (run by run ≤ 2 mm). A pad-quantized measurement with a uniformly illuminated pad cannot produce such a peak, so the
  B on/off TOF rφ tension of the least-squares fit is a change of the in-pad residual shape, not a geometry difference. Its
  origin (TOF clusterization or matching at B ≠ 0) is open.
- `dyPre` is not an ITS-only residual: its seed is the inward KF refit over all the points of the track (`refitTracks` →
  `fitTrack`), so TOF and TRD enter it, through the curvature in particular. At TOF it is biased towards zero by the point
  itself (median −0.1 mm for ITS+TRD+TOF at B0, against +8.5 for `dyUnb`) and is 2–3 times wider (68 % half-width
  19–32 mm against 9–10).
- `doc/ITSEnvTRDTOFch3_RB.json`, for the recording with `extraClsErrYTOF=0.5` (TOF σY 0.75 → 0.90 cm). It is
  `ITSEnvTRDTOFch2_RB.json` plus a `curl` equation over the TOF sectors, because per-sector priors cannot hold the common
  TOF rotation: with the full B ≠ 0 statistics its data error is ≈ 25 µm against 0.24 mm from 18 × 1 mm priors. The TOF
  global rotation thus stays at IBEdge, which the B = 0 data accept within ≈ 0.7 mm. TOF sector priors relative to
  IBEdge: rφ 1 mm, z and radial 2 mm, rotations 1 mrad. TRD sectors stay free without priors, to absorb the pad shift.
- Cosmics replayed run by run (123 runs, `solB1A2`, joint TRD/TOF rotation + radial). The TRD and TOF rotations are stable
  (χ²/ndf 82/83 and 123/83). The TRD sector radial is +4.6 mm in every period (LHC25aa +4.9, LHC26aa +4.4), with
  run-to-run scatter beyond the errors (χ²/ndf 336/83, about 1 mm extra), while the single-run collision datasets give
  B0 +1.3, B+ +0.1, B− −2.3 mm. A constant TRD drift-time offset δt moves every tracklet by vd·δt along the drift (radial)
  direction, which a sector radial shift reproduces exactly (4.6 mm ≈ 0.3 µs at vd ≈ 1.5 cm/µs). The radial requests are
  therefore read as the time-zero of each run (cosmic trigger vs collisions), not geometry; the earlier "polarity-odd" TRD
  radial is equally explained by B0, B+ and B− being three different runs. The TOF radial agrees in all datasets (joint
  fit: B0 −0.2 ± 0.2, B+ +0.1 ± 0.2, B− +0.0 ± 0.2, cosmics +0.2 ± 0.04 mm).
- `ITSEnvTRDTOFch3_RB.json`: TRD sector radial (`TZ`) fixed at IBEdge; cosmics not up-weighted.
- Chamber radial (sensor `TX`, fixed in `solB1A2`), replayed with all the collision files and the cosmics split by year:
  - per layer, after removing the common mean: B+ and B− show a layer-linear pattern (L0 −1.0/−1.5, L5 +3.0/+3.6 mm) that
    B0 (L0 +0.6, L5 +0.7, middle layers −0.2…−0.5) and the cosmics (2025 and 2026 alike, +0.4 → −0.4 mm) do not. Being
    polarity-even and absent at B = 0, it is read as a track-model effect at B ≠ 0 (curvature inside TRD: field map or
    energy loss), not geometry;
  - per chamber, after removing the sector and layer means, all dataset pairs correlate. The common part is 0.5–0.7 mm
    rms between collisions and cosmics, and cosmics 2025 vs 2026 (independent samples) correlate by 0.73 (1.0 mm). This is
    a chamber-level radial geometry pattern.
- `ITSEnvTRDTOFch3_RB.json`: chamber radial `TX` freed with a 1 mm prior. The sector radial stays fixed and its chamber
  mean is pinned (`pinChildrenMean` `TZ`), against the time-zero offsets. The per-layer means of L0–L4 are pinned by a
  `groupModes` `radialScaling` per layer, against the B ≠ 0 layer pattern; L5 follows from the sector pins. This needs
  `AlignParams.groupModesRelToRoots=false`, since the group volumes are the chambers themselves.
- Radial sensitivity of the TRD chambers per track (TX/TY information ratio ≈ ⟨tan²ψ⟩): B0 and cosmics 0.010, B± 0.065; the
  chamber error with all the files is 0.22–0.28 mm at B± against 0.43 mm at B0 (3× more tracks). In the debug trees the
  per-layer radial slope d(dy)/d(tanψ) at B± depends on pT. For |q/pT| < 0.7 (rms tanψ 0.12–0.18) the layer deviations go
  from +0.5…+0.8 (L0) to −0.5…−0.9 mm (L5), like the cosmics. For |q/pT| 1.0–1.5 (rms tanψ 0.37) they go from −1.4…−2.0 to
  +3.3…+4.1 mm. The B± layer pattern of the replay is therefore the low-pT tracks, i.e. a momentum-dependent track-model
  effect, while the B± high-pT tracks see the same small layer pattern as the cosmics.
- Summary of the TRD/TOF conclusions so far: [TRDTOFobservations.md](TRDTOFobservations.md).

## 2026-11-20: `solB1A3` (`ITSEnvTRDTOFch3_RB.json`, recording with minPt 1.2, LUT corrections)

- Recording `Mille_IBEdgeS4R/collHPLUT`, `cosmHPLUT`: `corrType=2`, `minPt=1.2` (cosmics 0.6), `meanPtB0=0.45`,
  `meanPtB0Cosm=6`, `extraClsErrY/ZTOF=0.5`, `extraClsErrY/ZTRD=0.05/1.5`. Pede χ²/ndf 0.96.
- Rejections at `solB1A3` (replay, single / PV records): B0 6.5 / 8.7 % (was 22 / 56 %), B± 9 / 24 % (was 5 / 16 %),
  cosmics 12 % (was 3.7–4.6 % with `meanPtB0Cosm=2`).
- Chamber radial: 6.8 % of the chambers beyond 3 mm and 0.9 % beyond 10 mm, up to −38 mm, against the 1 mm prior.
  For these chambers every dataset asks to undo the shift (e.g. sm12/st1/pl4 at −37.9 mm, cosmics +32.5 ± 1.3 mm back;
  sm17/st2/pl5 at −26 mm, all five datasets +19…+33 mm back). The coefficients of the sector `TZ` pins (1/30) and of the
  per-layer `radialScaling` equations (≈ 1) are uniform. With a hard equal-weight mean, every chamber of the equation gets
  the same Lagrange force, and the shift is that force over the chamber's curvature. When the data want a common layer
  shift, the well-measured chambers follow it and the weakly measured ones (dead in collisions, cosmics only) restore the
  mean by moving far the other way. `groupModes` have no statistics weighting.
- The layer radial pattern requested at `solB1A3` (joint fit with the TRD sector TX/RY/RX, layer TY and TOF modes): B0
  +1.2 (L0) → −1.9 mm (L5), cosmics +0.6…+1.6 → −1.7…−2.7, B± −2.1/−2.5 → +3.2/+3.4. B = 0 and B ≠ 0 ask for opposite
  layer patterns, also with `minPt=1.2`. Per chamber (sector and layer means removed), the cosmics of the two years
  correlate (0.75, 1.6 mm rms) and B0 with the 2026 cosmics (0.38), while B± anti-correlate with the cosmics (−0.2…−0.26):
  the B ≠ 0 chamber-radial information is affected too.
- The TRD sector RY (rotation about the beam) went to −3.0 ± 2.8 mrad, compensated by a layer-linear chamber TY
  (L0 −0.88 → L5 +1.21 mm). There is no mean-of-children equation for the sector RX/RY (no matching chamber DOF), so the
  pair is constrained only by the chamber tilt it leaves and by the priors. In the joint fit B0 and cosmics ask for
  +4.8…+7.1 mrad back, B± for −1.2…+1.2.
- Accepted by all datasets: the TRD rotation (pad shift; sector TX mean +6.6 mm; joint requests 0.0…+0.37 mm), the
  sector z offsets (rms 8.7 mm; requests ≤ 0.17 mm), the in-plane chamber RX (≤ 0.16 mrad). TOF: the `curl` holds the
  rotation; requests B0 −0.08, cosmics +0.24/+0.44, B± +2.05/+1.69 mm, the same tension as before. The common TRD radial
  (fixed) is requested as B0 +1.0, B+ −1.4, B− −4.1, cosmics +6.5/+5.9 mm, consistent with time-zero offsets.
- An ITS rφ ∝ r² mode is requested by every dataset with the same sign: d₆ = −13 (B0), −9/−11 (cosmics), −15/−17 µm (B±).
- The ITS rφ ∝ r² mode does not couple to the TRD layer radial pattern. In the joint fit at `solB1A3`, the correlation of
  d₆ with the layer-linear chamber radial mode (L5−L0) is 0.00 at B = 0 and 0.05 at B ≠ 0. Imposing d₆ = −13 µm (the B = 0
  request) changes the B± layer request (+5.3 / +6.3 mm) by 0.02 mm, and −50 µm by 0.1 mm. A charge-independent curvature
  gives a charge-independent residual, while the radial derivative (snp) is odd in q·B at B ≠ 0, so the two cancel over
  the charges. A polarity-even radial request at B ≠ 0 needs a residual odd in q·B, i.e. a curvature *scale* mismatch
  between ITS and TRD (Δκ ∝ κ). Possible sources: ITS radial scaling, the field map between the ITS and TRD radii, or the
  energy loss (in particular with the less precise material LUT, `corrType=2`). The pT dependence measured at IBEdge
  (flat for |q/pT| < 0.7, growing fast above) favours the energy loss, whose Δκ/κ grows at low p (∝ 1/β² for heavier
  particles). The two others would give a pT-independent radial request.
- Per sector at `solB1A3` (joint fit of the sector TX/RY/RX, layer TX and TY of every sector; sectors 2 and 3 have no
  complete layer set in collisions):
  - Layer-linear chamber radial (L5−L0): B0 −3.5, cosmics −3.5, B+ +5.2, B− +5.9 mm on average. B+ − B0 is
    +8.6 ± 2.1 mm (rms over sectors), B− − B0 +9.6 ± 2.1. The B on/off difference is common to all sectors. The
    sector-to-sector variation is shared by the datasets (correlation B0–cosmics 0.66, B0–B+ 0.61, B+–B− 0.87) and is
    geometry.
  - Common radial per sector: correlation 0.93–0.98 between all datasets, the differences being offsets (B+ − B0 −2.6,
    B− − B0 −5.6, cosmics − B0 +4.8 mm, with 0.7–2.9 mm rms over sectors). Only the global offset is the run time-zero;
    the sector-to-sector radial pattern (rms ≈ 4 mm, sectors 10–15 at +5…+10 mm) is shared and must be freed, with only
    the global mean pinned.
- ITS modes at `solB1A3`, from a replay with stave-level and per-half-barrel ITS groups (HB0 = top, HB1 = bottom):
  - ITS transverse radial scaling (each stave moved radially by ε·R): requested ε = −0.35·10⁻⁴ (B0), +0.04…+0.10·10⁻⁴
    (B±, cosmics; top/bottom alike), with errors below 10⁻⁶. It is uncorrelated with the TRD layer-linear radial mode. An
    imposed ε shifts that request by ≈ 0.04 mm per 10⁻⁴ equally at B0 and at B ≠ 0, so it cannot produce the B on/off
    difference (≈ 8 mm would need ε ≈ −2 %).
  - ITS rφ ∝ r² mode per half barrel: B0 −13.3 (top) / −14.9 µm (bottom), cosmics −10.3 / −7.1 µm. At B ≠ 0 the mode is
    degenerate with the track curvature and takes the TOF artefact (−145 / −66 µm); not usable.
  - Free layer rotations per half barrel (rφ displacement, uniform rotation removed), L0…L6: B0 +2.5 +2.6 +2.6 +3.7 +2.7
    −1.0 −3.4 µm (top), the same within 0.7 µm at the bottom; cosmics OB +2.3 +0.9 −0.3 −1.3 (top), +2.0 +0.7 −0.7 −1.3
    (bottom), IB −1.2 (top) and +1.6 (bottom). The pattern is an OB inner-vs-outer relative rotation, common to both halves
    and seen by both B = 0 samples (B0 about twice as large), plus an IB-vs-OB term that differs between B0 (rotation) and
    cosmics (opposite rφ in top and bottom, i.e. a shift in x). It is not a pure r² shape.
- The B ≠ 0 layer-linear TRD radial request grew with the material LUT (`corrType=2`): +5.2/+5.9 mm (L5−L0) at `solB1A3`
  with `minPt=1.2`, against ≈ +4 with TGeo corrections and `minPt=0.4` at `solB1A2`. This points further to an energy-loss
  (material) mis-modelling across TRD.
- The per-layer slope method on the debug trees (regression of `dyUnb` on the seed tanψ) is not valid at B ≠ 0. The seed
  includes the TRD points, so the bending angle is correlated with the residual. On the HPLUT debug trees it shows no B ≠ 0
  excess, while the replay does. The pT dependence derived from it at IBEdge (the argument for energy loss) is withdrawn.
- Replays at the IBEdge reference (mean vertex from the fits): layer-linear chamber radial (L5−L0) B0 +0.3, cosmics −0.1,
  B+ +9.2, B− +9.9 mm for the HPLUT recording (`corrType=2`, `minPt=1.2`); B0 +0.8, cosmics −0.2, B+ +4.0, B− +4.9 mm for
  the first recording (`corrType=1`, `minPt=0.8`). The B = 0 data accept the IBEdge layer radii. The B ≠ 0 excess doubled
  between the recordings, although the higher pT cut should have reduced an energy-loss effect (∝ 1/p).
- Energy-loss toy (ITS + TRD points, no scattering, constant-curvature model fitted to a track whose curvature grows): an
  underestimated energy loss gives outer layers requested outward, i.e. the sign of the B ≠ 0 request. The size is
  independent of pT for a given fractional curvature error: +1.4 mm (L5−L0) per 1 % curvature step between ITS and TRD,
  +1.0 mm per 1 % per TRD layer. The 4–10 mm excess would need 3–7 % of the momentum lost unaccounted (tens to ≈ 100 MeV at
  1.5 GeV), against a few MeV actually lost in the TPC and TRD material. Energy loss is excluded as the main cause.
- `doc/ITSEnvTRDTOFB0only_RB.json`, for B0 + cosmics only:
  - ITS half-barrel RZ free, with one `curl` group equation;
  - TRD sector TX, TY, TZ, RZ free, with one `globalModes` `radialScaling` equation over the sectors only (against the
    time-zero offset), sector RX/RY fixed;
  - chambers TX, TY, RX with priors, no layer pins;
  - TOF free with priors, no `curl`.
  It needs `groupModesRelToRoots=false`.

## 2026-11-20: the B ≠ 0 TRD radial excess sits in the composed vertex records

- `o2-dev-mille-replay` got two record selections: `-pv no|only` (records without / with vertex labels) and `-pt min,max`
  (all tracks of the record in the range of a pT estimated from the precision of their ITS multiple-scattering kinks; the
  records store no track parameters). The estimate reproduces the `minPt` edges of the recordings (0.8 and 1.2 GeV) with
  about 20 % resolution per track.
- Layer-linear TRD chamber radial (L5−L0, mm) at the IBEdge reference, B+ / B−:

  | recording | all records | single-track (no vertex) | composed vertex records |
  |---|---|---|---|
  | `corrType=1`, `minPt=0.8`, `maxTracksPerMultiTrackPV=10`, `extraClsErrYTRD=0.1` | +4.0 / +4.9 | +0.0 / +0.5 | +5.0 / +5.9 |
  | `corrType=1`, `minPt=1.2`, `maxTracksPerMultiTrackPV=4`, `extraClsErrYTRD=0.05` (coll_TGEOPt12, small sample) | +6.7 / +7.0 | +2.8 / +3.1 (±0.8) | +7.6 / +7.9 (±0.4) |
  | `corrType=2`, `minPt=1.2`, `extraClsErrYTRD=0.05` (HPLUT) | +9.2 / +9.9 | +4.1 / +4.1 | +10.0 / +10.8 |

  At B0 (first recording) single-track and vertex records agree (+0.6 / +0.9). The excess is therefore mostly a property of
  the composed vertex trajectories at B ≠ 0, the same for both polarities. The single-track records of the first recording
  agree with B0 and show only a soft-track term: +1.0 / +1.9 mm below an estimated pT of 1 GeV, −0.1…−1.6 above.
- Between the TGeo and LUT recordings at the same `minPt` the vertex records go from +7.6 to +10.0 (B+): the material
  treatment contributes, but the main part is common.
- The coll_TGEOPt12 recording has χ²/ndf 1.28 and rejects 18 % of the single-track and 49 % of the vertex records (first
  recording: 0.93, 5 %, 16 %), presumably from the tighter `extraClsErrYTRD`. Its mean-vertex reference is already at the beam,
  so it must be replayed with the HPLUT mean-vertex solution (`mvA3.res`), not the first recording's one (+0.3, −1.8 mm),
  which rejects 84 % of its vertex records.
- Debug trees, TRD `dyUnb` regressed per layer on the geometric part of the local angle, u = y/x (fit-independent), and on
  the bending part v = tanψ − u. A radial misplacement, or any chamber-local effect (drift, Lorentz angle), depends on the
  total tanψ and gives the same slope in u and v. L5−L0 of the slopes (mm), B+:

  | sample | single-track u / v | composed u / v |
  |---|---|---|
  | first recording (10 % debug) | −4.0 / +0.8 | −5.8 / +4.9 |
  | coll_TGEOPt12 B+ | −4.0 / +5.8 | −4.7 / +10.5 |
  | coll_TGEOPt12 B− | −3.3 / +6.6 | −4.3 / +11.2 |

  The replay radial request follows the v slope, which dominates its derivative (rms v ≈ 0.25 against rms u ≈ 0.08). The u
  slopes agree between the record types. The B ≠ 0 excess is therefore a residual proportional to the bending, with an
  antisymmetric layer pattern (about −5 at L0 to +6 mm at L5 per unit v). It is a track-model effect, not TRD geometry.
  It is stronger in the composed vertex trajectories, and in coll_TGEOPt12 also in the single tracks. The B0 v slopes
  (rms v ≈ 0.01) are dominated by the correlation of the seed direction with the residual and are meaningless.
- Reading the composed-trajectory code (`buildGBLVertex`, `computeVertexTransformation`, GBL `prepare` for a 2-row
  inner transformation) found nothing wrong. With a geometric constraint GBL gives every track its own curvature, the vertex
  fixes only the offsets at the vertex point, and d(Y,Z)/d(vertex) is correct to first order.
- coll_TGEOPt12NoMV (`useMultiTrackPVConstraint=0`, otherwise as coll_TGEOPt12). Its single-track vertex records must be
  replayed at the recorded reference: a few-µm mean-vertex correction (`mvA3.res`) rejects 99 % of them, the per-track
  vertex being tied tightly to the mean vertex. Layer-linear TRD radial (L5−L0, mm), B+ / B−: single-track with vertex
  +7.1 / +7.5, without vertex +3.2 / +3.3, against composed +7.6 / +7.9 (coll_TGEOPt12). The composition is not the cause,
  the vertex constraint enhances the effect. The debug trees agree (v slopes +9.8 / +10.4 with vertex).
- Per-layer ITS radial shifts, fitted jointly in the replay, are constrained to ≤ 0.3 µm even at B0, request a few µm,
  and leave the TRD mode unchanged: not the cause.
- In estimated-pT bins (vertex records of coll_TGEOPt12NoMV) the request is +25, +17, +7.5, +0.5, −1 mm from the lowest
  to the highest bin. Single tracks without vertex show the same.
- The whole effect comes from tracks with a large v (|v| > 0.3 at the TRD, i.e. a bending only a track below ≈ 0.6 GeV can
  have). With |v| < 0.3 the u and v slopes agree (L5−L0 −1 to −3.6 mm) in every sample and record type, and composed,
  vertex and single tracks no longer differ. Those tracks have a wrong GBL seed: their reference q/pT corresponds to
  ⟨pT⟩ ≈ 0.9 GeV although `minPt` = 1.2 cut the input tracks. The seed misses the outer ITS layers by up to 1 cm and the
  TRD by 30–60 cm. GBL recovers the positions, but it linearizes around a direction at the TRD that is wrong by ≈ 0.2.
  The radial derivative (tanψ of the seed) and the residuals inherit this.
- Fraction of accepted tracks whose seed misses the outermost TRD point by more than 25 cm: B0 0 %. At B ≠ 0, q·bz > 0
  (the same bending sense in the lab for both polarities): 8 % (single) / 20 % (composed) in the first recording,
  17–18 % / 28–31 % in coll_TGEOPt12(NoMV). For q·bz < 0: 2–7 %. The bad seeds come from the inward KF refit
  (`refitTracks` → `Track::fitTrack`, started from the input global track with `useStableRef`), and the vertex update
  makes them more frequent.
- To find which step produces the wrong seeds, the debug tree (`DbgTrack` v2) now stores the input track (`inputPar`), the
  state after the inward KF refit before the vertex update (`kfInnerPar`) and the GBL correction to the reference state
  (`gblRefCorr`, q/pt first).
- coll_TGEOPt12NoMV_LinRef0 / LinRef1 (`useStableRef=false` / `true`, otherwise as coll_TGEOPt12NoMV, with the new debug
  fields). The q/pT chain of the bad-seed tracks (seed off > 25 cm at the last TRD point), relative to the input track:
  inward KF refit +150 % (median, 10–90 %: +80…+200 %), reference after the vertex update the same, GBL fit +12 %. The
  wrong curvature comes from the inward KF refit, not from the vertex update. GBL brings it back only partly (linearization
  far from the solution). Good tracks: KF and GBL within a few % of the input.
- `useStableRef=false` cuts the bad seeds from 28 % to 6.5 % (q·bz > 0) and from 6–7 % to 2–3 % (q·bz < 0). From the same
  input it gives 33 % more accepted vertex records, rejections of 13 % instead of 21 % and ⟨χ²/ndf⟩ 1.20 instead of
  1.22–1.24. The layer-linear TRD radial (L5−L0) drops from +6.9 / +7.5 to +1.4 / +1.2 mm (B+ / B−), the debug v slopes from
  +10 to +1.5–2 mm. The B ≠ 0 radial excess is essentially an artefact of the KF seed with the stable linearization
  reference. The remaining bad seeds (still mostly q·bz > 0) point to the TRD/TOF-dominated start of the inward refit.
- `AlignParams.CollGBLSeedITSOnly` (default true): for collisions the GBL seed is the inward KF refit over the ITS points
  only, started from the outward refit state at the outermost ITS point. The outer detectors, misaligned at the start,
  then cannot bias the seed curvature, and `useStableRef` keeps its role of a stable linearization reference within the
  ITS.
- coll_TGEOPt12NoMV_LinRef1_ITSSEED / LinRef0_ITSSEED (`CollGBLSeedITSOnly=true`, `useStableRef=true` / `false`). Bad
  seeds (off > 25 cm at the last TRD point) 0.2–0.4 % / 0.0 %, with no q·bz asymmetry left (before 28 % / 6.6 % for q·bz > 0).
  The KF seed curvature agrees with the input within a few % (10–90 % range ±8 %). Layer-linear TRD radial (L5−L0), B+ / B−:
  −0.5 / −0.9 (stable reference) and −0.5 / −1.0 mm (without), the same for records with and without vertex. The debug u and
  v slopes now agree (−1.2 to −1.5 mm for the vertex records). The B ≠ 0 radial excess was entirely an artefact of the seed
  fitted through the misaligned TRD/TOF. With the stable reference there are still 8 % fewer accepted records
  (61k against 66k, B+) and a higher ⟨χ²/ndf⟩ (1.233 against 1.212).
- The GBL-fitted curvature differs from the input track by +1.1…+2.2 % for q·bz > 0 and −1.8…−1.9 % for q·bz < 0. This
  is the sign pattern of the one-pad TRD rotation still present at the IBEdge reference (about 1.5 % of the curvature
  expected at 1.2 GeV), and it should vanish once the TRD rφ is aligned.
- B0 added to coll_TGEOPt12NoMV_LinRef1_ITSSEED (same settings, `meanPtB0=0.45`). Replay at the recorded reference
  (L5−L0 layer-linear TRD radial / common TRD rφ shift at the sector origin): B0 +0.04 ± 0.13 / +4.7, B+ −0.49 ± 0.24 / +5.3,
  B− −0.90 ± 0.19 / +5.6 mm. Records with and without vertex agree. The B on/off radial difference dropped from 6–9 mm to
  0.5–0.9 mm (2–4 σ). The debug u slopes (L5−L0) agree between B0 (−1.9) and B± (−1.3 to −2.0). The B0 sample rejects
  32 % (single) / 23 % (vertex) of its records at ⟨χ²/ndf⟩ 0.98.
- With the ITS seed the rotation modes also agree between the fields (replay at the recorded reference, B0 / B+ / B−):
  TRD sector RY +3.0 / +3.4 / +3.4 mrad (old seed +2.5 / −2.9 / −1.6), TOF common rotation −0.8 / +0.3 / +0.5 mm at the
  TOF radius (old −0.6 / +2.1 / +1.6). The TOF B on/off tension is halved but not gone.
- `doc/ITSEnvTRDTOFch4_RB.json` for the recordings with the ITS seed. It is ch3 without the per-layer `radialScaling`
  group modes, which were there for the seed artefact and pushed the weakly measured chambers out. The sector `TZ` pin
  (time-zero) is statistics-weighted (`weightedMean`), and TOF keeps `curl` because of the remaining TOF rotation tension.

## 2026-11-20: axis tilts of the ITS layers and of TRD/TOF

- Fit of the displacement field (aligned − ideal) of the chip / chamber / strip centres:
  `d_T = t + s·z + r̂(a + b z) + φ̂(c + e z)`. The axis tilt is `s = (dx/dz, dy/dz)`. In µrad:

  | geometry | IB | OB L3–4 | OB L5–6 | TRD | TOF |
  |---|---|---|---|---|---|
  | V10N | −70, +55 | −60, −698 | −128, −415 | −48, +873 | +155, +597 |
  | IBEdge | −93, +105 | −53, −694 | −119, −372 | −35, +754 | +160, +511 |
  | solB1A3 (all fields) | = IBEdge | = IBEdge | = IBEdge | +161, +101 | +123, +495 |
  | solB1A4 (B0 + cosmics) | = IBEdge | = IBEdge | = IBEdge | +28, +381 | +101, +545 |

  - The ITS tilts come from V10N; IBEdge changed them by ≤ 50 µrad. The IB per-layer values scatter (L1 dy/dz +675).
  - The OB is tilted by −0.37…−0.70 mrad in dy/dz relative to the IB axis and to Z: −0.27 mm in y at |z| = 72 cm. This is
    a shear, not a rotation: its z part `dz/dy` is −0.2 mrad, where a rigid rotation would give +0.4…+0.7.
- The TRD/TOF fits did not take over the OB tilt.
  - TRD dy/dz went from +0.75 mrad toward zero (+0.0…+0.4), i.e. toward the IB axis and the beam line, not to the OB value
    (−0.40). TRD dx/dz went from −0.04 to +0.03…+0.16 mrad, away from the ITS (−0.1).
  - TOF is unchanged: it is held by priors and `curl`. Relative to the OB it stays +0.9 mrad in dy/dz, about 3 mm at the
    TOF ends.
- The fitted TRD z ∝ y (`dz/dy` −0.35 → −3.0 mrad, the sector "TZ −8.9 mm·sinφ" of `solB0A1`) is not an ITS tilt: the ITS
  has −0.2 mrad. It matches the −2.1 mrad already in the TOF CCDB alignment.
- The mean-vertex slopes fitted at IBEdge (−0.07, +0.03…−0.07 mrad relative to the recorded ones) are measured by the ITS.
  They follow the ITS (IB) axis and say nothing about the ITS tilt relative to the true beam.

## 2026-11-21: TGeo material with per-thread navigators loses OB scattering in the workflow

- `solB1A4` (`ITSEnvTRDTOFch4_RB.json`, `coll_TGEOPt08ITSSEED` + `cosmHP4LUT`): 32 % pede rejects. `o2-dev-mille-replay` gives the same
  rejection at the recording reference (no solution) as with the solution (B± 34–36 % vs 32–35 %), so it is in the records.
- Debug trees, single accepted tracks, B+: every recording with `corrType=1` made after e06a236 (`setMaxThreads` with `nthreads>1`)
  has ITS χ² ≈ 0.8–1.0 per coordinate and kink χ² ≈ 1.1–1.5 per kink coordinate, against 0.25 / 0.55–0.8 for the LUT
  (`collHPLUT`), and an OB unbiased track error of 15–18 µm against 24–29 µm. This includes the `coll_TGEOPt12*` test samples.
- Same 400 B+ TFs recorded locally three times: `corrType=1 --nthreads 4` → OB track error 15–18 µm, 36–44 % of the tracks above
  the 3σ χ² quantile; `corrType=1 --nthreads 1` → 24–29 µm, 7–10 %; `corrType=2 --nthreads 4` → identical to the serial TGeo.
  A standalone `GeometryManager::meanMaterialBudget` test (aligned IBEdge geometry, 4 std::threads after `setMaxThreads`) agrees
  with the serial result and with the LUT within 2 % per ITS layer, so the loss happens only in the workflow context.
- The old `coll` recording (11-16, `corrType=1`, 4 threads) predates the O2 mutex-free change and is not affected in the same way.
- Cure: `setMaxThreads` moved to the end of `AlignmentSpec::initOnFirstTF`, after `buildHierarchy()` and the TPC/TRD/vertexer
  initialisation (TGeo creates its per-thread volume data only once, in `SetMaxThreads`). The same 400 B+ TFs with
  `corrType=1 --nthreads 4` then reproduce the 1-thread recording exactly, with every OpenMP thread on its own TGeo slot and
  navigator. Which of the later initialisations invalidated the per-thread data was not isolated.
- `solB1A4` redone on the LUT re-recording `coll_LutPt08ITSSEED` (+ `cosmHP4LUT`): pede rejects 5.0 % (590k / 11.8M) instead of 32 %.
  Replay, reference → solution: B0 7.3 → 7.2 %, B+ 1.9 → 1.7 %, B− 2.1 → 2.0 %, cosmics 9.7/9.4 → 8.7/8.3 %. The residual
  B0 noPV rejection (6.4 %) is the known soft-track multiple scattering. The TRD/TOF geometry is practically that of the TGeo
  solution (layer rφ −0.15…−0.23 mrad, TOF −0.88 mrad, TRD tilt dx/dz +0.45, dy/dz −0.30 mrad, stacks within ~0.5 mm, sector
  DOFs correlated 0.97–1.00), so the lost OB material had inflated χ² and rejections without moving the outer detectors.
- The largest remaining pull is the TRD radial scale (`TRD_sect_TZ`, pinned via `pinChildrenMean`, and the per-stack chamber
  TX): +1.9 mm B+, −1.0 mm B−, +1.1 mm B0, +5.7…+6.2 mm cosmics. Its field-odd part points to a TRD calibration (drift
  velocity/t0, ExB) rather than to geometry.

## 2026-11-21: half-barrel-level change between pp LHC25ab 562837 and PbPb LHC25as 568664 (B+)

- In checkResid, late − early is the same for best/IBEdge and V10n (double difference < 0.5 µm), flat in pT, and absent
  for 564644 − 562837. So the ITS changed after 564644, independently of the alignment.
- IB: the top/bottom halves moved relative to the vertex (vertex ty ∓3.3/+3.6 µm), with Lr1 moving against Lr0/Lr2. The
  resulting ~1 µm IB sagitta changes the IB sub-track q/pT by ±0.05 GeV⁻¹, which explains the IB-OB degradation on p.9.
  About 85 % of it is rigid per HB; the rest (~1 µm) is stave-level.
- OB: Lr4 HB0 shows an apparent ty of −13 µm. The v21-type configs pin the OB HBs, so they cannot absorb it.
- Added `doc/ITS_RB_HBonly.json` (per-layer HB TX/TY/TZ/RZ plus the mean vertex) for a realignment with the late data.
  Details and tests are in [HBLevel_adjustment.md](HBLevel_adjustment.md).

## 2026-11-21: MilleRes ignored the rigid-body corrections of volumes without a configured DOF set

- `--output MilleRes,Regularize` on `solB1A4` without `dofConfigJson` produced the initial IBEdge TRD/TOF geometry,
  regularized: `result.json` held only the mean vertex. `Volume::getFittedLocalDelta` (and the `result.json` writing)
  required a `RigidBodyDOFSet`, which only the DOF configuration creates, so every rigid-body value of `millepede.res`
  was silently dropped.
- Fix: the values are looked up by label for every rigid-body alignable volume, configured or not. The same run then
  places the TRD chambers and TOF strips as the unregularized `solB1A4` objects do on the ideal geometry
  (≤ 4·10⁻¹⁰ cm TRD, 2·10⁻⁸ cm TOF), against up to 2.3 cm and 10 mrad before.
- The internal closure check reports 5·10⁻⁵ cm for TOF, before and after the fix. Applied to the ideal geometry, the
  output agrees with the input to 2·10⁻⁸ cm, so the deviation comes from the ideal TOF matrices derived from the loaded
  IBEdge geometry (G^id = P^-1 G), which match the IBEdge TOF objects only to that level.

## 2026-11-22: cosmic preparation and junction Jacobian performance

- Profiling ITS/TRD/TOF cosmics with material LUT corrections found numerical transport Jacobians dominated GBL construction;
  the leg junction consumed 57% of Jacobian time in 16% of calls. Parallelizing its five independent Ridder columns with OpenMP
  gave about 2.9x faster junction calculations while preserving the recorded matrices, error estimates and Mille binaries bit for bit.
  Parallel calls are disabled inside an existing OpenMP region; `--nthreads` sets the team size.
- At `ridderEps=1e-8`, raw-hit availability checks and early exits after detector selection reduced point collection from 2.032 s
  to 0.678 s (mean summed worker time over two runs each), and cosmic processing by about 22%. The 60,900-track sample used
  3 lanes and 5 junction threads per lane. TOF preparation calls fell from 121,716 to 3,336; all 851 accepted tracks, recorded
  fits, Jacobians and Mille binaries were unchanged. Removing collection restarts alone saved little: most restarts dropped ITS
  before processing the other detectors.
- Builds and benchmark installs were isolated from the shared installation used by concurrent jobs. The measurements cover
  B=0 ITS/TRD/TOF cosmics with material LUT corrections; they do not establish performance for TPC or TGeo material corrections.
- The production implementation was rebuilt in isolation and checked with 1 and 5 threads at `ridderEps=1e-8`: both LUT runs
  reproduced the original serial control's Mille binaries on 60,900 tracks (851 accepted). A bounded TGeo check on 2,031 tracks
  (38 accepted) also produced byte-identical Mille binaries with 1 and 5 threads.


## 2026-11-22: rφ ∝ r² sagitta of the best/IBEdge ITS in B=0 collisions and cosmics

- The TRD/TOF `pre` residuals of `solB1A4` keep B0 apart from B±. Each B=0 track was refitted with ITS points only as
  v = a + b·u + k·u·r (u along the track, v transverse, r the point radius), i.e. an rφ displacement k·r² on its ray;
  k = 10⁻⁶ cm⁻¹ is a 1.7 µm sagitta over the ITS. B± tracks absorb k into q/pT (about 1.3·10⁻³ (GeV/c)⁻¹ at 0.5 T), so
  only B0 and cosmics see it. Its B0 ITS seed extrapolation error, against the curved B± seed, predicts
  dy(B0) − dy(B±) ≈ 1–1.5 mm·k at the TRD/TOF radii. The per-TRD-sector `pre` difference B0 − B+ of `solB1A4` follows the
  10 % k(φ) map with correlation 0.92 and slope 0.9 mm per unit of k. B0 − B− behaves alike (correlation 0.91, slope 0.63)
  with an extra +0.33 mm offset; B− was recorded later. The `unb` B0 − B+ difference (±0.1–0.2 mm) does not follow k.
- 100 % debug of `Mille_IBEdgeS4R/coll_LutPt08ITSSEED` (B0 + B0Bis, 28.3M tracks with |tgl| < 1):
  - k = 1.38 ± 0.02, or 1.18 ± 0.03 with outer-barrel points only.
  - φ dependence: 2nd harmonic amplitude 1.9 with maxima near 170° and 350°; 1st harmonic 0.8 (1.65 for the outer barrel)
    peaking near −15°. The maximum, 7.3, is at 350–360°.
  - tgl dependence: k grows towards |tgl| → 1, with an odd part larger at positive tgl (3.1 at tgl 0.75…1 against 1.5 at
    −1…−0.75). Tracks with |tgl| > 1 have larger k still: the mean is 1.92 with no tgl cut.
  - B0 and B0Bis agree bin by bin.
- Leave-one-layer-out residuals per stave and z, measured against a straight line through the other layers, show the r²
  pattern. At φ ≈ 350° they are Lr0 +9, Lr1 +6, Lr2 +4, Lr3 −14, Lr4 −14, Lr5 0 and Lr6 +21 µm. The outer barrel staves
  also show U shapes along z. Against a line plus k·u·r, every stave and z bin is within 1–2 µm, with a few Lr6 staves up
  to 4.5 µm. So no single stave is displaced: the pattern is one smooth field over all layers, whose amplitude varies with
  φ and z.
- Cosmics, with the same cosmic data processed in both geometries:
  - Fitting each half as a collision track (≥ 2 inner and ≥ 3 outer barrel hits; about 2k halves) reproduces the collision
    k(φ) where cosmics are populated. For example, 140–160° gives 4.8/5.1 (2025/2026) against 4.2 in collisions.
    280–320° is in tension.
  - Over the whole cosmic, with a free top/bottom step, |d0| < 5 cm and the same runs: IBEdge k = −0.09 ± 0.16 (2025) and
    −0.24 ± 0.15 (2026); v10n 1.10 ± 0.28 and −0.42 ± 0.32.
  - The k spread between cosmics is 13 for IBEdge against 22–24 for v10n.
  - In both geometries the step depends on the cosmic direction, from about −20 µm at 20–60° to +30…+40 µm at 140–160°.
    It is the same in 2025 and 2026.
- What tracks cannot fix, per ray direction (φ, tgl). An rφ displacement a(φ, tgl) + b(φ, tgl)·r that is linear in r leaves
  every straight track from the origin straight.
  - The free vertex of the composed multi-track PV records constrains a, up to a translation of the whole ITS.
  - b(φ) is an azimuthal remap φ → φ + b(φ), the same at all radii. Nothing in the ITS stops it: collision tracks of any
    field are close to radial, and sensor overlaps move together with their neighbours.
  - A cosmic with d0 ≈ 0 joins the rays φ and φ+π. Its 2 line parameters, against 4 for the two rays, determine only
    a(φ) + a(φ+π) (the step) and b(φ) − b(φ+π) (a kink at the centre). An even-harmonic b, with b(φ) = b(φ+π), stays free.
  - Only cosmics with large d0, which cross the layers at different φ, or detectors outside the ITS are sensitive to b(φ).
  - So a correction of the measured k(φ, tgl) by δrφ(L) = −k·r_L² + b·r_L is equivalent for every b. Which layers "carry"
    the sagitta has to be fixed by a gauge choice (smallest displacement, or the hierarchical mean constraints), not by
    the data.
  - B± do not help either. Along a curved track the azimuth changes by Δφ ≈ r/(2R), at most about 3° at 0.5 GeV/c in
    the ITS. So b(φ)·r reads as b(φ₀)·r, which goes into the direction, plus b′(φ₀)·r²/(2R), which goes into q/pT as a
    scale factor 1 + b′(φ). Both charges and both field polarities scale the same way, so B0 + B+ + B− fix k (B0 sees it
    as curvature, B± as a q/pT offset) but not b. Cosmics in field compare q/pT between φ and φ+π, where an even b has
    the same b′. Only an external q/pT or direction reference fixes it: a mass peak versus φ, or TPC/TRD/TOF.
- The top/bottom step of the cosmic fit (line + k·u·r + step) is not a compromise with B−:
  - The ITS-only line of a collision track, extrapolated to r = 0 with the fitted beam subtracted (xb, yb ≈ −0.2, −1.9 mm),
    measures the collision side of the same quantity: (d0(φ) + d0(φ+π))/2 = a(φ) + a(φ+π), the same combination the step
    gives for a cosmic through the axis.
  - This is within ±1 µm in B0 (562862) and B+ (562837). B− (564644) differs by up to −2 µm at 100–180°.
  - The step depends on |d0|. Below 2 cm it is −0.3 ± 0.9 µm overall, with at most +12 ± 5 µm at 140–160°. At
    |d0| = 5–20 cm it runs from −27 to +42 µm.
  - A z-dependent rφ field c(φ)·z², which would make k rise with tgl² in collisions, predicts a step growing as z0².
    The data reject it.
  - Large-d0 cosmics cross the layers at different φ and are not radial, so the per-ray decomposition does not apply to
    them. They are also the only ITS data sensitive to b(φ). The step at the same |d0| differs between the IBEdge and v10n
    geometries by 15–26 µm in several φ bins, so it reacts to geometry differences that are invisible to collisions.
- Geometry-only decomposition (scratch macro `geoK.C`):
  - Method: the same chip-local points are placed with two ITS snapshots; the rφ part of the displacement is averaged per
    fine (φ, t = z/r) cell and layer and fitted per cell as a + b·r + k·r².
  - v10n → C1v21w3 accounts for all of it: Δk = +1.06, 2nd harmonic 1.41 at −15°, 4th harmonic 1.21. C1v21w3 → S4R and
    S4R → IBEdge change k by ≤ 0.05 and the ray rotation by ≤ 2 µrad.
  - Against the measured IBEdge map, Δk reproduces the 340–360° peak (7.1 of 7.1), most of the 140–170° peak, 1.4 of
    the 2.05 2nd harmonic, the 4th harmonic and the rise at positive tgl. Cell by cell the correlation is 0.60, because the
    geometry difference has stave-level structure that the track average smooths.
  - Measured minus Δk predicts for v10n k = 0.4 on average, 1st harmonic 1.5 at −32° (peak at 300–320°), 2nd harmonic 0.7.
    So v10n has its own, smaller sagitta. Still to be checked with v10n B0 data.
  - The displacement behind Δk: the outer barrel moved in rφ relative to the inner barrel. Its 2nd harmonic is 10–14 µm
    (phase 47–61°) and its 1st harmonic grows from 4 µm (Lr3) to 13 µm (Lr6), with up to −34 µm at Lr3, 340°. The inner
    barrel moved by ≤ 2.6 µm.
  - In the `rayModes` sense (mean of dv/r, inner and outer barrel weighted equally), the same step moved the ray rotation
    by +40 µrad on average (+70 µrad inner barrel, 0 outer barrel), 2nd harmonic 30–50 µrad and 3rd harmonic up to
    70 µrad. Over tgl it varies from 10 to 60 µrad, which is 30–180 µm at the TRD/TOF radii. That is the size of the b
    ambiguity between the two geometries.
  - The MC `o2sim_geometry.root` (no physical nodes) differs from both aligned geometries by mm and mrad, so it is not
    usable as an ideal reference.
- Replay of the IBEdge ITS-only Mille data (`Mille_IBEdgeS4R/*_ITS`, recorded at IBEdge, no solution applied):
  - Setup: one group per stave on TX, 192 groups. The modes were built with `_groupH.txt`:
    - sagitta modes δrφ = −k·r·(r − 16.3 cm), harmonics n ≤ 2;
    - ray rotations rφ ∝ r, harmonics n ≤ 2.
  - Hu of a mode is Σ J²·Hu, since every hit row belongs to one stave.
  - In collisions the sagitta modes are as degenerate as the ray rotations: Hp/Hu ≈ 1–5·10⁻⁴ for both. The B0 records
    carry multiple-scattering kinks with the errors of `meanPtB0` = 0.45 GeV/c (about 2–3 mrad per layer). These absorb
    a µm-level sagitta, which needs kinks of about 0.1 mrad. Cosmics reach Hp/Hu ≈ 2.5–6·10⁻³ for k, against 8·10⁻⁵ for
    the global rotation.
  - The requested sagitta removal (K0, joint K + rotation fit) differs between datasets:
    - B0 0.82 ± 0.02 (B0Bis 0.83). The debug-tree straight-line medians give 1.38.
    - B+ 2.35 ± 0.03 and B− 0.01 ± 0.03, opposite for the two polarities.
    - Cosmics 0.27 / 0.34 (2025 / 2026).
    - The harmonics and the single-mode fit of the B0 map do not agree with each other either.
  - In the Mille model the sagitta is therefore a quasi-weak mode. Systematic pulls far below their statistical errors
    set it, among them the field-odd B+/B− split. That is how C1v21w3 could build it. A correction fit must not leave
    it to pede.
- Direct sagitta correction relative to the IBEdge recording (scratch `mkSagCorr.py`):
  - Correction: δrφ = −k(φ, t)·r·(r − β), from the B0 map at |tgl| < 2, 36 × 16 cells, 1-2-1 smoothed. Applied as inner
    barrel stave TX and outer barrel module TX, since the `_ITS` recording has no chip DOFs.
  - β = ⟨r⟩_OB = 33.2 cm keeps the outer barrel ray rotation, and so the B0 straight-line directions; corrections stay
    within −34…+36 µm. The combined gauge, β = 18.3 cm, reaches −115 µm at the Lr6 ends.
  - Self-check over the |tgl| < 1 cells: Δk = −1.40 against −1.45, correlation 0.96.
  - Replay on the ITS-only Mille data: the correction raises ⟨χ²/ndf⟩ for B0 (+0.25 %), B± (+0.6 %) and cosmics (+1.1 %);
    the combined gauge raises it about half as much.
  - Scanning the B0 replay with the correction scaled by λ = −1, 0, 0.5, 1 gives a χ² that is nearly quadratic in λ,
    with its minimum at λ ≈ −0.07. The cost therefore comes from parts of the correction the tracks do see. A likely
    source is the module granularity: an Lr3/4 module spans Δt ≈ 1, an Lr5/6 module Δt ≈ 0.5, so the layers sample k(t)
    differently. The weak sagitta part itself is cheap. The Mille χ² cannot validate the correction; a B0 debug pass with
    the corrected geometry can.
  - Doubled cosmics: the half-cosmic k(φ) follows the collision map except at 300–340°, where cosmics are about 2.5 lower,
    and at 120–140° (3.1 against 1.6).
- B0 debug with the corrected geometries (`/data0/devAlg/tests/sag/v0` from `millepede_sagCorr.res`, `v1` from the
  combined-gauge `_comb.res`; 5.9M B0 tracks, compared with the same run in IBEdge):
  - The sagitta is removed. k changes from 1.36 ± 0.04 to 0.07 ± 0.05, and every harmonic is ≤ 0.25 (2nd harmonic 1.90 →
    0.21). The tgl columns are flat within ±0.4. Residual structure is left at stave scale: +2.3 at 350–360° next to −1.5 at
    0–10°, and +1.4 / −1.1 at 170–190°. With outer-barrel points only, a 1st harmonic of 1.0 remains: the curvature within
    the outer barrel is not purely the r² shape of the full ITS.
  - Leave-one-layer-out against a straight line, rms of the stave medians: Lr0 3.5 → 0.8, Lr3 5.4 → 1.2, Lr6 6.5 → 2.0 µm.
  - Against line + k, the outer barrel z profiles get worse: z rms 0.6 → 2.3 / 2.0 µm (Lr3/4) in v0, 0.9/1.3 → 3.8/4.6 µm
    (Lr5/6) in v1. The median per-track χ² rises by 0.2–0.4 % (B0) and 0.4–0.8 % (B±). This is the module-granularity
    cost the replay predicted.
- Reduced-fit inputs on top of v0 in `/data0/devAlg/tests/sag/v0R` (see its README):
  - Free: outer barrel module TX and RY, and the mean vertex. The inner barrel is fixed at v0.
  - Pins: 192 outer-barrel-only `rdphi` ray-mode cell equations, re-centred on v0. With the inner barrel fixed, they hold
    both b and k per cell.
  - Priors: module TX 10 µm centred on v0, module RY 1 mrad.
  - In a 2-file pede test, without the TX priors the Lr5/6 end modules drifted by 13–16 µm rms; with them, by 3–5 µm.
- checkResid of `sag/v0` against IBEdge (runs 562862 B0, 562837 B+, 564644 B−; per-sensor mean dy from the res trees,
  grouped into inner barrel staves and outer barrel modules):
  - B0 improves where the sagitta was: in Lr4 the module rms drops from 1.6 to 1.0 µm, and the change is −0.89 correlated
    with the applied correction.
  - B± get worse:
    - inner barrel stave rms 0.2–0.3 → 0.3–1.6 µm, most in Lr2;
    - Lr4 module rms 1.0 → 2.2 µm, with +5…+7 µm on stave 12 of both half-barrels;
    - outer barrel chips: the change grows with |z|/r, from about 1 µm rms at |t| < 0.5 to 2–3.6 µm at |t| > 1.
  - B+ and B− changes correlate at 0.93–1.00 per volume, so this is geometry, not statistics or field. A pure k r² + b r
    field would be invisible to B± tracks, so the changes are the parts of the correction that are neither.
  - Two sources:
    - Rigid inner barrel staves received k averaged over their chip positions, mostly at |t| > 1 where k is large, through
      the large gauge lever (r − β ≈ −30 cm, 1.15 µm per unit k in Lr2). The tracks through them mostly have |tgl| < 1.
    - Outer barrel modules got module-averaged k and the noisy high-|t| cells (errors 0.5–1.5), amplified by
      r(r − β): 2.7 µm per unit k in Lr3, 2.1 in Lr4.
  - The `v0R` reduced fit keeps the inner barrel fixed at v0 and pins the cell-level outer barrel values, so it would inherit
    both.
