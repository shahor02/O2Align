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
