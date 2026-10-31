# ITS weak modes, round 2: polar families, ray constraints, IB chip DOFs (handover)

Follow-up to [AlignmentWeakModes_handover.md](AlignmentWeakModes_handover.md) (sections 1-5: the
`ITSBarFix`/`solA1` investigation, `pinChildrenMeanNoCross`, harmonic modes, `ITS_RB_v2.json`).
Section 6 of that file holds the detailed numbers for this round; this file is the self-contained
summary and the starting point for the next step. Nothing described here is committed.

## Status in one paragraph

v2 is better than V10n in checkResid; the apparent OB degradation was a binning artefact. To
reach the large-scale Z/tgλ structures, v3 freed DOFs distributed in z (IB chips, OB modules) and
added the polar mode families. That opened **ray modes**: per-direction displacements ∝ r,
invisible to tracks from the vertex. v4 constrains them per (φ, tgλ) cell (`rayModes`) and adds
priors. Its corrections are 2.5× smaller than v3's at almost the same χ², and the OB improves in
checkResid. **The IB and vertex panels get worse than v2**, though. The cause is the new IB chip
DOFs, mainly a common per-layer chip in-plane rotation (`RY`) that the fit wants and the
reconstruction does not see. The root cause, a fit-vs-reconstruction discrepancy at IB chip level,
is still open.

## Where things are

| run dir (`/data0/devAlg/tests/`) | config | content |
|---|---|---|
| `resOp/solA1Iv1GC` | `ITS_RB_v1.json` | earlier reference (IB chip TX/TZ free, failed similarly) |
| `solA1Iv2GC` | `ITS_RB_v2.json` | reference: IB chips fixed; checkResid trees + `checkResid_hman.pdf` |
| `solA1Iv3` | `ITS_RB_v3.json` | polar families, z-distributed DOFs; `pede.log`, `millepede.res`, `result.json`, `weakModes.*` |
| `solA1Iv4` | `ITS_RB_v4.json` | ray constraints + priors; checkResid (`checkResid_hman.pdf`, hman roots), `algdump.root` |
| `solA1Iv4/runA` | v4 with IB stave `RY` | first v4 pass, superseded |
| `solA1Iv5` | `ITS_RB_v5.json` | v4 + `rMin: 10`; fit + MilleRes only, no checkResid |

All fits use the all-free Mille data of the `ITS_RB_all.json` production. The steering reads only
`mp2data_0.bin` of the 15 binaries per field (see open points).

## Code changes (uncommitted)

- **Polar mode families** (`Volume::GlobalMode`: `RadialPolar`, `TangentialPolar`,
  `LongitudinalPolar`). Field u = r^p·t^k (t = z/r) along r̂, φ̂ or ẑ, with p ∈ {0,1}, k ≤ 2
  (`GlobalPolarKMax`). Members that duplicate existing modes are excluded (`isPolarDuplicate`).
  The rotation field ω = ½∇×u was checked by finite differences. `NGlobalEq` = 39; names like
  `radialPolar_r0t1`.
- **`rayModes` clause** (`Volume::RayModesConfig`, `applyRayModesConfig`,
  `writeRayModeConstraints`). On the flagged volume it bins the sensors in (φ, t = z/r), optionally
  by r range. Per cell and component (rδφ, δz, δr) it writes Σ field·d/r = 0 over all free
  rigid-body DOFs on the sensor's path to the root. The displacement per DOF comes from a finite
  difference of the `AlignParam` matrix. Cells with too small a radial span or fewer than 2 terms
  are skipped. Keys: `nPhi, phi0, nT, tMax, rSpanMin, rEdges, rMin, rMax, components`. Documented in
  [DOFConfig_rules.md](DOFConfig_rules.md#constraining-ray-modes), including the rule not to
  combine it with in-span global modes.
- **pede grouping workaround** (`addPedeGroupingAnchor`, `minFreeRigidBodyLabel`). Every global-mode
  and ray equation gets the lowest free rigid-body label as an extra term with coefficient 1e-12.
  See "pede crash" below.
- **`WeakModes.cxx`:** `fitModes` uses `JacobiSVD` instead of `BDCSVD`. `modeScales` handles polar
  modes.
- **`macro/drawWeakModes.C`:** polar page (dr, rδφ, dz vs z/r per family).
- **Docs:** `DOFConfig_rules.md` (polar rows, ray section) and section 6 of the old handover.

Build: `ninja -C build -j 3 install`; run with `install/lib` prepended to `LD_LIBRARY_PATH`.

## Findings, in order

### 1. v2 vs V10n: the binning artefact; the residual structure sits below the stave

- **Binning:** the V10n hman files had 30 φ bins, the v2 ones 100. Compared at equal binning, v2 is
  better in the OB, e.g. Lr5 dy-vs-φ 28.1 → 19.7 µm and dz 23.7 → 8.9 µm (RMS over bins of
  truncated means, pp B0). The IB is at its ~1 µm floor.
- **Inside the OB modules:** the remaining OB dy structure is an offset between the two chip rows
  of each module (rows 7-13 vs 0-6: −10…−22 µm mean). It reproduces between pp B0 and OO B+ (per-chip
  correlation 0.73-0.84). No v2 DOF can touch it; OB chip `TX` would be needed.
- **IB Z structure:** the IB dy-vs-Z sawtooth (3 cm period) is a per-chip in-plane slope seen in
  V10n *and* v2. At Lr0 it goes from −3.5 to +2.5 µm/cm (chips 0 → 8), while its mean per layer
  is ≈ 0. Compare item 5, where the fit adds a non-zero common slope.

### 2. Polar families (v3)

For tracks from the vertex z ≈ r·t. A field r·F(φ,t) is absorbed by the track's φ, tgλ or
curvature, and a field F(t) is broken only by the vertex. With every DOF at z ≈ 0 (stave level)
this family was inaccessible. Freeing chips and modules opens it, so v3 constrains the polar
members.

### 3. v3 opened ray modes; the polar families are too few

v3 gained 2.8 % in χ² over v2, but its corrections were 2.5× larger (Lr5 module TX up to 864 µm).
Of the per-(φ, t)-cell dz/r variance, 87 % is common to all seven layers; rφ/r is common within
the IB and within the OB. The null space is a function of (φ, t) per component, not a short list of
modes. Hence `rayModes`.

### 4. v4: ray constraints + priors

| | v2 | v3 | v4 |
|---|---|---|---|
| objective at iteration 3 | 3.424·10⁷ | 3.328·10⁷ | 3.348·10⁷ |
| χ²/ndf | 1.12 | 1.025 | 1.030 |
| displacement RMS, all ITS | 59 µm | 104 µm | 42 µm |
| OB module TX RMS (max) | fixed | 73 (864) µm | 22 (91) µm |

- **Constraints:** 2655, full rank: 1974 hierarchical, 9 global, 720 ray.
- **Priors (assumed, not survey numbers):**
  - IB chips: TX, TZ 10 µm, TY 5 µm;
  - OB modules: TX, TY, TZ 30 µm.
- **Ray modes removed:** the ray-common dz/r is gone.
- **IB stave `RY` dropped.** With it, the fit had realised a chip rotation as stave RY plus a
  compensating chip-TX ramp (correlation 0.88, chip TX RMS 32 µm). This is another case of the rule
  of one representation per deformation.
- **Left over:** a small sagitta-like Lr3 vs Lr6 anticorrelated pattern (passes the ray constraint).

### 5. checkResid on v4: the IB and vertex get worse than v2

Same in pp B0 and OO B+. RMS over bins of truncated means, pp, µm:

| | V10n | v2 | v4 |
|---|---|---|---|
| Vtx dz vs φ | 4.3 | 1.8 | 5.4 |
| Vtx dz vs tgλ | 2.5 | 3.0 | 7.8 |
| IB dz vs φ, Lr0/1/2 | 2.0/1.3/2.5 | 1.4/0.8/1.2 | 2.1/1.4/2.0 |
| Lr0 dy vs Z | 1.8 | 1.7 | 4.6 |
| Lr2 dy vs Z | | 1.2 | 2.8 |
| Lr6 dz vs Z | 9.0 | 8.7 | 15.8 |

Improved: OB dy vs Z (Lr3-5), Lr5 dz vs φ, Vtx dy vs Z (3.9 → 1.3 µm).

- **The Lr0 dy-vs-Z sawtooth is chip `RY`.** The fit gives a common rotation per IB layer:
  +4.97/+1.15/−2.26·10⁻⁴ rad (Lr0/1/2, ±0.4·10⁻⁴). checkResid's V10n in-chip slopes are only
  +0.09/+0.02/−0.53 µm/cm. The chip-to-chip part of RY does follow the V10n slopes (correlation
  −0.46/−0.48 at Lr1/2), only the common part is unsupported.
- **Ruled out:**
  - geometry application: `algdump.root` shows exactly the fitted rotations, with no compensating
    chip-centre motion;
  - rotation convention: `AlignParam` is right-handed and active, as in `fillDerivatives`;
  - pivot transport: `computeJacobianL2T` is correct;
  - the IB ray cells. They are not null at the IB, since σ(z_v) ≈ 5 cm smears t by ±2, but v5
    (`rMin` 10) gives the same χ² and the same IB chip parameters.

## Bugs found on the way

- **pede `grpcon` (Millepede-II 4fd638f67f).** Merging of constraint groups is non-transitive and can
  leave empty groups. The redundancy check then indexes `globalParCons(0)` (`pede.f90:1902`), causing
  heap corruption and `free(): invalid next size` in LOOP1. This was localised with a bounds-checked
  pede build. v2 already had one empty group and survived by luck. The 1e-12 anchor term above makes
  pede group the long equations first. The coefficient cannot be 0, because pede's `addItem` drops
  zero factors. Worth reporting upstream.
- **`WeakModesDiag` fit B was void** (rank 0 or NaN conditioning, also in v2): `BDCSVD` fails on the
  exactly rank-deficient design matrix. Fixed with `JacobiSVD` (v3: rank 42/45, the 3 missing being
  the known duplicates). Disregard metric B in older `weakModes.json` files.

## Open points and next steps

1. **Decision pending, IB chip discrepancy:**
   - **(a) Workaround, v6:** v5 with the IB chips fixed (IB as in v2), keeping the OB module DOFs,
     priors and ray constraints. Run checkResid on it: the IB should recover while the OB gains are
     kept.
   - **(b) Diagnosis:** add a per-point residual dump to the alignment workflow (the `Debug` stream
     has only `MisRes`), run it on V10n for a few TFs, and compare the per-chip in-chip dy-vs-z slopes
     with checkResid. Suspects:
     - the IB measurement frame: `SensorITS::defineMatrixL2G` (chip node + ½ epitaxial shift) and
       `getFlatChipAlpha`, against the reconstruction's sensor matrix and tracking frame;
     - how cluster coordinates are carried into that frame;
     - different runs in the Mille data vs the checkResid samples.

   (a) and (b) can run in parallel. Chip-level IB DOFs have now failed twice (v1, v4), so (b) should
   be done before freeing them again.
2. **Lr6 dz vs Z got worse in v4** (8.7 → 15.8 µm), Lr3 too. The OB module DOFs may need the same
   per-chip/module scrutiny as the IB.
3. **OB chip `TX`** for the chip-row offsets in the modules: not yet freed. Expect to need `inversion`
   to handle the size.
4. **Lr5 dz** is anomalously broad (RMS 350 µm vs 220-245 µm elsewhere): not understood.
5. **Sagitta-like Lr3/Lr6 pattern** left by v4: watch it in B=0 vs B≠0.
6. **Priors** are assumed tolerances; replace them with survey numbers.
7. **Steering uses 1/15 of the data** (`mp2data_0.bin` only per field).

## Practical notes

- The DPL driver hangs on stdin when run inside a chained bash command; use `</dev/null`.
- `pgrep -f`/`pkill -f` on a substring of the command line matches the waiting shell itself.
- The O2 Python has numpy, but `python3 -I` breaks it (PYTHONHOME). The system python3 has no numpy.
- Analysis helpers used here are in the session scratchpad
  (`/tmp/claude-1000/-home-shahoian-alice-O2Align/29fd9c2a-.../scratchpad`, temporary). They include:
  - `perSens.C` / `chipFit.C`: per-sensor and per-chip residual fits from the checkResid trees;
  - `hm.C`: RMS over bins from the hman files;
  - `lab.py`: `millepede.res` → labelled JSON;
  - `rank.py`: constraint rank;
  - `grpsim.py`: replica of pede's grouping;
  - `ray.py`: ray-coherence test;
  - `gmrot.C`: applied chip rotations from `algdump.root`;
  - `mp2src/`: the bounds-checked pede.

  Move whichever are wanted into `macro/` before the scratchpad goes.
