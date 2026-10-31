# Weak Modes

## What a weak mode is

A *weak mode* is a coherent deformation of a whole detector branch that moves the sensors but
(almost) does not change the χ² of the tracks. A track fitted from the vertex absorbs such a
deformation into its own parameters (1/p<sub>T</sub>, φ, tgλ, d0, z0), so the data cannot tell the deformed
geometry from the nominal one. Millepede then has a (near-)null direction in the global matrix: a large
amplitude costs a negligible χ², and the value it takes is set by noise and systematics (material,
cluster-shape bias, B map), not by geometry.

Examples: a uniform radial scaling Δr = ε·r (absorbed by the momentum scale), a rigid tilt of the whole
barrel, an elliptical squashing of the layers, a twist along the beam.

The hierarchical mean constraint (`writeChildrenMeanConstraints`, `pinChildrenMean*`) does **not** remove
them: it removes only the net translation and rotation of the children of one parent, 6 equations per
parent. A deformation mode has no net rigid component of any single parent and passes through. It needs
its own operator, the `globalModes` / `rayModes` clauses of the DOF config
([DOFConfig_rules.md](DOFConfig_rules.md#constraining-global-deformation-modes)).

The diagnostic `WeakModesDiag` (`--output MilleRes,WeakModesDiag`) measures how much of every mode is
present in a fitted result, with or without the constraints, and `macro/drawWeakModes.C` plots it.

## Implemented modes

Cylindrical coordinates (r, φ, z), t = z/r (the tgλ of a track from the vertex), ε the amplitude.
The definitions are in `Volume::GlobalMode` ([Volume.h](../include/O2Align/Volume.h)); there are
`Volume::NGlobalEq` = 39 equations in total.

| mode (equations) | displacement field | meaning |
|---|---|---|
| `radialScaling` (1) | Δr = ε·r | uniform transverse dilation; absorbed by 1/p<sub>T</sub> (B≠0) and the impact parameter. The dominant weak mode of ITS |
| `zScaling` (1) | Δz = ε·z | longitudinal scale; absorbed by tgλ |
| `curl` (1) | rΔφ = ε·r | rigid rotation about the beam |
| `twist` (1) | rΔφ = ε·r·z | torsion: the two ends of the barrel rotated oppositely |
| `elliptical` (2) | Δr = ε·cos2φ, ε·sin2φ | layers squashed into ellipses (two axis orientations); r-independent amplitude |
| `bowing` (1) | Δr = ε·z² | radial sagitta along the beam |
| `radialHarmonic<ℓ>` (2 per ℓ) | Δr = ε·r·cos ℓφ, ε·r·sin ℓφ, ℓ = 1…4 | azimuthal shape change, amplitude ∝ r. ℓ=0 is `radialScaling` |
| `tangentialHarmonic<ℓ>` (2 per ℓ) | rΔφ = ε·r·cos ℓφ, ε·r·sin ℓφ, ℓ = 1…4 | ℓ=0 is `curl` |
| `globalTilt` (`_x`, `_y`) | rigid rotation of the branch about global X, Y | transverse partners of `curl` |
| `radialPolar_r<p>t<k>` (5) | Δr = ε·r<sup>p</sup>·t<sup>k</sup>, (p,k) = (0,0), (0,1), (0,2), (1,1), (1,2) | depends on the polar angle |
| `tangentialPolar_r<p>t<k>` (5) | rΔφ = ε·r<sup>p</sup>·t<sup>k</sup>, same (p,k) | |
| `longitudinalPolar_r<p>t<k>` (4) | Δz = ε·r<sup>p</sup>·t<sup>k</sup>, (0,1), (0,2), (1,0), (1,2) | (1,0) is the telescope mode Δz = ε·r |

The polar families: a track from the vertex is at z ≈ r·t at radius r, so a field r·f(t) (p=1) is absorbed
by the track's φ, tgλ or curvature whatever the vertex does, while f(t) (p=0) shifts d0/z0 and is broken
only by the vertex constraint. Members that coincide with another mode or with a rigid shift are left out.

Beyond this list the clause `rayModes` constrains the whole function space r·F(φ, t) of which the
harmonic and polar families are low-order members (per (φ, t) cell, the common angular displacement of
the sensors of all layers). It is a constraint only, not a column of the diagnostics.

### To which volumes it is applied

- **Mode constraints** are written over the free rigid-body DOFs of the volumes *flagged* by the
  `globalModes` clause. Flagging nested levels is legitimate: the equation constrains the total amplitude,
  whichever level carries it. In the ITS the flag typically sits on the **stave level and below**
  (`"match": "ITS/.+"`), i.e. over one whole detector branch, so that the equation spans the inner and
  the outer barrel together. The root has no DOFs, so each
  detector branch (ITS, TRD, TOF, TPC, …) is a separate, independent set.
- **Modes with a z-dependent coefficient** (`zScaling`, `bowing`, the polar families) project on a
  volume's local Z with a coefficient ∝ its z; they are meaningful only once DOFs *distributed in z*
  (modules, chips) are free. At stave level (centred at z = 0) they are skipped with a warning.
- `rayModes` is flagged on the detector envelope (`ITS`).
- **The diagnostics** are computed per detector branch (child of the root) with at least `NPar` = 6 +
  `NGlobalEq` sensors.

## What the diagnostic computes

Two numbers per branch and mode, deliberately independent:

- **Metric A — what the constraint sees.** `Σ cᵢ·pᵢ / Σ|cᵢ|`, over the free rigid-body DOFs of the
  branch: cᵢ is the coefficient the DOF would have in the mode equation of
  `writeGlobalModeConstraints`, pᵢ its fitted value (rotations enter with the lever arm of the mode
  field). In µm. It is exactly the left-hand side of the constraint, normalised: it is ≈0 for a mode that
  was constrained in the fit, and it depends on which DOFs are free.
- **Metric B — what the geometry did.** The fitted corrections of the whole hierarchy are composed
  into a displacement of the **centre of every sensor**. These displacements are fitted by a global shift,
  a global rotation and the mode fields (SVD, normalised columns). The amplitude is quoted as the
  displacement the mode produces at its largest lever arm on the branch (µm). It does not depend on which
  hierarchy level carries the correction. Also given: the fit error, the fraction of the
  displacement power explained by each mode (after the rigid part), the RMS before/after removing the
  rigid part and the modes, and the mode-to-mode correlation matrix.

Caveats: the columns are not orthogonal (`curl` and `globalTilt` duplicate the rigid rotation, several
modes overlap on a limited set of sensors), so the fit is rank-deficient by construction (the rank is
logged and stored in the JSON) and neighbouring amplitudes can be correlated: check the correlation
matrix before reading a single bar. Values are absolute: an injected misalignment (`misAlgJson`) is
*not* subtracted.

## Standalone diagnostic of the initial (CCDB) alignment

`--output WeakModesDiagIni` diagnoses the weak modes already present in the *initial* alignment
itself — the `o2::detectors::AlignParam` objects already baked into the loaded geometry via
`GRPGeomHelper::getAlignment(detID)` (CCDB `<Det>/Calib/Align`) — relative to the ideal geometry.
It needs no millepede.res and runs independently of `MilleRes`/`WeakModesDiag`: it rectifies
`ProcessData` off like `MilleRes` does, so it only has to reach the first `run()` call of the DPL
workflow (the same single dummy-triggered pass the standalone `MilleSteer` mode uses) for the CCDB
alignment and the hierarchy to become available.

The displacement composed per sensor is exactly the one `Detector::MP2AlignParams` calls `P^p_j`
(see `Detector.cxx`), not the millepede fit's `P^n_j`: it is the cumulative global delta of the
initial `AlignParam` objects alone, applied to the ideal sensor-centre position. Only **Metric B**
is produced — there is no fit, hence no free-DOF values for Metric A, and `weakModes.json` simply
omits `metricA_um`/`nTermsA` for every mode (the ROOT `modes` tree keeps the branches for schema
compatibility with `macro/drawWeakModes.C`, filled with 0, which must not be read as "constrained").
Everything else (RMS, rank, rigid shift/rotation, explained fractions, correlations, the per-sensor
`sensors` tree) is produced exactly as for `diagnoseWeakModes`, written to the same
`AlignParams.weakModesOutBase` files — point it elsewhere via `--configKeyValues` to avoid
overwriting a fit-based `WeakModesDiag` run from the same output directory.

## How to read the histograms

`macro/drawWeakModes.C` produces, per detector branch, a PDF page set. The first page is the summary:

- **Left, blue bars with error bars: metric B**, amplitude of each mode in µm. Red line: metric A.
  A mode that was constrained in the fit has A ≈ 0 (the constraint holds) — B can still be nonzero
  because the same displacement field can be produced by DOFs that were not part of the equation.
- **Right, green bars: the fraction of the displacement power explained by each mode** (after removing
  the rigid shift/rotation). Bars summing to ≈1 mean that the displacement is made of these modes.
  The modes are correlated, so the individual bars are not strictly additive.

The following pages show, versus the variable the mode depends on, two sets of points per sensor
centre, in µm:

- **black**: the **total** displacement of the sensor centre (the composed effect of all fitted
  corrections, not only of the weak modes);
- **red**: the model of **this panel's mode** only, see below.

### How red is built

The fitted corrections enter only once, to build the black points: they are composed down the hierarchy
into one displacement vector **d**ᵢ = (dx, dy, dz) per sensor centre i, at position **x**ᵢ. Red is then
computed from a separate least-squares fit of these vectors, which uses only the sensor positions and the
black displacements:

**d**ᵢ ≈ **t** + **ω** × **x**ᵢ + Σₘ aₘ · **u**ₘ(**x**ᵢ)

with 3 + 3 + `NGlobalEq` (39) = 45 unknowns:

- **t**: a rigid shift of the whole branch (3 numbers).
- **ω**: a rigid rotation of the whole branch (3 numbers). The rotation moves a sensor at **x**ᵢ by
  **ω** × **x**ᵢ.
- aₘ: the amplitude of mode (equation) m, one free number per mode (39 numbers).
- **u**ₘ(**x**): the fixed displacement field of the mode of unit amplitude, given by the formulas of the
  table above and the same as in the constraint equations. For `radialScaling`, **u** = (x, y, 0); for
  `curl`, **u** = (−y, x, 0); for `elliptical_cos2phi`, **u** = cos2φ · r̂.

All the sensors of the branch and all three components enter one fit, solved by SVD with normalised
columns. Its result is the 45 fitted numbers; the aₘ are the amplitudes of metric B on the summary page.

The sum runs over all the modes because they are fitted together; fitting each mode alone would give
different amplitudes, since the modes are correlated. Once the fit is done, **red** of a panel keeps only
the mode m named in the panel title (for a family, the sum over its members):

**red**ᵢ = **t** + **ω** × **x**ᵢ + aₘ · **u**ₘ(**x**ᵢ)

with **t**, **ω** and aₘ taken from the full fit, projected on the component of the panel (Δr, rΔφ or
Δz). The amplitudes of all the other modes are left out of this expression, but not out of the fit that
determined aₘ. For example, on the
`radialScaling` panel red is Δr = (**t** + **ω** × **x**ᵢ)·r̂ + a · r: the small rigid offsets plus a
straight line in r whose slope is the fitted ε.

The pads of one page show different projections and different x variables of the same black cloud; only
red changes from pad to pad.

Red is one term of one simultaneous fit, not a comparison of that mode with the data, so it can be far
from black even when the data are well described. The modes are correlated and share the amplitude
between them. A typical case is `elliptical` (Δr = ε·cos2φ with an r-independent amplitude): if the true
shape is `radialHarmonic2` (Δr ∝ r·cos2φ), that mode takes most of the amplitude and red shows only the
small leftover of `elliptical`. Red also includes the rigid part, so it can follow black on a pad where the
mode itself contributes little. To rank the modes use the summary page (amplitudes and explained
fractions), not the closeness of red to black. The harmonic and polar families come next, and finally the
(z, φ) maps described below.

### The last page: (z, φ) maps

Two colour maps per branch: the **mean Δr** (left) and the **mean rΔφ** (right) of the sensor centres,
in µm (colour scale), in bins of z (x axis, cm) and φ (y axis, rad); white bins contain no sensors. Each
bin averages the displacement of all the sensors in it, over all layers. They show the displacement as
a function of position, with no mode model overlaid: a uniform colour is a rigid or axisymmetric
displacement, stripes along φ at fixed z are z-dependent (scaling, twist, bowing, polar) patterns, and
stripes along z at fixed φ are the azimuthal harmonics (ℓ = 1, 2, …, tilt). Any structure that is
not smooth, e.g. a step at the boundary of the two half-barrels (φ ≈ 0, ±π) or patches limited in z,
is local and is not described by any of the global modes. Since the bins mix layers of different r, the
displacement of a mode ∝ r is averaged over them.

**Do these histograms show the final sensor corrections?** The black points do: they are the final
fitted corrections expressed as the displacement of each sensor centre (sensor rotations are not
shown). The summary bars (metric B) and the red curves do not: they are only the projection of those
displacements onto a few smooth global fields. The rest (the local, sensor-to-sensor part) is the
residual of the fit, `rmsAfterModes_um`. The corrections themselves are in `result.json`, the
displacements in the `sensors` tree. What the numbers mean:

- A large amplitude with a small χ² cost is the signature of a *weak-mode artefact*, not of a real
  misalignment: the data cannot distinguish it, so its value carries no information on the geometry
  (a real deformation would not be coherent as ε = amplitude/r across the independently mounted inner
  and outer barrel).
- A large amplitude is a *harm* to the corrections even if the residuals look fine: the sensors are
  moved coherently, with no benefit in χ², and distortions of this kind (momentum scale, a vertex position
  shift) appear in the physics. Leave-one-out residuals (checkResid) are blind to them by construction.
- Metric A ≈ 0 and B ≈ 0 for a mode: the mode is constrained and nothing leaks through other DOFs.
  A ≈ 0 and B large: the mode is constrained in the DOFs flagged, but it is realised by DOFs outside
  of the equation (e.g. a different level, or a rotation not flagged).

Typical use: run MilleRes once for a fit without the `globalModes`/`rayModes` clauses to see which modes
the data leave free (these are the ones to constrain), and again with the clauses on to verify that A
vanishes and that B does not reappear in an unconstrained neighbour.

## Files

| | |
|---|---|
| `weakModes.json` | per branch: RMS, rank, rigid part, per-mode B ± error, explained fraction, lever arm, A, correlation matrix |
| `weakModes.root` | trees `sensors` (per-sensor displacement, rigid part, `uMode[3*m+k]` contribution of mode m) and `modes` |
| `macro/drawWeakModes.C` | `root -l -b -q 'drawWeakModes.C("weakModes.root", "weakModes.pdf")'`; the mode column layout is read from the file |

The prefix is `AlignParams.weakModesOutBase`.
