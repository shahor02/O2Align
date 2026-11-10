# Handover: detecting the IB/OB kink and building a geometric cure for it

Written for: whoever continues the ITS alignment (a person or Claude), familiar with O2Align, Millepede-II and checkResid.
Date: 2026-10-31. State: v19 (`/data/devAlg/solA2v19`, config `doc/ITS_RB_v19.json`) and its cure `cure012`
(`/data/devAlg/solA2v19cure/cure012`), with checkResid on pp B0 (LHC25ab 562862) and OO B+ (LHC25ae) for both.
Sections 6 (pede evaluation of the cure) and 6a (PV test) were added on 2026-11-01, 6b (B0 vs B+) on 2026-11-03, 6c (B−) and 6d (IB/half-barrel equivalents) on 2026-11-04, 6e (v18, field-dependent response) on 2026-11-05.

## 1. What the kink is

checkResid fits each 7-cluster track twice and compares, at R = 12 cm (between Lr2 and Lr3), the IB sub-track
(vertex + Lr0..Lr2, `trIBOut`) with the OB sub-track (Lr6..Lr3, `trOBInw`), both in the same frame:

- dPar0 = y_IB − y_OB [µm], dPar2 = snp_IB − snp_OB [µrad] (`trIBOut.mP[i] − trOBInw.mP[i]`).
- Cells: 8 tgλ bins (edges −2, −1.2, −0.8, −0.4, 0, 0.4, 0.8, 1.2, 2; tgλ of the track at the vertex, `track.mP[3]`)
  × 24 φ bins (φ of the IB sub-track position at 12 cm). In each cell the Gaussian-fitted mean is taken.
- Per tgλ bin, the φ harmonics n = 0, |n = 1|, n = 2 (cos, sin) of these means.

The "kink" is a dPar2 pattern with dPar0 ≈ 0: the two sub-tracks cross at R0 = 12 − dPar0/dPar2 ≈ 12 cm, so it is
purely directional, i.e. a relative IB/OB rotation about R = 12 cm along the ray. In v18 it was an n = 2 pattern; in v19
it is n = 1 and n = 2, both growing with |tgλ| (pp, µrad):

| tgλ bin | v19 n=0 | v19 \|n=1\| | v19 n=2 (cos, sin) | cure012 n=0 | cure012 \|n=1\| | cure012 n=2 (cos, sin) |
|---|---|---|---|---|---|---|
| [−2, −1.2] | −12 | 141 | (−118, 138) | 8 | 28 | (−18, 25) |
| [−1.2, −0.8] | −29 | 111 | (−69, 45) | 1 | 11 | (7, −17) |
| [−0.8, −0.4] | −8 | 84 | (−34, 3) | −5 | 14 | (4, −5) |
| [−0.4, 0] | 15 | 39 | (−14, −40) | 2 | 6 | (8, −4) |
| [0, 0.4] | 15 | 65 | (−18, −57) | −3 | 17 | (0, −4) |
| [0.4, 0.8] | −17 | 88 | (−42, −5) | −11 | 14 | (2, −2) |
| [0.8, 1.2] | −49 | 147 | (−88, 66) | −9 | 11 | (−6, −9) |
| [1.2, 2] | −119 | 223 | (−142, 181) | −7 | 34 | (−3, 32) |

dPar0 is 1–2 µm (n = 0) and 3–7 µm (|n = 1|) in both.

Why it matters and why it is not trivially cured by pede: single-layer residuals are almost blind to it (it is a
relative rotation of two sub-detectors about a point between them), and the pede χ² prefers the configuration that
checkResid calls a kink. v18 → v19 (IB `tangentialHarmonic` gauge) removed the v18 n = 2 kink at the price of pede χ²
+51 M on the same data, and built the v19 one. See [history.md](../history.md). (The +51 M is the χ² definition of
the unfixed pede 05.01.06 used for v19 vs the fixed one used for v18; on one binary v18 and v19 are equal, see 6g.)

## 2. Detecting the kink: `macro/kinkMap.C`

Input: the checkResid `res` tree (`checkResid_<dataset>_AlgRS_PV1.root`, not the `_hman` file).

```bash
o2d    # O2PDPSuite environment
cd <a scratch dir>   # ACLiC writes .so/.d/.pcm next to a copy of the macros; copy macro/kink*.{C,h} macro/makeKinkCure.C there
root -l -b -q 'kinkMap.C+("v19=/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root,cure=/data/devAlg/solA2v19cure/cure012/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
root -l -b -q 'kinkMap.C+("v19OO=/data/devAlg/solA2v19/checkResid_LHC25ae_AlgRS_PV1.root", true)'   # B != 0: also dPar2 - k*dPar4
```

The argument is a comma-separated list of `label=file`. Optional arguments: `curvCorr` (bool) and `nmax` (number of
tracks per tree; the 1.5 GB pp tree takes 2–3 min in full).

How to read it:
- **Use pp B0 (562862) as the kink reference.** With B = 0 all tracks have the same assigned pT, so dPar2 is free of the
  curvature terms.
- **OO B+ dPar2 is dominated by the curvature difference** dPar4 (about −8 µrad per 10⁻³ (GeV/c)⁻¹); `curvCorr = true`
  removes the linear part. OO and pp do not agree on the kink: in v19 they differ by about (−170, +220) µrad in n = 2
  at large |tgλ|, and a geometry change moves both by the same amount (section 5), so no single OB geometry zeroes both.
- **The edge bins |tgλ| > 1.2 have less statistics but are not noise** here: the kink is largest there.
- The quick look without macros is the checkResid PDF page "IB-OB par2 vs tgl" / "vs phi", but its binning does not
  separate the φ harmonics per tgλ.

## 3. Response model: `macro/kinkResponse.C`

To design a cure one needs to know how dPar2 reacts to a geometry change. Model, per (tgλ, φ) cell:

  Δ dPar2 = Σ_L c_L d_L,

d_L [µm] being the mean apparent rφ shift of the layer-L hits of that cell when the geometry is changed (straight ray from
the origin through the hit, local-x crossing point in the new vs the old chip). `kinkResponse()` fits the 7 c_L over the
|tgλ| < 1.2 cells using several geometries against one reference, each with a checkResid tree on the same data:

```bash
root -l -b -q 'kinkResponse.C+("/data/devAlg/v10n/algdump.root", "/data/devAlg/v10n/checkResid_LHC25ab_562862_AlgV10n_PV1.root",
  "/data/devAlg/solA2v18/algdump.root,/data/devAlg/solA2v19/algdump.root,/data/devAlg/solA2v19rmc/algdump.root",
  "/data/devAlg/solA2v18/checkResid_LHC25ab_562862_AlgRS_PV1.root,/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root,/data/devAlg/solA2v19rmc/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
```

Geometry arguments are aligned geometry files (`o2sim_geometry-aligned.root`) or algDump outputs (`v10n` has only the
latter). Result for pp B0 (V10n reference; v18, v19, v19rmc):

- dPar2: c_L = −33.9, −42.3, 75.7, 20.7, −21.6, 1.2, 0.03 µrad/µm (Lr0..Lr6); χ² over 432 cells 4887 → 93.
- dPar0: c_L = −2.6, −3.5, 7.0, −2.6, 1.8, −0.2, 0.07 µm/µm; χ² 5796 → 225.
- For an OB rφ shift ω·(r − 12 cm) the model gives Δ dPar2 = −0.857 ω (the cure gain).

`kinkPredict(solGeom, newGeom, solRes)` (same file, `.L kinkResponse.C+` first) applies these c_L to predict the dPar2
harmonics of a modified geometry from those measured on the solution, without running checkResid:

```bash
root -l -b -q -e '.L kinkResponse.C+' -e 'kinkPredict("/data/devAlg/solA2v19/o2sim_geometry-aligned.root", "/data/devAlg/solA2v19cure/cure012/o2sim_geometry-aligned.root", "/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root")'
```

For cure012 it predicted |n = 2| 182/82/35/42/60/42/110/230 → 32/3/2/5/11/2/6/41 µrad; measured (section 1)
31/18/6/9/4/3/11/32. The model is fitted on pp B0 only; it does not describe OO.

## 4. The cure and how it was produced

**Model of the cure.** An OB-only rotation about R = 12 cm along each ray: every OB chip (Lr3..Lr6) gets the rφ shift

  d(φ, t) = ω(φ, t) · (r − 12 cm),  ω = b(φ, t) / 0.855,  t = z/r of the chip centre,

b being the v19 pp dPar2: its φ harmonics n = nmin..nmax per tgλ bin, linear in t between the bin centres and constant
beyond the outermost ones. The shift is written as a change of the chip TX (local x): TX += d / (x̂·φ̂). Only free OB chip
TX labels are touched; everything else in `millepede.res` is copied unchanged. Variants:
- `cure2`: n = 2 only; rms shift Lr3..6 = 9, 12, 26, 28 µm, max 74 µm.
- `cure012`: n = 0, 1, 2; rms 13, 19, 39, 44 µm, max 125 µm.

**Steps** (all in the O2 environment):

1. Build the modified `millepede.res` with `macro/makeKinkCure.C` from the solution's pp tree, aligned geometry,
   `mp2param.txt` (to map the TX labels to chips through their `! V TX <symname>` comments) and `millepede.res`:
   ```bash
   root -l -b -q 'makeKinkCure.C+("/data/devAlg/solA2v19/checkResid_LHC25ab_562862_AlgRS_PV1.root", "/data/devAlg/solA2v19/o2sim_geometry-aligned.root",
     "/data/devAlg/solA2v19/mp2param.txt", "/data/devAlg/solA2v19/millepede.res", "millepede_cure012.res", 0, 2)'
   ```
   (`nmin, nmax` = 2, 2 for cure2; last argument = gain, default 0.855). It prints the harmonics used, the number of
   labels modified (23 688) and the rms shift per layer.
2. Put it as `millepede.res` in a new directory and run [doc/resToCCDB.sh](../resToCCDB.sh) there with the DOF config
   of the solution:
   ```bash
   mkdir cure012 && cp millepede_cure012.res cure012/millepede.res && cd cure012
   Config=~/alice/O2Align/doc/ITS_RB_v19.json ~/alice/O2Align/doc/resToCCDB.sh
   ```
   It runs `--output MilleRes` (via `doc/proc0.sh`; `result.json`, `ITS_alignment.root`), creates the aligned geometry
   with `o2-create-aligned-geometry-workflow` and fills the local CCDB layout `ITS/Calib/Align`, `GLO/Config/GeometryAligned`,
   `ITS/Config/Geometry`. `MVJSON` (mean-vertex time slots, default `/data/devAlg/MV_B0BPBM.json`) and `TimeStamp`
   (default 1746914681240) must be those of the solution.
3. Run checkResid on the new directory as for any solution, i.e. with
   `CCDB="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned"`.

**Checks done for cure012.**
- MilleRes on the unmodified v19 `millepede.res` (`/data/devAlg/solA2v19cure/v19check`) reproduces
  `solA2v19/ITS_alignment.root` exactly, so the chain adds nothing of its own.
- Compared with v19 only OB chip translations change; IB and all rotations are identical.
- `makeKinkCure.C` reproduces the original scratch-macro output (`/data/devAlg/solA2v19cure/millepede_cure012.res`)
  exactly.

## 5. Result of the cure (checkResid)

- pp B0: dPar2 n = 2 drops from 39–223 to 6–34 µrad, n = 0 to ≤ 11, |n = 1| to 6–34 (table of section 1). dPar0 is
  unchanged. Single-layer residuals improve slightly (rms of the dy vs φ means over all layers 2.5 → 2.2 µm, dy vs Z
  1.7 → 1.4 µm); dz is unchanged in the core (the rise of the dz vs tgλ rms comes only from outermost-tgλ bins whose fit
  failed in v19).
- OO B+: dPar2 n = 2 gets worse by about (+140, −140) µrad at |tgλ| > 0.8 (e.g. [−2, −1.2]: (49, −84) → (190, −223)),
  i.e. the same change as in pp. Single-layer residuals unchanged.
- Hence the pp − OO dPar2 difference is not geometric: either a field-dependent effect (curvature, field map, material)
  in the B ≠ 0 data, or a real difference of the detector between the two periods.

## 6. Does the Mille fit reject the cure? (pede evaluation)

pede on the v19 Mille data with the ITS labels fixed, once at the v19 values (test 2, `runTest2.sh`) and once at the
cure012 values (test 1, `runTest1.sh`), each on B0, B+, B− separately. Inputs and outputs are in
`/data/devAlg/solA2v19cure/pedeEval` (`run_<sol>_fixedMVfree_<set>/`).

- Parameter files `param_{v19,cure012}_fixed.txt`: every label `label value -1`; they differ only in the 23 688 OB chip TX.
  Mean-vertex slot 6 (no data in v19, absent from `millepede.res`) is fixed at 0. pede refuses these
  (`STOP ... no variable global parameters`), so the runs use `param_*_fixedMVfree.txt` (`./runTestN.sh fixedMVfree`),
  where the 35 mean-vertex labels are free and refitted in both tests.
- Steering `steer_<param>_<set>.txt`: the v19 options (`fullMINRES-QLP 4 0.10`, `chisqcut 20 5`, `outlierdownweighting 4`,
  `dwfractioncut 0.15`, `entries 55`, `skipemptycons`) and v19's `mp2con.txt`. The iterations move only the mean vertex;
  they reproduce v19's outlier treatment. The value to compare is the last `Sum(Chi^2)/Sum(Ndf)` line of `millepede.log`,
  which is the χ² sum.

Result (2026-11-01):

| set | records | χ² v19 | χ² cure012 | Δχ² | Δχ²/record | rejects v19 → cure |
|---|---|---|---|---|---|---|
| B0 | 4 058 920 | 99 885 863 | 100 089 713 | +203 850 (+0.20 %) | 0.050 | 88 148 → 88 185 |
| B+ | 4 873 246 | 185 250 781 | 185 513 432 | +262 651 (+0.14 %) | 0.054 | 137 978 → 138 079 |
| B− | 5 073 458 | 195 322 229 | 195 585 543 | +263 313 (+0.13 %) | 0.052 | 148 447 → 148 659 |
| sum | 14 005 624 | 480 458 873 | 481 188 688 | +729 815 (+0.15 %) | 0.052 | |

- The v19 sum is within 0.02 % of the χ² of the v19 fit (4.8036E+08), as expected.
- The Mille data disfavour the cure, but by 0.05 per record (the +51 M v18 → v19 gap once quoted here is a pede-binary
  effect, 6g). The value is an
  upper bound: all other labels stay at v19, none is allowed to readjust to the cure.
- The cost per record is the same in B0, B+ and B−, so the Mille preference for v19 is not field-related. Hence pp B0
  data prefer cure012 in checkResid and v19 in the Mille fit.
- The χ² is however almost blind to the kink. The χ² the kink itself carries is Σ (A/σ_A)² of its harmonics: for v19 pp B0
  562862 (2.9 M tracks) the n = 2 terms give ≈ 1.3k, i.e. ~4·10⁻⁴ per track. Removing a kink seen equally by the Mille
  tracks would gain O(1k) on the 4 M B0 records, 100–200 times less than the +204k measured. The pede Δχ² is therefore set
  by other properties of the cure (it is applied per chip, from the chip-centre phi and tgl), and does not test the kink.

## 6a. Is the checkResid kink set by the vertex? No

The checkResid IB sub-track contains the refitted PV as a point (`checkresid.addPVAsCluster = true`, the default), which
could set the IB sub-track direction. In the Mille data the vertex enters only through the composed free-vertex fit of PVs
with ≤ 15 tracks (`useMultiTrackPVConstraint = 2`). Test (2026-11-01): checkResid on pp B0 562862 with and without the PV
point (`_PV1` / `_PV0` files in `/data/devAlg/v10n`, `solA2v19`, `solA2v19cure/cure012`). `kinkMap.C` gives the same
dPar0 and dPar2 for PV0 and PV1 in every tgl bin, within 0.5 µm and 3 µrad, for all three alignments, e.g. v19
[+1.2,+2.0]: dPar2 n2 (−134, 186) with the PV, (−135, 184) without; cure012 stays at ≤ 33 µrad. The refitted PV is
made of the same tracks, so on average it lies on their IB extrapolation and does not bias the IB direction. The kink
is therefore in the cluster positions, and the checkResid/Mille disagreement of section 6 is not caused by the vertex.
The test is sensitive: matching tracks by gid, the PV changes the IB snp at 12 cm by 167 µrad rms per track (OB: 0.9 µrad),
so any vertex bias relative to the IB would move the PV1 means away from the PV0 ones; in PV0 the IB sub-track is the three
IB clusters only. The mean vertex dy in the same cells has |n2| < 0.8 µm.

## 6b. B0 vs B+ of the same period: the kink is B0-only

checkResid on pp B+ 562837 (LHC25ab, next to B0 562862, both used for the alignment; 0.9 M IB/OB tracks) for V10n, v19 and
cure012: `checkResid_LHC25ab_562837_*_PV1.root` in the same directories (2026-11-03). In B ≠ 0 the plain dPar2 is
dominated by the IB sub-track curvature error (cell errors ~250 µrad vs 25 µrad in B0); `kinkMap.C` with `curvCorr = true`
also gives dPar2 with both sub-tracks constrained to their common q/pT (`kink::DPar2QC`), cell errors 75–140 µrad.
dPar2 n = 2 (cos, sin) [µrad], q/pT-constrained for B+:

| tgl bin | v19 B0 (±6–10) | v19 B+ | cure012 B+ |
|---|---|---|---|
| [−2.0,−1.2] | (−125, 138) | (9, 119) ± 41 | (73, 80) ± 42 |
| [−1.2,−0.8] | (−62, 29) | (41, −27) ± 30 | (77, −26) ± 30 |
| [−0.8,−0.4] | (−32, −2) | (−20, 46) ± 25 | (22, 37) ± 24 |
| [−0.4, 0.0] | (−12, −42) | (−22, 23) ± 21 | (9, 82) ± 21 |
| [ 0.0,+0.4] | (−19, −53) | (17, 35) ± 21 | (46, 82) ± 21 |
| [+0.4,+0.8] | (−44, −5) | (29, 34) ± 25 | (62, 41) ± 25 |
| [+0.8,+1.2] | (−88, 56) | (−25, 72) ± 31 | (86, 12) ± 31 |
| [+1.2,+2.0] | (−134, 186) | (28, 65) ± 38 | (110, −64) ± 40 |

- v19 B+ shows no B0-like kink: the B0 pattern (negative cos everywhere, large at |tgl| > 0.8) is absent, the two differ by
  3–4σ in several bins ([+1.2,+2.0] cos, central sin).
- The cure moves B+ in the same direction as B0 (cos up by 30–110 µrad), as a geometric change must, and makes B+ worse.
  So in B+ checkResid prefers v19, like the Mille fit.
- The geometry is the same and the runs are adjacent, so the B0 kink is not an alignment error of v19; the cure should not
  be used. A straight-line fit of unbiased clusters gives no mean kink, and the vertex is excluded (section 6a), so the
  IB-only and OB-only cluster sets of B0 tracks disagree. Candidates: an ITS movement between magnet off and on (v19 is
  fitted on 10 M B± vs 4 M B0 records, and pede is blind to the kink), or a field-dependent cluster position effect
  (Lorentz-like shift, odd in B; incidence angle of radial vs bent tracks). A B− run of the same period separates them:
  B− like B+ means field on vs off, B− − B0 opposite to B+ − B0 means an effect odd in B.
- Errors: B+ statistics are ~3× lower and the per-track error ~3× larger than B0; ~4× more B+ tracks would make the
  B0/B+ difference > 6σ per bin. A B− run of the same period would tell whether the B ≠ 0 residual pattern changes with
  the field sign.

## 6c. B+ and B−, more statistics: v19 is nearly flat at B ≠ 0

2026-11-04: B+ 562837 with ~3× the statistics of 6b, and pp B− LHC25ah 564644 (another period), for V10n, v19, cure012
(`checkResid_LHC25ah_564644_*_PV1.root`). (The B+ and B− columns of the first version of this table were swapped; corrected
2026-11-05.) v19 dPar2 n = 2 (cos, sin) [µrad], B± q/pT-constrained:

| tgl bin | B0 562862 (±6–10) | B+ 562837 | B− 564644 |
|---|---|---|---|
| [−2.0,−1.2] | (−125, 138) | (19, 97) ± 23 | (34, 24) ± 14 |
| [−1.2,−0.8] | (−62, 29) | (40, 27) ± 18 | (4, −8) ± 11 |
| [−0.8,−0.4] | (−32, −2) | (−8, 30) ± 14 | (1, −2) ± 9 |
| [−0.4, 0.0] | (−12, −42) | (−21, 22) ± 12 | (0, −1) ± 8 |
| [ 0.0,+0.4] | (−19, −53) | (7, 24) ± 12 | (13, 3) ± 8 |
| [+0.4,+0.8] | (−44, −5) | (31, 18) ± 14 | (−5, 25) ± 9 |
| [+0.8,+1.2] | (−88, 56) | (8, 48) ± 18 | (5, 26) ± 11 |
| [+1.2,+2.0] | (−134, 186) | (−9, 71) ± 23 | (37, 43) ± 15 |

- v19 is flat in B− (≤ 3σ) and nearly so in B+ (sin +20…+100); only B0 has the large kink.
- cure − v19 is the same in B+ and B− (outer bins ≈ (+105, −100) and (+120, −145)) and in B0: the cure is linear in r,
  and such a pattern has the same dPar2 response in all fields (section 6e).
- This does not make the B0 kink non-geometric: section 6e shows that geometric patterns that are not linear in r
  respond differently at B = 0 and B ≠ 0.

## 6d. Can the cure be expressed as IB (half-barrel) motions?

Written as a phase, b ∝ cos 2(φ − ψ) per tgl bin.

- **Rigid half-barrel motions do not reproduce it.**
  - A translation of a half-barrel (IB or OB) moves all its layers equally in rφ. dPar2 hardly responds: Σc_L = −0.5
    (IB) and +0.3 (OB) µrad/µm. The shift goes to dPar0, while the kink is directional (dPar0 ≈ 0, R0 ≈ 12 cm).
  - Opposite rotations of the top and bottom halves about z give rΔφ ∝ r·sign(sin φ), which has only odd harmonics.
  - Half-barrel tilts (RX, RY) are linear in z, i.e. odd in tgl, while the kink is mostly even in tgl at large |tgl|.
- **Phase versus the half-barrel structure.** The halves meet at φ = 0, π. A deformation that is the same in the top
  and bottom halves, mirrored, has rΔφ odd under φ → −φ: pure sin 2φ (2ψ = ±90°). cos 2φ is antisymmetric between the
  halves; any n = 2 term is invariant under a 180° rotation. v19 B0, n = 2 (cos, sin) [µrad] and ψ:

  | tgl bin | (cos, sin) | ψ |
  |---|---|---|
  | [−2.0,−1.2] | (−125, 138) | 66° |
  | [−1.2,−0.8] | (−62, 29) | 77° |
  | [−0.4, 0.0] | (−12, −42) | 127° |
  | [ 0.0,+0.4] | (−19, −53) | 125° |
  | [+0.8,+1.2] | (−88, 56) | 74° |
  | [+1.2,+2.0] | (−134, 186) | 63° |

  - Central bins: ψ ≈ 125°, close to a mirror-symmetric deformation.
  - |tgl| > 0.8, where most of the kink is: ψ ≈ 63–77°, mixed symmetric and antisymmetric, not aligned with the split.
  - The field-independent B0 − B+ difference has the same phase at the outer bins (≈ 61° V10n, 63° v19).
  - |n1| is as large as |n2| (90–240 µrad).
  - So the kink is not a clean half-barrel deformation, as a field-on/off movement of the half-barrels would suggest.
- **IB-side equivalent and the role of r.** The IB spans 2.2–3.8 cm, so an IB displacement acts through two numbers per
  (φ, tgl) cell: its mean, which moves the IB sub-track at 12 cm (dPar0), and its slope across Lr0..Lr2 (dPar2, gain
  Σc_L·r_L·10⁻² = 0.83 µrad per µrad).
  - Since dPar0 n2 ≈ 0, the IB equivalent of the OB cure is a rotation about R = 12 cm: rΔφ = −ω(r − 12).
  - For ω = 150 µrad: about 13 µm common to Lr0..Lr2 (±1 µm), and only 2.3 µm differential between Lr0 and Lr2. The
    differential cures dPar2; the common part keeps dPar0 unchanged.
  - `tangentialHarmonic2` (rΔφ = ε·r, pivot at the origin) is not this mode: it moves the IB sub-track at 12 cm by 12ε
    (≈ 18 µm), i.e. it creates a dPar0 n2 that is not observed.
- **Half-barrel content of the IB equivalent.**
  - The common 13 µm is an r-independent n = 2 tangential shift. Relative IB/OB half-barrel translations produce it:
    top/bottom separation in y gives ≈ 0.85·Δy·sin 2φ; opposite sliding in x gives cos 2φ plus an n = 0 term.
  - The differential needs per-layer IB half-barrel shifts in the ratio (12 − r_L) ≈ 9.8 : 9.0 : 8.2, about 1 µm per
    layer step (per-layer half-barrel DOFs).
  - The kink sits at |tgl| > 0.8, i.e. |z| ≈ 3–8 cm on the IB, so the shifts would also have to vary along z: a
    half-barrel bending, not a rigid motion.
- The kink appears only at B = 0 (section 6c), so this geometric reading is given only for completeness.

## 6e. v18 at B0, B+, B−: the dPar2 response depends on the field

2026-11-05: checkResid of v18 on the same three runs (`/data/devAlg/solA2v18/checkResid_{LHC25ab_562862,LHC25ab_562837,
LHC25ah_564644}_AlgRS_PV1.root`). dPar2 n = 2 (cos, sin) [µrad], B± q/pT-constrained:

| tgl bin | v18 B0 (±6–10) | v18 B+ | v18 B− |
|---|---|---|---|
| [−2.0,−1.2] | (−123, 161) | (5, 106) ± 24 | (10, 31) ± 14 |
| [−0.4, 0.0] | (−102, 91) | (4, −5) ± 12 | (18, −33) ± 8 |
| [ 0.0,+0.4] | (−115, 67) | (30, 1) ± 12 | (30, −21) ± 8 |
| [+1.2,+2.0] | (−204, 73) | (−7, 122) ± 23 | (26, 78) ± 15 |

v18 has a large B0 kink in every tgl bin, but at B± it is as flat as v19. The v18 − v19 change is (−90, +133) µrad at
B0 in the central bins and only ≈ (+20, −30) at B±.

**The response differs by field.** `kinkResponse.C` with `qConstr = true` on B+ and B− (reference V10n; v18, v19,
cure012), and on B0 with the same set (cells |tgl| < 1.2; B+: χ² 3583 → 236 over 432 cells). c_L [µrad/µm]:

| field | Lr0 | Lr1 | Lr2 | Lr3 | Lr4 | Lr5 | Lr6 | OB gain |
|---|---|---|---|---|---|---|---|---|
| B0 | −33.2 | −42.4 | 75.7 | 20.7 | −21.5 | 1.2 | 0.0 | 0.845 |
| B+ | −51.3 | −28.6 | 79.9 | 45.0 | −43.0 | −31.5 | 29.6 | 0.944 |
| B− | −51.3 | −28.7 | 81.3 | 44.4 | −44.2 | −29.0 | 28.3 | 0.894 |

- B+ and B− agree. In B0 the OB sub-track direction at 12 cm is set by Lr3 − Lr4 only (pT is fixed at 0.6 GeV/c, so the
  multiple scattering weights favour the innermost points). At B ≠ 0 Lr5 − Lr6 also enters, and the common q/pT,
  dominated by the OB, absorbs curvature-like OB patterns.
- A pattern linear in r (rotation, as the cure) has almost the same response in all fields. A pattern that is not linear
  in r across the OB, such as a bow of Lr4/Lr5 against the Lr3–Lr6 line, gives a B0 kink but little B± dPar2.
- Folding the v18 − v19 apparent layer shifts (`kink::cellHitShifts`) with each c_L reproduces the measured changes:
  central bins predicted (−90, +133) at B0 and (+14, −32) at B±. The OB part of v18 − v19 reaches ~50 µm at Lr6 and is
  not linear in r (Lr4/Lr5 ~3 µm off the Lr3–Lr6 line); the IB part is −1.5 … −2.4 µm, growing with r.
- **dPar4 at B±** ((q/pT)_IB − (q/pT)_OB, n = 2 sin): v18 − v19 is +3…+7 ·10⁻³ (GeV/c)⁻¹ in B+ and −2…−7 ·10⁻³ in B−.
  It flips sign with the field, as a geometric curvature-like change must. cure − v19 is 0, as for a pattern linear in r.
  1 µm of Lr4/Lr5 bow gives ≈ 1.8 ·10⁻³ (GeV/c)⁻¹ at 0.5 T.

Consequences:
- The B0 kink can be geometric: an OB rφ pattern that is not linear in r, which B± tracks see mostly as a q/pT bias
  (weak mode for B ≠ 0) and B0 tracks as an IB/OB kink. The conclusions "field-off only, non-geometric" of the first
  versions of 6b/6c are withdrawn.
- cure012, linear in r, is the wrong shape: it removes the B0 kink by creating one at B±. A cure that is consistent with
  all fields must lie in the B±-blind subspace: per cell, OB layer shifts with c_B0·d = −b, c_B±·d = 0, and no dPar0
  change, then checked against dPar4 at B±.
- The Mille χ² hardly constrains this (section 6). B0 data are the lever on it; the v18 → v19 choice was made on B0
  checkResid, where v19 is better than v18 at |tgl| < 0.8.

## 6f. B = 0 cosmics as an independent test (pede, labels fixed)

2026-11-06: B = 0 cosmics of LHC25aa (`/data0/devAlg/tests/cosmNoGold`, 104 runs, `cosmiGoldOnly=false`, v10n geometry,
`doc/ITS_RB_all.json`, assigned pT 2 GeV/c; ~272k accepted records, ndf 2.07 M, no mean-vertex labels). Setup in
`/data0/devAlg/tests/cosmEval`: `mkCosmEval.sh` writes one parameter file per solution with every ITS label fixed at its
millepede.res value (v10n = 0), except ITS TX (to which straight cosmics are blind), free with presigma 1e-5 cm so that
pede runs (it stays at 1e-12); no constraints. `runCosmEval.sh` runs the variants `cut` (chisqcut 20 5,
outlierdownweighting 4, as in production) and `nocut` (one pass); `sumCosmEval.sh` writes `summary.txt`.

| solution | χ² (cut) | rejected | Δχ² vs v10n (cut) | Δχ² vs v10n (nocut) |
|---|---|---|---|---|
| v10n | 3 292 086 | 58 979 | 0 | 0 |
| v11 | 3 165 976 | 56 060 | −126 110 | −64 797 |
| v16 | 3 209 408 | 59 009 | −82 678 | +280 280 |
| v17N | 3 160 662 | 55 874 | −131 424 | −192 464 |
| v18 | 3 160 806 | 55 806 | −131 280 | −193 834 |
| v19 | 3 164 866 | 55 928 | −127 219 | −202 434 |
| v19rmc | 3 161 995 | 56 323 | −130 091 | −102 481 |
| cure012 | 3 151 035 | 55 386 | −141 051 | −234 760 |

- All Mille solutions improve the cosmics by ~4 % of χ², an independent confirmation, since no cosmics entered the fits. v16 (module RY
  fixed, OB per-layer bow) is the worst of them and is worse than v10n without the cut.
- **cure012 is better than v19 on the cosmics**: Δχ² = −13.8k (cut, with 540 fewer records rejected) and −32.3k (nocut).
  The collision pede gives the opposite sign (+0.73 M, section 6), but per record the cosmic preference is 4–5 times
  larger: about −0.05 (cut) and −0.12 (nocut) vs +0.05 for collisions. B = 0 non-radial tracks are the data most sensitive to a
  ray-wise OB rotation.
- v17N, v18, v19 and v19rmc are within 4k (cut); without the cut the order changes (v19 −8.6k vs v18, v19rmc +91k vs v18), so the
  tails dominate the nocut numbers.
- This does not settle the field dependence of section 6e (cure012 creates a B± kink). It says the cosmics carry a
  B = 0 lever on the OB ray-wise patterns that pede should get in the combined fit.

## 6g. Combined fit B0 + B+ + B− + cosmics, v18 configuration (`solC1v18`)

2026-11-06: `/data0/devAlg/tests/solC1v18`: v18 DOF configuration, the 72 collision files of `solA2v18` plus the 104
cosmic files of 6f, all at weight 1; pede on epn224 with the outlier fix reverted (banner 05.01.06). Checks: checkResid of pp B0 562862, B+ 562837, B− 564644,
and the cosmic evaluation of 6f with the C1v18 labels.

**pede χ² and the pede binaries.** Fit χ²: v18 428 995 802 (ndf 477.44 M, 0.8985), v19 480 360 488 (1.0061), C1v18
483 866 510 (479.51 M, 1.0091; cosmic part ≈ 3.0 M). v18 was fitted on epn000 with a pede carrying the outlier
down-weighting fix (b91d586; its `$Id: cb3f0945…$` is that patched pede.f90, not a public blob). v19, v19rmc, C1v18
and the pedeEval runs of section 6 used the epn224 pede (banner 05.01.06), on which the fix had been reverted.

Fixed-label evaluation (ITS fixed, mean vertex free, production cuts), local pede 05.01.07, file `mp2data_0.bin` of each
field (934k records, 6.7 % of the data), `solC1v18/collEval`:

| file | v18 | v19 | C1v18 |
|---|---|---|---|
| B0 (270 539 rec.) | 5 994 810 | 5 998 356 | 6 000 259 |
| B+ (325 691) | 11 099 518 | 11 099 629 | 11 106 145 |
| B− (337 891) | 11 693 808 | 11 693 580 | 11 700 610 |
| sum | 28 788 136 | 28 791 566 (+3.4k) | 28 807 014 (+18.9k) |

- **v18 and v19 are equal on the collisions** (+0.004 per record, ≈ +51k if scaled to all data). The "v18 → v19 gap
  (+51 M)" of sections 1 and 6 is not a cost of the v19 IB gauge. It is the χ² definition of the two binaries: v19 χ²/ndf
  with 05.01.07 here (0.613, 1.026, 1.028) is 0.894 ± 0.003 times the epn224 05.01.06 evaluations of section 6 (0.683,
  1.148, 1.150), the same ratio as the v18 / v19 fit χ² (0.893). With `outlierdownweighting 4` the fixed pede sums Cauchy
  terms, which lie below the Huber terms of the unfixed pede at every residual.
- C1v18 costs the collisions +0.02 per record (≈ +0.28 M scaled) and gains 127k on the cosmics, so on the combined
  objective of the fixed pede it is *worse* than v18. The C1v18 fit (05.01.06, with the inconsistent-weight bug) did not
  find the minimum of the corrected objective.
- **Consequence: v19, v19rmc, C1v18 and the section 6 pedeEval were all run with the unfixed pede.** Their χ² are on the
  unfixed scale, and the v19/v19rmc/C1v18 solutions may carry the bias that the fix removes. The ±0.05 per record
  comparisons within section 6 are internally consistent (same binary), but the fits should be redone with ≥ 05.01.07.
- The `measurement` priors give 45k (v18), 50k (C1v18); constraints are satisfied to 4e-7.

Cosmic evaluation (6f setup, same table):

| solution | χ² (cut) | rejected | Δχ² vs v10n (cut) | Δχ² vs v10n (nocut) |
|---|---|---|---|---|
| v18 | 3 160 806 | 55 806 | −131 280 | −193 834 |
| cure012 | 3 151 035 | 55 386 | −141 051 | −234 760 |
| C1v18 | 3 034 119 | 49 486 | −257 967 | −646 318 |

The cosmics gain 127k (cut) / 453k (nocut) over v18, and the collisions lose nothing visible. Geometry change C1v18 − v18
(aligned chip centres):
- IB: a rigid transverse shift of all three layers by (Δx, Δy) ≈ (−2.2, +3.7) µm, with the mean vertex following
  (ΔX, ΔY rms 2.4, 4.3 µm, ΔZ 22 µm). Otherwise < 0.2 µm.
- OB: no harmonic (n ≤ 2, × {1, z}) pattern above ~1–2 µm. The median |Δrφ| is 0.7–2.2 µm. The rms (Lr3 16, Lr4 8,
  Lr5 18, Lr6 7 µm) comes from tails: 184 chips move by > 60 µm, up to 370 µm. These are stave-end chips at |z| ≈ 37–41
  cm (Lr3, modules 0/3) and 63–72 cm (Lr5), almost outside the collision acceptance (5 of 400k pp B0 tracks hit Lr3 at
  |z| > 36 cm). They are already large in v18 (76 µm rms) and only weakly determined by either sample.

Kink (kinkMap, 600k tracks): the C1v18 pp B0 dPar2 n = 2 is the v18 one within errors, e.g. (−105, +55) vs (−111, +60)
µrad at tgλ ∈ [0, 0.4], R0 ≈ 11–12 cm. At B± (q/pT-constrained dPar2) both are consistent with 0, as in 6e. **The
cosmics in the fit do not remove the B0 kink.** Their gain comes from directions orthogonal to it, consistent with
the small cure012 − v18 difference on the cosmics (−10k).

## 6h. Refits with the fixed pede (05.01.07): v19 and C1v18

2026-11-07: v19 and C1v18 refitted on epn224 with Millepede 05.01.07, which has the fix. The inputs are those of the
first fits. Copies are in `/data/devAlg/solA2v19` and `solC1v18` (alternative PC); the earlier results there are
overwritten. `solC1v18/millepede_labeled.res` is still the one of the unfixed fit.

| fit | χ² | χ²/ndf | rejected |
|---|---|---|---|
| v18 (epn000, fixed) | 428 995 802 | 0.8985 | 375 052 |
| v19 (05.01.07) | 428 841 347 | 0.8982 | 374 591 |
| C1v18 (05.01.07, + cosmics) | 432 067 126 | 0.9011 | 424 867 |

- **v19 is 155k below v18 on the same data.** The IB `tangentialHarmonic` gauge costs nothing; the "+51 M" was the
  unfixed pede.
- Cosmic evaluation (6f setup) of the refits: C1v18 3 037 443 (cut) / 6 108 709 (nocut), v19 3 164 892 / 6 538 745. The
  unfixed C1v18 gave 3 034 119 and the unfixed v19 3 164 866, so the fix barely changes the cosmic side.
- Split of C1v18, against v18:
  - cosmics 3.037 M vs 3.161 M: **−123.4k**;
  - collisions 432.067 − 3.037 = 429.030 M vs 428.996 M: **+33.9k** (+0.002 per record);
  - combined 432.067 M vs 432.157 M: **−89.5k**.
- At weight 1 the two samples are compatible: the cosmics gain 123k in directions the collisions hardly constrain, at a
  collision cost of 34k.
- v19 against v18: −155k on collisions, +4k on cosmics. The cosmics are neutral on the v19 OB z-odd patterns.
- Next: a combined fit in the v19 configuration (C1v19, same 72 + 104 files, weight 1, 05.01.07).
- Geometry relative to v18 (aligned chip centres):
  - **C1v18:** IB rigid shift of 1.8 µm. OB medians 0.6–2 µm, plus the same stave-end tails as in 6g, no harmonic
    pattern ≥ 1 µm. The fixed-pede C1v18 is again essentially v18.
  - **v19:** IB about 5 µm (n = 1 and n = 2). The OB carries large z-odd rφ patterns: c1·z +16 / +18 / +52 / +53 µm
    and s2·z −10 / −12 / −25 / −33 µm at Lr3–Lr6 (amplitudes at the stave ends).
- Kink (kinkMap, pp B0 562862, 600k tracks), dPar2 n = 2 (cos, sin) in µrad:

| tgλ | v18 | C1v18 | v19 |
|---|---|---|---|
| [−0.8, −0.4] | (−82, 101) | (−75, 100) | (−14, −9) |
| [0, 0.4] | (−110, 58) | (−108, 55) | (−6, −64) |
| [0.8, 1.2] | (−185, 62) | (−186, 52) | (−101, 79) |
| [1.2, 2] | (−201, 79) | (−213, 84) | (−106, 204) |

  v19 keeps the behaviour of the unfixed v19: kink removed at |tgλ| < 0.8, not at the ends. C1v18 is v18. At B± (q/pT
  constrained) all three are consistent with 0 within 2–3σ.

## 6i. Combined fit in the v19 configuration (`solC1v19`)

2026-11-07: v19 configuration, the same inputs as C1v18 (B0 + B+ + B− + cosmics, weight 1), epn224 Millepede 05.01.07.
A copy is in `/data/devAlg/solC1v19` (alternative PC).

| fit | χ² | χ²/ndf | rejected |
|---|---|---|---|
| C1v18 (05.01.07) | 432 067 126 | 0.9011 | 424 867 |
| C1v19 (05.01.07) | 431 908 691 | 0.9007 | 423 206 |

- **C1v19 is 158.4k below C1v18.** NREC and ndf are identical (14 277 030 records), so the v19 gain on collisions
  (155k, 6h) survives with the cosmics added.
- Cosmic evaluation (6f setup): C1v19 3 013 496 (cut) / 6 025 601 (nocut), the lowest of all solutions; v19 alone has
  3 164 892 / 6 538 745.
- Split of C1v19 (collisions = total − cosmics = 428 895 195):

| against | cosmics | collisions | combined |
|---|---|---|---|
| v19 | −151.4k | +53.8k | −97.5k |
| v18 | −147.3k | −100.6k | −247.9k |
| C1v18 | −23.9k | −134.5k | −158.4k |

- The cosmics cost the collisions 54k in the v19 configuration, against 34k in the v18 one: about +0.004 per record, still
  compatible at weight 1. The cost is most likely the halved OB s1·z term (below).
- Geometry (aligned chip centres):
  - **IB:** within 0.3 µm of v19.
  - **OB, relative to v18:** keeps the v19 z-odd rφ pattern. The c1·z term is +14 / +17 / +48 / +50 µm at Lr3–Lr6
    (v19: +16 / +18 / +52 / +53) and s2·z is −10 / −11 / −29 / −31 µm.
  - **Cosmics vs v19:** they reduce the s1·z term about by half: −7 / −7 / −18 µm at Lr4–Lr6, against −14 / −16 / −38
    for v19. Near-vertical cosmics constrain exactly this pattern.
  - **Tails:** 432 chips move > 60 µm relative to v19 (284 on Lr5, 117 on Lr6), mostly stave-end modules, as in 6g.
- Kink (kinkMap, pp B0 562862, 600k tracks), dPar2 n = 2 (cos, sin) in µrad:

| tgλ | v19 | C1v19 |
|---|---|---|
| [−2, −1.2] | (−119, 123) | (−87, 107) |
| [−1.2, −0.8] | (−50, 34) | (−55, 45) |
| [−0.8, −0.4] | (−14, −9) | (−14, 6) |
| [0, 0.4] | (−6, −64) | (−9, −39) |
| [0.8, 1.2] | (−101, 79) | (−96, 90) |
| [1.2, 2] | (−106, 204) | (−115, 188) |

  Errors are 13–22 µrad. C1v19 equals v19 within about 1.5σ, so the kink at |tgλ| > 0.8 remains. The dPar2 n0
  term moves by up to +22 µrad at |tgλ| < 1.2. At B± (q/pT constrained, 562837 and 564644) all n = 2 terms are
  within 3σ of 0, as for v18, v19 and C1v18.

### Cosmics weighted ×3 (`solC1v19w3`)

Same input and configuration as C1v19, cosmic files with weight 3 (steering `mpSteerV12CCw3.txt`).
- Fit: converged, Sum(W·χ²) = 437 815 712, rejected 420 330 (C1v19: 423 206).
- cosmEval: cosmics 2 924 129 (cut, 44 855 rejected) / 5 790 055 (nocut). Collisions = 437 815 712 − 3 × 2 924 129 =
  429 043 325. Relative to C1v19 (both at weight 1):

| | cosmics | collisions | combined |
|---|---|---|---|
| C1v19w3 − C1v19 | −89.4k | +148.1k | +58.8k |

  The extra 89k gained on cosmics cost 148k on collisions (+0.012 per collision record, against +0.004 for C1v19 vs
  v19), so the two samples are in tension beyond weight 1. C1v19 stays the preferred solution.
- Geometry: IB within 0.2 µm of C1v19. OB rφ z-odd terms relative to v18, Lr3–Lr6:

| term | v19 | C1v19 | C1v19w3 |
|---|---|---|---|
| s1·z | −3 / −14 / −16 / −38 | 0 / −7 / −7 / −18 | +2 / −1 / +6 / 0 |
| c1·z | +16 / +18 / +52 / +53 | +14 / +17 / +48 / +50 | +12 / +15 / +41 / +43 |
| s2·z | −10 / −12 / −25 / −33 | −10 / −11 / −29 / −31 | −8 / −9 / −27 / −27 |

  The weight 3 removes the s1·z term completely and reduces c1·z by 15 %. The chip-level scatter grows: relative to v18
  the OB dz rms is 23 / 22 µm at Lr5 / Lr6 (C1v19: 13 / 13), dr 16 / 18 µm (9 / 11). 301 chips move > 60 µm
  relative to C1v19.
- Kink, pp B0 562862, dPar2 n = 2 (cos, sin) in µrad (errors as above):

| tgλ | C1v19 | C1v19w3 |
|---|---|---|
| [−2, −1.2] | (−87, 107) | (−60, 101) |
| [−1.2, −0.8] | (−55, 45) | (−53, 51) |
| [−0.8, −0.4] | (−14, 6) | (−20, 25) |
| [0, 0.4] | (−9, −39) | (−13, −28) |
| [0.8, 1.2] | (−96, 90) | (−99, 87) |
| [1.2, 2] | (−115, 188) | (−125, 163) |

  Within 1.3σ of C1v19: removing the s1·z term does not change the kink, so it is not what drives it. The |n1| amplitude
  drops at |tgλ| < 1.2 (e.g. 65 → 23 µrad at [−1.2, −0.8]). B± (q/pT constrained): all n = 2 terms within 3σ of 0
  except B− [1.2, 2] sin, 177 ± 56 µrad (C1v19: 87 ± 51).

### v20, cosmics ×3 (`solC1v20w3`)

`doc/ITS_RB_v20.json` = v19 plus OB chip TY, TZ (priors 30 µm; the end-chip rules repeat them) and a stat-weighted
module mean (`weightedMean` TX, TY, TZ, RY). Cosmics weighted 3 as in C1v19w3. The steering files were regenerated
without `AlignParams.statInFile`, so the `weightedMean` rule fell back to the equal split: the module means are
unweighted, as in v19 (`mp2con.txt` is identical to that of A2v20).
- Fit: Sum(W·χ²) = 432 973 990, rejected 388 056; vs C1v19w3 −4.84 M with 32k fewer records rejected.
- cosmEval: cosmics 1 887 592 (cut, 19 422 rejected) / 3 384 762 (nocut); C1v19w3: 2 924 129 (44 855) / 5 790 055.
  Collisions = 432 973 990 − 3 × 1 887 592 = 427 311 214. At weight 1:

| | cosmics | collisions | combined |
|---|---|---|---|
| C1v20w3 − C1v19w3 | −1 036.5k | −1 732.1k | −2 768.6k |
| C1v20w3 − C1v19 | −1 125.9k | −1 584.0k | −2 709.9k |

  The cosmics gain 35 % of their χ² (nocut −2.41 M) and reject 57 % fewer records. The ~47k new chip parameters are far
  too few to explain this as fitting freedom, so the cosmics support the chip moves, including the IP-projective
  component that collisions cannot see. The cosmics are in the fit, though: an independent check needs v20 without
  cosmics (A2v20) evaluated with cosmEval.
- OB dz vs tgλ (checkResid, fitted means): the sigmoid shrinks and sets in later in all three runs. pp B0 562862,
  tgλ ≈ +1.5: Lr3 +24 → +14, Lr5 +35 → +18 µm; tgλ ≈ −1.9: Lr5 −84 → −64. The Lr6 one (−25 at +1.6, leakage of Lr5
  into the Lr6 extrapolation) is gone.
- Chip motion vs C1v19w3 (chips with ≥ 200 tracks in the three checkResid runs):
  - Live chips move by 30–45 µm rms within their module, and chips with little or no data by 50–120 µm. Module means
    move by 8–18 µm.
  - The moves lie along the projective direction ΔZ = t0·ΔR, where t0 is the chip's mean tgλ (corr 0.86–0.94 at
    |t0| > 0.8). Tracks from the IP do not see this direction; only the z_PV spread, the cosmics and the priors fix it.
  - The combination the collision tracks measure, c = ΔZ − t0·ΔR, has an rms of 13–25 µm. Its stave-to-stave rms grows
    from 3 µm at t0 ≈ 0 to 25–40 µm at |t0| > 1.4. Only 10–27 % of its variance is common to all staves at the same z.
  - The common part is odd in t0, and is significant only at Lr5 and Lr6: +20 / +24 µm at t0 = −1.8, −32 / −36 µm at
    +1.8.
  - A single incidence bias f(tgλ) would be common to all staves. The stave-specific part therefore indicates geometry.
    The common odd part, of the same sign at Lr5 and Lr6, is what a cluster-z incidence bias would look like. Chip
    TY + TZ absorb any f(tgλ) locally, so they cannot tell the two origins apart.
- Geometry vs v18, OB rφ: c1·z +7 / +8 / +23 / +28, s1·z +10 / +9 / +26 / +26 (Lr3–Lr6). The s1·z term removed by
  C1v19w3 comes back with the sign opposite to v19, and c1·z halves. Chip dz rms vs v18 is 54 / 59 µm at Lr5 / Lr6
  (C1v19w3: 23 / 22).
- Kink, pp B0 562862, dPar2 (cos, sin) in µrad, C1v19w3 → C1v20w3:

| tgλ | C1v19w3 | C1v20w3 |
|---|---|---|
| [−2, −1.2] | (−60, 101) | (−71, 140) |
| [−1.2, −0.8] | (−53, 51) | (−57, 62) |
| [−0.8, −0.4] | (−20, 25) | (−15, 40) |
| [0, 0.4] | (−13, −28) | (−15, −33) |
| [0.8, 1.2] | (−99, 87) | (−112, 77) |
| [1.2, 2] | (−125, 163) | (−137, 140) |

  The n = 2 kink is unchanged, and |n1| roughly halves at |tgλ| > 0.8 (93 → 46 and 175 → 141 µrad). B± (q/pT
  constrained): all terms within 2.7σ of 0; B− [1.2, 2] sin 108 ± 49 µrad (C1v19w3: 177 ± 56).

### v20 without cosmics (`solA2v20`)

Same steering files as C1v20w3 (equal-split module means, see above), collisions only.
- Fit: Sum(χ²) = 427 056 594, rejected 367 783; vs A2v19 (428 841 347, 374 591) −1.78 M, close to the collision gain
  of C1v20w3 over C1v19w3 (−1.73 M).
- Chip moves vs v18, live chips (≥ 200 tracks), A2v20 / C1v20w3:
  - c = ΔZ − t0·ΔR: rms 13 / 14, 12 / 13, 16 / 19, 22 / 25 µm (Lr3–Lr6), correlation 0.88–0.92 between the two fits.
    The collisions alone fix it.
  - Along the ray, a = (ΔR + t0·ΔZ)/√(1 + t0²): rms 27 / 51, 30 / 47, 34 / 58, 37 / 63 µm, correlation 0.53–0.79.
    The difference C1v20w3 − A2v20 lies on the ray (corr(ΔZ, t0·ΔR) 0.96–0.98 at |t0| > 0.8). The cosmics roughly double
    the projective component, which the collisions fix only weakly (z_PV spread and the priors).
  - rφ: rms 8–38 µm, correlation 0.82–0.84.
- cosmEval: 2 819 555 (cut, 43 538 rejected) / 5 530 044 (nocut). At weight 1:

| | cosmics | collisions | combined |
|---|---|---|---|
| A2v20 − A2v19 | −345.3k | −1 784.8k | −2 130.1k |
| C1v20w3 − A2v20 | −932.0k | +254.6k | −677.3k |
| C1v19w3 − A2v19 (for comparison) | −240.8k | +202.0k | −38.8k |

  - Without seeing them, v20 improves the cosmics by 345k (11 %; nocut −1.01 M): an independent confirmation of the
    collision-determined part c.
  - This is 27 % of the cosmic gain of C1v20w3 over A2v19. The rest, 932k, comes from the cosmic-driven projective
    component. It costs the collisions 255k at cosmic weight 3, a ratio of 0.27, against 0.84 in v19: the chip TY/TZ
    freedom makes the two samples much more compatible.
  - The cosmics are needed to fix the projective component. Whether it is geometry or a cosmic-specific bias, e.g. of
    the cluster position at the cosmic incidence, is not decided by these numbers.
- Weak modes (`weakModes.json`, metric B, µm), A2v19 / C1v20w3 / A2v20: `elliptical_cos2phi` 9.2 / 15.5 / 21.6 and
  `radialHarmonic2_cos2phi` −11.9 / −19.4 / −27.8, `curl` −0.1 / −1.8 / 4.1 and `tangentialPolar_r0t0` 0.0 / 2.9 / −6.8.
  Each pair is correlated at −0.97 / −0.96 and cancels largely, so the bars cannot be read singly. Per-chip radial
  cos2φ amplitude, A2v20 − A2v19 (Lr3–Lr6): +4.9 / +5.1 / −0.8 / −4.3 µm, an ellipse of Lr3/4 against Lr6. It lies on
  the IP ray: along-ray part +4.2 / +5.8 / −1.1 / −3.9, collision-measured c ≤ 2.4 µm. C1v20w3 − C1v19w3 has the same
  Lr3/4 cos2φ (+6.2 / +5.8), so the cosmics support that part. A2v20 − C1v20w3 is again on the ray (a up to 6 µm, mostly
  sin2φ and Lr3/Lr6; c ≤ 2.5 µm). The layer rφ offsets change by ≤ 2 µm (the curl pair).

### v21: v20 + OB chip RX/RY/RZ, cosmics ×3 (`solC1v21w3`)

Steering also without `statInFile` (equal split). Results in `/data0/devAlg/tests/solC1v21w3`.
- Fit: Sum(W·χ²) = 429 696 760, rejected 377 195; vs C1v20w3 (432 973 990, 388 056) −3.28 M.
- cosmEval (cut): 1 767 358 (18 699 rejected). At weight 1:

| | cosmics | collisions | combined |
|---|---|---|---|
| C1v21w3 − C1v20w3 | −120.2k | −2 916.5k | −3 036.8k |

  The collision gain is larger than that of the whole v20 step (−1.73 M), and the collision part (424.39 M) is below
  the collisions-only A2v20 (427.06 M).
- Chip rotations, live chips, rms (max) in mrad: RX 1.0–1.2 (≤ 8), RY 0.15–0.8 (≤ 7.5), RZ 2.4–3.6 (≤ 17).
- checkResid (Gaussian core of the 2D histograms; C1v19w3 / C1v20w3 / C1v21w3, all three runs alike):
  - OB dz core σ: −1.0…−3.4 µm vs C1v20w3, e.g. Lr5 201.9 / 200.9 / 197.5 (562837). The rms over φ of the
    per-bin dz means roughly halves (Lr5 4.5 / 4.0 / 2.1, Lr3 2.2 / 1.9 / 1.2 µm).
  - dy: σ −0.05…−0.1 µm, means vs φ unchanged within ±0.3 µm. IB and vertex unchanged.
  - dz vs tgλ: unchanged at |tgλ| < 1.5. At the edges Lr5 improves slightly; Lr3 (tgλ → +2) and Lr6 (tgλ → −2,
    +42 vs +11 µm in 564644) partly return towards C1v19w3. Chip-level rotations do not change the stave-common sigmoid.
- Weak modes: rmsDisplacement 80.4 → 72.7 µm, no mode amplitude grows beyond the v20 range.

## 6j. IOB dz vs tgλ: OB vs (IB+PV) extrapolation, and the IB/OB rz-slope mismatch

Merged from `doc/obsoleteNotes/handover_OB_dz_tgl_IOB_07_11_26.md` (2026-11-07, solA2v19 trees), which has the full
tables. The IOB residual is the OB cluster vs the IB+PV sub-track extrapolation (and the IB cluster vs the OB sub-track),
as in `checkResid_PV1_iob_hman.pdf` p.4 (dz vs Z) and p.8 (dz vs tgλ).

**Established on 2026-11-07:**
- **B≠0 IOB OB panels are an estimator artefact.** The IB+PV sub-track barely measures curvature: median
  σ(q/pT)/|q/pT| is 0.37–0.60, against 0.10 for the OB sub-track. The extrapolated z is quadratic in C, so
  dz = −tgλ·σ_C²R³/24. At fixed pT, the Lr6 slope is about −50 µm per unit tgλ for |C_IB/C_OB − 1| < 0.3 and about −230
  above. That accounts for the 100–250 µm curves and for the deeper tgλ ≈ 1.2 dip in B± than in B0. Using the IOB residual
  at B≠0 needs the IB sub-track q/pT fixed to the full-track value in CheckResidSpec (not done).
- **B0 IOB is clean.** In A2v19, Lr3/Lr4 (−9/−13) match the IB−OB par1 + par3·(R − 12). The IB+PV line is steeper than the
  OB line, Δtgλ ≈ +7·10⁻⁵·tgλ. par3 has no curvature bias and is field-independent. It differs between V10n and MP2, so it is
  geometry. A common IB radial shift (~2 µm) or an IB/OB radial scale (7·10⁻⁵) does not also fit the vertex IOB value.
- **Rigid Lr5/Lr6 offsets** are limited to a few µm: the B0 IOB excess over the Lr3–Lr4 line is −5/−11, and the standard
  Lr5/Lr6 slopes (−2.7/−1.1) don't confirm even that.

**Response to v19w3–v21 (2026-11-09).** These are linear slopes for |tgλ| < 1, in µm (par3 in 10⁻⁶) per unit tgλ, from
±2σ Gaussian core fits of the per-run hman 2D histograms. The par3 values come out higher than the tree values of 11-07
(A2v19 B0 87 vs 71), so compare only within this table.

| B0 (562862) | IOB Vtx | IOB Lr3 | IOB Lr4 | IOB Lr5 | IOB Lr6 | IB−OB par1 | IB−OB par3 | std Lr5 |
|---|---|---|---|---|---|---|---|---|
| A2v19 | −9.1 | −8.8 | −12.8 | −26.2 | −33.8 | 3.8 | 87 | −2.8 |
| C1v19w3 | −11.0 | −9.0 | −13.5 | −27.1 | −33.8 | 3.1 | 99 | −1.9 |
| A2v20 | −14.1 | −8.5 | −13.1 | −24.4 | −32.5 | 1.8 | 108 | 0.8 |
| C1v20w3 | −14.9 | −9.4 | −14.4 | −23.0 | −27.7 | 1.8 | 121 | 4.9 |
| C1v21w3 | −13.7 | −9.6 | −14.0 | −22.7 | −27.8 | 2.8 | 113 | 2.9 |

- **The IB/OB rz-slope mismatch is not cured by the OB chip DOFs. It grows:** par3 goes 87 → 108 with chip TY/TZ (A2v20)
  and → 121 with the cosmics on top (C1v20w3). B+ and B− show the same ordering (B+ 65 / 77 / 88 / 100 / 93, B− 60 / 72 / 84 /
  95 / 83 for A2v19 / C1v19w3 / A2v20 / C1v20w3 / C1v21w3). The IOB vertex slope follows (−9 → −15).
- **The Lr5/Lr6 IOB excess drops by about 5 µm** in C1v20w3/C1v21w3 (Lr6 −34 → −28), while the collisions-only A2v20 barely
  moves it. So this is the cosmic-driven projective component (6i), and the Lr5/Lr6 part of the user's IOB observation is
  partly real and partly fixed by it.
- So the remaining IOB signal is an IB-vs-OB slope difference that the collisions do not remove even with free OB chip
  positions, and which the cosmics push further. Candidates: an IB-side (radial or z-scale) term, or a non-geometric bias of
  the IB sub-track at the vertex. Next: par3 vs (z_PV, tgλ) per solution from the trees, and the IB half-barrel/stave
  radial and z-scale freedoms (6d) in a fit with cosmics.

## 6k. Per-chip residual maps (`rmap_*`, `macro/drawChipResidualMaps.C`) and the stave-end sigmoid

The maps are per-chip Gaussian means of the checkResid residuals, per half barrel, for C1v21w3 (2026-11-09).
- **Chip-to-chip structure at large |Z| is mostly statistical.** Group the chip columns by distance from the stave centre
  and compare the three runs (B0 562862, B+ 562837, B− 564644). In every group the median |value| ≈ the mean error, and
  χ²/n ≈ 2 at all |Z|, so the edges are not more significant than the centre. The error grows from 3–6 µm at the centre to
  20–26 µm at the Lr5/Lr6 ends. The run-to-run correlation of the chip values is 0.1–0.35, about the same at all |Z|. The
  reproducible part is √(corr·rms₁·rms₂) ≈ 5–8 µm per chip. With the median ± 4·MAD colour scale the noisy edge chips
  saturate.
- **The stave-common part at the edges is significant.** It is the sigmoid. Standard dz core means for
  tgλ ∈ [−2, −1.6] / [1.6, 2] (B0, µm):

| | A2v19 | C1v19w3 | A2v20 | C1v20w3 | C1v21w3 |
|---|---|---|---|---|---|
| Lr3 | −9.8 / 28.5 | −11.6 / 29.9 | −3.1 / 16.2 | −3.5 / 15.7 | −5.9 / 17.2 |
| Lr5 | −30.2 / 52.5 | −34.8 / 55.3 | −16.0 / 32.5 | −14.3 / 32.3 | −10.7 / 29.2 |
| Lr6 | 0.5 / −18.1 | 1.7 / −19.6 | −5.5 / −2.3 | −3.7 / 0.8 | 4.2 / −5.5 |

  (B− 564644 is the same within 1–10 µm.) A2v20 = C1v20w3, so the cosmics neither cause nor remove what is left.
- **Candidate reason why free chip DOFs stop at about half.** An odd-in-tgλ dz common to the staves is a radial displacement
  even in z (bow-like, dz = t·Δr) or a Δz odd in z. These shapes are what the v17+ gauge pins, as sensor-projected
  equations that chip DOFs cannot evade: the per-layer `bowing` groupModes (Δr = ε·z²), and the OB-wide `radialPolar_r0t0`,
  `r0t2` (Δr = ε, ε·t²) and `longitudinalPolar_r0t1` (Δz = ε·t). The stave-specific part is absorbed, the common part
  stays. The IB−OB tgλ mismatch ∝ tgλ (6j) is also of the pinned r0t0 / longitudinal r0t1 type. These polar p = 0 modes
  are broken by the vertex and by the cosmics, so they need not be pinned in a fit with cosmics.
- Test (not run): v21 without the per-layer `bowing` groupModes and the OB polar p = 0 modes, as C1 w3 and as A2. Look at
  the edge sigmoid, at IB−OB par3 and at the χ², and check that the v15-era per-layer bow is supported by the cosmics.
  Fast version prepared in `/data0/devAlg/tests/reducedFit` (README there). Everything is fixed at C1v21w3 except the OB
  stave/half-stave/module DOFs (plus module RX) and the mean vertex, 10 223 free. The constraints are reduced with the fixed
  labels moved to the right-hand side. S1 releases the five gauge equations above; S2 also releases `Ray mode dr`.
  - **Result (11-09, ~1/7 of the collision statistics, MINRES-QLP, full-rank constraints):** S1 reached Sum(W·χ²) 61.706 M,
    S2 61.685 M. No run with everything fixed at v21 on the same sample exists yet, so the absolute gain is unknown.
  - **S1 changes nothing:** the common module shifts are ≤ 4–6 µm, module RX ≤ 0.1 mrad. Edge dz (B0, tgλ ∈ [−2, −1.6] /
    [1.6, 2]): Lr3 −5.0 / 16.9, Lr5 −9.9 / 28.2 µm, against −5.9 / 17.2 and −10.7 / 29.2 for C1v21w3. Released, the gauge is
    not used: the collisions + cosmics χ² does not want the stave-common edge shapes at module level.
  - **S2 moves the sigmoid rather than curing it:** Lr5 at +2 goes 29 → 18 µm, but Lr6 worsens (−2 edge +4 → +11, +2
    edge −5 → −11). The common module TY oscillates along the stave (Lr5 +12 / −13, Lr6 +16 / −10 µm at modules 1 / 3),
    which looks like the released projective direction rather than a bow.
  - **Conclusion:** the gauge-pinning hypothesis is not confirmed. With the freedom to remove the stave-common edge dz, pede
    does not do it. Most likely its track model sees no such residual at the edge chips, while checkResid does. Next: compare the
    per-chip GBL residuals (`debugTracksFraction` → `gbltrack`, MilleData stage run with the C1v21w3 aligned geometry on
    a small sample) with the checkResid maps for the same geometry.
  - **GBL vs checkResid, same tracks (11-09).** Inputs: debug trees `/data0/devAlg/tests/Mille_C1V21/mille{B0,BPlus,BMinus}`
    (geometry of reducedFit/S1, PV mode 2, 2000 TFs) and checkResid of reducedFit/S1. 300 701 7-cluster B0 tracks matched by
    the sensor and z of their Lr3 and Lr6 clusters (tgλ difference rms 6·10⁻³). Core means for tgλ ∈ [−2, −1.6] / [1.6, 2], µm:

| | checkResid dz | GBL dzUnb | per-track difference at [1.6, 2] (corr) |
|---|---|---|---|
| Lr3 | −6.5 / 16.5 | −0.5 / 1.2 | 13.4 (0.75) |
| Lr4 | 1.5 / 8.5 | 2.9 / −5.9 | 14.2 (0.19) |
| Lr5 | −12.5 / 30.9 | −1.6 / 7.1 | 23.7 (0.58) |
| Lr6 | 3.1 / −3.8 | 1.9 / 1.6 | −2.5 (0.92) |

    The sigmoid is a property of the residual estimator, not of the geometry or of the track sample. PV contributors only
    (`t.composed`) give the same flat GBL values. The checkResid Lr5 residual is in effect an outward extrapolation (rms
    1050 µm vs 345 µm for GBL); at Lr6 both are one-sided and agree.
  - **Where the GBL fit puts it.** The fitted kinks have systematic means, odd in tgλ, much larger at large |tgλ|
    (`kinkZ`, B0, core means, µrad, tgλ = −1.9 … +1.7):
    - Lr0: +256 … −452
    - Lr1: +86 … −88
    - Lr3: about −100 at tgλ 0.9–1.3
    - Lr4: +374 at 1.7

    The straight B0 lines through the aligned geometry are not consistent at large |tgλ|. GBL absorbs this as scattering,
    which is cheap with the p = 0.6 GeV assumed for all B0 tracks (σ_MS ≈ mrad). So pede pays little for it, and
    checkResid shows it as the inner/outer prediction mismatch.
  - **B±, `debug_algPM.root`.** With pT > 2 the OB kinks are as large as with pT < 1, i.e. they do not scale with 1/p. They
    are a fixed-angle misfit, not scattering:
    - kinkZ at tgλ +1.8: Lr3 +286, Lr4 +300, Lr5 +133 µrad;
    - at −1.8: −75 / −109 / +12 µrad;
    - pulls 0.2–0.3.

    pede does see it in the B≠0 data. At fixed cluster Z (B±, pT > 1), the Lr3 kink at z ∈ [24, 32] cm changes sign with
    tgλ: −131 at 1.1, −24 at 1.4, +139 at 1.6, +468 µrad at 1.9. Lr4 at z ∈ [30, 40] cm goes −27 … +559 µrad. Tracks of
    different tgλ through the same chip come from different z_PV and cross other chips in the other layers, so this alone
    does not exclude geometry. checkResid dz at fixed tgλ shows no significant cluster-size dependence (±8 µm, at the
    statistical level).
  - **S3: edge chips free together with the released gauge** (`reducedFit/S3`, S1 plus TY, TZ, RX of the 8 064 chips of
    the end modules, priors 30 µm / 1 mrad; 34 415 free, 3 410 constraints, full rank; reduced statistics as S1/S2).
    - pede Sum(W·χ²) 61 324 771 vs S1 61 706 226 (−381k) and S2 61 685 055, with 63 667 rejected records vs 66 268.
    - The end chips move a lot: median |Δ| vs C1v21w3 TY 16–31 µm, TZ 6–41 µm (Lr4 smallest, Lr5 largest), RX 0.4–0.9 mrad;
      p90 TZ 90–130 µm at Lr3/5/6. All Lr5/6 end chips have < 2k entries in this sample, so the moves are prior-limited.
    - checkResid dz edge bins (µm, tgλ ∈ [−2,−1.6] | [−1.6,−1.3] | [1.3,1.6] | [1.6,2]), pp B0 562862:

      | | S1 | S3 |
      |---|---|---|
      | Lr3 | −5.0, 0.2, 2.8, 16.9 | −1.4, 2.1, 2.1, 14.6 |
      | Lr4 | 1.6, 3.9, −5.0, 6.5 | 7.4, 7.3, −8.0, 0.2 |
      | Lr5 | −9.9, −3.5, 6.5, 28.2 | 7.2, 5.2, 3.0, 19.0 |
      | Lr6 | 1.4, −0.6, 1.1, −4.9 | 5.1, 1.9, −0.7, −6.0 |

      B+ 562837 and B− 564644 behave the same (Lr5 B−: −16.2 … 28.9 → −2.4 … 20.1). The odd part at Lr5 drops from 19 to
      6 µm, but an even part of +5…+13 µm appears at Lr4/Lr5, so the dz vs tgλ curves are not flatter overall.
    - GBL kinks do not shrink, they are redistributed. B0 kinkZ at tgλ +1.8: Lr3 165 → 61, Lr4 374 → 272 µrad, but Lr5 at
      −1.8 14 → 113. B±, pT > 2: Lr3 at −1.8 −75 → +81, Lr5 at −1.8 12 → 181 µrad. GBL dzUnb (B0) was flat already in S1
      and stays flat.
    - **χ² weight of the misfit.** Averaged over all tracks, the mean OB kink pull per tgλ bin is ≤ 0.06 (B0 and B±, S3
      trees, debugTracksFraction = 1, 2000 TFs). Σ N·⟨pull⟩² over the tgλ bins is 0.1–1.8k per layer, against ~2M of
      Σ pull². Removing the mean kinks entirely would gain pede only a few k per 2000 TFs. The sigmoid is nearly free in
      χ² terms, so the fit has almost no lever to remove it, and moves of the end chips by tens of µm buy much more χ²
      elsewhere than they cost or gain at the edges.
- **IB staves at the half-barrel boundaries (dY).** Weighted per-stave means over the 9 chips, C1v21w3.
  - The standard dY is flat everywhere, within ±1 µm.
  - IOB dY in B0 (562862, IB cluster vs OB sub-track extrapolation): last stave of each half barrel +17…+21 µm in all
    three layers (L0 S5/S11 18/21, L1 S7/S15 17/19, L2 S9/S19 17/17). First stave of HB1 +11, first of HB0 +1…+3,
    middle staves 0…−3. It is cos2φ-like with a positive mean, and decreases from L0 to L2 by up to ×1.2, as expected
    from a direction kink about R ≈ 12 cm ((12 − r) lever).
  - B+ (562837) and B− (564644) give HB0 +6…+7 / −5…−7 and HB1 −7…−8 / −10…−11 at the same staves.
  - This is the B0 IB/OB kink of §1–6 in map form. It depends on the field, so it is not visible in the full-track
    residuals and not fixed by a single geometry. B0 cosmics also preferred the kink cure (§6f).

## 7. Notes for Claude (environment)

- Run ROOT through `bash -ic 'o2d >/dev/null 2>&1; root -l -b -q ...' < /dev/null`. `o2d` is an interactive alias; the
  "O2PDPSuite/latest was not found" message it prints is harmless. Without `< /dev/null` backgrounded DPL workflows hang.
- Compile the macros (`.C+`) from a copy in the scratchpad, not in `macro/`.
- If the shell environment points to `sw_old` again, load the new stack in a clean environment:
  `env -i HOME=$HOME PATH=/usr/bin:/bin bash -c 'ALIBUILD_WORK_DIR=$HOME/alice/sw alienv printenv O2PDPSuite/latest' > newenv.sh`,
  then run ROOT as `env -i HOME=$HOME PATH=/usr/bin:/bin bash -c '. newenv.sh; root ...'`.
- Since 2026-11-03 the new stack (`~/alice/sw/ubuntu2604_x86-64`, ROOT in `$ROOTSYS`) compiles the macros:
  `env -u PYTHONHOME $ROOTSYS/bin/root -l -b -q 'kinkMap.C+(...)'`; the libssl.so.1.1 / HistoManagerS messages of the
  logon are harmless.
- After the 2026-11 OS upgrade, ACLiC of the old (`sw_old`, ubuntu2404) ROOT fails on the system pthread headers
  (`__GTHREAD_COND_INIT`), and cling crashes on the O2/TGeo headers of `kinkUtils.h`. Until the new stack is built,
  `kinkMap.C` runs interpreted (no `+`) from a scratch copy whose `kinkUtils.h` drops `Frame`/`loadFrames` and the
  O2/TGeo includes; the geometry-reading macros need the new stack.
- The O2 Python breaks `python3`; use `env -u PYTHONHOME -u PYTHONPATH /usr/bin/python3 -I`.
- `TTree::Draw` with vector leaves needs `SetEstimate(nTracks*9+1)` or it overflows (done in `kink::cellHitShifts`).
- The `makeCDB.sh` scripts found in `/data/devAlg/*` query the CCDB at the current time and fail; `resToCCDB.sh` passes the
  timestamp explicitly.
- The `algDump.C` installed with O2 writes only the `gm` tree; the version with `gmx/gmy/gmz` is
  `~/alice/mymacrosO2/algDump.C`. The kink macros read the aligned geometry directly and do not need it.
- The original scratch macros of this study (`kmap*.C`, `emu*.C`, `cure.C`, `cmpalg.C`) are superseded by `macro/kinkMap.C`,
  `macro/kinkResponse.C` and `macro/makeKinkCure.C`.
