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
+51 M on the same data, and built the v19 one. See [history.md](../history.md).

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
- The Mille data disfavour the cure, but by 0.05 per record, 70 times less than the v18 → v19 χ² gap (+51 M). The value is an
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
