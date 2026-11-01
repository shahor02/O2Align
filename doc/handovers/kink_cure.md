# Handover: detecting the IB/OB kink and building a geometric cure for it

Written for: whoever continues the ITS alignment (a person or Claude), familiar with O2Align, Millepede-II and checkResid.
Date: 2026-10-31. State: v19 (`/data/devAlg/solA2v19`, config `doc/ITS_RB_v19.json`) and its cure `cure012`
(`/data/devAlg/solA2v19cure/cure012`), with checkResid on pp B0 (LHC25ab 562862) and OO B+ (LHC25ae) for both.
Section 6 (pede evaluation of the cure) was added on 2026-11-01.

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
  data prefer cure012 in checkResid and v19 in the Mille fit: the two estimators disagree on this mode, as they did for
  v18 vs v19.

## 6a. Next: is the checkResid kink set by the vertex?

The checkResid IB sub-track contains the refitted PV as a point (`checkresid.addPVAsCluster = true`, the default). With
σ(PV) of a few µm and the IB at 2.2–3.9 cm, the PV dominates the IB sub-track direction: a 5 µm PV–IB transverse offset
alone gives ~150 µrad in dPar2, the size of the kink. In the Mille data the vertex enters only through the composed
free-vertex fit of PVs with ≤ 15 tracks (`useMultiTrackPVConstraint = 2`), so it cannot impose such a direction.
Test: rerun checkResid on pp B0 562862 for v19 (and cure012) with `checkresid.addPVAsCluster=false`, so the IB sub-track
is Lr0..Lr2 only, and compare with `kinkMap.C`. If the v19 kink disappears, it is a PV/IB effect and the OB cure
compensates it in the wrong place; if it stays, the IB/OB disagreement is in the tracks themselves.

Other cheap tests not yet done: checkResid on the second B0 run 562837; a pair of runs of the same period with B = 0 and
B ≠ 0.

## 7. Notes for Claude (environment)

- Run ROOT through `bash -ic 'o2d >/dev/null 2>&1; root -l -b -q ...' < /dev/null`. `o2d` is an interactive alias; the
  "O2PDPSuite/latest was not found" message it prints is harmless. Without `< /dev/null` backgrounded DPL workflows hang.
- Compile the macros (`.C+`) from a copy in the scratchpad, not in `macro/`.
- The O2 Python breaks `python3`; use `env -u PYTHONHOME -u PYTHONPATH /usr/bin/python3 -I`.
- `TTree::Draw` with vector leaves needs `SetEstimate(nTracks*9+1)` or it overflows (done in `kink::cellHitShifts`).
- The `makeCDB.sh` scripts found in `/data/devAlg/*` query the CCDB at the current time and fail; `resToCCDB.sh` passes the
  timestamp explicitly.
- The `algDump.C` installed with O2 writes only the `gm` tree; the version with `gmx/gmy/gmz` is
  `~/alice/mymacrosO2/algDump.C`. The kink macros read the aligned geometry directly and do not need it.
- The original scratch macros of this study (`kmap*.C`, `emu*.C`, `cure.C`, `cmpalg.C`) are superseded by `macro/kinkMap.C`,
  `macro/kinkResponse.C` and `macro/makeKinkCure.C`.
