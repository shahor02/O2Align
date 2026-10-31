# Closure test with an injected misalignment

`AlignParams.misAlgJson` makes the workflow shift the measurements of the tracks by a known misalignment
while the Millepede data are filled. The fit is then expected to return that misalignment, and
the stage that converts `millepede.res` subtracts it, so what is left in `result.json` and `closure.json` is
the **error of the alignment**, to be compared with zero.

What is injected, in which volumes and with which magnitude is defined by a file of rules, in the style
of the DOF configuration ([DOFConfig_rules.md](DOFConfig_rules.md)). Nothing is injected if `misAlgJson` is empty.
There are no other switches.

## What the test does and does not test

The values are injected **linearly through the derivatives of the fit itself**: a measured point of a
track gets `res += d(prediction)/d(DOF) * value` for all the DOFs of its sensor, of the volumes above it
in the hierarchy (rigid body) and of the volume right above the sensor (calibration DOFs), see
`AlignmentSpec::buildPointGlobals`. Therefore the test checks
the constraints, the weak modes, the statistics, `pede` and the translation of its result, **not the
derivatives**: an error of a derivative would be injected and recovered consistently. The derivatives are
tested by injecting a misalignment into the geometry, see [closure_chipRY.sh](closure_chipRY.sh) and
`macro/injectChipRY.C`.

An exception is the opt-in `"model": "exact"` of the ITS3 `legendre` and `inextensional` calibration DOFs: the
non-linear deformation of the half barrel is evaluated by its own code (`evaluateLegendreShift`,
`evaluateInextensionalShift`), independent of the linear derivatives of the fit, which then
tests them.

## The misalignment file

A JSON object with the array `"rules"`, or the bare array of rules. Comments: any object member with a key starting with
`#` or `_comment` is dropped, at any depth, and so is an array element left empty by that (the same rule as in the
DOF config). There is no `defaults`: giving it is fatal.

```json
{
  "#": "L3 HB0 and two of its staves",
  "rules": [
    { "match": "ITS/ITSULayer3/ITSUHalfBarrel0",
      "rigidBody": { "TX": 30e-4, "TY": -20e-4, "TZ": 50e-4, "RZ": 5e-5 } },
    { "match": "ITS/ITSULayer3/ITSUHalfBarrel0/ITSUStave2", "rigidBody": { "TX": 40e-4, "RZ": 1e-4 } },
    { "match": "ITS3Layer0/ITS3CarbonForm0",
      "calib": { "type": "legendre", "model": "exact", "values": { "L(1,0)": 8e-4 } } }
  ]
}
```

| Key | Meaning |
|---|---|
| `match` | regex of the symbolic volume names, matched exactly like in the DOF config: against the full name, or with an implicit leading `.*`. Pseudo volumes (the ITS3 tiles) and the root are not matched. Give the full name ending at the volume (`ITS/ITSULayer3/ITSUHalfBarrel0/ITSUStave2`): a pattern ending with `.*` would also match the staves 20..., and a bare `ITSUStave2` the stave 2 of every half barrel and layer. |
| `rigidBody` | object `DOF -> value`, DOF in `TX,TY,TZ,RX,RY,RZ`, in the local frame of the volume, **cm and rad**: exactly the DOFs of the fit. Unlisted DOFs get 0. An array of the 6 values `[TX,TY,TZ,RX,RY,RZ]` is also accepted (its zeros are not checked against the DOF config). |
| `calib` | `{ "type": ..., "model": ..., "values": {...} }`. `type`: `legendre`, `inextensional` or `tpcvdrift`. `model`: `linear` (default) or `exact`. `values`: object `DOF name -> value`, the names being those of the configured DOF set (the ones in `mp2param.txt`): `L(i,j)`; `a_n,b_n,c_n,d_n` (n = 2..), `alpha`, `beta`; `VDRIFT`, `DRIFTOFF`. The value is injected in every calibration time slot of the volume. |

At least one of `rigidBody` and `calib` is required, no other keys are allowed.

**Rules are applied in order, and a later matching rule replaces the earlier one, it does not merge with it.**
If two rules with a `rigidBody` clause match a volume, only the values of the later one are injected; the
same for `calib`. (A `calib` clause does not replace a `rigidBody` one.) This is the semantics of the DOF config.

### Checks against the DOF configuration

The file is validated when the hierarchy is built, after the DOF config and the config macro have been
applied (so a standalone `MilleSteer` stage validates it as well). Any mismatch is `LOGP(fatal, ...)`
with a message naming the rule and the volume. Fatal are:

- the file cannot be opened, `defaults` is present, a rule is not an object, has no string `match`, has neither
  `rigidBody` nor `calib`, has an unknown key (also in `calib`), or has an invalid regex;
- a rule matching no volume (a typo);
- `rigidBody` on a volume which is not rigid-body alignable (the root, the envelopes without geometry), an unknown
  DOF name, a value which is not a finite number, an array of other size than 6;
- a listed rigid-body DOF (even with the value 0) which is not free in the DOF config, or the volume has no rigid-body
  DOF set at all: the fit could not recover it;
- `calib` on a volume without a calibration DOF set, of another `type` than the configured one, with an unknown DOF name
  or with a DOF which is fixed;
- `calib` of the type `meanvertex` (not supported), or on a volume none of whose children is a measurement leaf: the
  calibration derivatives are collected only from the direct parent of the leaf, other would not enter the fit;
- `"model": "exact"` for a type other than `legendre`/`inextensional`, on a volume which is not an ITS3 half barrel
  (the parent of the pseudo tiles, sensor id 0..5), or exact `legendre` when processing data without MC (it uses the
  MC truth of the tracks; the `inextensional` one does not).

At the start the table of the injected volumes, DOFs and values is logged (`Injected misalignment`).

## Stages

The same executable is used as for an alignment, see the [README](../README.md); `AlignParams.*` are `--configKeyValues` entries:

1. `--output MilleData,MilleSteer` on the data, with `AlignParams.dofConfigJson` and `AlignParams.misAlgJson`.
   The injection happens when the GBL points are filled.
2. `pede` on the produced `mp2data.bin` and the steering files.
3. `--output MilleRes` (reads no data): `millepede.res` → `result.json`, `<DET>_alignment.root`,
   with the injection **subtracted** from `result.json` and `closure.json`, if `misAlgJson` is given. Use the same
   `misAlgJson` and `dofConfigJson` as at stage 1.

### Closure on Monte Carlo

On MC with the ideal geometry the track data have no misalignment on their own, so the fit has to return
the injected values: `result.json` and `closure.json` are compatible with 0 within the statistical error of the fit.

### Closure on real data

The data are misaligned on their own, the fit returns that **plus** the injection. Hence a **control run** is needed:
the same data, the same DOF config and `pede` steering without `misAlgJson`. At stage 3 of the injected run, give the
`millepede.res` of the control run as `AlignParams.misAlgRefMilleRes=<path>`. The result is then
`fit - control - injected` per label. (`misAlgRefMilleRes` without `misAlgJson` is ignored with a warning.)
The control run is also the one to compare with to see which part of the residual is the statistical noise of the particular
data selection and DOF set.

## Reading the result

- **`result.json`** (`AlignParams.milleResOutJson`): per volume with free DOFs, `fitted - control - injected` of the DOFs, i.e. zeros (within
  errors) if the injected DOF is recovered. The ITS3 `matrix`/`inextensional` entries are in the same layout as before.
  The other output of the stage, `<DET>_alignment.root` (`AlignParams.algParamsOutFile`), is the alignment objects of the
  absolute fitted values, with nothing subtracted.
- **`closure.json`** (`AlignParams.closureOutJson`): the same at the **sensor level**. The sensors are displaced by the composition of
  the corrections of all the levels above them (as done for the weak modes diagnostics, `computeSensorDisplacements`);
  for each detector branch and layer the file (and the log) gives the number of sensors, the RMS and the maximum of the 3D
  displacement of the sensor centres of the residual in µm (the group is the whole layer of the branch, with the sensors
  which are not misaligned, so the injected RMS is diluted by them), the RMS of the displacement of the injection itself, and the RMS
  of the rotation of the residual and of the injection (µrad). The residual should be much smaller than the injection, at the level
  of the statistical precision. For an ITS3 half barrel (a sensor made of the pseudo tiles) the displacement of the
  local origin and the rotation are given, the non-rigid deformation DOFs only in `result.json`.

### Why the per-volume numbers can be non-zero

When a DOF is free at two nested levels, as in the example (b) (the half barrel and all of its staves), the split of a
displacement between the levels is not defined by the data: the automatic hierarchical constraint (the weighted mean of the movement
of the children of a free parent is zero) moves the transported mean of the staves' injection into the half barrel. The fit
returns then, for the half barrel, its injected value **plus** the mean of the stave injections, and for the staves their injected values
**minus** that mean. In `result.json` this appears as non-zero opposite residuals of the half barrel and of
every one of its staves (also the non-misaligned ones!) while the sensors are at the right places: this is what
`closure.json` shows, with a residual about 0.

For the same reason all the 12 staves of the half barrel are free in the example (b), although 2 of them are misaligned: with
only the 2 staves free, the constraint would force their mean to 0, which biases the result as nothing else absorbs it.

## Example (a): ITS3 layer 0 half barrel 0

ITS3 MC is needed (`--with-its3`, and no
`--disable-mc`) because of the exact Legendre model.

Files: [misalign_ITS3_L0HB0.json](misalign_ITS3_L0HB0.json), [misalign_dof_ITS3_L0HB0.json](misalign_dof_ITS3_L0HB0.json). It injects
into `ITS3Layer0/ITS3CarbonForm0` the rigid body `TX=10 µm, TY=5 µm, RY=100 µrad` and the exact Legendre mode `L(1,0)=8 µm`, the DOFs being
free in `misalign_dof_ITS3_L0HB0.json`. Dropping `"model": "exact"` from the file injects `L(1,0)` linearly (and
then no MC is needed).

```bash
D=/data/devAlg/tests/closureITS3; mkdir -p $D && cd $D
source /home/shahoian/alice/O2Align/install/init.sh      # if O2's libO2Align.so gets picked, prepend install/lib to LD_LIBRARY_PATH
DOC=/home/shahoian/alice/O2Align/doc
OPT="AlignParams.dofConfigJson=$DOC/misalign_dof_ITS3_L0HB0.json;AlignParams.misAlgJson=$DOC/misalign_ITS3_L0HB0.json;AlignParams.minPt=0.1"

# 1. Millepede data and steering files from the ITS3 MC reconstruction (the input of the workflow as for any ITS3 alignment)
o2-dev-alignment-workflow --with-its3 --track-sources ITS --output MilleData,MilleSteer \
  --configKeyValues "$OPT" -b --run | tee stage1.log
# 2. pede: mpSteer.txt = your usual steering file (method, options) listing mp2con.txt, mp2param.txt and mp2tree.txt, see pede_hints.md
echo mp2data.bin >> mpSteer.txt
pede mpSteer.txt | tee pede.log
# 3. result.json and closure.json, the injection subtracted
o2-dev-alignment-workflow --with-its3 --track-sources ITS --output MilleRes \
  --configKeyValues "$OPT" -b --run | tee stage3.log
```

The expected result: the entries of `ITS3Layer0/ITS3CarbonForm0` in `result.json` (`rigidBody` and `matrix`) and the layer 0 line of
`closure.json` are consistent with 0. The Legendre mode is injected non-linearly, its recovery is
approximate, to the level of the linearisation.

## Example (b): ITS layer 3 half barrel 0 and 2 of its staves

Files: [misalign_ITS_L3HB0.json](misalign_ITS_L3HB0.json), [misalign_dof_ITS_L3HB0.json](misalign_dof_ITS_L3HB0.json). The
rigid body misalignment is injected into `ITS/ITSULayer3/ITSUHalfBarrel0` (`TX=30 µm, TY=-20 µm, TZ=50 µm, RZ=50 µrad`), `ITSUStave2`
(`TX=40 µm, RZ=100 µrad`) and `ITSUStave7` (`TZ=-60 µm, RY=50 µrad`), all in the local frames. In the DOF config the 6 DOFs of
the half barrel and of all its 12 staves, and the mean vertex, are free, everything else is fixed and defines the frame: no global-mode
constraints are needed.

Run it on real data, with a control run, using the processing script [proc0.sh](proc0.sh): `DataDir` and `Config` (the DOF
config) are required, `OutType` selects the stage, `ConfString` replaces the default `AlignParams` string (`EXTRACONF` is appended to it)
and the readers are started only for the data-processing stage, any non-empty `DataDir` (e.g. `N/A`) does for the others. The default `ConfString`
has the production options (`useMultiTrackPVConstraint`, `MVTimeSlotsJson` of the `/data0` machine, ...): set what you need, but the same in
both cases. `pede` needs `libopenblas` in `LD_LIBRARY_PATH`. The test area of this machine is `/data/devAlg`;
[closure_chipRY.sh](closure_chipRY.sh) still has the `/data0/devAlg` paths of another one.

```bash
DOC=/home/shahoian/alice/O2Align/doc
DATA=/data/devAlg/data/_alice_data_2025_LHC25ab_562862_apass1_alignment_its_0010_o2_ctf_run00562862_orbit0261432416_tf0000451373_epn313
W=/data/devAlg/tests/closureL3HB0
export Config=$DOC/misalign_dof_ITS_L3HB0.json
BASE="AlignParams.dofConfigJson=$Config;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20"
INJ="$BASE;AlignParams.misAlgJson=$DOC/misalign_ITS_L3HB0.json"
mkdir -p $W   # with mpSteer.txt: the usual pede steering file listing mp2con.txt, mp2param.txt, mp2tree.txt (see pede_hints.md); for a single TF set "entries 5"
for CASE in ctl inj; do
  mkdir -p $W/$CASE && cd $W/$CASE
  CONF=$BASE; [[ $CASE == inj ]] && CONF=$INJ
  # 1. Millepede data and steering files (the same data for both cases)
  DataDir=$DATA ConfString="$CONF" OutType=MilleData,MilleSteer Nthreads=5 MAXTF=1 runOPTL="-b --run" $DOC/proc0.sh | tee stage1.log
  # 2. pede
  cp ../mpSteer.txt . && echo mp2data.bin >> mpSteer.txt
  pede mpSteer.txt | tee pede.log
done
# 3. the result of the injected case, relative to the control one and to the injection (no data read)
cd $W/inj
DataDir=N/A ConfString="$INJ;AlignParams.misAlgRefMilleRes=../ctl/millepede.res" OutType=MilleRes runOPTL="-b --run" $DOC/proc0.sh | tee stage3.log
```

The initial alignment (CCDB or a remapped `ITS_alignment.root`, see `--condition-remap` in `closure_chipRY.sh`) has to be the same in both cases.

Expected: `closure.json` has, for `ITS` layer 3, a residual RMS much below the injected one (about 40 µm RMS of the displacement,
the exact number is in the `injectedRmsUm` field), at the level of the statistical precision of the data set used,
which is what the control run shows as well; the other layers have 0 residual, all their DOFs being fixed. In `result.json` the half barrel and
its 12 staves have the non-zero residuals described [above](#why-the-per-volume-numbers-can-be-non-zero).

## No-data check of a configuration

`--output MilleSteer` without `MilleData` reads no data but builds the hierarchy, applies the DOF config and validates the
misalignment file against it (a fatal on any mismatch), so it is a quick check of a pair of files:

```bash
o2-dev-alignment-workflow --disable-mc --hbfutils-config none --output MilleSteer --detectors ITS --track-sources ITS \
  --configKeyValues "AlignParams.dofConfigJson=misalign_dof_ITS_L3HB0.json;AlignParams.misAlgJson=misalign_ITS_L3HB0.json" -b --run
```

The geometry comes from the CCDB (`ALICEO2_CCDB_LOCALCACHE` can point to a local copy).

## Debug output

With `--output MisRes` and `misAlgJson` set, a tree `mis` of the points with a non-zero injected shift
is written to `debug_alg[_<lane>].root` (see [DebugInfo.md](DebugInfo.md)): the shifted residual `dy,dz`, the injected `ddy,ddz`, the detector `det`, the sensor `sens`, `z`, `phi`.
