# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A standalone package built against ALICE O2, producing `libO2Align` plus one DPL workflow
executable, `o2-its3-alignment-workflow` ([src/alignment-workflow.cxx](src/alignment-workflow.cxx)).
It aligns ALICE barrel detectors by refitting tracks with GBL (General Broken Lines) and writing
Millepede-II input. Originally cloned from `O2/Detectors/Upgrades/ITS3/alignment`, so ITS/ITS3 is
still the most developed part.

## Build and run

```bash
./build.sh                      # loads the O2 env (alienv O2PDPSuite + Millepede-II), cmake + ninja install
ninja -C build -j 3 install     # incremental; build.sh uses -j 20, which can get the job OOM-killed
source install/init.sh          # puts o2-its3-alignment-workflow on PATH / LD_LIBRARY_PATH
```

There is no test suite and no linter config. Verification means running the workflow — see the
closure test in [README.md](README.md), which injects a known misalignment (`misAlgJson`) and checks
that the fit recovers it.

`Params` is an O2 `ConfigurableParam` whose prefix is **`AlignParams`** (`O2ParamDef(Params,
"AlignParams")` in [Params.h](include/O2Align/Params.h)). The README examples still use
`AlignmentParams.*`, which no longer matches.

Three-stage workflow, the same executable each time:
1. `--output MilleData,MilleSteer` over reconstructed data — writes `mp2data.bin` and the steering
   files `mp2tree.txt` / `mp2con.txt` / `mp2param.txt`. Note that `process()` constructs a
   `gbl::MilleBinary` per TF and GBL opens that file with `ios::out`, i.e. truncating, so a
   multi-TF run currently leaves only the records of the last TF in the binary.
2. `pede` (external, from Millepede-II) → `millepede.res`.
3. `--output MilleRes` — reads no data at all (the workflow inserts `NoInpDummyOutSpec`), only
   rebuilds the hierarchy and translates `millepede.res` into `result.json`, subtracting the
   injected misalignment if `misAlgJson` is given.

`--output` is a comma-separated `EnumFlags` list: `VerboseGBL, MilleData, MilleSteer, MilleRes,
MisRes, Debug` ([AlignmentSpec.h](include/O2Align/AlignmentSpec.h)).

Which `Detector` objects exist is decided by `--cluster-sources`, not `--track-sources`: the spec
derives `detMask` from the cluster sources, and the workflow allows only `ITS` there. So TPC/TRD/TOF
contribute measurements to tracks but currently get no `Detector` instance and no hierarchy branch.

## Architecture

### The volume hierarchy is the central abstraction

All detectors live in **one** tree ([Volume](include/O2Align/Volume.h)), rooted in a fictitious
volume from `Volume::makeRoot()` whose children are the detector top volumes. The root owns no
rigid-body DOFs, so unrelated detector branches are never mutually constrained. Every operation —
`applyDOFConfig`, `finalise`, `writeTree`/`writeRigidBodyConstraints`/`writeParameters`,
`writeMillepedeResults` — is applied once, to the root.

`Volume::finalise()` is the ordering-sensitive heart of the class: it assigns levels, makes sensors
and fictitious volumes define their own `L2G`/`T2L`, derives each volume's local-to-parent matrix
`mL2P` and the jacobians `mJL2P`/`mJP2L` from the (possibly pre-aligned) geometry, and auto-disables
a parent's DOFs when no child is active. Call it only on the root, after the DOF configuration.

Conventions that are easy to violate:
- **`isLeaf()` is used as the proxy for "is a sensor"** in the matrix code: a leaf defines its own
  matrices, a non-leaf takes them from `TGeoPhysicalNode`. A childless volume that is neither a
  sensor nor virtual would therefore silently get an identity `L2G`; `finalise` makes that fatal.
- **`virt = true`** means no alignable entry in the geometry (envelopes, TPC sectors, the mean
  vertex). Such a volume must define `L2G` itself, or leave it at the identity deliberately.
- **`setRigidBodyAllowed(false)`** marks volumes that cannot be shifted (the root, the detector
  envelopes w/o geometry): a matching `rigidBody` rule is ignored and no constraint over their
  children is generated. Calibration DOFs are unaffected.
- **Millepede labels** ([Label.h](include/O2Align/Label.h)) are bit-packed
  `DOF|CALIB|ID|SENS|DET`. Uniqueness across detectors relies purely on the `DET` field, since every
  detector restarts its own volume counter from 0. `Label::DET_GLOBAL` is reserved for the root.

### Detector plugins

Each detector is a [Detector](include/O2Align/Detector.h) subclass with three jobs: `prepareData`
(once per TF, e.g. transforming clusters/tracklets that do not depend on the track), `prepareTrack`
(per track, appending `FrameInfoExt` + `Measurement` to the `Track`), and the protected
`buildHierarchy` (its own stand-alone branch). `Detector::attachTo` grafts that branch under the
common root, discards a branch containing no sensor, and remembers the top as `getTopVolume()`.

`Sensor*` classes (`SensorITS`, `SensorIT3`, `SensorTPC`, `SensorTRD`, `SensorTOF`) exist only to
override `defineMatrixL2G` / `defineMatrixT2L` / `computeJacobianL2T` for volumes whose measurement
plane differs from the geometry volume or has no geometry at all.

`DetectorPVT` is a virtual detector with a single sensor-like volume, the mean interaction vertex.
It provides no measurement: frames with `isVertex()` bypass the hierarchy derivative chain, and the
vertex position enters the fit through `mMeanVtxLabels`. Its DOFs are configured on the volume built
for slot 0, while the emitted labels carry the current mean-vertex time slot
([TimeSlotsSet](include/O2Align/TimeSlotsSet.h)), so each slot is aligned independently.

### From a measurement to global derivatives

In [AlignmentSpec.cxx](src/AlignmentSpec.cxx) the per-point globals are built by walking from the
measurement leaf up to the root: the base derivative in the tracking frame is converted to the
sensor's local frame via `computeJacobianL2T`, then transported level by level with each child's
`getJL2P()`, each volume with a `RigidBodyDOFSet` contributing its columns. The top volume of the
detector participates; the root, having no DOFs, ends the chain.

This block exists **twice**, in `process()` and in `fillGBLPoints()`, in near-identical form — keep
them in sync. Calibration derivatives are collected from one volume only, the direct parent of the
measurement leaf (for ITS3 that is the half-barrel sensor above the pseudo tiles), so a calibration
DOF set placed on a leaf or on an envelope is written to the parameter file but does not yet enter
the fit.

`writeRigidBodyConstraints` implements the hierarchical constraint: for each free DOF of a parent,
the weighted mean of the corresponding child DOFs (transported with `getJP2L()`) is required to
vanish.

### DOF sets and JSON configuration

[DOFSet](include/O2Align/DOFSet.h) is the DOF-block interface (`nDOFs`, free/fixed mask,
`fillDerivatives` from a `DerivativeContext`): `RigidBodyDOFSet` (TX,TY,TZ,RX,RY,RZ in the local
frame), `LegendreDOFSet` and `InextensionalDOFSet` (deformation modes of ITS3 half-cylinders).

Two JSON files drive a run, both documented by example in the README:
- `dofConfigJson` — `defaults` plus `rules`, each with an fnmatch `match` pattern (also retried with
  an implicit leading `*`), a `rigidBody` and/or `calib` clause. Rules are applied in order, a later
  match replaces an earlier DOF set. This is the only place where DOFs are freed or fixed; the
  hierarchy construction itself fixes nothing.
- `misAlgJson` — misalignment injected per sensor ID ([MisalignmentUtils](include/O2Align/MisalignmentUtils.h)),
  reused at stage 3 as the reference to subtract from the fitted values. Injection forces
  single-threaded processing (the `nthreads` option is otherwise honoured, OpenMP).

## Repository notes

- `sav/` is a tracked snapshot of the pre-refactor code (`AlignableVolume`, `AlignmentDOF`, …). It
  is not built and not maintained — read it for history only, never edit it.
- `macro/` holds ROOT macros installed to `share/macro`, e.g. `CreateTimeSlots.C`, which produces
  the mean-vertex time-slot JSON consumed by `AlignParams.MVTimeSlotsJson`.
- `build/`, `install/`, `*.root` and `*.so` are git-ignored; stray `core_dump_*` and `*.json` run
  products in the working tree are not, so do not `git add -A`.
