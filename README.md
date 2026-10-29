# Simulate ITS/ITS3 misalignment and re-alignment

A known misalignment, defined by a file of rules (volumes, DOFs, values), is injected into the data
by setting `AlignParams.misAlgJson`; the result of the alignment is then given relative to it. See
[doc/MisalignmentTest.md](doc/MisalignmentTest.md) for the format of the file, the checks against the DOF
configuration, two examples (ITS3 layer 0 half barrel 0; ITS layer 3 half barrel 0 and two of its staves) and
the instructions of the closure test.

```bash
o2-dev-alignment-workflow --with-its3 --track-sources ITS --output MilleData,MilleSteer --configKeyValues "AlignParams.minPt=0.1;AlignParams.dofConfigJson=doc/misalign_dof_ITS3_L0HB0.json;AlignParams.misAlgJson=doc/misalign_ITS3_L0HB0.json" -b --run
```

## Running workflow

`o2-dev-alignment-workflow` is a DPL workflow; the same executable is used at all three stages
(Millepede data/steering, `pede`, results). If its output is piped to other shell commands, add
`--run` before the first non-workflow element: `o2-dev-alignment-workflow ... --run | tee my.log`.
Parameters of the alignment itself are set via `--configKeyValues "AlignParams.<name>=<value>;..."`.

### Options defined by the workflow

| Option | Default | Meaning |
|---|---|---|
| `--disable-mc` | off | MC labels are read and propagated unless this is given |
| `--track-sources` | all | comma-separated list of track sources, e.g. `ITS,ITS-TPC,ITS-TPC-TRD`; `cosmic-ext` switches to the `TrackCosmicsExtended` input. Only `ITS, TPC, ITS-TPC, ITS-TPC-TRD, ITS-TPC-TOF, ITS-TPC-TRD-TOF` are accepted |
| `--detectors` | `ITS` | detectors participating in the alignment and the accepted combinations of them, see below |
| `--with-its3` | off | ITS3 alignment mode |
| `--output` | empty | comma-separated list of `VerboseGBL, MilleData, MilleSteer, MilleRes, MisRes, WeakModesDiag, WeakModesDiagIni, Regularize`. `MilleData` and `VerboseGBL` imply reading the data; `MilleRes`, `WeakModesDiag`, `WeakModesDiagIni` and `Regularize` read no data (a dummy source is inserted). `Regularize` cannot be combined with `VerboseGBL`/`MilleData` |
| `--disable-root-input` | off | do not add the root-files input reader |
| `--configKeyValues` | empty | semicolon-separated `key=value` list of `ConfigurableParam` settings (`AlignParams.*`, `pvertexer.*`, ...); applied after the built-in primary-vertex defaults, so it overrides them |

The alignment device itself adds `--nthreads` (number of OpenMP threads of the track refit, default 1)
and `--config-macro` (see [User configuration macro](#user-configuration-macro)).

Options added by shared O2 helpers, not described here (use `--help full`): the HBF utils options
(`--hbfutils-config`, ...), the TPC correction-map options (needed when `TPC` is among the
detectors) and the ITS configuration option.

### Selecting detectors and their combinations

`--detectors` takes a comma-separated list whose items are either a single detector or several
detectors joined by hyphens. Allowed detectors are `ITS`, `TPC`, `TRD`, `TOF`. The list is
interpreted in two ways:

1. **Participating detectors**: the union of all items. Only these have their clusters/tracklets
   loaded, their volumes in the hierarchy and their measurements added to the tracks.
2. **Accepted combinations**: every item gives one combination. After the refit, a track is kept
   only if the set of detectors that actually provided its points is *exactly* one of the listed
   combinations; other tracks are discarded (the track as a whole, not just the detector's data).
   This makes it possible to reject a track on the basis of one detector while still keeping
   tracks built with the others.

Example: `--detectors ITS,ITS-TRD,TOF-ITS` uses ITS, TRD and TOF (like `ITS,TRD,TOF` did before), and
accepts tracks seen by ITS alone, by ITS+TRD, or by ITS+TOF. An ITS-TPC-TRD-TOF track is rejected,
since ITS+TRD+TOF is not in the list, while an ITS-TPC-TRD track passes: TPC is not participating,
so it counts as ITS+TRD. The track source must of course be among `--track-sources`.

Plain lists (`--detectors ITS,TPC`) thus accept only tracks having exactly the ITS, or exactly the
TPC, points; to also accept the tracks with both, add `ITS-TPC`.

## In-existensional modes
```json
{
  "defaults": { "rigidBody": "fixed" },
  "rules": [
    {
      "match": "ITS3Layer1/ITS3CarbonForm0",
      "calib": {
        "type": "inextensional",
        "order": 2,
        "free": ["a_2", "b_2", "c_2", "d_2", "alpha", "beta"]
      }
    }
  ]
}
```

```json
[
  {
    "id": 2,
    "inextensional": {
      "modes": {
        "2": [0.0008, -0.0005, 0.0006, -0.0007]
      },
      "alpha": 0.0004,
      "beta": -0.0003
    }
  }
]
```


## TPC drift calibration

The TPC envelope owns the drift calibration DOFs: `VDRIFT`, a relative correction to the drift
velocity, and `DRIFTOFF`, an offset of the drift length in cm (the equivalent of a t0 shift). The
cluster Z being reconstructed as `z = s * (Lz - l)` with `l = (t - T0 - tOffset) * vDrift` and
`s = +1 (-1)` on the A (C) side, the two displace the measurement along Z by
`dz = -s * l * VDRIFT - s * DRIFTOFF`, and are separated by the lever arm in `l`.

```json
{
  "defaults": { "rigidBody": "fixed" },
  "rules": [
    { "match": "TPC_envelope", "calib": { "type": "tpcvdrift" } }
  ]
}
```

`"free": ["VDRIFT"]` (or `"fix": ["DRIFTOFF"]`) fits the drift velocity alone.

The drift is calibrated independently in consecutive time intervals, declared in
`AlignParams.VDTimeSlotsJson` in the same format as `AlignParams.MVTimeSlotsJson` (see
`macro/CreateTimeSlots.C`): every slot gets its own set of global parameters, the slot ID being
encoded in the volume ID of the Millepede label. Without the file the whole run is a single slot.
Stage 3 writes one record per slot under `tpcDrift` in `result.json`.

Note that `DRIFTOFF` is degenerate with a common Z shift of the TPC sectors, and would also become
degenerate with a `TZ` DOF of the envelope should the latter be made rigid-body alignable.

The fitted correction is relative to the drift calibration already accounted for by the correction
maps: whenever TPC clusters are requested, the workflow asks for the per-TF `TPC/TPCCORRMAP` input
(produced upstream by `o2-tpc-scaler`, as for any other TPC cluster consumer) and for the drift
CCDB objects, and builds the `GPUParam` used for the cluster errors from the current field.


## User configuration macro

`--config-macro mySetup.C` compiles the macro with ACLiC in `init` and executes the function named
after the file:

```cpp
int mySetup(std::vector<o2::alignrs::Detector*>* dets, int)
{
  for (auto* det : *dets) { /* ... */ }
  return 0; // anything else is fatal
}
```

The call happens at the end of `buildHierarchy`, right after `dofConfigJson` has been applied and
before `finalise()`, so the macro sees the volumes of every detector (`Detector::getTopVolume`) and
has the last word on their DOFs. A non-zero return value aborts the workflow.
