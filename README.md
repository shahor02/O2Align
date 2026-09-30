# Simulate ITS3 misalignment and re-alignment


```bash
o2-its3-alignment-workflow --track-sources ITS --output MilleData,MilleSteer --configKeyValues "AlignmentParams.minPt=0.1;AlignmentParams.doMisalignmentLeg=true;AlignmentParams.doMisalignmentRB=true;AlignmentParams.misAlgJson=test_closure.json;AlignmentParams.extraClsErrZ[0]=10e-4;AlignmentParams.extraClsErrY[0]=10e-4;AlignmentParams.extraClsErrZ[3]=10e-4;AlignmentParams.extraClsErrY[3]=10e-4;AlignmentParams.dofConfigJson=dofSet.json" -b --run
```

test_closure.json:
```json
[
  {
    "id": 0,
    "rigidBody": [0.001, 0.0005, 0.0, 0.0, 0.0001, 0.0],
    "matrix": [[0.0], [0.0008, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
  }
]
```

dofSet.json:
```json
{
  "defaults": { "rigidBody": "fixed" },
  "rules": [
    {
      "match": "ITS3Layer0/ITS3CarbonForm0",
      "rigidBody": ["TX", "TY", "RY"],
      "calib": { "type": "legendre", "order": 1, "fix": [0, 2] }
    }
  ]
}
```


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
