# DOF Configuration Rules

`Volume::applyDOFConfig` traverses the hierarchy and applies matching rules to non-pseudo volumes. A `defaults` object, if present, is applied first to every volume; rules then run in order, and a later matching clause replaces the earlier DOF set of the same kind. Patterns are matched against symbolic volume names, with an implicit `*` prepended. Since `*` can match `/`, use exact-depth patterns when targeting a hierarchy level.

## Mean Vertex

The mean vertex volume is `PVT/meanVertex`. Its X, Y, and Z calibration DOFs default to free. Enable all three with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex" } }
```

Fix all three with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex", "fixed": true } }
```

Keep X, Y, and Z either all free or all fixed: the vertex fit only uses the position labels when all three are free.

## TPC Drift

The TPC drift calibration belongs to `TPC_envelope` and is not enabled by default. Enable both `VDRIFT` and `DRIFTOFF` with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift" } }
```

Fix both with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift", "fixed": true } }
```

## ITS Inner Layer

For standard ITS (not ITS3), Layer 0 symbolic names have the form `ITS/ITSULayer0/ITSUHalfBarrelN` and `ITS/ITSULayer0/ITSUHalfBarrelN/ITSUStaveN`. These rules free translations only on the Layer 0 half-barrels and staves:

```json
[
  {
    "match": "ITSULayer0/ITSUHalfBarrel?",
    "rigidBody": ["TX", "TY", "TZ"]
  },
  {
    "match": "ITSULayer0/ITSUHalfBarrel?/ITSUStave?",
    "rigidBody": ["TX", "TY", "TZ"]
  }
]
```

The single-character `?` matches the half-barrel or stave index without matching their descendants. ITS3 uses different symbolic names, such as `ITS3CarbonForm`, so these patterns do not select ITS3 volumes.