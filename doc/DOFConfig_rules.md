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

## ITS Alignable Match Strings

The patterns below are based on the output from the `macro/printAlignableVolumes.C` running over Run3 TGeometry root file. 
The ITS envelope is the `ITS` alignable. The layer name is a path component, not a separate alignable entry. In the match patterns, layer IDs use `[0-2]` for the inner barrel and `[3-6]` for the outer barrel. Half-barrel, half-stave, module, and chip IDs in these paths are single digits, matched with `?`. Inner-barrel stave IDs are single digits; outer-barrel stave IDs can be one or two digits, so both exact-width forms are listed there. Avoid a trailing `*` when matching a specific level: `fnmatch` allows it to match `/` and descendants too.

### Inner Barrel (Layers 0–2)

The inner-barrel hierarchy has half-barrels, staves, and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?` |
| Stave | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?` |
| Chip | `ITS/ITSULayer[0-2]/ITSUHalfBarrel?/ITSUStave?/ITSUChip?` |

### Outer Barrel (Layers 3–6)

The outer-barrel hierarchy adds half-staves and modules between staves and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?` |
| Stave, one-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?` |
| Stave, two-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]` |
| Half-stave, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?` |
| Half-stave, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?` |
| Module, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?/ITSUModule?` |
| Module, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?/ITSUModule?` |
| Chip, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave?/ITSUHalfStave?/ITSUModule?/ITSUChip?` |
| Chip, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel?/ITSUStave[0-9][0-9]/ITSUHalfStave?/ITSUModule?/ITSUChip?` |

For example, replace a `match` value with one of the strings above in a rule such as `{ "match": "...", "rigidBody": ["TX", "TY", "TZ"] }`. ITS3 uses different names, such as `ITS3CarbonForm`, and is not selected by these patterns.