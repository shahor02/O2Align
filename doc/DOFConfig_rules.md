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

## Pinning the Common Mode of Children

`Volume::writeRigidBodyConstraints` forces, for every *free* rigid-body DOF of a parent, the
weighted mean of its active children's movement in that DOF to vanish. This breaks a real
degeneracy: moving the parent and moving all its children together in the same direction are
otherwise indistinguishable. When the parent DOF is fixed instead, that degeneracy doesn't arise,
so no such constraint is written by default — the children are then free to collectively drift in
that DOF, which may be a genuinely unconstrained (or only weakly constrained) mode of the fit.

The `pinChildrenMean` clause asks for the same mean-zero constraint on a DOF of the parent even
though that DOF is fixed there, to remove such a mode explicitly. It takes a boolean (all DOFs), the
string `"all"`, or an array of DOF names (`TX,TY,TZ,RX,RY,RZ`), and — like `rigidBody`/`calib` — a
later matching rule replaces the flags of an earlier one rather than merging with them:

```json
{ "match": "TRD/sm[0-9][0-9]", "rigidBody": "fixed", "pinChildrenMean": ["TY", "RZ"] }
```

With the chambers free, this pins the mean `TY` shift and `RZ` rotation of the chambers of every
supermodule to the (fixed) supermodule position, without making the supermodule itself a fit
parameter. The flags don't depend on the volume's `rigidBody` clause: they also work on a volume with
no rigid-body DOF set at all. A pinned DOF that no free child DOF contributes to (no active children,
or none of their free DOFs feeds it through the child-to-parent jacobian) is skipped with a warning.
`pinChildrenMean` is ignored, with a warning, on a volume that is not rigid-body alignable (the root,
or an envelope without its own geometry such as `TRD_envelope`), so the common mode of the
supermodules of a detector can't be pinned this way.

## Aligning a Subtree as One Rigid Body

The opposite situation to `pinChildrenMean` is freeing a parent's rigid-body DOF while giving none
of its descendants a `rigidBody` rule at all. This is not auto-disabled and needs no mean-of-children
constraint: the parent-to-child jacobian is applied unconditionally down to every leaf measurement
regardless of whether a descendant carries a DOFSet, so the parent's free DOF is fully observable on
its own. A constraint is only needed to resolve the degeneracy of a DOF being free at *both* a parent
and a child level simultaneously (the case this section opened with); with nothing free below the
parent, there is no such degeneracy. Use this to align, say, a half-barrel as a single rigid body
while keeping its staves/modules/chips fixed to it, without needing to give them any DOF config.

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

## TRD Alignable Match Strings

TRD has a dummy `TRD_envelope` (not a geometry alignable), followed by alignable supermodules and chambers. The geometry log uses two digits for supermodule IDs and one digit each for stack and plane IDs:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TRD_envelope` |
| Supermodule | `TRD/sm[0-9][0-9]` |
| Chamber | `TRD/sm[0-9][0-9]/st?/pl?` |

## TOF Alignable Match Strings

TOF has a dummy `TOF_envelope` (not a geometry alignable), followed by alignable supermodules and strips. Both IDs in the geometry names are formatted as two digits:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TOF_envelope` |
| Supermodule | `TOF/sm[0-9][0-9]` |
| Strip | `TOF/sm[0-9][0-9]/strip[0-9][0-9]` |

The envelopes are not rigid-body alignable; rigid-body clauses on them are ignored. Calibration clauses can be configured on these synthetic volumes.