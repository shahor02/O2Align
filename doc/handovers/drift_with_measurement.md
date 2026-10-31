# Handover: ITS envelope drifts despite a `Measurement` on its TX/TY/TZ

Date: 2026-10-20.

## Problem

Config `doc/ITS_RB_from_staves.json`: the ITS envelope is free and carries `Measurement 0 0.01` on
TX/TY/TZ (labels 268435712..714, written to `mp2con.txt` by the new `measurement` clause).
Pede output (/data0/devAlg/testNVMultiSingleNoPV/millepede.*) nevertheless gives corrections of about
-1.07, -1.02, -42.5 cm for ITS TX/TY/TZ and -1.03, -1.20, -42.6 cm for the mean vertex (MV, labels
134218112..114).

**Caveat found afterwards:** in the config used, the measurement was in the same rule as
`"rigidBody": ["TX","TY","TZ"]`, which also fixed ITS RX/RY/RZ. The runs below were done with that
config, so the ITS rotations were fixed. This has not been rechecked with rotations free. The mix is now
rejected (fatal) and the JSON split into a DOF rule and a measurement rule.

## Implementation (uncommitted)

- `Constraint(name, value, sigma = -1)`: with sigma > 0 `write()` emits `Measurement value sigma`
  (`Constraint.h/.cxx`, ClassDefNV 3).
- `Volume::applyMeasurementConfig`, `writeRigidBodyDOFMeasurements`, member `mDOFMeasurement`; called from
  `applyDOFConfig` and `writeRigidBodyConstraints` (leaves included). Documented in `doc/DOFConfig_rules.md`.
- `applyDOFConfig` rejects a rule having both `rigidBody` and `measurement`.

## Findings

- Pede does read the 3 measurements (pede.f90:3706 adds w(sum c p - value)^2, w = 1/sigma^2).
  The syntax and the output are fine.
- Without the measurements (run B): "Rank defect 2", ITS TX/TY/TZ undetermined. The ITS+MV common
  translation is a flat mode of the ITS-only fit, so the measurement is needed.
- With them (run A, all 1155 parameters, inversion): ITS (1.556, -1.599, -42.58) cm, MV (1.59, -1.78, -42.71);
  ITS error is exactly 0.0100, i.e. the data add no curvature along this mode. MV - ITS equals run B
  (0.035, -0.18, -0.127), so only the common mode runs away. fcn goes 4.37e7 -> 1.07e7 -> 2.2e7 -> 2.85e7 -> 2.9e7,
  not converging cleanly. Not a MINRES issue (inversion behaves the same).
- Run C (only the 6 ITS+MV translations free, the rest fixed): the same runaway.
- Sigma scan (C, C_1e-3, C_1e-5): the ITS shift scales as sigma^2 (p/sigma^2 is constant, about 7.7e4 TX,
  3.4e5 TY, -6.2e5 TZ). That constant is the net data gradient along the common mode:
  sigma=1e-3 gives (0.078, 0.337, -0.620); sigma=1e-5 gives (7.7e-6, 3.4e-5, -6.2e-5).
- A constant force with zero curvature is impossible for a consistent linear least squares
  (Cauchy-Schwarz), so gradient and Hessian must be inconsistent, or the mode is not really flat.

## Root cause: NOT found. Candidates

1. A vertex measurement that does not move with the MV, e.g. the single-track PV constraint
   (commit ae20205 "Make single-track PV constraint to follow the MV prior").
2. Sign mismatch between the ITS and MV derivatives: `addMeanVertexDerivatives` uses `-trans` because the
   measurement moves with the MV (src/AlignmentSpec.cxx ~1335-1370).
3. Dependence on pede outlier rejection (about 2.7% records rejected): fcn is not comparable between iterations.
4. Float derivatives in the Mille files.
5. The fixed ITS rotations (see the caveat above).

## Suggested next steps

- Rerun with the corrected JSON (ITS RX/RY/RZ free) and compare.
- Dump global derivatives of a few records: the sum over the vertex point and the ITS points of d/dITS-TX
  plus d/dMV-TX must be 0 for a flat mode (the same for TY, TZ).
- Test with the single-track PV constraint disabled.
- Workarounds meanwhile: sigma <= 1e-4 cm, a hard `Constraint`, or fixing the ITS TX/TY/TZ.

## Scratch material

/tmp/claude-1000/-home-shahoian-alice-O2Align/8e105e27-7f84-41dc-9007-bc5cb43630ac/scratchpad/{A,B,C,C1,C_1e-3,C_1e-5}
(symlinks to mp2data_*.bin, mp2param.txt, mp2con.txt, mpSteer.txt, plus pede.log and millepede.res).
Pede: `eval "$(alienv printenv Millepede-II/latest)"`, source in /home/shahoian/alice/Millepede-II/pede.f90.
