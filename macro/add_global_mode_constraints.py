#!/usr/bin/env python3
"""Append global ITS deformation-mode constraints to an mp2con.txt written by MilleSteer.

Stop-gap equivalent of the "globalModes" clause of the DOF config (Volume::writeGlobalModeConstraints),
which derives the coefficients from each volume's own getL2G() and also covers the rotation DOFs; prefer
that. This script is useful for trying a mode on an mp2con.txt that already exists, with no new MilleSteer
pass. Its two approximations: it uses nominal layer radii and infers phi from the stave index, so the
azimuthal phase is offset from the real geometry, and the IB stave tilt (which gives radial modes a
projection on the stave TX) is missing.

The hierarchical "mean of children = parent" constraint (Volume::writeRigidBodyConstraints, pinChildrenMean)
is a *rigid-body* constraint: 6 scalar equations per parent, removing only the net translation/rotation of
that parent's set of free descendants. Deformation modes of the barrel as a whole -- uniform radial
scaling, elliptical (2-phi) squashing, twist -- carry no net rigid component per parent and pass through it.
They are also (very nearly) null modes of a helix fit, so they need explicit constraints.

This script writes these constraint equations directly into the steering file.

Coefficients use the local-frame convention of RigidBodyDOFSet for ITS sensors:
local X = r-phi (in plane), local Y = sensor normal (~ radial), local Z = beam.

Usage:
    add_global_mode_constraints.py mp2param.txt mp2con.txt [-m radial,ellipse] [-o out.txt]
"""

import argparse
import math
import re
import sys

# nominal ITS2 mid-layer radii [cm]; only the ratios matter for a constraint
LAYER_RADIUS = {0: 2.34, 1: 3.15, 2: 3.93, 3: 19.6, 4: 24.55, 5: 34.39, 6: 39.34}
# staves per half-barrel, used to turn a stave index into an azimuth
STAVES_PER_HALFBARREL = {0: 6, 1: 8, 2: 10, 3: 12, 4: 15, 5: 21, 6: 24}
# highest azimuthal order of the harmonic families; keep in step with Volume::GlobalHarmonicLMax
LMAX = 4

PARAM_RE = re.compile(
    r"^\s*(\d+)\s+\S+\s+\S+\s+!\s+(?P<fv>[VF])\s+(?P<dof>\w+)\s+(?P<path>\S+)\s*$"
)
STAVE_RE = re.compile(
    r"^ITS/ITSULayer(?P<lay>\d+)/ITSUHalfBarrel(?P<hb>\d+)/ITSUStave(?P<st>\d+)$"
)


def read_free_stave_dofs(path):
    """label -> (layer, half-barrel, stave index, dof name), for free stave-level DOFs only."""
    out = {}
    with open(path) as f:
        for line in f:
            m = PARAM_RE.match(line)
            if not m or m.group("fv") != "V":
                continue
            s = STAVE_RE.match(m.group("path"))
            if not s:
                continue
            out[int(m.group(1))] = (
                int(s.group("lay")),
                int(s.group("hb")),
                int(s.group("st")),
                m.group("dof"),
            )
    return out


def stave_phi(lay, hb, idx):
    """Azimuth of a stave centre, assuming staves are numbered consecutively in phi."""
    return math.pi * (idx + 0.5) / STAVES_PER_HALFBARREL[lay] + hb * math.pi


def build(dofs, modes):
    """Return [(title, [(label, coefficient), ...]), ...]."""
    cons = []

    def add(title, dof, weight):
        terms = [
            (lab, weight(lay, hb, idx))
            for lab, (lay, hb, idx, d) in sorted(dofs.items())
            if d == dof
        ]
        terms = [(lab, c) for lab, c in terms if abs(c) > 1e-12]
        if terms:
            cons.append((title, terms))

    if "radial" in modes:
        # uniform radial scaling dr = eps * r: a pure 2D dilation maps both lines and circles
        # onto lines and circles, so it is absorbed entirely by the 1/pt scale and is invisible
        # to the track chi2 -- the dominant unconstrained mode.
        add("global radial scaling (dr = eps*r)", "TY", lambda lay, hb, i: LAYER_RADIUS[lay])
    if "curl" in modes:
        # rigid rotation about the beam axis, d(r*phi) = eps * r
        add("global curl (r*dphi = eps*r)", "TX", lambda lay, hb, i: LAYER_RADIUS[lay])
    if "ellipse" in modes:
        # elliptical (2-phi) squashing of each layer, dr = eps * cos/sin(2 phi)
        add(
            "global elliptical, cos2phi (dr = eps*cos2phi)",
            "TY",
            lambda lay, hb, i: math.cos(2 * stave_phi(lay, hb, i)),
        )
        add(
            "global elliptical, sin2phi (dr = eps*sin2phi)",
            "TY",
            lambda lay, hb, i: math.sin(2 * stave_phi(lay, hb, i)),
        )
    if "zscale" in modes:
        add("global z scaling (dz = eps*z)", "TZ", lambda lay, hb, i: LAYER_RADIUS[lay])
    # The generic azimuthal-harmonic families, amplitude proportional to r, matching
    # Volume::RadialHarmonic / Volume::TangentialHarmonic. The radial component projects on the stave
    # TY (sensor normal), the tangential one on TX (r-phi). These are the modes left over once the
    # radial scaling, the curl and the radial l=2 are constrained away.
    for name, dof in (("radialHarmonic", "TY"), ("tangentialHarmonic", "TX")):
        for l in range(1, LMAX + 1):
            if name in modes or f"{name}{l}" in modes:
                for trig, fn in (("cos", math.cos), ("sin", math.sin)):
                    add(
                        f"global {name} l={l} {trig} (eps*r*{trig}({l}*phi))",
                        dof,
                        lambda lay, hb, i, l=l, fn=fn: LAYER_RADIUS[lay] * fn(l * stave_phi(lay, hb, i)),
                    )
    return cons


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("param", help="mp2param.txt from MilleSteer")
    ap.add_argument("con", help="mp2con.txt from MilleSteer")
    ap.add_argument("-m", "--modes", default="radial,ellipse",
                    help="comma-separated subset of radial,curl,ellipse,zscale,radialHarmonic[l],"
                         "tangentialHarmonic[l] (default: radial,ellipse)")
    ap.add_argument("-o", "--output", help="output file (default: overwrite <con>)")
    args = ap.parse_args()

    modes = {m.strip() for m in args.modes.split(",") if m.strip()}
    known = {"radial", "curl", "ellipse", "zscale"}
    for name in ("radialHarmonic", "tangentialHarmonic"):
        known |= {name} | {f"{name}{l}" for l in range(1, LMAX + 1)}
    if not modes <= known:
        sys.exit(f"unknown mode(s) {sorted(modes - known)}, known are {sorted(known)}")

    dofs = read_free_stave_dofs(args.param)
    if not dofs:
        sys.exit(f"no free stave-level rigid-body DOFs found in {args.param}")
    cons = build(dofs, modes)

    with open(args.con) as f:
        body = f.read()
    out = args.output or args.con
    with open(out, "w") as f:
        f.write(body)
        if not body.endswith("\n"):
            f.write("\n")
        for title, terms in cons:
            f.write(f"\n!!! {title}\nConstraint 0\n")
            for lab, coef in terms:
                f.write(f"{lab} {coef:.12g}\n")

    nst = len({(v[0], v[1], v[2]) for v in dofs.values()})
    print(f"{len(dofs)} free stave DOFs over {nst} staves; appended {len(cons)} constraints to {out}:")
    for title, terms in cons:
        print(f"   {len(terms):>4} terms  {title}")


if __name__ == "__main__":
    main()
