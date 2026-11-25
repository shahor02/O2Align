# Studying the ITS sagitta with the q/pT spectrum

[macro/getSagitta.C](../macro/getSagitta.C) estimates the curvature bias of ITS-only tracks from the `res` tree of the
CheckResid study (`o2-check-resid-workflow`, `checkResid_<period>_<run>_*.root`).

The bias is caused by a curl/sagitta-like deformation, r·Δφ ∝ r². The leave-one-out residuals of checkResid are blind
to it: a deformation that keeps every track on a helix only changes the fitted q/pT.

## Method

A geometric curvature bias shifts every measured q/pT by the same amount δ: `q/pT_meas = q/pT_true + δ`.

The q/pT spectrum has a deep hole around 0, because the spectrum falls steeply with pT. Apart from the charge ratio,
it is mirror symmetric around the true zero. For each trial δ:

1. **Count tracks in mirrored windows.** For every window k in |q/pT − δ|, from 0.06 to 0.85 (GeV/c)⁻¹ (pT ≈ 1.2–17 GeV/c):
   - N+(k) = N(δ + [a_k, b_k]),
   - N−(k) = N(δ − [b_k, a_k]).

   The counts come from linear interpolation of the cumulative spectrum, so δ is not quantised by the histogram bins.
2. **Fit the ratio.** Fit N+(k)/N−(k) with a constant, which is the fitted charge ratio.
3. **Pick δ.** The estimate is the δ that minimises the χ² of that constant fit. The error is the half width of the
   χ²_min + 1 interval.

The scan goes over ±0.025 in steps of 1e-4 (GeV/c)⁻¹. The same estimate is made for all φ together and for `nPhi`
equal sectors of the azimuth of the track direction at the vertex.

The track q/pT is that of `tr.track`: the ITS-only inward fit at the vertex, before the vertex update, so the vertex
does not constrain the curvature.

**Sensitivity.** δ changes the mirror ratio by about 2nδ/|q/pT|, where n ≈ 2.6 is the local slope of the spectrum.
The effect is largest at high pT. One unit of the curl mode k (r·Δφ = k·10⁻⁶ cm⁻¹·r²) corresponds to about
13e-4 (GeV/c)⁻¹ at 0.5 T.

### Charge abundances and other non-geometric effects

The charge ratio is a free parameter, so a pT-independent charge asymmetry is absorbed. Only a change of N+/N− with pT
inside the windows biases δ: a 1 % drift across the window fakes about 1.7e-4 (GeV/c)⁻¹. Three checks:

- **Field polarity.** A geometric bias flips sign with the field; a physics charge asymmetry does not. For the same
  geometry, (δ(B+) + δ(B−))/2 bounds the non-geometric part.
- **Window restriction.** Restricting the windows to the low-pT or the high-pT part must give compatible values.
- **Same tracks, two geometries.** The difference between two geometries refitted on the same data is free of any
  abundance effect.

The per-sector values are less reliable as absolute numbers. Charge-dependent acceptance holes (dead chips, the
cluster requirement, TPC sector edges) remove positive and negative tracks at different φ, because the two charges
bend in opposite directions. This cancels when integrating over φ, but not per sector, and it flips with the field like
a geometric bias. Sector-by-sector differences between geometries on the same data remain valid.

## Usage

```bash
root -l -b -q '$O2Align_ROOT/share/macro/getSagitta.C+("checkResid_LHC25as_568664_AlgRS_PV1.root", 6, "hq2pt.root")'
```

Arguments:
- the checkResid tree file;
- the number of φ sectors (it must divide 36);
- the output file for the `hq` histogram of q/pT vs φ (empty string: not written).

Output:
```
delta(q/pT) [1e-4 (GeV/c)^-1] (+- chi2+1 half width), N+/N- for all phi; file ...
all phi                          [  0, 60)      [ 60,120)      [120,180)      [180,240)      [240,300)      [300,360)
 -28.0 +-  3.0  N+/N- 1.006  -20.0 +-  6.0   -0.0 +-  6.5  -53.0 +-  9.0  -18.0 +-  7.0  -14.0 +- 12.0  -67.0 +- 10.0
```

Run it on several geometries refitted on the same data, and on B+ and B− data with the same geometry, to separate the
geometric part from the rest.

### Example results (November 2026)

| geometry | data | δ [1e-4 (GeV/c)⁻¹] |
|---|---|---|
| IBEdge | pp B+ 562837 | −24 ± 4 |
| IBEdge | pp B− 564644 | +30 ± 3 |
| v10n | PbPb 568664 | −14 ± 3 |
| IBEdge | PbPb 568664 | −28 ± 2.5 |
| solPbHBA0v2b | PbPb 568664 | −28 ± 3 |
| solPbHBA0v2 | PbPb 568664 | −46 ± 2.5 |
| solPbHBA0 | PbPb 568664 | +116 ± 3 |

- The sign flip between B+ and B− shows that IBEdge has a geometric sagitta, k ≈ 2.
- HBonly_v2 adds about Δk ≈ +1.4 on top of it; HBonly_v2b leaves it unchanged.
- The HB-only fit solPbHBA0 fell into the helix family, with k ≈ −11.
