# Rigid-body consistency test of the IB chips from the Mille data (`macro/milleSplit*`)

Written for: whoever continues the ITS alignment tuning (knows O2Align, Millepede-II/GBL and checkResid).

## 1. What the test asks

A real displacement of a chip moves its measurement plane, so it shows up in both measured coordinates of
every hit. A normal shift δX, for example, gives dy = δX·tan φ_inc and dz = δX·tan λ_inc. The same holds for
every rigid-body DOF: the y rows and the z rows of a chip must ask for the same parameter values.

The test reads the Mille records written by `--output MilleData` and profiles out each track's local
parameters. It then fits every ITS IB chip (432 sensors) with two independent 6-DOF sets (TX, TY, TZ, RX,
RY, RZ), one seen only by the chip's y rows and one only by its z rows. Comparing this split model with the
rigid one (both copies tied) measures how much of the data no rigid-body motion can describe, independently
of the DOF configuration actually used in a pede run.

The test needs no pede, constraint file or run configuration. It does need the Mille data to carry all six
derivatives of every IB chip (record with `doc/ITS_RB_all.json` or another wide-open `dofConfigJson`), and a
`mp2tree.txt` from the same hierarchy, so that the labels resolve.

## 2. Files

| file | role |
|---|---|
| `macro/milleSplit.cxx` | stand-alone C++ (no ROOT/O2). Reads a Mille binary and writes the reduced normal equations of the 432×12 split parameters |
| `macro/milleSplitModels.py` | `labels`: writes the chip label → sensor-ID map from `mp2tree.txt`. `fit`: solves the models from one or more accumulations |
| `macro/sensDump.C` | ROOT macro (O2 environment): writes the per-sensor table (layer, half barrel, stave, chip in stave, ref. α, ref. x, centre) that `fit` uses for grouping |

## 3. Running it

```bash
cd macro
g++ -O2 -std=c++17 -o milleSplit milleSplit.cxx
python3 -I milleSplitModels.py labels /path/to/mp2tree.txt iblabels.txt
./milleSplit /data/devAlg/milleB0/mp2data_0.bin iblabels.txt ms_B0.bin        # optional 4th arg: max records
root -l -b -q 'sensDump.C+("o2sim_geometry-aligned.root","sens.txt")'         # in the O2 environment
python3 -I milleSplitModels.py fit sens.txt ms_B0.bin [ms_BPlus.bin ...]       # several files are summed
```

- `labels` should report 432 labels. Any other count means the tree doesn't belong to an IB hierarchy with
  per-chip volumes.
- `milleSplit` takes about 20 min and 1 GB of memory per 6 GB file, and writes a ~215 MB file. Run one
  process per file in parallel.
- `fit` holds several dense 5184×5184 matrices: allow ~10–20 min and a few GB per call. Set
  `OMP_NUM_THREADS` when several run in parallel.

### `milleSplit` details

- **Records:** each Mille record is one GBL trajectory: a single track, or a multi-track vertex. Its local
  parameters (track parameters, kinks, vertex) are profiled out exactly: N_gg −= C·A_ll⁻¹·Cᵀ, as pede does
  internally.
- **Skipped records:** more than 600 local parameters, a singular local system, or χ²/ndf > 6 after the
  local fit (an outlier guard). The last log line gives the used and rejected record counts and ⟨χ²/ndf⟩.
- **y/z rows:** a row is a "z row" of a chip when |∂/∂TZ| > |∂/∂TX|, otherwise a "y row".
- **Mean vertex:** its calibration globals (X, Y, Z, SlopeX, SlopeY of up to 8 time slots) are accumulated as
  free parameters, after the chip ones, and `fit` frees them in every model. The vertex records measure the
  mean vertex relative to its reference at (0, 0), while the 2025 beam sits at about (+0.3, −1.8) mm. Keeping
  it fixed pulls every vertex by ~40σ, transversely only, and fakes a large rφ-only chip term.
- **Other rows:** OB, TPC, vertex and every other row help determine the local parameters, but carry no
  parameters of their own: these globals stay at their recorded reference values. Check that this is
  acceptable (compare with `millepede_labeled.res`) before trusting a result.

## 4. Output of `fit`

1. The χ² decrease (exact, pseudo-inverse) for these models:
   - rigid: 6 DOFs per chip, shared by the chip's y and z rows;
   - rigid plus a y-only TY term β (derivative tan φ_inc on the y rows only), grouped per sensor, per layer
     × chip index, per layer, and per layer × (end chip, middle chip);
   - the full split: 12 DOFs per chip.

   Compare each extra gain with its number of extra parameters (the expected gain if the rigid model were
   enough).
2. β per layer × chip index, averaged over the staves, in µm.
3. For chips 0, 4 and 8 of each layer, the stave-to-stave mean, rms and mean error of the per-sensor β.

Every chip is free in all six DOFs, which leaves global weak modes. Parameter values are therefore quoted
with a 50 µm / 5 mrad prior; the gains do not use it. Only y/z differences and chip-to-chip patterns are
meaningful: the absolute values are not.

## 5. Interpreting the result

- If the split model gains far more than its number of extra parameters, the data are not described by rigid
  chip motion. A y-only term β ∝ tan φ_inc (no dz counterpart) that carries most of the gap is not a
  displacement: it points to the cluster positions (dictionary, centroid) or to a track-frame/incidence effect.
- Per stave, a β that is a first harmonic in the stave azimuth ∝ R·cos α is the signature of a fixed mean
  vertex (see the mean-vertex bullet of section 3), not of the chips.
- The per-sensor β is correlated with the chip RZ (free in the rigid model); check by removing RZ from the
  rigid set before reading it.
- A pooled RY-shaped end-chip tilt cannot represent a y-only ∝ tan φ_inc term.

## 6. Next steps

- **Cluster origin:** split the checkResid residuals by cluster size or pattern (this needs checkResid to
  store them), and compare with MC. If the effect is cluster-level, the remedy is in the cluster positions
  (dictionary or centroid), not in the alignment.
- **Fit parameter:** if one is wanted anyway, it is a calibration DOF on the IB chips with derivative
  tan φ_inc on the y measurement only.
- **Extending the tool:** to another detector, layer or DOF grouping, change `iblabels.txt` (label →
  index) and `NS` in both files.
