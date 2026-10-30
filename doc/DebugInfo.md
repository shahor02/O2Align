# Debug output of tracks and points

The data-processing stage can write, for a random fraction of the tracks, everything needed to look at a track fit in detail:
the track parameters at the reference point and, for each point of the GBL trajectory, the measurement, the residuals before and
after the GBL fit and the kinks. The records are `DbgTrack` and `DbgPoint` ([AlignmentDbg.h](../include/O2Align/AlignmentDbg.h)).


## Activation

```bash
... --configKeyValues "AlignParams.debugTracksFraction=0.05;..."
```

| Parameter | Default | Meaning |
|---|---|---|
| `AlignParams.debugTracksFraction` | `-1` | fraction of the tracks written to the `gbltrack` tree. `<=0`: none, `>=1`: all |

- The selection is random, one draw per track, with a generator seeded by the lane number, so that a run is reproducible and the lanes
  do not select the same tracks.
- The output file is `debug_alg.root`, or `debug_alg_<lane>.root` if the workflow runs with more than one lane (`--pipeline
  barrel-alignment:N`), in the working directory. Nothing is created when there is nothing to write: the file exists only if
  `debugTracksFraction>0` or if a misalignment is injected with `--output MisRes` (it then holds the `mis` tree, see
  [MisalignmentTest.md](MisalignmentTest.md#debug-output)). The `MilleRes` stage never creates it.
- The cost is paid only by the selected tracks. The tracks are selected before the GBL points are built, whether or not they are
  eventually accepted.
- Tracks rejected by the `maxChi2Ndf` cut of the GBL fit are written too, with `accepted=false`: they are often the ones to look at.
  Tracks whose GBL trajectory could not be built or fitted at all (too few points, a failed propagation, a failed fit) are not written.

Reading it needs the library dictionary:

```cpp
gSystem->Load("libO2Align");
TFile f("debug_alg.root");
auto t = (TTree*)f.Get("gbltrack");   // branch "t" of class o2::alignrs::DbgTrack, one entry per track
t->Draw("points.dyFit*1e4", "(points.flags&1) && t.accepted");  // fitted dy residuals of the measured points, in um
```

## `gbltrack` tree: `DbgTrack`

`DbgTrack` derives from `o2::track::TrackPar` (float): its parameters are the track state at the **reference point**, the innermost
point of the trajectory (the vertex point if the track has one, otherwise the innermost cluster), as refitted by the KF before the
GBL fit. It is the seed of the GBL linearization; the GBL correction is not applied to it.

| Member | Meaning |
|---|---|
| (`TrackPar`) | `X`, `alpha`, `Y`, `Z`, `snp`, `tgl`, `q/pt` at the reference point, e.g. `t.getX()`, `t.getPt()`, `t.getSnp()`; the parameters are in the tracking frame `alpha` |
| `gid` | `GlobalTrackID` of the track |
| `kfFit` | KF refit: `chi2`, `ndf`, `chi2Ndf`. The KF itself accumulates only the chi2; the ndf is computed here as the number of measured coordinates (2 per point, 1 for a point with the Z taken from the track, the vertex point included) minus 5 |
| `gblFit` | GBL fit: `chi2`, `ndf`, `chi2Ndf`. Filled also for the rejected tracks. For a track of a composed trajectory (see below) this is the fit of the **whole** trajectory, shared by all its tracks |
| `accepted` | the GBL fit passed the `maxChi2Ndf` cut and the trajectory went to the Millepede output |
| `composed` | the track was fitted in a composed trajectory with the other tracks of its vertex (`useMultiTrackPVConstraint` 1 to 4); otherwise it is a single-track trajectory |
| `bz` | nominal field, kG (0 for B=0 data) |
| `tfCounter`, `firstTForbit`, `runNumber` | timing info of the TF |
| `points` | `std::vector<DbgPoint>`, the points of the GBL trajectory in the order GBL sees them (from the innermost outward) |

## `DbgTrack::points`: `DbgPoint`

One `DbgPoint` per **GBL point**, which is not the same as one per cluster:

- a point with a measurement (`HasMeas`) is a cluster, or the vertex point if it is measured;
- the vertex point of a track in a composed trajectory (`IsVertex` w/o `HasMeas`) carries no measurement of its own, the common
  vertex position enters through the inner transformation;
- when a step to a point spans a wide radial gap with a lot of material (the ITS-TPC transition, see `splitMSThreshold` and
  `splitMSStepMinDX`), an extra scatterer-only point is inserted at the midpoint. It has `det=-1`, only `x`, `alpha`, `snp` and
  the kink fields are meaningful.

| Member | Meaning |
|---|---|
| `det` | detector index of [Detector.h](../include/O2Align/Detector.h): 0 = PVT (the vertex), 1 = ITS, 2 = TPC, 3 = TRD, 4 = TOF; -1 for the scatterer-only points |
| `sens` | sensor ID within the detector (the volume ID of the Millepede label), -1 for the vertex |
| `lr` | detector-specific layer-like index (ITS layer, TPC pad row group, ...) |
| `flags` | `HasMeas` (1): has a measurement; `HasKink` (2): has a scatterer, the kinks are defined; `IsVertex` (4): the primary vertex; `ITSOverlap` (8): the cluster is on an ITS overlap |
| `x`, `alpha` | position of the point in its tracking frame (`x` along the local radial direction) and the rotation angle of that frame |
| `y`, `z` | the measured coordinates in the tracking frame (0 if no measurement) |
| `sigY`, `sigZ` | cluster errors, including the extra errors (`extraClsErr*`) |
| `snp` | sin of the track azimuthal angle in the tracking frame at the point, from the linearization track |
| `dyPre`, `dzPre` | residual that was fed to GBL: measurement minus the linearization track (see below) |
| `dyFit`, `dzFit` | residual after the GBL fit: measurement minus the fitted track |
| `kinkY`, `kinkZ` | fitted kink of the scatterer at the point, along the tracking-frame `snp` and `tgl` |
| `pullKinkY`, `pullKinkZ` | kink divided by its fitted error |

The methods `getLabX()`, `getLabY()`, `getLabR()` give the lab-frame coordinates `x cosα - y sinα`, `x sinα + y cosα` and `sqrt(x²+y²)`
of the measured point; the lab `Z` is `z`. `hasMeas()` and `hasKink()` test the flags.

### How to interpret the values

**Residuals.** All residuals are *measurement minus track* in cm, in the tracking frame of the point.

- `dyPre`, `dzPre` are relative to the *linearization track* of the GBL, which is propagated from point to point and, when
  `useStableRef` is set (the default), displaced to each cluster after the point is accounted. They are therefore the residuals with
  respect to the broken-line reference, they can be as large as mm (the RMS in the first tests was 0.9 mm in Y and 1.4 mm in Z) and they are
  not track residuals in the usual sense. They are the input to the fit: if a misalignment is injected, it is included in them.
- `dyFit`, `dzFit` are the residuals with respect to the GBL-fitted track. They are *biased*: the point itself pulls the fit, so their
  RMS is below the cluster error (about 0.5 σ for ITS tracks with ~7 points), as for any residual to a fit that includes the point.
  Their mean as a function of position (sensor, phi, z, ...) is what shows the misalignment that the fit could not absorb. The unbiased
  version is not stored.
- Residuals of the same point in a composed trajectory belong to a fit with the common vertex of all its tracks.

**Kinks.** The kink is the change of direction at the scatterer, which stands for the material crossed since the previous point
(thin scatterer at the point). `kinkY` and `kinkZ` are the fitted values (not GBL's "residual", which has the opposite sign) along the
`snp` and `tgl` of the tracking frame, in rad. They are expected to be zero on average and to have the width given by the
multiple scattering of the crossed material. The pull should have unit RMS if the material and the cluster errors are right: a
pull RMS above 1 means that the scattering is underestimated or that the points around it are not well described. The pulls are
given in the diagonal basis of the scatterer, which is the `snp`/`tgl` basis here.

**Fit quality.** `kfFit.chi2Ndf` is the chi2 per ndf of the KF refit; `gblFit.chi2Ndf` should be close to 1 for good tracks. For the tracks of a
composed trajectory all tracks have the same `gblFit`.

## Examples

Distribution of the fitted residuals of the accepted tracks, in um, versus the ITS layer:

```cpp
t->Draw("points.dyFit*1e4:points.lr", "(points.flags&1) && t.accepted && points.det==1", "colz");
```

Mean fitted residual in the ITS sensors (a candidate misalignment) and its dependence on the sensor ID:

```cpp
t->Draw("points.dyFit*1e4:points.sens", "(points.flags&1) && t.accepted && points.det==1 && points.lr==3", "prof");
```

Kink pulls of the accepted tracks, the points with an available kink only:

```cpp
t->Draw("points.pullKinkY", "(points.flags&2) && points.pullKinkY!=0 && t.accepted");
```

The `Scan` of one track with its points:

```cpp
t->Scan("points.det:points.sens:points.x:points.dyPre:points.dyFit:points.kinkY", "Entry$==0", "", 1);
```

Tracks of poor fit and their worst point:

```cpp
t->Scan("t.gid.asString():t.gblFit.chi2Ndf:Max$(abs(points.dyFit/points.sigY))", "!t.accepted");
```

## The `mis` tree

With `--output MisRes` and `AlignParams.misAlgJson`, one entry per point to which a misalignment was injected is written to the same
file: the shifted residual `dy,dz`, the injected shift `ddy,ddz`, the detector `det`, sensor `sens`, `z` and `phi`. See
[MisalignmentTest.md](MisalignmentTest.md).
