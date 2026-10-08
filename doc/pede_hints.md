# Reading pede output

`millepede.log` describes the fit and its numerical progress; `millepede.his`
contains diagnostic histograms **and graphs**, in text format. Neither alone
establishes that an alignment is physically correct. Keep the steering files,
solver version, `millepede.res`, and captured standard output/error with each run:
some final warnings are printed to the terminal rather than the log.

This guide concerns Millepede II. Titles, available plots, and stop codes depend
on the version, method, and steering options. The histogram IDs below were checked
against the local Millepede-II `pede.f90`; prefer the titles and the explanations
printed by your executable.

## 1. Read `millepede.log`

Start at the end for completion, rejection statistics, and the solver's
**“Explanation of iteration table”**, then check the initialization and iterations.

| Check | What deserves attention |
|---|---|
| Inputs and steering | Correct binary files, record counts, solver method, cuts, priors, and constraints. Unexpectedly low counts can mean missing files or read errors; some read errors are treated as EOF. |
| Global parameters | Expected number of variable/fixed parameters and labels excluded by the `entries` threshold. Entries count connected measurements, not necessarily distinct tracks. Check coverage again after rejection: initially sufficient statistics may disappear. |
| Constraints and rank | Empty, dependent, or inconsistent constraints; rank defects of the global matrix or parameter groups. Investigate missing reference freedoms, disconnected detector parts, redundant parameterizations, and poor illumination. A solver returning a solution to a singular system does not make all parameters identifiable. |
| Local-fit failures | Separate failures of the track fit/conditioning from rejection by residual, chi-square, or downweight-fraction cuts. Large fractions suggest bad tracks, wrong errors, derivatives, units, or an inadequate model. |
| Final fit quality | Compare the reported chi-square with its reported NDF and read any downweighting correction. A ratio near one is useful only with an appropriate error model and understood selection. |

### Iteration table

| Column | Meaning and practical reading |
|---|---|
| `it`, `fc` | Global iteration and function-evaluation count. Iteration 0 is the initial evaluation; line searches can require several evaluations per iteration. |
| `fcn_value` | Twice the objective function. In ordinary least squares this is chi-square; with robust weights/penalties its interpretation changes. Look for improvement and eventual stability. |
| `dfcn_exp` | Expected reduction of the objective. Compare it with the printed convergence limit (`Delta F`); a small observed change alone does not establish convergence. |
| `step`, `ls` | Line-search step multiplier and status. A full Newton step is usually 1. Repeated tiny steps, exhausted searches, rounding limits, or non-descending directions warrant investigation. `ls=1` denotes successful line-search conditions in the checked version. |
| `iit`, `st` | Inner iterative-solver steps and stop code. Decode `st` using the log's method-specific legend: MINRES and MINRES-QLP use different codes. Iteration limits, ill-conditioning, or symmetry/preconditioner failures require follow-up. |
| `cutf`, `rejects` | Chi-square cut factor and number of rejected local-fit objects (usually tracks). In the checked version, `cutf=1` uses the NDF-dependent three-sigma chi-square threshold, and 0 disables this cut. Compare fractions using the appropriate record count. |

Changing cuts or robust weights changes the effective sample/objective, so a drop
in `fcn_value` need not be an alignment improvement. Check that the accepted sample
and rejection rate settle. The solver traditionally warns above 1% rejection;
this is a diagnostic convention, not a universal acceptance criterion.

Distinguish convergence from reaching an iteration limit or stopping in
`subito`/`-s` mode. Internal pede iterations reuse the supplied linearization;
large geometry corrections can require applying the alignment, reconstructing
tracks, and producing new Mille records.

## 2. Read and plot `millepede.his`

Use the converter supplied with your Millepede installation. Current versions
provide:

```sh
readPedeHists -i millepede.his -r millepedeHistos.root -p millepedeHistos.pdf
```

Older installations provide `tools/readPedeHists.C`. In ROOT, load that macro
using its installed path, then call:

```cpp
.L /path/to/Millepede-II/tools/readPedeHists.C+
readPedeHists("print write", "millepede.his");
```

The older macro produces PostScript and ROOT output; the newer executable supports
PDF. `histprint` in the steering file additionally prints ASCII histograms.

**Check the title, snapshot/version, entries, range, and underflow/overflow.**
The same histogram ID may occur repeatedly for different data passes; compare
initial and later snapshots rather than merging them. Histograms are filled at
different stages of rejection/downweighting, so do not assume every plot describes
only the final accepted tracks. Missing pulls can reflect disabled covariance
calculation or a limited number of pull-producing passes.

| Histogram title / usual ID | Interpretation and warning signs |
|---|---|
| Entries per label (2) | Measurement coverage of global labels. A long low-statistics tail motivates checking which detector elements/DOFs are poorly covered and which fail `entries`. This distribution does not identify individual labels. |
| Normalized residuals, global/local (3/12) | Post-local-fit residual divided by the input measurement sigma. “Global” means the measurement has global derivatives; “local” means it has none, e.g. a local scattering constraint. Check shifts, asymmetry, broad tails, and clipping. These are not global-parameter residuals. |
| Pulls, global/local (13/14) | Residual divided by its residual uncertainty, accounting for the local fit covariance. In a well-modelled Gaussian sample, expect approximately zero mean and unit width. Ordinary normalized residuals can be narrower than one because fitting consumes freedom; do not force their width to one. Pulls here do not include the full fitted global-parameter covariance. |
| Chi-square/NDF after local fit (4) | Broad tails indicate poor track fits, outliers, or underestimated errors. Width depends strongly on track NDF; low-NDF tracks naturally have broad distributions. Selection and robust weighting distort the nominal chi-square distribution. |
| NDF / number of local parameters (5/11) | Verify expected track structure and useful redundancy. Unexpected populations can reveal incomplete records or an incorrect local model. |
| Downweight fraction (6) | A large fraction of downweighted measurements means the fit relies strongly on robust treatment. Check affected track/file populations instead of merely relaxing the rejection cut. |
| Local-fit log10 condition (16) | Large values indicate sensitivity to rounding and a nearly degenerate local problem. Examine parameter scales, track lever arms, and local-model redundancies. |
| Positive eigenvalues (7), low-end eigenvalue graph | Available for suitable methods. Very small eigenvalues indicate weakly determined combinations. Inspect eigenvectors to identify modes; positive-only histograms cannot reveal all zero/negative modes. Numerical scales depend on parameter units. |

Also inspect the **function-value and inner-solver iteration graphs**, residual/
pull location and dispersion versus **record number**, and, where available,
rejection fraction and mean chi-square/NDF versus **input file**. Abrupt changes
can expose a bad file, run period, or mixed error models that pooled distributions
hide. Record number is input order, not automatically time or detector position.

## 3. Before using the alignment

- Resolve numerical failures, unexplained rank defects, and empty/redundant constraints.
- Verify convergence, adequate accepted statistics, and understood rejection/downweighting.
- Investigate chi-square/pull discrepancies before rescaling uncertainties: model,
  material, correlations, derivatives, and selection can all be responsible.
- Inspect `millepede.res` for implausible corrections and uncertainties, respecting
  the application's units/sign conventions. Error availability depends on method;
  a missing or zero error is not evidence of perfect precision.
- Validate on independent tracks and detector/kinematic subdivisions. Weak modes
  can preserve small residuals while biasing momentum or geometry; use closure
  tests, external references, or physics observables where appropriate.

## References

- [Official Millepede II manual](https://millepede.pages.desy.de/millepede-ii/draftman_page.html):
  steering, iterations, outlier treatment, constraints, and output files.
- [pede source and iteration-table legend](https://www.desy.de/~kleinwrt/MP2/doc/html/pede_8f90_source.html):
  histogram booking/filling, rejection logic, and method-dependent stop codes.
- [Current histogram converter usage](https://millepede.pages.desy.de/millepede-ii/readPedeHists_8cpp.html).
