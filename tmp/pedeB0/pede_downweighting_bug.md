Subject: pede: outlier down-weighting breaks the gradient along directions absorbed by local parameters

Version: V04-17-06 (HEAD 2f0f8a2, main); also seen in 4fd638f. The relevant code in pede.f90 is unchanged between the two.

Summary
-------
With `outlierdownweighting n` (n >= 1), the global normal equations stop being consistent with the
local fits from the 2nd iteration on: the local fit and the global update use different weights for
the same measurements. For any global direction v that the local parameters can fully compensate
(A_g v in span(A_l)), pede then sees a non-zero gradient with zero curvature. If that direction is
held only by a soft `Measurement` with sigma s, the parameters run away by about F*s^2, where F is
the spurious force, until the measurement term balances it. fcn increases from one iteration to the
next.

Where
-----
- Local fit, loopbf, outlier loop (pede.f90:4718-5048): in pass `iter` the weights are computed from
  the residuals of the previous pass (`resid=rmeas-localCorrections`, l.4747). They are Huber for
  iter <= 3 and Cauchy for iter > 3 (l.4754). The local solution blvec of the last pass (l.4778)
  therefore minimises chi^2 with the weights W_(nter-1), taken from the residuals of pass nter-1.
- Global update, "fourth loop" (l.5112 ff.): the weights are recomputed from the *final* residuals,
  and always with Huber (l.5126-5131), even when the last local pass used Cauchy. The global vector
  (l.5143), the global-global block and the global-local block (used for G^T Gamma^-1 G) all use
  these weights W_4. Gamma^-1 itself comes from the local pass, with W_(nter-1).

Then g.v = sum (A_l d)^T W_4 r. This is not zero, because r is orthogonal to A_l with respect to
W_(nter-1), not W_4. The two weights differ for every measurement beyond chuber*sigma: Huber weights
come from different residuals, and in the last pass the formula itself differs (Cauchy vs Huber).
The curvature computed along v is ~0 in both, so the step along v is limited only by the external
measurement.

Reproduction (ALICE ITS alignment, ~150k GBL records, 1155 fit parameters)
--------------------------------------------------------------------------
The ITS envelope translation and the mean vertex position (from a vertex prior in composed GBL
trajectories) have a common translation that is exactly flat in the data. We checked this directly
on the records: A_g v is reproduced by the local parameters; curvature ~1e-19, gradient ~1e-9 along
Z. The ITS translation carries `Measurement 0 sigma`.

| steering                      | ITS TX, TY, TZ result (cm) | fcn over iterations         |
|-------------------------------|----------------------------|-----------------------------|
| outlierdownweighting 4, s=1e-2| 1.06, -3.37, -15.3         | 4.26e6 -> 5.63e6 -> 6.67e6  |
| same, method inversion        | 1.12, -3.52, -15.0         | same pattern                |
| same, mrestol 1e-12           | 1.06, -3.37, -15.3         | same pattern                |
| outlierdownweighting 3 (Huber)| 2.49, -2.02, +17.3         | 4.26e6 -> 6.34e6 -> 7.31e6  |
| no outlierdownweighting       | -4e-5, -1.5e-4, -9e-9      | 4.92e6 -> 4.67e6, monotonic |
| outlierdownweighting 4, s=1e-4| 2.5e-4, -4.0e-4, -1.7e-3   | stable                      |

Result/s^2 is the same at s=1e-2 and s=1e-4, i.e. a fixed force F of about (1e4, -4e4, -1.6e5).
At iteration 2, dfcn_exp = 1.4e6 is predicted, but fcn grows by 1.4e6. In the final chi^2 the
increase over the run without the problem is entirely the measurement term sum(p/s)^2: the data
chi^2 is unchanged.

Possible fix
------------
Accumulate the global vector and matrices (fourth loop) with exactly the weights used in the last
local-fit pass (store them per equation), or redo the local solution with the final weights. Then
the global gradient is again orthogonal to the local-parameter directions.
