# Algorithm notes

## Model parameterization

DiceKriging parameterizes its built-in correlations in terms of range values.
The Fortran kernels follow the package's C implementation, including the exact
internal scalings for Gaussian and Matern correlations. This matters: using a
textbook Matern formula with a different range convention gives different
maximum-likelihood estimates even though the covariance family has the same
name.

For a design matrix `F`, covariance matrix `C`, and observations `y`, the model
uses Cholesky factors and whitened least squares. When trend coefficients are
unknown, they are profiled by generalized least squares. Depending on the model
case, process variance and/or nugget are profiled or optimized consistently
with the upstream likelihood routines.

## MLE and PMLE

`km_estimate` reproduces DiceKriging's profile-likelihood cases:

1. no nugget/no observation noise: covariance ranges/shapes are optimized and
   process variance is profiled;
2. estimated homogeneous nugget: the signal fraction is optimized and total
   variance is profiled;
3. known nugget or heteroskedastic noise: signal variance is an explicit
   positive optimization parameter.

SCAD penalized likelihood is implemented for the PMLE path.

The R package delegates bounded optimization to `stats::optim` (L-BFGS-B) and
optionally `rgenoud`. This translation uses a self-contained bounded BFGS
method in logistic coordinates with randomized multistarts. The objective and
analytic gradients are the translated DiceKriging mathematics; only the
optimizer implementation differs.

## LOO

The LOO criterion and gradient use the inverse-covariance/Dubrule identities
rather than refitting `n` models. `leave_one_out` similarly uses the analytic
shortcut for the common UK-with-reestimated-trend and SK-with-fixed-trend
cases, with an exact leave-one-out reconstruction fallback for other cases.

## Prediction

`km_predict` returns simple- or universal-kriging means. Universal-kriging
variance includes the trend-estimation correction. `bias_correct=.true.`
multiplies UK conditional covariance by `n/(n-p)`, matching DiceKriging.

For standard-error output, 95% intervals follow upstream behavior:

- SK: standard normal 0.975 quantile;
- UK: Student-t 0.975 quantile with `n-p` degrees of freedom.

The exact common t critical values for 1--30 degrees of freedom are embedded;
for larger degrees of freedom a third-order asymptotic expansion is used.

## Covariance re-estimation on update

`km_update(..., cov_reestimate=.true.)` is intentionally a full refit after
appending the new observations. This is the computational behavior behind the
upstream `update.km(..., cov.reestim=TRUE)` path and is the functionality needed
for KrigInv's `CovReEstimate=TRUE` mode.

With `cov_reestimate=.false.`, the existing covariance parameters are retained
and only the factorization/trend state is rebuilt. `km_update_response` changes
already-stored responses without changing the covariance factorization.

## Scaling covariance

The nonlinear scaling model integrates a piecewise-linear field of inverse
range parameters over each input coordinate. The translated implementation
includes the value, input derivative, parameter derivative, and multidimensional
application. The optimizer samples starting `eta` values on the inverse-range
scale, matching DiceKriging's initialization logic; naive uniform sampling in
the raw very-wide eta bounds is numerically poor and does not reproduce the
upstream fit.

## Numerical linear algebra

The package is deliberately self-contained. SPD solves, Cholesky
factorization, normal-equation trend solves, covariance inverses where required
by diagnostics, and Gaussian simulation are implemented in modern Fortran.
There is no BLAS/LAPACK dependency. For large production designs, replacing the
internal dense kernels with optimized BLAS/LAPACK would be a performance
improvement without changing the public API.

## Translation boundaries

`covUser` is a holder for an arbitrary R function. Its numerical value cannot
be translated without knowing that user function, so arbitrary R callbacks are
not part of the self-contained Fortran covariance type. The built-in covariance
families used by DiceKriging and KrigInv are implemented.

`covStruct.create` also mentions `matern5_2add0`/`covAdditive0`, but DiceKriging
1.6.1 does not provide that class implementation in its source tree (the class
union even comments it out). It is therefore not reproduced as if it were an
available upstream algorithm.
