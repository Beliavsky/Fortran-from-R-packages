# Porting notes

## Scope

The computational content of the four exported R functions was translated:

| R entry point | Fortran entry point |
| --- | --- |
| `meanEstimation` | `mean_estimation` |
| `covEstimation` | `covariance_estimation` |
| `semidevEstimation` | `semideviation_estimation` |
| `optimalPortfolio` | `optimal_portfolio` |

The package contains 4 mean methods, 10 executable covariance methods, 2
semideviation methods, 7 portfolio methods, and 4 constraint modes. All are
represented in the Fortran API.

## R-specific features not translated

- R list parsing and partial matching
- R error and warning object behavior
- The bundled `Industry_10.rda` dataset as an executable Fortran data object
- R documentation objects and package namespace machinery
- `quadprog`, `nloptr`, `MASS::factanal`, and other R dependency interfaces

Original R source and metadata are retained under `original/`.

## Numerical replacements

The R package uses `quadprog::solve.QP` for bounded quadratic programs and
`nloptr::slsqp` for nonlinear constrained problems. The Fortran port uses:

- Exact linear-system solutions for unconstrained minimum variance,
  mean-variance, and maximum decorrelation cases
- Projected gradient optimization with Barzilai-Borwein step estimates and
  backtracking for bounded and nonlinear cases
- Projection onto budget plus box constraints by scalar bisection
- Dykstra projection for the intersection with an L1 gross-exposure ball

The resulting objectives and constraints match the R routines, but iteration
paths and final floating-point values need not be bit-for-bit identical.

## Factor covariance

The original calls `MASS::factanal`, which performs maximum-likelihood factor
analysis through R's optimization stack. The port uses a principal-component
factor reconstruction of the correlation matrix, retaining the largest K
eigencomponents and diagonal uniquenesses. This preserves the K-factor
covariance model and positive diagonal residual variance without adding a
second nonlinear optimizer. It is the one estimator that is intentionally not
a bit-for-bit translation of the R dependency result.

## Corrected source issues

The following source-level inconsistencies were handled explicitly:

1. Bayes-Stein mean sample size

   The R code sets `T <- length(mu)`, which equals the number of assets rather
   than the number of observations. The Fortran routine uses `size(rets,1)`,
   consistent with the Bayes-Stein formula and the package's covariance code.

2. Maximum decorrelation without constraints

   The R unconstrained branch solves with `Sigma`, while constrained branches
   optimize the correlation matrix `Rho`. The Fortran implementation uses the
   correlation matrix in every branch, consistent with the method name and
   constrained implementation.

3. Large-dimensional shrinkage bounds

   The R `large` estimator does not clamp its computed shrinkage coefficient.
   The Fortran port clamps it to [0,1], as the other Ledoit-Wolf routines do,
   preventing unintended extrapolation outside the sample-target segment.

4. Risk-efficient decile construction

   The R code stores ten decile groups in an N by N matrix and can create empty
   groups for small or tied cross sections. The Fortran port computes the same
   decile-median score directly and leaves an asset's semideviation as a safe
   fallback when a bin is empty.

5. Internal `rtm` label

   `.ctrCov` lists `rtm`, but `covEstimation` has no dispatch branch or
   implementation for it. No `rtm` method is exposed by this port.

## Preserved source behavior

The unconstrained mean-variance R branch normalizes `Sigma^-1 mu`; its later
normalization makes the supplied risk-aversion value cancel. The Fortran port
preserves this executable behavior. Risk aversion affects constrained
mean-variance optimization, matching the R objective.

The EWMA powers and observation order follow the R source exactly: rows are
oldest to newest, and weights use powers T through 1.
