# Porting notes

## R object system

The R package returns lists and data frames and uses plotting side effects.
The Fortran port uses derived types and explicit status codes. Plotting is
omitted.

## Mixture fitting

The upstream code calls `mixtools::normalmixEM` with many random restarts. The
Fortran implementation uses four deterministic quantile initializations and
selects the highest-likelihood two-normal EM fit. Components are ordered by
log-scale mean. This avoids dependence on R's global RNG while retaining the
same mixture model.

The upstream component-label test evaluates a lognormal density at
`min(log(y))`, which can be nonpositive and therefore produces zero density.
The port uses the mathematically intended ordering by component log mean.

The change point is obtained by solving the equality of the two weighted
lognormal densities. A crossing between the component log means is preferred.

## `infoset`

The original loop permits at most two change points and rejects change points
at or above the sample median. The port preserves that behavior without
raising an R exception when no split exists.

## `LR_cp`

The R source overwrites user-supplied `FT` and `ov` with 1290 and 125. The
Fortran routine honors the supplied arguments. The histogram calculation uses
equal-width Sturges bins; R's default `hist` uses aesthetically rounded breaks,
so small numerical differences are expected.

If a full-sample split cannot be identified, the port uses the 10th percentile
of gross returns as a fallback and returns status `infoset_no_split`.

## Portfolio construction

The original `quadprog::solve.QP` constraints and linear terms are preserved.
`Matrix::nearPD` is replaced by symmetric eigendecomposition with eigenvalue
flooring. The EDC construction retains the original 5% threshold and centering
rule.

The R code rounds weights to six decimals before out-of-sample evaluation.
The Fortran implementation retains full solver precision so equality
constraints remain accurate.

The upstream out-of-sample index range contains `ov + 1` observations, despite
`summary_ptf` later reading only 125 fixed observations. The Fortran result
stores all `overlap + 1` values and summarizes the entire supplied matrix.
