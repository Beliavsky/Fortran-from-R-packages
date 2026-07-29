# Porting notes

## Scope

The port covers every exported computational routine in FRAPO 0.4-2:

`capser`, `cr`, `dr`, `mrc`, `PAveDD`, `PCDaR`, `PERC`, `PGMV`, `PMaxDD`,
`PMD`, `PMinCDaR`, `PMTD`, `returnconvert`, `returnseries`, `rhow`, `sqrm`,
`tdc`, `trdbilson`, `trdbinary`, `trdes`, `trdhp`, `trdsma`, and `trdwma`.

The S4 classes are represented by the typed Fortran `portfolio_result` object.
Simple R accessors such as `Weights` and `DrawDowns` correspond directly to its
components.

## Deliberately omitted R infrastructure

The following are not numerical algorithms and were not reimplemented:

- S4 registration and method dispatch
- plotting and formatted `show` methods
- `timeSeries`, `ts`, `xts`, and `zoo` metadata propagation
- example-file discovery, copying, editing, and execution
- bundled dataset loading as R objects

The original data and book examples remain in `original/FRAPO-master/`.

## Optimization substitutions

FRAPO calls `cccp` for long-only convex quadratic programs and risk parity, and
`Rglpk` for drawdown linear programs. The Fortran port instead provides:

- a dense infeasible-start, predictor-corrector primal-dual interior-point
  solver for convex quadratic programs with equality and inequality constraints;
- cyclic coordinate descent for equal-risk-contribution portfolios;
- the same interior-point solver with a negligible diagonal regularization for
  linear drawdown programs.

These substitutions remove all non-BLAS/LAPACK dependencies. In problems with
multiple optimal LP solutions, auxiliary slack variables can differ from GLPK's
chosen basic solution while weights and objective values remain optimal.

## Drawdown reporting correction

The R functions report drawdowns from LP high-water-mark variables. Those
variables are constrained correctly but are not always uniquely minimized, so
a solver may return inflated slack values that still satisfy all constraints.
The Fortran port reports the economically meaningful path instead:

1. calculate the optimized portfolio equity path;
2. calculate its running high-water mark, starting at zero;
3. report high-water mark minus equity.

For CDaR portfolios, the displayed threshold is the R type-7 empirical
`alpha` quantile of that realized drawdown path, and `risk_value` is the mean of
observations at or above the threshold. The optimization formulation itself is
unchanged.

## Numerical details

- Sample covariance divides by `n-1`, matching R's default `cov` behavior.
- Tail dependence uses average ranks for ties.
- The untrimmed first return is represented by IEEE NaN.
- `sqrm` is restricted to real symmetric matrices. Negative eigenvalues are
  rejected instead of producing complex results.
- Weighted moving averages preserve the coefficient orientation used by
  `stats::filter(..., sides=1)`.
- HP filtering is solved as `(I + lambda D'D) trend = data` with LAPACK.

## License

The attached FRAPO package is GPL version 3 or later. The translated source has
`SPDX-License-Identifier: GPL-3.0-or-later` headers and includes the complete
GPL version 3 text.
