# sharper-fortran

`sharper-fortran` is a modern Fortran 2018 translation of the computational
core of the R package SharpeR 1.4.0. It provides estimators, distributions,
inference procedures, and hypothesis tests for ordinary and optimal Sharpe
ratios without requiring R or an external numerical library.

The public API is exposed by the umbrella module `sharper` and is organized
into typed modules for distributions, estimation, inference, tests, linear
algebra, and numerical support.

## License

SharpeR is copyright 2012-2025 Steven E. Pav and is distributed under the GNU
Lesser General Public License, version 3 or later. This translation preserves
that license as `LGPL-3.0-or-later`. See `LICENSE`, `COPYING.GPL`, and `NOTICE`.
An unmodified copy of the supplied R package is retained under
`original/SharpeR-1.4.0`.

## Build

With Fortran Package Manager:

```text
fpm build
fpm test
fpm run sharper_demo
fpm run --example basic_sharpe
fpm run --example optimal_sharpe
```

A direct GNU Fortran validation build is also available:

```text
sh scripts/build_validate.sh
```

### Windows covariance-symmetry fix

The inverse-second-moment covariance returned by `ism_vcov` is explicitly
symmetrized after the delta-method matrix products. This prevents tiny
compiler-dependent differences between opposite matrix entries that caused the
original `test_unified` assertion to fail with GNU Fortran on Windows.

The sources use free-form Fortran 2018, `implicit none`, allocatable arrays,
derived result types, generic interfaces, procedure-independent numerical
routines, and `real(dp)` with `dp = kind(1.0d0)`.

## Main capabilities

- Sample Sharpe ratios for vectors and matrices of returns.
- Annualization, higher-moment bias corrections, variance estimates,
  standard errors, confidence intervals, and prediction intervals.
- Exact normal-return Sharpe-ratio density, CDF, quantile, and simulation.
- Optimal Sharpe-ratio and Hotelling T-squared distributions.
- Noncentral t, noncentral F, and noncentral chi-square calculations.
- Markowitz weights, sample optimal Sharpe ratio, SRIC, and spanning deltas.
- One-sample, paired, unpaired, equality, maximum, conditional, and optimal
  Sharpe-ratio tests.
- Power and required-sample-size calculations.
- Second-moment covariance and inverse-second-moment covariance propagation.
- F and T-squared noncentrality inference, including KRS, MLE, and unbiased
  estimators where defined.

See `COVERAGE.md` for the R-to-Fortran mapping and `PORTING_NOTES.md` for
semantic differences and documented corrections.

## Minimal example

```fortran
program basic_sharpe
   use sharper, only: dp, sr_result, test_result, fit_sr, sr_test
   implicit none
   real(dp) :: returns(8)
   type(sr_result) :: estimate
   type(test_result) :: test

   returns = [0.012_dp, -0.004_dp, 0.009_dp, 0.006_dp, &
              0.015_dp, -0.002_dp, 0.011_dp, 0.005_dp]
   estimate = fit_sr(returns, ope=252.0_dp, higher_order=.true.)
   test = sr_test(estimate, zeta=0.0_dp, alternative='greater')

   print '(a,f10.5)', 'annualized Sharpe ratio: ', estimate%value(1)
   print '(a,f10.6)', 'one-sided p-value:      ', test%p_value
end program basic_sharpe
```

## Design notes

R S3 objects are represented by explicit derived types such as `sr_result`,
`sropt_result`, `del_sropt_result`, and `test_result`. R generic dispatch is
replaced by typed procedures and generic Fortran interfaces. Time-series class
adapters, printing methods, and plotting are intentionally not part of the
compiled library.
