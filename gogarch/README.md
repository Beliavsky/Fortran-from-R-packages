# gogarch-modern-fortran

A modern Fortran translation of the computational core of the R package
`gogarch` 0.7-6 by Bernhard Pfaff.

This release targets numerical GO-GARCH functionality while excluding
plotting, S4 classes, R method dispatch, time-series labels, and other
R-specific infrastructure.

## Tested computational coverage

The following features are implemented and exercised by the test suite:

- covariance initialization and whitening corresponding to `goinit`
- two-dimensional rotations (`Rd2`)
- products of Euler/Givens rotations (`UprodR`)
- orthogonal-column matching (`Umatch`)
- symmetric `vech` and `unvech` conversion
- matrix-process autocovariance/autocorrelation calculation (`cora`)
- univariate GARCH(p,q), including orders greater than one
- APARCH(p,o,q) powered-scale recursions with leverage terms
- conditional filtering, log likelihood, constrained fitting, simulation,
  and multi-step forecasting
- standardized conditional distributions:
  - Normal (`norm`)
  - Fernandez-Steel skew-Normal (`snorm`)
  - Student-t (`std`)
  - Fernandez-Steel skew-Student (`sstd`)
  - generalized error distribution (`ged`)
  - Fernandez-Steel skew-GED (`sged`)
- shape estimation for Student-t and GED families
- skew estimation for skewed distribution families
- APARCH power-delta estimation
- symmetric FastICA estimation
- GO-GARCH estimation by FastICA, methods of moments, nonlinear least
  squares, and Euler-angle maximum likelihood
- direct construction and likelihood evaluation from Euler angles
- conditional covariance, variance, and correlation paths
- covariance-standardized residuals with fitted factor means removed
- legacy first-lag coefficient extraction and full higher-order parameter
  extraction
- multi-step asset mean and covariance forecasts
- simulation from fitted higher-order and APARCH GO-GARCH models
- a CSV-fitting command-line program supporting model orders and conditional
  distributions

## Matrix convention

For observations stored by rows, the implementation uses

```text
X = F A^T
H_t = A diag(h_t) A^T
```

where `F` contains independent factor returns, `A` is the asset-by-factor
mixing matrix, and `h_t` contains factor conditional variances. Tests verify
both exact data reconstruction under this convention and covariance
construction.

## Univariate model definitions

GARCH(p,q) uses

```text
h_t = omega
    + sum_i alpha_i epsilon_(t-i)^2
    + sum_j beta_j h_(t-j)
```

APARCH(p,o,q) uses

```text
sigma_t^delta = omega
              + sum_i alpha_i (abs(epsilon_(t-i))
                - gamma_i epsilon_(t-i))^delta
              + sum_j beta_j sigma_(t-j)^delta
```

with `gamma_i = 0` for ARCH lags beyond the requested leverage order `o`.
Conditional innovations are standardized to mean zero and variance one.

The optimizer enforces positive `omega`, nonnegative ARCH/GARCH coefficients,
`abs(gamma) < 1`, bounded positive delta, valid shape/skew parameters, and a
coefficient-sum stability bound below 0.995. APARCH simulation additionally
computes the innovation-power persistence and rejects a supplied parameter set
when that persistence is not below one. This is not claimed to reproduce every
constraint or initialization choice made internally by `fGarch`.

## Distribution parameterization

`std` is a variance-one Student-t distribution with shape greater than two.
`ged` is a variance-one generalized error distribution. The skewed families
use the Fernandez-Steel two-piece transformation and are subsequently centered
and scaled to mean zero and variance one.

The tests numerically check density normalization and random-number mean and
variance for all six distributions. Likelihood fitting is exercised for all
six. Separate tests fit shape, skew, APARCH delta, and leverage parameters.

## Public specification type

```fortran
type(univariate_spec) :: spec

spec%model = "aparch"       ! "garch" or "aparch"
spec%distribution = "sstd" ! norm, snorm, std, sstd, ged, sged
spec%p = 1
spec%o = 1
spec%q = 1
spec%delta = 1.5_dp
spec%shape = 8.0_dp
spec%skew = 1.1_dp
spec%fit_delta = .true.
spec%fit_shape = .true.
spec%fit_skew = .true.

fit = fit_gogarch_ica(data, factor_spec=spec)
```

The same `factor_spec` optional argument is accepted by `fit_gogarch`,
`fit_gogarch_ica`, `fit_gogarch_mm`, `fit_gogarch_nls`,
`fit_gogarch_ml`, `gogarch_from_angles`, and `gogarch_negloglik`.

## Important limitations

This is not a complete replacement for the R package or for `fGarch`.
The following are not included or not claimed:

- R formulas and arbitrary formula parsing
- ARMA conditional-mean models
- distributions beyond `norm`, `snorm`, `std`, `sstd`, `ged`, and `sged`
- exact `fGarch` recursion initialization or optimizer equivalence
- analytical Hessians, parameter covariance matrices, and inference
- exact numerical equivalence to R's `fastICA`, `optim`, or `nlminb`
- separate model specifications for different factors; one specification is
  applied to all independent factors in a GO-GARCH fit
- S4 classes and methods
- `ts` attributes, row names, formulas, calls, update methods, and validation
  methods tied to R objects
- plotting, printing, and summary presentation methods
- bundled `.rda` datasets and PDF talks

## Requirements

- a Fortran 2018 compiler
- BLAS
- LAPACK

Validation used GNU Fortran 14.2.0 with the system BLAS and LAPACK libraries.

## Build and test

On a Unix-like system:

```sh
make check
```

or:

```sh
./run_checks.sh
```

This performs source-license checks, a runtime-checked warnings-as-errors
build, an optimized warnings-as-errors build, four numerical test suites, the
demo, and CSV runs for default GARCH, higher-order Student GARCH, and
skew-Student APARCH.

A Windows batch file, `run_checks.bat`, is included for GNU Fortran setups with
BLAS and LAPACK available to the linker. It was not executed in the validation
environment.

An `fpm.toml` file is included. `fpm` was not installed in the validation
environment, so an `fpm` build is not claimed as tested.

## CSV example

```sh
build/debug/fit_csv data/sample_returns.csv ica
build/debug/fit_csv data/sample_returns.csv ica garch std 2 0 1
build/debug/fit_csv data/sample_returns.csv ica aparch sstd 1 1 1
```

Arguments are:

```text
FILE METHOD MODEL DISTRIBUTION P O Q
```

The first CSV column is treated as a date or label and the remaining columns
as asset returns.

## Main public API

```fortran
use gogarch
```

Primary univariate procedures include `filter_garchpq`, `fit_garchpq`,
`simulate_garchpq`, `filter_aparch`, `fit_univariate`, `simulate_aparch`, and
`forecast_univariate`.

Primary multivariate fitting functions include `fit_gogarch`,
`fit_gogarch_ica`, `fit_gogarch_mm`, `fit_gogarch_nls`, and `fit_gogarch_ml`.
Use `factor_coefficients_full` to retrieve all lag vectors and distribution
parameters.

See `API_MAP.md` for the mapping from R functions and methods to Fortran
procedures.

## License

The original package declares `GPL (>= 2)`. This translation preserves that as
`GPL-2.0-or-later` in `fpm.toml` and in every Fortran source file. `LICENSE`
contains the GNU General Public License version 2 text.
