# apt-fortran

Modern Fortran/FPM translation of the computational algorithms in the R package **apt 4.0**, *Asymmetric Price Transmission*.

The library estimates threshold cointegration and two-equation error-correction models for a pair of price series. It is array based: callers pass equally sized `real(dp)` vectors rather than R `ts` objects.

## Implemented algorithms

- TAR and momentum-TAR threshold cointegration
- Consistent TAR/MTAR with nonzero thresholds
- Trimmed least-squares threshold search
- AIC/BIC lag selection with optional common estimation windows
- Symmetric two-equation error-correction models
- Asymmetric ECMs with TAR, MTAR, or zero-threshold linear adjustment
- Optional positive/negative splitting of distributed price changes
- Equilibrium-path, Granger-causality, distributed-lag, and cumulative-asymmetry F-tests
- OLS coefficient inference and general linear restrictions
- Durbin-Watson statistics and approximate two-sided p-values
- Ljung-Box tests at arbitrary lags
- AIC, BIC, R-squared, adjusted R-squared, and equation F statistics

## Build

```sh
fpm test
fpm run
fpm run --example threshold_search
fpm run --example asymmetric_ecm
```

FPM links LAPACK and BLAS as specified in `fpm.toml`.

Without FPM:

```sh
./run_gfortran_tests.sh
```

On Windows with GNU Fortran and compatible BLAS/LAPACK libraries:

```bat
run_gfortran_tests.bat
```

## Minimal use

```fortran
use apt, only : dp, apt_mtar, ci_tar_fit_result, ci_tar_fit
real(dp), allocatable :: y(:), x(:)
type(ci_tar_fit_result) :: fit

call ci_tar_fit(y, x, fit, model=apt_mtar, lag=2, threshold=0.0_dp)
print *, fit%threshold_regression%coefficients
print *, fit%no_cointegration_test%f_statistic
```

See `API.md`, `PORTING.md`, and the programs under `example/`.

## Scope

R plotting, S3 print/summary methods, `ts` metadata, and the bundled `daVich` dataset are not converted into Fortran runtime features. The original package tree is retained under `original/` for provenance, including the dataset and documentation.

## License

The original package declares `GPL (>= 2)`. This translation is distributed under **GPL-2.0-or-later**. Complete GPL-2.0 and GPL-3.0 texts are included under `licenses/`.
