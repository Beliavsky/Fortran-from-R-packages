# imputeFin-fortran

Modern Fortran translation of the computational core of the R package
`imputeFin` 0.1.2.9000. The project is an FPM library and has no external
runtime dependencies.

## Implemented API

- `fit_ar1_gaussian`
- `impute_ar1_gaussian`
- `impute_rolling_ar1_gaussian`
- `fit_ar1_t`
- `impute_ar1_t`
- `fit_var_t`
- `impute_ohlc`
- `impute_vol`

The AR routines accept IEEE NaNs as missing values. Vector and matrix overloads
are supplied where the R package accepted one or several series. Leading and
trailing NaNs are retained; only missing values between observed endpoints are
imputed.

## Build

```text
fpm build
fpm test
fpm run --example gaussian_ar1_example
fpm run demo_imputefin
```

GNU Fortran users can also run `tools/test_gfortran.sh` or
`tools/test_gfortran.bat`.

## Numerical scope

The Gaussian AR(1) estimator uses the original conditional-moment EM formulas.
Student-t AR fitting includes the fast consecutive-pair IRLS method and a
portable Gibbs/SAEM path. Student-t imputation uses latent-scale Gibbs updates.
The VAR estimator uses multivariate Student-t IRLS and sequential conditional
mean filling for missing contemporaneous observations.

See `PORTING.md` for differences from the R implementation.

## License

GPL-3.0-only. The original package sources and metadata are retained under
`original/` for license attribution and comparison.
