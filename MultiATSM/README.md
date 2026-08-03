# MultiATSM Fortran

A modern Fortran 2018 translation of the reusable computational core of the R
package **MultiATSM 1.5.1-2**. The library is organized as an FPM project and
uses variables-by-time matrices throughout.

The port covers the numerical workflow used by multi-country Gaussian affine
term-structure models:

- principal-component/spanned yield factors;
- unrestricted and restricted VAR(1), VARX, GVAR, and JLL dynamics;
- trade-weight transition matrices and foreign/star factors;
- affine yield-loading recursions and latent-to-observed rotations;
- Gaussian state and measurement likelihoods;
- BFGS and Nelder-Mead optimization with numerical derivatives;
- stationarity and positive-semidefinite covariance transforms;
- forecasts, fitted yields, IRFs/GIRFs, FEVD/GFEVD, forward rates, expected
  short-rate components, and term premia;
- residual, wild, and moving-block bootstrap utilities;
- stochastic-approximation VAR bias correction and eigenvalue shrinkage.

The R S3/list interface, Excel/data-frame/date processing, packaged datasets,
folder management, plotting, and report-generation code are not translated.
See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for the detailed mapping and
known differences.

## Requirements

- A Fortran 2018 compiler
- BLAS and LAPACK
- FPM, or GNU Fortran for the supplied scripts

Tested with GNU Fortran 14.2.0 and the system BLAS/LAPACK libraries.

## Build with FPM

```text
fpm build
fpm test
fpm run demo_multiatsm
fpm run --example affine_pricing_example
```

## Build without FPM

On Unix-like systems:

```text
./tools/test_all.sh checked
./tools/test_all.sh optimized
```

On Windows with `gfortran`, run:

```text
tools\test_all.bat checked
tools\test_all.bat optimized
```

## Main modules

| Module | Purpose |
|---|---|
| `multiatsm` | Convenience umbrella module |
| `multiatsm_pca` | PCA weights and spanned factors |
| `multiatsm_var` | VAR, restricted OLS, VARX, GVAR, transition weights |
| `multiatsm_jll` | JLL orthogonalization and restricted dynamics |
| `multiatsm_affine` | Affine A/B loading recursions and rotations |
| `multiatsm_likelihood` | Gaussian ATSM likelihood and measurement variance |
| `multiatsm_outputs` | Forecasts, responses, decompositions, premia |
| `multiatsm_optimization` | Optimizers, derivatives, parameter transforms |
| `multiatsm_bootstrap` | Bootstrap resampling and confidence bounds |
| `multiatsm_bias` | VAR bias correction and stability shrinkage |

## Array conventions

- Time series are generally `n_variables x n_observations`.
- Yield matrices are `n_yields x n_observations`.
- Country panels used by GVAR routines are
  `n_countries x n_country_variables x n_observations`.
- Impulse responses are `response x shock x horizon`.
- Maturity inputs are integer numbers of model periods.

## License

The upstream package declares `GPL-2 | GPL-3`. This translation is therefore
licensed under **GPL-2.0-only OR GPL-3.0-only**. Complete license texts and the
unmodified upstream source snapshot are included.
