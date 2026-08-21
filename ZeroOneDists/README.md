# ZeroOneDists-fortran

Modern Fortran 2018 translation of the computational code in the R package
`ZeroOneDists` 1.0.0.

The port preserves the upstream MIT license and intentionally excludes R class,
formula, and plotting/presentation machinery. The numerical family callbacks
used by GAMLSS are retained as ordinary Fortran procedures, and a generic
array/design-matrix maximum-likelihood regression layer is provided so the
families can be fitted without an R runtime.

## Implemented distributions

- Beta Rectangular (`BER`): `dber`, `pber`, `qber`, `rber`
- Mean-parameterized Beta Rectangular (`BER2`): `dber2`, `pber2`, `qber2`, `rber2`
- Unit Half Logistic-Geometric (`UHLG`): `duhlg`, `puhlg`, `quhlg`, `ruhlg`
- Unit Maxwell-Boltzmann (`UMB`): `dumb`, `pumb`, `qumb`, `rumb`
- Unit-Power Half-Normal (`UPHN`): `duphn`, `puphn`, `quphn`, `ruphn`

The d/p/q functions are scalar-elemental where mathematically practical, so
array arguments work naturally. RNG routines accept either scalar parameters
or same-length parameter arrays.

## GAMLSS-family computational layer

`zero_one_families` provides:

- family IDs and parameter counts
- default link functions and inverse links
- link derivatives (`mu_eta`)
- log density, CDF and quantile dispatch
- first-derivative scores
- the upstream working second-derivative convention
- deviance increments
- upstream initial parameter values
- parameter/response validity checks
- BER and BER2 mean functions

`zero_one_fit` provides `fit_zero_one`, a generic design-matrix regression
fitter with separate model matrices for `mu`, `sigma`, and `nu`. It uses BFGS
on linked coefficients and returns fitted parameters, coefficient covariance,
log likelihood, AIC, convergence status, and coefficient blocks.

## Build

With FPM:

```sh
fpm test
fpm run --example basic
```

A compiler-only validation script is also included:

```sh
./run_tests.sh
```

The validation build used:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## Dependencies

The active numerical library is self-contained and has no external Fortran
package dependency. The supplied `gamlss-fortran-v0.9.0` and
`gamlss.dist-fortran-v0.3.0` translations were inspected as compatibility
references for link and GAMLSS semantics. They are not linked into the default
build; both supplied FPM manifests declare `GPL-3.0-only`, whereas
ZeroOneDists is MIT licensed.

See `docs/API_MAP.md`, `PORTING_NOTES.md`, and `LICENSES.md` for details.
