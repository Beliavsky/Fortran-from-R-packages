# Porting map

## Directly translated

| FatTailsR area | Fortran implementation |
|---|---|
| `kashp`, `ashp`, `dkashp_dx` | `fattailsr_math` |
| `logit`, `invlogit` | `fattailsr_math` |
| standardized logistic `d/p/q/r`, derivatives, tail means, ES | `fattailsr_distributions` |
| K1/K2/K3/K4/K7 `d/p/q/r`, derivatives and logit transforms | `fattailsr_distributions` |
| Kiener VaR, LTM, RTM, ES, `dtmq`, `c`, `h` | `fattailsr_distributions` |
| parameter conversion functions | `fattailsr_params` |
| `kmoment`, `kcmoment`, moment summaries | `fattailsr_moments` |
| `fiveprobs`, `sevenprobs`, `elevenprobs` | `fattailsr_estimation` |
| `estimkiener5`, `estimkiener7`, `estimkiener11` | `fattailsr_estimation` |
| K4 quantile regression | `fit_kiener_k4` |
| `elevate`, `replaceNA`, `fatreturns`, `logreturns` numerical core | `fattailsr_returns` |
| Laplace-Gauss normal numerical primitives | normal functions in `fattailsr_math` |

## Representation differences

The seven values `(m,g,a,k,w,d,e)` are represented by the derived type
`kiener_parameters`. Constructors `make_k1`, `make_k2`, `make_k3`, and
`make_k4` populate all equivalent fields consistently. This replaces R's
vectors, matrices, arrays, lists, names, and class attributes.

Scalar procedures are elemental where useful, so they can operate directly on
Fortran arrays. Operations that require inversion or fitting are ordinary pure
functions or subroutines.

## Not ported

- Plotting and graphics.
- R object classes and printing methods.
- Dynamic dispatch across lists, data frames, and third-party time-series
  classes.
- Bundled example datasets and extraction/download helpers.
- Parallel-cluster wrappers.
- Formatting helpers whose only purpose is R display names or rounded tables.

These omissions do not remove the underlying distribution, risk, moment, or
estimation mathematics.
