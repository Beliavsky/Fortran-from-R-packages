# Porting notes

## R to Fortran API mapping

| R export | Fortran procedure |
|---|---|
| `UnifiedEst` | `unified_est` |
| `RealizedEst` | `realized_est` |
| `RealizedEst_Option` | `realized_est_option` |

R lists are represented by `garchito_result`. Optional R arguments become
Fortran optional arguments, and an explicit `garchito_control` replaces
implicit optimizer defaults.

## Option function argument correction

The upstream R signature is

```text
RealizedEst_Option(RV, JV=NULL, NV=NULL, homogeneous=TRUE)
```

but its vignette also shows `RealizedEst_Option(RV, NV)`, which passes the
second positional argument to `JV` and therefore cannot work as written. The
Fortran interface removes this ambiguity by taking required `rv, nv` first and
making `jv` optional.

## Optimization

The original R code delegates optimization to `Rsolnp::solnp`. The previously
translated Rsolnp project is GPL-2.0-only, while this package is GPL-3.0-only.
Copying or statically combining those source files would make the distributed
work license-incompatible.

This translation therefore uses an independent optimizer implementation:

- projected bounded Nelder-Mead simplex iterations;
- exact projection onto the model's nonnegative stationarity triangle;
- three deterministic starts for the dynamic coefficients;
- bounded coordinate polishing after simplex termination;
- finite-objective guards and explicit convergence status.

It preserves the R bounds, starting values, stationarity restriction, and
likelihood functions, but optimizer iteration paths and last digits need not
match Rsolnp.

## Numerical safeguards

The R bounds permit exactly zero error standard deviation and exactly unit
persistence, although either value makes the likelihood singular. The
Fortran implementation uses a tiny positive lower bound for scale parameters
and projects persistence below one by `1e-8`.

Unlike the R source, all three estimators consistently reject negative realized
volatility and jump variation before optimization.

## Data and plotting

The upstream `sample_data.rda` file is retained under `original/data/`. It is an
R serialization and is not parsed by the Fortran library. Examples use
self-contained synthetic series. Vignette plotting is presentation code and is
not translated.
