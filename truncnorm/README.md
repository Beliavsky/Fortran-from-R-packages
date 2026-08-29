# truncnorm-fortran

Modern free-form Fortran/FPM translation of the computational core of the R
package **truncnorm 1.0-9**.

The upstream package provides density, CDF, quantile, random generation,
expectation and variance for the univariate truncated normal distribution.
This port preserves those numerical functions and the upstream accept/reject
sampler while replacing R `.Call` registration and vector-object handling with
typed Fortran interfaces.

## Public API

```fortran
use truncnorm
```

The umbrella module exports:

- `dtruncnorm`
- `ptruncnorm`
- `qtruncnorm`
- `rtruncnorm`
- `etruncnorm`
- `vtruncnorm`
- explicit `*_recycle` array forms implementing R-style recycling
- `dp` and `set_seed_int` from the supplied `r_mod`

Scalar arguments use the same defaults as R: `a=-Inf`, `b=Inf`, `mean=0`,
`sd=1`. `dtruncnorm`, `ptruncnorm` and `qtruncnorm` also accept a vector first
argument with scalar truncation parameters. `rtruncnorm(n,...)` returns `n`
draws for scalar parameters. The `*_recycle` functions accept array arguments
and recycle each array to the maximum input length, matching the upstream R
interface.

## Algorithms retained

- the Geweke-style accept/reject truncated-normal sampler and the original
  branch thresholds from `src/rtruncnorm.c`;
- the R/NETLIB Brent `zeroin` algorithm used by `qtruncnorm`;
- the Foulley (2000) decomposition used by `vtruncnorm` for two-sided
  truncation;
- the upstream extreme-tail midpoint/uniform approximations when the interval
  is far beyond `mean +/- 6*sd`.

The supplied MIT-licensed `r_mod.f90` is used for Normal density/CDF/quantile,
RNG and seeding helpers. Only tail-stable log-Normal probability and
`log1p`/log-space subtraction helpers that are genuinely absent from `r_mod`
were added locally.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The manifest links BLAS/LAPACK because the supplied general-purpose `r_mod`
contains routines that reference them.

## Scope

R registration (`exports.c`), SEXP recycling macros and documentation/UI glue
are not translated. All exported computational functionality is represented.
See `API_MAPPING.md` and `PORTING_NOTES.md` for details.
