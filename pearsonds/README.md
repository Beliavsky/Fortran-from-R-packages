# PearsonDS Fortran

Modern Fortran/FPM translation of the computational core of the R package
**PearsonDS 1.3.2** (2025-03-24).

The project implements the Pearson distribution system types 0 through VII,
including density, CDF, quantile, random generation, analytical moments,
method-of-moments classification, maximum-likelihood fitting, and the model
selection criteria computed by `pearsonMSC`.

## Scope

Translated computational functionality:

- Pearson type 0: normal
- Pearson type I: shifted/scaled beta
- Pearson type II: symmetric shifted/scaled beta
- Pearson type III: shifted/reflected gamma
- Pearson type IV: native normalization, density, CDF, quantile, and RNG
- Pearson type V: shifted/reflected inverse gamma
- Pearson type VI: shifted/scaled beta-prime via the F distribution
- Pearson type VII: shifted/scaled Student t
- generic `dpearson`, `ppearson`, `qpearson`, `rpearson` dispatch
- `pearson*_moments` and generic `pearson_moments`
- `emp_moments`
- `pearson_fit_m` (method-of-moments Pearson-family identification)
- `match_moments`
- family-specific and all-family maximum-likelihood fitting
- `pearson_msc` (ML, AIC, AICc, BIC, HQC)

Intentionally skipped:

- `pearsonDiagram` and all plotting code
- R namespace/loading/registration glue
- display-oriented list/data-frame construction

The upstream package's double-double and quad-double C implementations used by
one optional hypergeometric acceleration path are not required by the public
Pearson-IV CDF when GSL is unavailable. This translation therefore follows the
upstream package's no-GSL numerical-integration path for public `ppearsoniv`.
A compact double-complex `hypergeom_2f1` series is included, but it is not a
replacement for the upstream extended-precision F21 machinery near difficult
analytic-continuation regions.

## r_mod.f90 reuse

The supplied MIT-licensed `r_mod.f90` is included in `src/` and is used for
R-compatible standard distributions, RNGs, numerical integration, and
optimization. The translation adds no duplicate implementations for helpers
already available there.

`r_mod.f90` contains C-preprocessor conditionals, so the FPM manifest enables
CPP preprocessing for `.f90` files.

## Build

```text
fpm build
fpm test
```

The project links LAPACK and BLAS because the supplied monolithic `r_mod.f90`
contains routines that reference those libraries even though PearsonDS itself
does not directly need linear algebra.

A direct gfortran build equivalent is:

```text
gfortran -cpp -std=f2018 -O2 -c src/r_mod.f90 -J build -o build/r_mod.o
gfortran -cpp -std=f2018 -O2 -I build -c src/pearsonds_mod.f90 -J build -o build/pearsonds_mod.o
gfortran -cpp -std=f2018 -O2 -I build test/test_pearsonds.f90 \
  build/pearsonds_mod.o build/r_mod.o -llapack -lblas -o build/test_pearsonds
```

## Parameter convention

`pearson_params_t%family` uses the same numeric family codes as PearsonDS:

| code | family | `par` order |
|---:|---|---|
| 0 | Pearson 0 | mean, sd |
| 1 | Pearson I | a, b, location, scale |
| 2 | Pearson II | a, location, scale |
| 3 | Pearson III | shape, location, scale |
| 4 | Pearson IV | m, nu, location, scale |
| 5 | Pearson V | shape, location, scale |
| 6 | Pearson VI | a, b, location, scale |
| 7 | Pearson VII | df, location, scale |

For bounded/one-sided MLEs, the Fortran optimizer uses smooth transformed
parameters to enforce positivity and sample-support constraints. This keeps the
same likelihood model while avoiding invalid finite-difference steps that can
occur in a literal direct-parameter translation of R's `nlminb` calls.

## Example

See `example/demo.f90`.

## Licensing

The PearsonDS-derived translation is GPL-2.0-or-later, matching upstream
`License: GPL (>= 2)`. The supplied `r_mod.f90` remains MIT-licensed. See
`LICENSE`, `NOTICE`, and `LICENSES/`.
