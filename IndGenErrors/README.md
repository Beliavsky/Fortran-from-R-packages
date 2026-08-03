# IndGenErrors-fortran

Modern Fortran translation of the computational code in R package
`IndGenErrors` 0.1.6, packaged for the Fortran Package Manager (FPM).

The library tests independence between innovations of two or three time
series using cross-correlations, rank-based cross-dependence measures, and
Cramer-von Mises Mobius statistics.

## Implemented routines

- `crosscor_2series` and `crosscor_3series`
- `crossdep_2series` and `crossdep_3series`
- `cvm_2series` and `cvm_3series`

The dependence routines include Spearman, van der Waerden, and Savage
Mobius scores adapted from the earlier `MixedIndTests` Fortran translation.
The finite-sample Cramer-von Mises calibration tables for sample sizes 50
and 100 are translated exactly from the original C sources.

## Build

```text
fpm build
fpm test
fpm run
```

With GNU Fortran but without FPM:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## API conventions

Lag arrays are returned in natural order from `-lag` through `+lag`.
Three-series lag pairs are ordered lexicographically from `(-lag,-lag)`
through `(lag,lag)`. These are the same final orders produced by the R
wrappers after reordering the native C output.

All input vectors must have equal lengths. Circular lagging follows the
original package's native implementation.

## Scope

`dependogram` and `CrossCorrelogram` are plotting-only R functions and are
not translated. The bundled R datasets are retained under `original/data`
but are not compiled into the Fortran library.

## License

The original `IndGenErrors` package is GPL-2.0-or-later. This self-contained
translation includes adapted code from `MixedIndTests`, which is
GPL-3.0-only, so the combined Fortran project is distributed under
GPL-3.0-only. See `NOTICE.md` and `LICENSE`.
