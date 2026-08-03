# INFOSET-fortran

Modern Fortran/FPM implementation of the computational core of the R package
**INFOSET 4.1.1**.

INFOSET estimates a left-tail informative set from gross asset returns by
repeatedly fitting two-component lognormal mixtures. It also computes a rolling
Left Risk measure and constructs Markowitz or extreme-downside-correlation
portfolios, optionally adjusted by Left Risk.

## Implemented

- `g_ret`: sorted gross returns from prices.
- `tail_mixture`: shifted two-component lognormal-mixture fit and change point.
- `infoset_estimate`: up to two adaptive left-tail change points.
- `create_overlapping_windows`: fixed-length rolling windows.
- `lr_cp`: rolling Left Risk based on the first full-sample change point.
- `ptf_construction`: `M`, `C_M`, `EDC`, and `C_EDC` portfolios.
- `summary_ptf`: six-number and mean summary of out-of-sample returns.

The library is self-contained. The normal-mixture EM routine is adapted from
the earlier `mixtools-fortran` work, and the Goldfarb-Idnani QP solver is
adapted from `quadprog-fortran`.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example tail_mixture_example
fpm run --example left_risk_example
fpm run --example portfolio_example
```

## Direct GNU Fortran validation

Unix-like systems:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

Windows:

```text
scripts\test_gfortran.bat
```

## Interface notes

R data frames, lists, S3 methods, and strings such as `plot="T"` are replaced
by explicit arrays, derived result types, and status codes. Plotting routines
are not included. IEEE NaNs are not accepted by this package because the
upstream routines require complete positive price series.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.
