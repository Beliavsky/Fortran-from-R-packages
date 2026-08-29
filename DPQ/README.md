# DPQ-fortran

Modern Fortran/FPM translation of the computational core of the R package
**DPQ 0.6-1** (Density, Probability, Quantile computations).

DPQ is primarily a numerical-research package: many exported functions are
alternative, asymptotic, historical, or diagnostic implementations of the
same probability calculations. This port therefore exposes numerical routines
rather than reproducing R's object, plotting, or MPFR interfaces.

## Implemented numerical areas

- R-style DPQ lower/upper-tail and log-probability helpers.
- Stable log-space arithmetic, `log1mexp`, `log1pexp`, `log1pmx`, `lgamma1p`,
  `logcf`, Chebyshev evaluation, `rexpm1`, `rlog1`, and related series.
- `bd0`, Stirling error, `lgammacor`, Poisson/Gamma/binomial/negative-binomial
  density kernels, and discrete quantile replications.
- Beta/Gamma/chi-square starting approximations and exact DPQR wrappers.
- Normal tail and quantile approximations.
- Beta and noncentral-beta approximations.
- Hypergeometric approximations and recurrence utilities, including the
  upstream `pdhyper` probability-to-density ratio recurrence.
- Central and noncentral chi-square densities, CDFs and quantiles, including
  centered Poisson-mixture and Bessel forms and the package's named
  approximation families.
- Central/noncentral t helpers, Guenther/R noncentral-t series, JKB density,
  approximations, and quantile inversion.
- Wiener-germ noncentral-chi-square approximations.
- TOMS-1006-facing log-Gamma/generalized incomplete-Gamma functionality.

The umbrella module is `dpq`; component modules are available separately for
users who want a smaller namespace.

## Build

```text
fpm build
fpm test
fpm run --example dpq_demo
```

The project links BLAS/LAPACK because the supplied `r_mod.F90` helper module
contains linear-algebra procedures even though most DPQ routines do not need
them directly.

## Dependencies

The only vendored helper dependency is the user-supplied MIT-licensed
`r_mod.f90`. `src/r_mod.F90` is a free-form-line-wrapped build copy; the exact
supplied file is retained as `upstream/r_mod-original.f90`.

## Scope boundary

Plotting (`pl2curves`, `plRpois`), R formatting/type predicates, S3/Rmpfr
integration, console diagnostics, and other presentation/glue code are not
ported. Scalar Fortran procedures naturally replace several R `.1` versus
vectorized pairs.

DPQ deliberately contains multiple historical ways to calculate the same
quantity. Where an R variant mainly exists to expose a historical internal R
implementation, this port may route multiple public variants through a common
stable numerical kernel. See `PORTING_NOTES.md` and `API_MAPPING.md`; original
R/C/Fortran sources are retained under `upstream/` for algorithm archaeology.

## Validation

The included tests cover exact identities, p/q inversions, discrete-density
kernels, hypergeometric recurrence calculations, noncentral chi-square and
noncentral-t values against independent SciPy references, Wiener-germ sanity,
and the TOMS-facing generalized incomplete-Gamma interface.

The validation build used GNU Fortran with Fortran 2018, runtime checking, and
implicit external interfaces promoted to errors.
