# MM-fortran

Modern Fortran translation of the computational code in the R package
`MM` 1.7-0 by Robin K. S. Hankin and P. M. E. Altham.

The upstream package implements the multiplicative multinomial distribution
and related multiplicative-binomial utilities.  This port keeps the numerical
algorithms and replaces R's S4/S3/GLM/optimizer machinery with typed Fortran
APIs and native numerical implementations.

## Implemented computational functionality

- `paras` parameter construction and `p` / `theta` extraction and replacement.
- `lmultinomial` and `multinomial`.
- `mm_single_log`, `mm_single`, `normc_log`, `normc`, and `dmm`.
- Equal-row-sum and differing-row-sum log likelihoods.
- Sufficient statistics, support-form likelihood, and expected sufficient
  statistics.
- Lindsey-Mersch Poisson log-linear fitting for the multiplicative
  multinomial.
- Bivariate `Lindsey_MB` Poisson fit for multiplicative-binomial data.
- Maximum-likelihood refinement in the upstream log-parameterization.
  The default smooth optimizer is a finite-difference BFGS analogue of R's
  `nlm`; `method="Nelder"` selects Nelder-Mead.
- `gunter` support-table aggregation for multinomial and MB observations.
- Metropolis-Hastings simulation via `rmm`.

Plotting, printing methods, R formula parsing, S3/S4 dispatch, and `Oarray`
container behavior are not translated.  Their computational equivalents are
represented with Fortran derived types where useful.

## Dependencies

The supplied Fortran translations are vendored under `src/deps/` and used
by the implementation:

- `partitions-fortran` supplies exact composition enumeration used by the
  multiplicative-multinomial normalizing constant and support tables.
- `quadform-fortran` supplies the quadratic-form operations used for the
  interaction term and its row-wise evaluation.

The exact dependency archives supplied with this translation are retained in
`provenance/dependencies/`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example demo_mm
```

No external C, R, BLAS, or LAPACK library is required.

## Main module

```fortran
use multiplicative_multinomial
```

A small example:

```fortran
use multiplicative_multinomial
implicit none

type(paras_type) :: par
real(dp) :: th(2,2)
integer :: y(2)

par = paras(2)
th = 1.0_dp
th(1,2) = 2.0_dp
call set_theta(par, th)
y = [1, 1]

print *, dmm(y, par)   ! 2/3
```

## Numerical notes

`NormC` in the R package performs direct summation over all weak
compositions.  The Fortran port enumerates the same exact support through the
translated `partitions` package but uses log-sum-exp accumulation to reduce
underflow/overflow risk.

The number of support points is combinatorial in total count and dimension,
so exact normalization can still become expensive for large problems.  This
is an inherent property of the upstream algorithm.

`MM_allsamesum_A()` is retained as `mm_allsamesum_a()` with the upstream
formula exactly, including its use of `NormC` rather than `log(NormC)` in the
first term.  The regular `mm_allsamesum()` routine is the recommended
likelihood implementation.

## License and provenance

The upstream `MM` DESCRIPTION declares `License: GPL-2`.  This translation is
provided under GPL-2 and retains the exact uploaded upstream archive, selected
metadata/source files, and the exact supplied dependency archives under
`provenance/`.
