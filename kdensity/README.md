# kdensity-fortran

Modern Fortran/FPM translation of the computational parts of the R package
`kdensity` 1.2.0.

The library implements parametrically guided univariate kernel density
estimation, including asymmetric kernels for positive and unit-interval data.
Plotting and R's S3 function-object machinery are intentionally omitted.

## Implemented kernels

- Gaussian/normal
- Epanechnikov
- Rectangular/uniform
- Triangular
- Biweight
- Cosine and optimal cosine
- Triweight and tricube
- Laplace
- Gaussian-copula
- Bias-corrected and biased gamma kernels
- Bias-corrected and biased beta kernels

## Parametric starts

Built-in starts and aliases include uniform/constant, normal, lognormal,
exponential, gamma, Weibull, beta, inverse Gaussian/Wald, Gumbel, logistic,
Cauchy, Laplace, Pareto, and Lomax. `fit_kdensity_custom` accepts user-defined
start and kernel procedure pointers, replacing the dynamic R registries.

## Bandwidth selectors

- `nrd0` and `nrd`
- Hallberg Szabadvary beta-kernel rule (`HS`)
- Jones-Henderson Gaussian-copula rule (`JH`)
- Hjort-Glad Hermite reference rule (`RHE`)
- Leave-one-out unbiased cross-validation (`ucv`)
- `bcv` and `SJ` compatibility selectors

## Basic use

```fortran
program demo
  use kdensity
  implicit none

  real(dp) :: x(8) = [0.08_dp, 0.12_dp, 0.18_dp, 0.25_dp, &
                       0.32_dp, 0.45_dp, 0.60_dp, 0.72_dp]
  type(kdensity_options) :: options
  type(kdensity_fit) :: fit

  options%kernel = 'beta'
  options%start = 'beta'
  options%bandwidth = 'HS'
  options%support = [0.0_dp, 1.0_dp]
  options%support_supplied = .true.

  fit = fit_kdensity(x, options)
  if (fit%status /= kd_ok) error stop trim(fit%message)

  print *, fit%pdf(0.30_dp)
end program demo
```

## Build

With FPM:

```text
fpm test
fpm run --example kdensity_example
```

With GNU Make:

```text
make check
make optimized
make example
```

The Makefile uses GNU Fortran's heap trampoline implementation so callback
procedures do not require an executable stack.

## Layout

- `src/`: library modules
- `test/`: five deterministic test programs
- `example/`: demonstration program
- `docs/`: additional documentation
- `upstream/`: preserved upstream source and original archive

See `PORTING_NOTES.md`, `API_MAP.md`, and `VALIDATION.md` for details.
