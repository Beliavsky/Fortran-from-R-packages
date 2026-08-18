# copula-fortran

A dependency-free modern Fortran/FPM translation of a substantial numerical
core from the R package `copula` 1.1-7.

The upstream package is a large R/S4 framework with about 45,000 lines of R and
C source. This project translates the most reusable array-based numerical
algorithms rather than claiming a mechanical translation of every S4 method,
plotting function, or specialized test.

## Features

- Independence, Gaussian, Student-t, Clayton, Gumbel, Frank, AMH, Joe, FGM,
  Plackett, Marshall-Olkin, lower/upper Frechet-Hoeffding, Galambos,
  Husler-Reiss, and asymmetric logistic/Tawn copulas.
- Bivariate rotations by 90, 180, and 270 degrees.
- Copula CDFs, densities, log densities, conditional distributions, and seeded
  simulation.
- Multivariate Gaussian and Student-t copulas with typed correlation matrices.
- Kendall's tau, Spearman's rho, tail-dependence coefficients, and parameter
  inversion from tau or rho.
- Pseudo-observations, empirical copulas, sample rank correlations, and
  permutation tests for independence, exchangeability, and radial symmetry.
- Rosenblatt and inverse Rosenblatt transforms for bivariate models.
- Inversion-of-tau and maximum-pseudo-likelihood fitting for supported
  one-parameter families, including numerical standard errors.
- Mixture copulas, Khoudraji asymmetrization, and a two-level nested-Clayton CDF.
- Pickands dependence functions for implemented extreme-value families.
- Stirling numbers, Eulerian numbers, Sibuya probabilities/simulation, and
  logarithmic-series probabilities/simulation.
- Native normal, Student-t, beta, gamma, matrix, random-number, and quadrature
  support. No BLAS, LAPACK, GSL, R, or external statistics library is required.

## Quick start

```fortran
program example
  use copula
  implicit none

  type(copula_model) :: model
  real(dp), allocatable :: sample(:,:)
  logical :: ok

  model = clayton_copula(2.0_dp)

  print *, pCopula([0.3_dp,0.7_dp],model)
  print *, dCopula([0.3_dp,0.7_dp],model)
  print *, tau(model)

  call rCopula(1000,model,sample,ok,12345_i8)
  if (.not. ok) error stop 'simulation failed'
end program example
```

## Build

```text
fpm build
fpm test
fpm run copula_demo
fpm run --example common_families
fpm run --example fit_and_rosenblatt
```

Direct GNU Fortran validation scripts are also included:

```text
scripts/validate.sh
scripts/validate_optimized.sh
```

On Windows with GNU Fortran:

```bat
scripts\validate.bat
```

## Numerical conventions

- Ordinary bivariate family densities that do not have a dedicated closed-form
  implementation are evaluated by stable central differences of the CDF.
- Bivariate Gaussian and Student probabilities use deterministic numerical
  integration. Higher-dimensional Gaussian probabilities use conditional
  shifted-Halton integration; higher-dimensional Student probabilities combine
  chi-square mixing with that Gaussian engine.
- Maximum-pseudo-likelihood fitting is currently one-dimensional. Student-t
  fitting estimates the correlation with degrees of freedom supplied by the
  caller.
- Singular copulas such as Frechet bounds and Marshall-Olkin contain singular
  mass. `dCopula` returns the ordinary absolutely continuous density only and
  therefore returns zero for the purely singular implementations.
- Seeds use an internal xorshift generator, providing reproducibility independent
  of the Fortran processor's intrinsic random-number generator.

## Scope

See `COVERAGE.md` for the exact translated and excluded functionality and
`PORTING_NOTES.md` for numerical and API differences from R.

## License

GNU GPL version 3 or later (`GPL-3.0-or-later`), preserving the upstream
license and attribution. See `LICENSE`, `NOTICE`, and the retained original
package.
