# mixtools-fortran

Modern Fortran translation of the computational core of the R package
`mixtools` 2.0.0, packaged for the Fortran Package Manager (FPM).

The library provides typed, array-based interfaces for finite mixture models,
regression mixtures, semiparametric mixtures, censored reliability mixtures,
simulation, bootstrap/model selection, and supporting numerical utilities.
Plotting, R formula handling, S3 objects, Plotly/Shiny integration, and bundled
R data objects are not part of the compiled library.

## Main capabilities

- Univariate and multivariate Gaussian mixture EM
- Gamma, multinomial, and repeated-measures normal mixtures
- Equality and linear constraints for normal mixture parameters
- Linear, logistic, Poisson, segmented, mixed, local-proportion, and
  mixture-of-experts regression mixtures
- Product-kernel `npEM`, `spEM`, `mvnpEM`, symmetric-location mixtures, and
  semiparametric regression mixtures
- Exponential and Weibull censored reliability mixtures
- Normal-mixture bootstrap comparison and standard errors
- AIC, BIC, CAIC, and ICL model selection
- Bayesian regression-mixture sampling and coefficient intervals
- Weighted KDE/quantiles, Dirichlet and multivariate-normal densities,
  matrix square roots, ellipse construction, component CDFs, FDR curves,
  depth, and random mixture generation

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The manifest uses a plain semantic version (`2.0.0`) accepted by FPM.

## R and Fortran multivariate-mixture comparison

From the `mixtools` directory, generate a reproducible mixture of two
bivariate normal distributions, fit it with the R package, and fit the same
observations with the Fortran translation:

```text
Rscript example\generate_mvnormal_mixture.R
Rscript example\fit_mvnormal_mixture.R
fpm run --example fit_mvnormal_mixture
```

On Windows, the same workflow can be run with:

```text
run_mvnormal_comparison.bat
```

On Linux, macOS, WSL, or Git Bash, use:

```text
bash run_mvnormal_comparison.sh
```

From Windows `cmd` with Git Bash installed, use:

```text
bash mixtools/run_mvnormal_comparison.sh
```

Both runners create the observation file only if it does not already exist.
Their `fpm run` command compiles or rebuilds the Fortran program when needed.

The generator writes `example/mvnormal_mixture_data.txt`. The fitting programs
write `example/mvnormal_mixture_fit_r.txt` and
`example/mvnormal_mixture_fit_fortran.txt`, respectively. These reproducible
data and report files are ignored by Git. Both R scripts require the CRAN
`mixtools` package. Each fit report includes elapsed wall-clock seconds for
reading the observations, fitting the mixture, and total startup-through-fit
time.

## Minimal example

```fortran
program example
  use mixtools
  implicit none
  type(rng_state) :: rng
  type(mixture_result) :: fit
  real(dp), allocatable :: x(:)

  call rng_seed(rng, 1234)
  allocate(x(400))
  call rnormmix(rng, 400, [0.4_dp, 0.6_dp], [-2.0_dp, 2.0_dp], &
    [0.5_dp, 0.8_dp], x)
  call normalmixEM(x, 2, fit)
  print *, fit%lambda
  print *, fit%mu
end program example
```

## Numerical adaptations

R's formula/list/S3 interfaces are represented by explicit arrays and derived
result types. Several specialized routines depended on external R packages or
stochastic control flow; their documented Fortran equivalents preserve the
model family but are not bit-for-bit ports. See `PORTING.md` and
`TRANSLATION_COVERAGE.md` before comparing exact output with R.

## Licensing

The original package is GPL (>= 2). This translation is distributed under
GPL-2.0-or-later. Original package metadata, R code, and C code are retained in
`original/` with the compiled binary removed.
