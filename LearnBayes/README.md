# LearnBayes Fortran

`learnbayes-fortran` is a modern Fortran/FPM translation of the computational
functionality in the R package **LearnBayes 2.15.2**.

The goal is to make the package's Bayesian teaching and simulation algorithms
usable directly from Fortran while keeping the original formulas recognizable.
R graphics, S3 presentation, formulas, and dynamic R function objects are not
recreated.

## Highlights

The public facade is:

```fortran
use learnbayes
```

and re-exports `dp`, the numerical/statistical routines, callback/result types,
and the deterministic native random-number generator.

Implemented areas include:

- beta and normal prior elicitation;
- beta-binomial, normal-normal, Poisson-gamma, and discrete predictive models;
- one- and two-parameter discrete Bayes calculations and HPD sets;
- contingency-table Bayes factors;
- normal, beta, gamma, Poisson, binomial, multivariate-normal and
  multivariate-t numerical infrastructure;
- Dirichlet, inverse-gamma, multivariate-normal, multivariate-t, and truncated
  normal simulation;
- generic user log-density callbacks;
- Laplace approximation using Nelder-Mead optimization plus a numerical
  Hessian;
- Metropolis-within-Gibbs, random-walk Metropolis, independence Metropolis,
  importance sampling, rejection sampling, SIR, and gridded contour sampling;
- beta-binomial, Poisson-gamma, normal-normal, logistic, grouped-normal,
  Bradley-Terry, Cauchy-error, g-prior regression, transplant, and Weibull
  posterior functions;
- Bayesian linear regression, expected/predictive responses, Bayesian residual
  diagnostics, probit Gibbs sampling, and exhaustive Bayesian model selection;
- normal posterior simulation/prediction, robust Student-t Gibbs sampling,
  hierarchical Gibbs sampling, order-restricted Gibbs sampling, independence
  Bayes-factor Monte Carlo, and influence diagnostics;
- `regroup`, career-trajectory data preparation, contour-grid values, and
  beta-prior/likelihood/posterior triplot data.

The plotting-only `predplot()` calculation is already exposed through `pbetap()`;
`summary.bayes()` corresponds to `summarize_discrete()`.

## Example

```fortran
program example
   use learnbayes
   implicit none

   type(rng_state) :: rng
   type(blinreg_result) :: fit
   real(dp) :: x(6, 2)
   real(dp) :: y(6)
   integer :: info

   x(:, 1) = 1.0_dp
   x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   y = [1.1_dp, 2.8_dp, 5.2_dp, 6.9_dp, 9.1_dp, 11.2_dp]

   call rng_seed(rng, 20260830_i8)
   call blinreg(rng, y, x, 4000, fit, info)
   if (info /= 0) error stop "blinreg failed"
end program example
```

See `example/learnbayes_example.f90` for a complete example.

## FPM

```text
fpm build
fpm test
fpm run --example learnbayes_example
```

The package has no external numerical-library dependency.

## Validation

`tools/run_strict_tests.sh` builds the maintained source with GNU Fortran using:

```text
-std=f2018
-Wall
-Wextra
-Werror
-Wimplicit-interface
-fimplicit-none
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The tests include deterministic formula fixtures and seeded simulation smoke
checks. `tools/check_source_rules.py` also enforces the maintained Fortran source
rules.

## Source conventions

- `dp = real64` is defined once in `learnbayes_kinds` and re-exported publicly.
- Maintained real variables use `real(dp)` and real constants use `_dp` suffixes.
- Every procedure dummy has explicit `INTENT` or `VALUE`.
- Each dummy is declared separately.
- Each dummy declaration has a meaningful trailing FORD `!!` documentation
  comment.
- Maintained source stays within the normal free-form 132-column limit and is
  written to be compatible with `fprettify`.

## Deliberate interface adaptations

LearnBayes is an educational R package and several functions accept arbitrary R
functions plus `...`. The Fortran translation replaces those with typed
`log_density_callback` and `likelihood_callback` objects carrying optional
numeric context.

The upstream generic `rtruncated(n, lo, hi, pf, qf, ...)` is not reproduced as
a variadic runtime-function interface. `rtruncated_normal` implements the
truncated-normal use required internally by LearnBayes; other distributions can
be built using the public CDF/quantile and RNG infrastructure or typed callback
APIs.

Graphical functions (`plot.bayes`, `plot.bayes2`, `predplot`, `mycontour`, and
`triplot`) are not rendered by Fortran. Their useful calculations are available
as ordinary arrays through existing probability routines, `contour_grid`, and
`triplot_data`.

See `TRANSLATION_NOTES.md` for detailed correspondence and qualifications.

## License

GPL-2.0-or-later. See `LICENSE` and `NOTICE.md`.
