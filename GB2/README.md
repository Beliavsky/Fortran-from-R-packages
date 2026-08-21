# GB2-fortran

Modern Fortran/FPM translation of the computational core of the R package
**GB2 2.1.2**. The port targets Fortran 2018 and keeps the upstream
`GPL (>= 2)` licensing as `GPL-2.0-or-later`.

The public facade is the `gb2` module. Lower-level modules are also available
for applications that want a smaller namespace.

## What is translated

The port covers the substantive numerical functionality of the upstream
package:

- GB2 density, CDF, quantile and random generation.
- Ordinary and incomplete moments; moments of `log(X)`.
- Poverty and inequality indicators, including the Thomae-transformed
  hypergeometric Gini formula.
- Pointwise log density, analytic scores and Hessians, Fisher information,
  Fisk starting values, full weighted/household-weighted pseudo-likelihood,
  profile likelihood and profile scores.
- Full and profile GB2 fitting by BFGS.
- Weighted empirical indicators, robust weights, high-level model-vs-empirical
  fitting summaries, and nonlinear indicator fitting.
- Compound GB2 right/left decompositions: component factors, component and
  mixture density/CDF, interval probabilities, moments/incomplete moments,
  quantiles, RNG, mixture-logit parameterization, score/Hessian, fitting, and
  poverty/inequality indicators.
- Compound GB2 with auxiliary covariates: softmax probabilities, starting
  coefficients, likelihood, score, Hessian, fitting and covariance helpers.
- Score-sandwich parameter covariance, numerical indicator Jacobians,
  delta-method indicator covariance, and survey-design score covariance using
  the translated `survey-fortran` dependency.

Plotting/presentation-only routines are intentionally not translated:
`plotsML.gb2`, `saveplot`, `contprof.gb2`, `contindic.gb2`, `dplot.cgb2`, and
`dplot.cavgb2`.

See `docs/TRANSLATION_COVERAGE.md` for a detailed R-to-Fortran mapping.

## Build with FPM

```sh
fpm build
fpm test
fpm run --example basic
```

The default FPM graph uses the included `vendor/survey-fortran` dependency.
No external R installation is required.

## Basic use

```fortran
program demo
  use gb2, only : dp, dgb2, pgb2, qgb2, moment_gb2
  implicit none
  real(dp), parameter :: a=2.3_dp, b=4.2_dp, p=1.7_dp, q=3.4_dp

  print *, dgb2(3.1_dp,a,b,p,q)
  print *, pgb2(3.1_dp,a,b,p,q)
  print *, qgb2(0.5_dp,a,b,p,q)
  print *, moment_gb2(1.0_dp,a,b,p,q)
end program demo
```

For these parameters the reference values are approximately
`0.297754045600`, `0.540341717870`, `2.967145869832`, and
`3.198039862238`.

## Numerical implementation notes

The upstream R package delegates one-dimensional integration to `cubature`
and the Gini calculation's generalized hypergeometric series to `hypergeo`.
The supplied Fortran translations of those two packages have incompatible
GPL-version declarations when linked together (`GPL-3.0-or-later` versus
`GPL-2.0-only`). To keep the default build legally coherent and self-contained,
GB2-fortran implements a 15-point adaptive Gauss-Kronrod integrator and the
convergent real `3F2(1)` series internally. The supplied translations remain in
`vendor/*-reference` but are not linked. See `LICENSES.md`.

Positive GB2 parameters are optimized on the log scale. This has the same
likelihood target as upstream BFGS while preventing invalid negative trial
parameters. Compound CDF integration callbacks and fitting callbacks use saved
module context, so simultaneous calls to the same fitting/CDF routines from
multiple threads are not currently re-entrant.

## Validation

The test suite includes independent reference values from standard special-
function calculations, finite-difference checks of analytic scores/Hessians,
full-vs-profile MLE agreement, compound right/left decomposition checks,
auxiliary-mixture derivative checks, survey sandwich covariance checks,
nonlinear-fit recovery, and RNG moment/range checks.

The project has been compiled and tested with GNU Fortran using:

```sh
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

No nonstandard unlimited-source-line option is required.

## License

GB2-fortran: `GPL-2.0-or-later`, following upstream `GB2` 2.1.2.
See `LICENSE`, `LICENSES.md`, and `UPSTREAM.md`.
