# good-fortran

`good-fortran` is a modern Fortran/FPM translation of the computational code in the R package **good 1.0.2** (Good Regression), originally by Jordi Tur, Amanda Fernandez-Fontelo, David Morina, Pere Puig, Argimiro Arratia, Alejandra Cabana, and David Agis.

The original package is licensed under **GPL (>= 2)**. This translation therefore remains **GPL-2.0-or-later**. See `LICENSE`, `NOTICE`, and `original/DESCRIPTION`.

## Implemented computational functionality

- Good distribution probability mass function: `dgood`
- CDF: `pgood`
- Quantile function: `qgood`
- Random generation: `rgood`
- Mean and variance: `goodmean`, `good_moments`
- Numerically stable Good/polylog normalizer using positive-term log-sum-exp summation
- Good regression by maximum likelihood: `glm_good`
  - log link
  - logit link
  - identity link
  - feasibility-preserving BFGS optimization
  - observed Hessian and covariance matrix
- Predictions and delta-method standard errors: `predict_good`
- Model summaries: `summary_good`
  - coefficient standard errors, z statistics, and p-values
  - AIC/BIC
  - likelihood-ratio comparison with an intercept-only Good model
  - intercept-only comparisons with logarithmic and geometric models

R-specific formula parsing, S3 dispatch/printing, data-frame machinery, and package dataset loading are not translated because they are interface/infrastructure rather than numerical algorithms. The source package contains no plotting routines.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The public umbrella module is `good`.

## Minimal example

```fortran
program demo
    use good, only : dp, dgood, goodmean
    implicit none

    print *, dgood(4, 0.6_dp, -3.0_dp)
    print *, goodmean(0.6_dp, -3.0_dp)
end program demo
```

For regression, pass the response as an integer vector and the already-constructed design matrix explicitly. Unlike R's `glm.good`, the Fortran routine does not parse formulas.
