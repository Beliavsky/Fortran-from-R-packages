# tolerance-fortran

Modern Fortran 2018/FPM translation of the computational portions of the R package `tolerance` 3.0.0 by Derek S. Young and Kedai Cheng.

The upstream package computes statistical tolerance intervals and regions for a broad collection of discrete, continuous, nonparametric, regression, multivariate, and fiducial models. This translation preserves those numerical/statistical kernels while replacing R-specific formulas, S3/data-frame output, optimizers, and graphics with explicit Fortran APIs and derived types.

## Build

With Fortran Package Manager (FPM):

```text
fpm build
fpm test
fpm run --example basic
```

The release was also validated directly with GNU Fortran 14.2 using Fortran 2018 mode, `-Wall -Wextra -Wimplicit-interface -fcheck=all`.

## Main modules

- `tolerance`: convenience module re-exporting the public API.
- `tolerance_normal`: normal K factors, classical/Bayesian/simulation tolerance intervals, and differences of normal means.
- `tolerance_discrete`: binomial, Poisson, negative-binomial, hypergeometric, negative-hypergeometric, UMA, and acceptance-sampling procedures.
- `tolerance_continuous`: exponential, two-parameter exponential, uniform, Laplace, gamma, logistic, Cauchy, Pareto/power, Gumbel, and Weibull tolerance intervals.
- `tolerance_custom`: fitted discrete Pareto, Poisson-Lindley, Zipf, Zipf-Mandelbrot, and Zeta models.
- `tolerance_nonparametric`: Wilks, Wald, Hahn-Meeker, Young-Mathew, beta/order-statistic, and distribution-free calculations.
- `tolerance_regression`: linear, nonlinear, nonparametric regression tolerance bands and ANOVA group limits.
- `tolerance_multivariate`: multivariate-normal and nonparametric multivariate tolerance regions and multivariate-regression regions.
- `tolerance_fiducial`: fiducial two-sample binomial/Poisson/negative-binomial procedures and semicontinuous models.
- `tolerance_distributions`: package-specific probability distributions and RNGs.
- `tolerance_support`: K tables, normal OC inversions, sample-size calculations, and Bonferroni combination.
- `tolerance_math` / `tolerance_optimize`: self-contained probability, integration, linear-algebra, RNG, root-finding, and optimization support.

## API design

R formula/model objects are intentionally not reproduced. Regression routines accept numerical design matrices, and nonlinear regression accepts a model callback. Higher-order R routines such as `bonftol.int` similarly have callback-based Fortran equivalents. Results use derived types such as `tolerance_interval`, `discrete_tolerance_interval`, `regression_band`, and `mv_tolerance_region`.

Random-number streams are not bit-for-bit compatible with R, but the same target distributions are used.

## Plotting and R infrastructure

The `plotly_*` functions, `plottol`, graphical portions of `norm.OC`, package startup code, data-frame formatting, and R formula/deparse helper `rFUN` are omitted. The numerical calculations behind `norm.OC` remain available as `norm_oc_content` and `norm_oc_alpha`.

## License and provenance

The upstream package declares `GPL (>= 2)`. This translation is therefore distributed under GPL-2.0-or-later. Upstream metadata is retained in `upstream/DESCRIPTION` and `UPSTREAM_DESCRIPTION`; translation details are in `TRANSLATION_NOTES.md`.
