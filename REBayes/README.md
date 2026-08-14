# REBayes-fortran v0.1.0

Modern free-form Fortran/FPM computational port of REBayes 2.60.

REBayes implements empirical-Bayes and density-estimation methods centered on
finite-grid Kiefer-Wolfowitz nonparametric maximum likelihood.  The upstream R
package uses the proprietary MOSEK optimizer through the optional Rmosek package.
This port replaces that dependency with a native finite-grid KW likelihood solver,
and uses native Fortran implementations for distributions, posterior summaries,
repeated-measures estimators, regularized logistic regression, and MEDDE.

## Implemented numerical areas

- Kiefer-Wolfowitz primal/dual likelihood fitting on a finite support grid.
- Gaussian location mixtures (`GLmix`).
- Binomial and bivariate-binomial mixtures (`Bmix`, `B2mix`).
- Binomial-Poisson and normal-Poisson mixtures (`BPmix`, `NPmix`).
- Poisson mixtures (`Pmix`), including exposure and truncated-support likelihoods.
- Gamma and Gaussian-variance mixtures (`Gammamix`, `GVmix`).
- Weibull and Gompertz mixtures.
- Student-t location and noncentral-t mixtures (`TLmix`, `Tncpmix`).
- Huber least-favourable mixtures (`HLmix`) and Huber density utilities.
- Uniform mixtures (`Umix`) and Cosslett binary-response mixtures.
- Joint Gaussian mean/variance mixtures (`GLVmix`).
- Weighted/repeated-measures `WGVmix`, `WGLVmix`, `WLVmix`, and `WTLVmix`.
- Posterior mean/quantile/mode calculations and local-FDR calculations.
- KW smoothing, bandwidths, random sampling, quantiles, bivariate marginal
  quantiles, trapezoidal integration, threshold-FDR, inverse-root helper, and
  step-distribution L1 distance.
- General-D regularized logistic regression (`RLR`) via ADMM with Newton
  inner solves.
- `medde` entropy/Renyi density estimation and numerical quantile/simulation
  helpers.

The umbrella module is `rebayes`.

## Deliberately deferred in v0.1.0

Three specialized numerical routines are not silently approximated:

- `BDGLmix` (natural-spline penalized Bayesian deconvolution)
- `HuberSpline` (rotated-quadratic-cone spline problem)
- `HodgesLehmann` (Huber/Mallows contamination SOCP)

`Rxiv` is also omitted because it is an R/LaTeX file-archiving utility rather
than numerical code.  R plotting, S3 classes, formula/data-frame handling, and
MOSEK-control objects are outside the Fortran API.

The experimental nonzero-`alpha` extension of upstream `KWDual` is not exposed
by the v0.1.0 KW solver; Renyi objectives used by `medde` are implemented in the
MEDDE module itself.

## Build

```text
fpm build
fpm test
fpm run rebayes_example
```

All source is `.f90` free form.  The manifest deliberately uses:

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

No BLAS, LAPACK, R, Rmosek, or MOSEK installation is required.

## Basic example

```fortran
use rebayes_kinds, only : dp
use rebayes_mixtures, only : mixture_fit, glmix

type(mixture_fit) :: fit
real(dp) :: x(10), grid(5), sigma(1)

x = [-1.2_dp,-1.0_dp,-0.8_dp,-1.1_dp,-0.9_dp, &
      0.8_dp, 1.0_dp, 1.2_dp, 0.9_dp, 1.1_dp]
grid = [-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp]
sigma = 0.25_dp
call glmix(x, grid, sigma, fit)
print *, fit%mass
```

## Licensing

The supplied REBayes 2.60 package declares `GPL (>= 2)`.  This translation is
distributed under GPL-2.0-or-later.  Complete GPL-2 and GPL-3 texts and the
upstream DESCRIPTION/CITATION are retained in the release.
