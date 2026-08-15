# Translation notes

## Upstream

- Package: VGAM
- Upstream version: 1.1-14
- Upstream date: 2025-12-02
- Principal author: Thomas W. Yee
- Upstream license field: GPL-3

The original `DESCRIPTION`, `NAMESPACE`, `NEWS`, `ChangeLog`, and
`LICENCE.note` are retained in `original/`.

## Design choices

### Arrays instead of R objects

R's `vglm`, `vgam`, and family objects depend heavily on formulas, model frames,
slots, environments, expression evaluation, and callbacks. The Fortran port
exposes those numerical operations through ordinary arrays and derived result
types. Predictors are passed as design matrices explicitly.

### Numerical family API

One-parameter exponential-family models use IRLS in `vgam_vglm`. Families with
coupled or nuisance parameters use a BFGS likelihood optimizer plus numerical
Hessian in `vgam_extended_models`, `vgam_count_models`, `vgam_qreg`, and
`vgam_timeseries`.

### Categorical models

`fit_multinomial` uses a joint baseline-category likelihood and the full block
information matrix. `fit_ordinal` parameterizes cutpoint increments on the log
scale so cutpoints remain strictly ordered during optimization.

### Constraints

VGAM's constraint-matrix idea is retained directly. A matrix `C` maps free
parameters to the full stacked coefficient vector. `parallel_constraint`
constructs a common VGAM pattern where selected covariate effects are shared
across linear predictors.

### Reduced-rank VGLM

Version 0.2 adds `fit_rrvglm`. It follows the central alternating reduced-rank
parameterization used by upstream RR-VGLMs:

```text
eta = X1 B1 + (X2 C) A^T + offset
```

`X1` is the unrestricted block (by default a constant first column), while the
coefficient block for `X2` has rank at most `Rank`. Initialization starts from
ordinary per-response VGLM fits and a symmetric-eigendecomposition low-rank
projection. Alternating working-weight least-squares updates then estimate the
loadings/unrestricted coefficients and latent coefficient matrix. Component
scales/signs are normalized only for identifiability; fitted predictors are
unchanged.

This is the reusable numerical RR-VGLM core, not a complete clone of all
upstream RR control paths. Version 0.4 adds a separate `fit_drrvglm` path. It
represents each loading column as `A(:,r) = H.A(:,:,r) alpha_r` and each row of
`C` as `C(k,:)^T = H.C(:,:,k) gamma_k`, with optional active-column counts so
Fortran 3-D arrays can represent the upstream list of differently sized
constraint matrices. The constrained likelihood is optimized directly. The
full R corner-mask normalization, labeling, QR-summary, and HDE-test machinery
is still outside the numerical API.

### Quadratic reduced-rank VGLM

Version 0.3 adds a QRR-VGLM computational core. For latent scores `z = X2 C`,
each response has

```text
eta_j = X1 b_j + A_j z + z^T Q_j z + offset_j
```

where `Q_j` is a full symmetric `Rank x Rank` matrix. Thus the implementation
has the upstream `Rank*(Rank+1)/2` quadratic degrees of freedom rather than an
axis-aligned squared-term approximation. Given `C`, the response-specific
linear and quadratic terms are estimated by the ordinary family IRLS engine;
`C` is then updated against the actual summed family deviance. Latent columns
are rescaled with compensating transformations of both `A` and `Q`, preserving
fitted predictors. A `dzero` mask implements the common upstream case in which
selected responses are linear, not quadratic, in the latent variables.

`qrrvglm_result_t%optima` solves `A_j + 2 Q_j z = 0` when `Q_j` is nonsingular
and returns the Hessian/curvature `2 Q_j`. Version 0.4 exposes this machinery
through `fit_cqo`, adds direct evaluation of latent response surfaces, and adds
reverse calibration that estimates latent scores from new multivariate
responses by minimizing the fitted family deviance. A minimum-norm inverse map
from latent scores to the reduced environmental predictor block is also
provided. This is the numerical CQO layer, not the full R control/S4 layer.

### Constrained additive ordination

Version 0.4 adds `fit_cao_rank1`. This follows the important current upstream
restriction that CAO is rank 1 with an intercept-only unrestricted component.
For a canonical score `nu = X C`, each response is fitted against `nu` with a
cubic penalized spline using the vendored splines backend. The algorithm
alternates between the response smooths and the identified canonical coefficient
vector. It supports the Gaussian, Poisson, and binomial families in the present
Fortran family core. Upstream smart prediction, control-object options, and
compiled derivative accelerators are not replicated.

### GAITD

The upstream GAITD machinery dynamically constructs large families with
separate alteration/inflation/deflation submodels. Version 0.2 ports a reusable
finite-support computational kernel rather than the R family generator.
`gaitd_transform_pmf` starts from a base PMF, applies truncation and multiplicative
deflation, reserves user-specified final altered masses and additive inflation
masses, and renormalizes the nonspecial mass. Poisson and negative-binomial
wrappers, moments, CDF, quantiles, and RNG are supplied.

Version 0.4 adds a regression layer for Poisson and negative-binomial bases.
The ordinary mean and the special masses have separate design matrices. Special
masses are parameterized jointly with a multinomial logit, with the remaining
category representing the normalized baseline distribution. Altered points are
removed from the baseline before their modeled mass is added; inflated points
retain their baseline contribution and receive modeled extra mass. Fixed
truncation is supported.

Version 0.5 adds the direct MLM interpretation used by upstream `dgaitd*`
routines. Altered points replace the parent mass with a supplied/modelled
probability, inflated points add probability to the scaled parent mass, and
deflated points subtract probability from it. The nonspecial parent mass is
scaled by the upstream normalizer

```text
Delta = (1 - sum(p_altered) - sum(p_inflated) + sum(p_deflated)) /
        (1 - parent_mass(truncated) - parent_mass(altered))
```

for unbounded Poisson/NB parents (and by the corresponding finite parent CDF in
the finite-support constructor). `fit_gaitd_mlm_*_regression` gives these direct
special probabilities covariate-dependent multinomial-logit predictors and
checks pointwise deflation/inflation validity during likelihood optimization.

Version 0.6 adds direct outer-distribution `a.mix`/`i.mix`/`d.mix` semantics for
Poisson and negative-binomial parents. Each supplied total mix mass is allocated
within its named point set in proportion to a separately parameterized outer
Poisson/NB distribution, matching the recursive restricted-support construction
in upstream `dgaitdpois`/`dgaitdnbinom`. These mix components may coexist with
direct MLM altered/inflated/deflated masses and truncation.

The complete generated GAITD parameterization still has gaps: especially
covariate-dependent regression for the outer-distribution parameters/mix totals,
generated equality/parallel constraint and EIM machinery, and every
distribution-specific wrapper.

### Bivariate copulas

Version 0.5 ports five important copula kernels from `family.bivariate.R`:
Clayton, Frank, FGM, Gaussian, and Plackett. Version 0.6 adds the
Ali-Mikhail-Haq (AMH) density/CDF/RNG and integrates AMH into the same
one-parameter regression API. The Fortran formulas retain VGAM's
parameter conventions; in particular Frank uses a positive `apar` with
`apar = 1` representing independence. Density, CDF, and simulation routines are
provided. The Gaussian copula CDF is evaluated by deterministic one-dimensional
quadrature of the conditional-normal representation.

`fit_copula_regression` fits a single copula dependence predictor against an
ordinary design matrix using maximum likelihood. Positive parameters use a log
link; FGM/Gaussian dependence use bounded hyperbolic-tangent links. Intercept
starts use sample dependence information, mirroring the purpose of VGAM's
correlation/Kendall initialization and avoiding boundary saturation in numerical
optimization.

Version 0.6 also ports the upstream bivariate Student-t density and a numerical
likelihood fit with design matrices for degrees of freedom and correlation. A
self-contained univariate Student-t CDF/quantile/RNG supports the upstream
`dbistudenttcop` density and Student-t copula simulation/fitting without Rmath.
The Student-t CDF implementation was cross-checked against SciPy on a grid with
maximum absolute error approximately 4.3e-16.

The remaining bivariate catalog (other specialized copulas, bivariate
exponential/gamma families, exchangeable binary families, and their full VGLM
expected-information machinery) remains future work.

### Altered/inflated beta helpers

Version 0.6 ports the direct `dzoabeta`, `pzoabeta`, `qzoabeta`, and `rzoabeta`
computations for a beta law with explicit masses at zero and one. It also ports
the zero/one-inflated beta-binomial probability/CDF/RNG kernel in the
shape-parameter form. R family-object constraints and EIM generation are not
replicated.

### Yeo-Johnson / quantile work

The upstream `yeo.johnson` transform, inverse, derivative with respect to
lambda, and Jacobian derivative are represented directly. `fit_yj_normal`
provides a likelihood-based regression in which transformed responses are
normal with linear mean and constant scale, including the transformation
Jacobian and quantile prediction on the original scale.

Version 0.3 adds `fit_lms_yj`, with separate design matrices for a varying
transformation parameter `lambda(x)`, transformed location `mu(x)`, and
`log(sigma(x))`. The full Jacobian likelihood is optimized jointly and
original-scale conditional quantiles are available. Passing B-spline or
natural-spline design matrices from `vgam_smoothing`/`splines` gives an
additive LMS-style computational fit without reproducing R's smoother/formula
objects. Upstream control paths and every `lms.yjn` smart-prediction option are
not yet mirrored one-for-one.

### Time series

Version 0.2 translated the computational core of `dAR1` and a Gaussian AR(1)
MLE with exact stationary or conditional likelihood.

Version 0.3 adds GARMA(p,0). This matches the current upstream restriction that
`q.ma.lag` must be zero. The implemented observation-driven predictor is

```text
g(mu_t) = x_t beta + sum_i phi_i [g(y_{t-i}) - x_{t-i} beta]
```

with stable boundary handling for binary and positive-response links, direct
likelihood/quasi-likelihood optimization, numerical Hessian covariance, and
recursive forecasting.

Version 0.4 adds nested reduced-rank autoregression (`fit_rrar`). Lag ranks must
be non-increasing. The lag matrices share an identified left subspace while each
lag has its own right factor at the requested rank; parameters are estimated by
concentrated Gaussian likelihood. The result includes lag matrices, innovation
covariance, optional numerical-Hessian covariance, and recursive forecasts. This
is the Ahn-Reinsel/VGAM numerical model without the VGLM object wrapper.

### Smoothing

The user-supplied modern Fortran translation of the R `splines` dependency is
vendored. `vgam_smoothing` constructs B-spline or natural-spline bases and
second-difference (or configurable-order) penalties, then feeds them to the
IRLS core.

### Special functions

The Fortran library avoids a dependency on Rmath. Special functions needed by
translated distributions are implemented in `vgam_special`, using stable
series, continued fractions, asymptotic expansions, and iterative inverses.

### Version 0.7 additions

Version 0.7 adds `vgam_gaitd_mix_regression`. The parent parameter, three
special-mass logits, and three outer-distribution means are explicit regression
blocks. Inactive mix classes receive exactly zero mass. The probability kernel
uses the GAITD parent renormalization together with restricted outer-law weights
on the named `a.mix`, `i.mix`, and `d.mix` supports. The negative-binomial
version currently shares one estimated size across the parent and outer NB laws;
separate outer-dispersion regressions remain future work.

`vgam_censored` adds exact/left/right/interval likelihood contributions matching
the upstream censoring conventions for normal, Poisson, exponential, and
Rayleigh families. Fitting is by BFGS with a numerical Hessian rather than by
reproducing each R family's generated expected-information expression.

`vgam_bivariate_extra` adds the bivariate normal model, the bivariate logistic
distribution used by `bilogistic`, and Freund's bivariate exponential model. The
bivariate logistic log-density is evaluated with log-sum-exp; this is important
because directly capping `exp(-z)` can create a false likelihood ridge when a
scale becomes very small.

`vgam_information` supplies reusable information-matrix operations needed by
translated VGLM families: numerical observed information, score outer-product
information, projection through a coefficient-constraint matrix, and covariance
lifting from the free parameterization to the full coefficient space. This is
not a claim that every family-specific generated VGAM EIM expression has been
ported.

### Version 0.8 additions

Version 0.8 adds `vgam_gaitd_nb_dispersion`.  The parent negative-binomial
size and the altered, inflated, and deflated outer-distribution sizes each have
an explicit log-dispersion regression block.  The probability calculation still
uses the same parent renormalization and restricted outer-support weights as the
v0.6/v0.7 GAITD mix code; only the dispersion parameterization is generalized.
When an outer mix class is active, its mean and size are therefore allowed to
vary independently with covariates.

`vgam_zero_altered` ports direct computational helpers for positive-normal and
positive-geometric laws and for common zero-altered count distributions.  The
zero-altered Poisson, negative-binomial, geometric, and binomial functions use an
explicit observed mass at zero and the corresponding zero-truncated parent law
on the positive support.  The zero-inflated/deflated binomial and geometric
helpers retain the upstream admissible negative `pstr0` range instead of silently
clamping deflation to zero.

`vgam_zero_altered_models` adds likelihood regression for those four
zero-altered count families, with separate parent-parameter and observed-zero
design matrices.  The negative-binomial version jointly estimates its size
parameter.  These are numerical model APIs rather than translations of the R
family-object initialization and generated EIM expressions.

`vgam_multivariate_extra` adds the upstream trivariate-normal computational
family and the FGM-exponential (`bifgmexp`) construction.  The trivariate normal
uses the full three-correlation matrix with positive-definiteness checks, direct
simulation, maximum likelihood, and numerical-Hessian covariance.  Its analytic
log-density was independently checked against SciPy's multivariate normal over a
random grid; the largest absolute log-density difference in that check was about
`1.1e-12`.  `bifgmexp` combines unit-exponential margins with VGAM's FGM
parameterization and provides density, CDF, simulation, and dependence
regression.  A direct Kendall-tau helper is also exposed.


### Version 0.9 additions

Version 0.9 adds `vgam_dirichlet`, including direct Dirichlet density/simulation,
joint log-shape regression, optional parallel slopes, and exact expected
information matrices.  For shape parameters `alpha`, the raw expected information
is `diag(trigamma(alpha)) - trigamma(sum(alpha))`, with the log-shape version
obtained by the corresponding diagonal Jacobian transformation.  This is the
first family-specific analytic EIM added beyond the generic numerical/outer-score
information tools.

`vgam_normal_special` ports the upstream Tobit point-mass semantics and the
generalized folded-normal helpers.  Tobit regression allows separate latent mean
and log-scale design matrices with fixed lower/upper censoring limits.  Folded
normal regression likewise supports separate mean/scale designs and arbitrary
positive `a1`/`a2` fold constants.

`vgam_positive_count` adds full d/p/q/r APIs for zero-truncated Poisson and
negative-binomial laws and an NB likelihood regression with estimated size.
`vgam_zoa_beta_models` adds a zero/one-altered beta likelihood in which beta mean,
precision, and both endpoint masses may have separate predictor blocks; a
three-category softmax enforces valid endpoint/interior probabilities.

## Not a semantic target

The following R behaviors are intentionally outside the port:

- formula parsing / `model.frame` / contrasts machinery;
- S3/S4 object identity and dispatch;
- `terms`, `update`, `anova`, print/summary methods whose work is primarily R
  object orchestration;
- plotting, preplot objects, color and graphics helpers;
- dynamic R callbacks and expression evaluation.

## Computational work remaining after v0.9

Major advanced areas not yet represented as complete ports include:

- remaining DRR corner-mask normalization, summary covariance partitioning, and
  R control/object paths beyond the direct H.A/H.C constrained numerical fit;
- higher-level CQO/CAO control, smart-prediction, and diagnostic object methods;
- the complete generated GAITD family system beyond the outer-mix and separate-
  dispersion regressions, especially exact family-specific EIM/control generation
  and the full wrapper catalog;
- the remaining copula/bivariate family catalog beyond
  Clayton/Frank/FGM/Gaussian/Plackett/AMH, Student-t, bivariate normal/logistic,
  Freund, trivariate normal, and FGM-exponential additions;
- remaining LMS smart-prediction/control variants and other quantile-regression
  smoothing families;
- specialized positive, altered, inflated, and censored variants beyond the
  positive count, Tobit/folded-normal, zero/one-altered beta, and existing
  censoring likelihoods now provided;
- every compiled legacy kernel reachable only from those still-unported areas.

This boundary is explicit so downstream users do not mistake v0.9 for a
bit-for-bit replacement of the entire R package.

## Validation

A clean GNU Fortran build is run with:

```text
gfortran -std=f2018 -Werror=implicit-interface -fcheck=all
```

The tests exercise special functions and distributions, link inversions,
actuarial/extreme-value inversions, ordinary IRLS, multinomial and ordinal
models, beta/negative-binomial/zero-inflated-Poisson regression, constraint
matrices, spline additive prediction, reduced-rank fitting, new count models,
GAITD transforms and regressions (including MLM deflation and mix semantics),
copula kernels/regression, Student-t computations, altered beta helpers,
Yeo-Johnson identities/fitting, AR(1) fitting,
full symmetric QRR/CQO surfaces and calibration, rank-1 CAO, constrained DRR,
three-predictor LMS/Yeo-Johnson fitting, GARMA(p,0), and nested RRAR.
