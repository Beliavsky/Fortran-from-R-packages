# Translation notes

## Upstream

- Package: `gamlss`
- Upstream version: 5.5-0
- Upstream license declaration: GPL-2 | GPL-3
- Upstream package date: 2025-08-19

Selected original computational R sources and the compiled C source are retained
under `upstream/` for provenance and algorithm comparison. v0.2 additionally
retains fractional-polynomial, LOESS, monotone-spline, categorical-fusion and
random-effect reference sources. v0.3 retains the upstream random-effect,
stepwise-GAIC and profile-likelihood sources used for the new numerical layers.

## Numerical design

R `gamlss()` creates up to four model frames and delegates each parameter update
to GLIM/backfitting code. The Fortran port receives numeric design matrices
directly:

    eta_mu    = X_mu    beta_mu    + offset_mu
    eta_sigma = X_sigma beta_sigma + offset_sigma
    eta_nu    = X_nu    beta_nu    + offset_nu
    eta_tau   = X_tau   beta_tau   + offset_tau

Links and likelihoods come from the vendored `gamlss.dist` translation.

### RS and automatic penalty estimation

`GAMLSS_METHOD_RS` performs parameter-wise Fisher/IRLS updates. v0.2 can also
estimate a scalar penalty multiplier during an RS update. For a quadratic
penalty `P`, the update uses the working residual variance and the penalized
coefficient variance implied by `beta' P beta`, with effective degrees of
freedom correcting for the null space of `P`. Geometric damping is used when
the variance ratio changes substantially.

This provides the numerical variance-component role needed by random effects
and ML-style penalized smoothers without introducing an R smoother object.

### Marginal prediction (`getMarginal`)

v0.9 translates upstream `getMarginal()` from `R/random.R`.  The R routine first
forms `etam = lp - rt`, then marginalizes the selected parameter inverse link
over `N(0, sigma_b^2)`.  The Fortran low-level API accepts `etam` directly; the
random-intercept adapter reconstructs it from the fitted parameter linear
predictor, group labels and stored random effects.

The deterministic quantile method intentionally uses the exact upstream grid
0.001 through 0.999 by 0.001.  The Monte Carlo method intentionally defaults to
10,000 draws.  `random_intercept_result_t%sigma_b` uses the same variance
component summary as upstream `random()`: the square root of the random
coefficient sum of squares divided by its effective degrees of freedom.

### Random effects

`fit_gamlss_random_intercept` remains the compact scalar-variance adapter from
v0.2. v0.3 adds `fit_gamlss_random_effects`, which accepts an arbitrary
within-group random-effect design matrix. The group blocks share either a
diagonal or full covariance. After each GAMLSS fit the covariance is updated
from the conditional group effects plus their posterior covariance blocks, then
inverted to form the next quadratic penalty. The update is geometrically damped.
For Gaussian `mu` models the supplied `nlme` `fit_lme` covariance can initialize
the iteration.

This covers common correlated random intercept/slope models. The full `nlme`
residual-correlation and variance-function catalog is not replicated inside the
GAMLSS working-response loop.

### Censoring

`fit_gamlss_censored` maximizes the exact censored likelihood directly. For a
continuous family the contributions are

- exact: `log f(y)`;
- left: `log F(u)`;
- right: `log(1-F(l))`;
- interval: `log(F(u)-F(l))`.

For discrete families the censoring bounds are converted to the corresponding
integer CDF cutoffs. `family_cdf` dispatches over all 62 generic backend family
constants. This is a numerical `Surv`/censored-family adapter rather than an R
class implementation. v0.3 additionally conditions every contribution on
`Y > entry` when delayed-entry times are supplied, and adds interval2 and
counting-process array adapters.

### Smoothers

P-splines continue to use the supplied spline translation. v0.2 adds:

- fractional-polynomial exhaustive search on the upstream power grid;
- local polynomial LOESS with tricube neighborhoods in one or two predictors;
- varying-coefficient spline designs formed by row-wise multiplication by the
  modifier variable;
- monotone P-splines using active first-difference inequality penalties.

The monotone solver preserves the computational idea of upstream `pbm()` but
uses direct penalized normal equations rather than reproducing its R QR/SVD
object structures.

### Penalized categorical fusion

`fit_pcat` uses the translated all-pairs difference matrix `D`, updates adaptive
weights

    omega_j = 1 / (|D_j beta|^(2-Lp) + kappa^2),

and refits the penalized weighted least-squares problem. The optional ML-style
lambda update follows the local variance-ratio calculation used by upstream
`gamlss.pcat`.

### Selection and profiling

The matrix stepwise and profiling routines preserve the numerical refitting
operations while intentionally omitting R formula manipulation. A profiled
coefficient is fixed by removing its design column and placing its contribution
in the parameter offset, then refitting all remaining coefficients.


### Additive and tensor smooths

v0.3 composes multiple persistent P-spline terms with explicit sum-to-zero
coefficient contrasts, avoiding the constant-space aliasing that would arise
from naively concatenating several spline bases plus an intercept. The tensor
2D builder uses the row-wise tensor basis and anisotropic penalty

    lambda_x (P_x kron I_y) + lambda_y (I_x kron P_y).

Both additive and tensor specifications retain knots for prediction.

### Bootstrap and confidence intervals

The case bootstrap resamples complete design/response rows and refits the same
matrix model, so it is available to all translated backend families without a
family-specific random generator. Profile confidence intervals use the usual
one-degree-of-freedom likelihood-ratio cutoff and linearly interpolate the
profile grid crossings.

### Residual correlation and variance functions

v0.4 adds `fit_gamlss_no_gls`, which delegates an exact Gaussian/NO location
model to the supplied `nlme` GLS implementation.  This exposes its complete
translated correlation and variance-function catalogs without approximating
them through independent working weights.  The full correlated Gaussian
log-likelihood, fitted covariance matrix and correlation/variance parameters
remain available in `result%gls`; the embedded GAMLSS-style `model` is a
convenient view.  Its per-case deviance field contains marginal standardized
residual squares, not an additive decomposition of the correlated joint
log-likelihood.

For non-Gaussian families, residual correlation is still not iterated through
every RS/CG working response.

### Multiple-parameter random intercepts

`fit_gamlss_multi_random_intercept` augments each selected parameter design with
the same compressed group-indicator matrix.  Each active parameter has its own
identity penalty block and RS estimates its scalar penalty multiplier.  The
random effects are independent across distribution parameters in this wrapper.
This is intentionally distinct from v0.3 `fit_gamlss_random_effects`, which can
fit an unstructured intercept/slope covariance for one selected parameter.

### Parameter-wise selection and diagnostics

`stepwise_gaic_parameter` repeats the matrix scope search while replacing only
the target parameter design.  The other parameter designs remain fixed.  Worm
plot support is numerical only: residuals are split into equal-count covariate
bins, sorted against normal scores, and a cubic polynomial is fitted to the
detrended Q-Q coordinates.

### Callback portability and reentrancy

The v0.4 build removes GNU nested-procedure trampolines from the optimizer and
integration paths that previously caused executable-stack linker warnings.
The replacement callback contexts are module-level saved state because the
existing optimizer/integrator interfaces accept only a bare procedure argument.
Consequently those specific callback-based fit/integration routines are not
thread-safe or reentrant.  A future context-aware optimizer interface could
remove that limitation without reintroducing executable-stack trampolines.

## Remaining differences from R gamlss

The following remain incomplete or intentionally excluded:

- formula/model-frame/terms/contrast parsing;
- S3/S4 object infrastructure and update methods;
- plotting, centile fans/plots, worm/bucket plots and graphical diagnostics;
- parallel and formula-scope stepwise wrappers;
- exact Gaussian/NO covariance likelihood remains available through `fit_gamlss_no_gls`;
  v0.5 additionally embeds correlation matrices in non-Gaussian RS working responses,
  but does not jointly optimize non-Gaussian variance-function parameters;
- exact replication of every smoother control heuristic and R prediction object;
- R profile/bootstrap presentation objects and graphical bootstrap diagnostics.

The numerical kernels implemented through v0.3 are exposed directly rather
than being represented as R-compatible objects.


### Correlated RS working responses (v0.5)

For parameter j, v0.5 first computes the ordinary GAMLSS Fisher working
response `z_j` and diagonal curvature weights `w_j`.  GLS is then applied with
working covariance `D_j V R V D_j`, where `D_j = diag(1/sqrt(w_j))`, `R` is an
`nlme` correlation matrix, and `V` contains the supplied base variance-function
standard deviations.  The location block can update a shared non-fixed
correlation parameter vector; later parameter blocks reuse that correlation.
This remains a correlated working-response/estimating-equation algorithm. v0.6 adds a separate exact
continuous-margin Gaussian-copula joint-likelihood API rather than changing the RS estimator itself.

### Multi-parameter random slopes (v0.5)

Each selected distribution parameter may receive the same q-dimensional
within-group design shape but with its own observed design values and its own
q by q covariance matrix.  The covariance update is the average posterior
second moment across groups.  The implementation supports correlation among random terms within a parameter. v0.6 adds
`fit_gamlss_joint_random_effects` when cross-parameter covariance is required.


### Continuous-margin Gaussian copula likelihood (v0.6)

`fit_gamlss_gaussian_copula` is intentionally separate from the v0.5
correlated-RS working-response estimator.  For continuous margins, observations
are transformed to `z_i = Phi^{-1}(F_i(y_i))`; within each supplied correlation
group the exact Gaussian-copula contribution is

    -0.5 * [log|R| + z' R^{-1} z - z' z].

This is added to the marginal GAMLSS log densities, and the marginal regression
coefficients and non-fixed correlation parameters are optimized together.
The routine supports the translated `nlme` temporal/spatial correlation
structures.  Discrete and mixed point-mass margins are deliberately rejected:
their exact Gaussian-copula likelihood is a multivariate normal rectangle
probability, not the continuous copula density above.

### Joint cross-parameter random effects (v0.6)

`fit_gamlss_joint_random_effects` uses one Gaussian random-effect covariance
`Sigma` over all active `(parameter, random-term)` components.  Conditional on
`Sigma`, all fixed and random coefficients are optimized jointly under the
GAMLSS log likelihood plus `0.5 sum_g b_g' Sigma^{-1} b_g`.  The covariance is
then updated from the average posterior second moment
`b_g b_g' + Var(b_g | y)` using the optimizer inverse-Hessian approximation.
This directly estimates cross-parameter terms such as `Cov(b_mu,b_sigma)` and
extends naturally to cross-covariances between random slopes.  It is a
penalized/Laplace-style covariance iteration rather than exact high-dimensional
quadrature of the marginal random-effects likelihood.


### Discrete and mixed Gaussian copulas (v0.7)

For an atomic observation the latent Gaussian coordinate is not observed at one
point.  It lies in the interval

    Phi^-1(F(y-)) < Z <= Phi^-1(F(y)).

For a mixed correlation group, v0.7 conditions the discrete latent coordinates
on the exactly transformed continuous coordinates.  If `C` denotes continuous
coordinates and `D` atomic coordinates, the group likelihood is

    prod_{i in C} f_i(y_i)/phi(z_i)
      * phi_Rcc(z_C)
      * Pr(a_D < Z_D <= b_D | Z_C=z_C).

Thus no continuous copula-density approximation is used for discrete margins.
The final rectangle probability is evaluated numerically with a deterministic
Genz sequential transform and antithetic Halton points.  It is the correct
rectangle likelihood formula, with numerical integration error controlled by
`n_qmc`.  `BEINF` endpoints are treated as atoms while observations strictly
inside `(0,1)` remain continuous.

### Direct marginal random-effects quadrature (v0.7)

`fit_gamlss_joint_random_effects_ghq` parameterizes the full joint random-effect
covariance by a Cholesky factor and maximizes the group-marginal likelihood
obtained by tensor Gauss-Hermite integration.  The integral includes all active
GAMLSS parameter/random-term combinations simultaneously, so covariance between
`mu` and `sigma` random effects is part of the integrated model rather than a
post-hoc block estimate.

Tensor quadrature is deliberately limited to four combined latent components.
At order 7, four dimensions already require 2401 likelihood evaluations per
group per objective evaluation.  Larger blocks should use the v0.6
Laplace/posterior-moment fitter.  This is a computational scaling boundary, not
a difference in the intended Gaussian random-effects model.


### Adaptive importance sampling (v0.8)

Tensor Gauss-Hermite quadrature scales as `order**r` for latent dimension `r`.
The v0.8 AIS path instead constructs, for every group, a Gaussian proposal at
the conditional posterior mode.  Its covariance is the inverse numerical
posterior Hessian, multiplied by a configurable proposal-scale factor.  The
marginal group likelihood is evaluated as the log mean of deterministic
importance weights from Halton/antithetic normal points.  Posterior means and
covariances use the same normalized weights; ESS is reported to diagnose poor
proposal overlap.

The scalable default uses the v0.6 Laplace fit as the parameter estimate and
uses AIS to improve marginal-likelihood and posterior-moment evaluation.  An
optional full refinement reoptimizes fixed effects and the Cholesky parameters
of the joint random-effect covariance under the AIS likelihood.  This remains
a numerical integration method: increasing `qmc_points` is the primary way to
check stability.
