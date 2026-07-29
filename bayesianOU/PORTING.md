# Porting notes

## Scope

The port translates the package's model equations, transformations, fitting
workflow, diagnostic calculations, multiple-imputation rules, and geometry
engine. R S3 classes, graphics, spreadsheet/text export formatting, cmdstanr and
rstan object plumbing, parallel R orchestration, and progress messages are not
reproduced.

## Inference backend substitution

The most important difference is the posterior engine.

### Original

The R package compiles `ou_nested.stan` and uses Stan NUTS. The Stan model
jointly samples:

- hierarchical Level-1 parameters;
- latent Level-1 stochastic-volatility states;
- latent Level-2 production prices;
- optional Level-2 stochastic volatility and Student-t tails;
- measurement uncertainty and optional capital reconstruction; and
- all hyperparameters.

### Fortran

The Fortran implementation is self-contained and uses:

- nonlinear least squares and Nelder-Mead for cubic Level-1 initialization;
- OLS and LAPACK for pooled effects and Level-2 initialization;
- training-weighted COM-based empirical-Bayes partial pooling;
- a conditional AR(1) log-squared-residual filter for the SV path;
- the standardized production-price anchor as the conditional latent path;
- blocked random-walk Metropolis updates for structural parameters; and
- the translated priors and likelihood for Metropolis acceptance.

This preserves the structural equations and provides actual posterior-oriented
simulation, but it is an approximation to the original joint Stan posterior.
It does not claim draw-for-draw equality, identical uncertainty, or NUTS-level
mixing. R-hat and ESS are therefore reported and must be inspected rather than
assumed satisfactory.

The generic `ou_geom_hmc` implementation is a separate exact static-HMC engine
for arbitrary callback-defined differentiable targets. It supports constant
Euclidean mass matrices and position-dependent SoftAbs metrics.

## Nested paths

For `n_levels >= 2`, the original Stan program samples `Phi_lat`. The native fit
uses the standardized production-price anchor as the conditional path and fits
the Level-2 transition around that path. `simulate_ou_nested` does simulate a
latent production path from the translated transition equation.

The Level-3 value-coupling term `m_v * V` is translated. Optional uncertain-K
reconstruction inputs are represented in `ou_input`, but the native conditional
fit does not sample the full `z_K` field. The original Stan implementation is
retained for provenance.

## Stochastic volatility

The stationary AR(1) transformation and scale convention are retained. The
native fit obtains a robust conditional path from log squared residuals and
estimates AR(1) parameters by conditional Gaussian regression. Posterior draws
for SV summaries are local conditional draws rather than a joint latent-state
NUTS sample.

## LOO

The original uses the `loo` package's PSIS implementation. The Fortran port
implements a self-contained PSIS-style calculation:

- raw leave-one-out importance ratios;
- generalized-Pareto method-of-moments tail-shape estimation;
- finite-sample weight truncation; and
- pointwise ELPD aggregation.

It produces the same concepts (`elpd_loo`, `p_loo`, `looic`, pointwise values,
and Pareto-k diagnostics), but can differ numerically from `loo::loo`.

## Geometry

The Euclidean static HMC path follows the original explicit leapfrog. The
position-dependent path follows the generalized implicit leapfrog with
fixed-point iterations. SoftAbs masses use an eigendecomposition of the
negative log-density Hessian and finite-difference metric derivatives.

No `ou_geom_bridge` is needed: Fortran callers directly attach procedure
pointers to an `ou_geom_target_type`.

## Multiple imputation

Rubin's rules, degrees of freedom, fraction of missing information, and t-based
intervals are translated. The Fortran entry point accepts a `[T,S,D]` array.
Checkpoint files and resumed R object serialization are presentation/runtime
features and are omitted.

## Naming and array order

Fortran uses column-major arrays, matching R and Stan storage. Time is the first
index and sector is the second. Public procedure names use the package's
snake-case names where they are valid Fortran identifiers.

## Random numbers

The port uses an explicit portable xorshift-based RNG with Box-Muller normal,
Marsaglia-Tsang gamma, and normal/chi-square Student-t generation. It does not
reproduce R or Stan random streams.
