# Porting notes

## Scope

`mev` 2.2 exports 229 R-level names and contains roughly 45,000 lines of R and
C++ code. Version 0.3.0 continues the numerical translation begun in v0.1.0.
Plotting, S3 print/summary methods, formula/model-frame handling, vignettes and
packaged-data presentation helpers remain intentionally outside the Fortran API.

## Version 0.3 additions

- GEV Cox-Snell bias, implicit bias correction and Firth-score correction.
- BAB Monte Carlo testing and simultaneous confidence envelope.
- GPD/GEV tangent-exponential-model (TEM) profile corrections.
- Multivariate generalized-Pareto likelihoods (`mgp.ll` numerical layer) for
  logistic, negative-logistic, Brown-Resnick and extremal-t models.
- Censored multivariate generalized-Pareto likelihoods (`mgp.cll` numerical
  layer) for the same four model families.
- Conditional Gaussian and Student-t probabilities for censored likelihoods via
  the native deterministic QMC routines added in v0.2.
- Extremal-t and covariance-form Brown-Resnick spectral simulation.
- `rmev` support for those two additional spectral models.
- General `rparp` and `rgparp` rejection sampling for the translated spectral
  models and common risk functionals.

## GEV Cox-Snell implementation

The upstream package contains a very large generated symbolic expression for
the GEV Cox-Snell bias. The Fortran port instead evaluates the defining
Cox-Snell cumulants deterministically: Gauss-Legendre quadrature is applied to
an exact GEV quantile transform, analytic score/information routines are reused,
and the required parameter derivatives are evaluated by centered differences.
This is the same first-order Cox-Snell correction mathematically, but not a
line-for-line transcription of the generated algebra. The code uses an explicit
near-zero-shape interpolation and is tested for location/scale equivariance,
shape invariance and 1/n scaling. Shapes extremely close to the regularity
boundary xi=-1/3 remain numerically more delicate than central parameter values.

## Censored multivariate likelihoods

The upstream Brown-Resnick/extremal-t censored likelihoods call `mvPot` for
quasi-Monte-Carlo multivariate Gaussian/t probabilities. Version 0.3 uses the
independent deterministic Halton/antithetic QMC implementation already present
in `mev_mgp`; therefore `mvPot` is no longer required for these likelihoods.
The results are deterministic for a fixed `nqmc`, but high-dimensional values
are approximate rather than bit-for-bit identical to `mvPot`.

The translated censored likelihoods support inferred censoring from `mthresh`
or an explicit logical censoring matrix, and mgp, Poisson and binomial
likelihood contributions.

## Upstream multivariate-likelihood corrections

Several source-level inconsistencies become visible when uncensored and
censored forms are tested against one another. Version 0.3 corrects them rather
than reproducing internally inconsistent formulas:

1. `intensBR` uses `+Lambda[1,-1]` although the conditional Brown-Resnick mean
   and the code comments require `+2*Lambda[1,-1]`.
2. `intensBR` uses `log(pi)` in the multivariate-normal normalizing constant;
   the standard Gaussian density and the censored implementation require
   `log(2*pi)`.
3. The uncensored negative-logistic mixed-derivative expression has gamma and
   logarithmic terms that do not reduce to the full mixed derivative used by
   the censored likelihood. The Fortran formulas are generated from the mixed
   derivative itself and censored/uncensored paths agree when no component is
   censored.
4. The uncensored Poisson branch of `mgp.ll` has `+ntot*exponentMeasure`, while
   `mgp.cll` and the Poisson-process likelihood require the negative integrated
   intensity. The Fortran port uses the negative sign in both paths.

Regression tests explicitly exercise the censored-to-uncensored reduction.

## Simulation choices

`rparp` and `rgparp` use direct rejection from the translated spectral measure.
This is exact for the represented proposal/conditioning construction but can be
less efficient than upstream specialized site samplers or composition samplers.
The latter use model-specific extremal-function/truncated-normal machinery and
are therefore retained as optional future performance/parity work rather than a
correctness dependency.

## License handling for supplied dependency translations

The supplied dependency archives were inspected before integration.
`nleqslv-fortran` is GPL-2.0-or-later and `expint-fortran` is GPL-3.0-or-later;
both can be combined with this GPL-3.0-only project. Their notices/licenses are
retained under `vendor/`.

The supplied `Rsolnp-fortran` and `mvtnorm-fortran` translations identify
themselves as GPL-2.0-only. GPL-2.0-only code cannot be distributed as part of
a GPL-3.0-only combined work, so those two translations are not vendored,
linked or copied into this archive. The affected `mev` routines use independent
native algorithms: transformed/constrained optimization and deterministic QMC.

## Numerical implementation choices

- GEV/GPD retain explicit near-zero shape handling.
- MLE/profile/constrained fits use deterministic transformed pattern
  optimization unless a compatible translated solver is explicitly used.
- GPD and GEV nonlinear bias-correction equations use `nleqslv-fortran`.
- BAB/lower-trimmed-Hill calculations use `expint-fortran`.
- MV-normal/MV-t upper probabilities use deterministic low-discrepancy QMC.
- Scores/information for very large symbolic upstream likelihood expressions
  may use stable finite differences of the translated likelihood.
- `rmev` implements upstream Algorithm 1 using translated spectral samplers.

## Remaining lower-priority gaps

The principal standalone numerical workflows are now covered. Remaining gaps
are narrower or dependency-heavy:

- specialized `rparpcs` composition samplers and model-specific direct
  extremal-function/site samplers, mainly performance alternatives to the new
  general rejection simulation;
- some higher-order profile/TEM variants for every reparameterized return-level,
  expected-shortfall and N-block quantity;
- optional GMM and Bayesian paths depending on `gmm`/`revdbayes`;
- less commonly used pairwise-beta/pairwise-exponential and weighted
  Ballani-Schlather spectral families;
- R-specific plotting, S3, formula/model-frame and presentation infrastructure.

A translation of `TruncatedNormal` would still be useful if exact parity with
the specialized Brown-Resnick/extremal-t composition samplers is desired, but
it is no longer required for censored likelihood inference or general R-Pareto
simulation.

## Provenance

Original R/C++ source used for translation is retained under `orig/`. Upstream
metadata identifies `mev` version 2.2 and `License: GPL-3`.
