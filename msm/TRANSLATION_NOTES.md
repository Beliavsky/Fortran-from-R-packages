# Translation notes

## Source mapping

| Original msm code | Fortran translation |
| --- | --- |
| `src/pijt.c` | `src/msm_linalg.f90`, `src/msm_ctmc.f90` |
| `src/analyticp.c` | numerically consolidated into `transition_matrix` |
| `src/hmm.c` | `src/msm_emissions.f90` |
| `src/hmmderiv.c` | `emission_derivative` in `src/msm_emissions.f90` |
| `src/lik.c` ordinary branch | `ctmc_minus2loglik`, aggregate likelihood and gradients |
| `src/lik.c` censor branch | `ctmc_censored_minus2loglik` |
| `src/lik.c` hidden/Viterbi branch | `src/msm_hmm.f90` |
| `R/simul.R` | `src/msm_simulation.f90` |
| `R/pexp.R` | piecewise-exponential routines in `msm_distributions.f90` |
| `R/phase.R` | two-phase routines in `msm_distributions.f90` |
| `R/tnorm.R` | truncated-normal routines in `msm_distributions.f90` |
| `R/totlos.R`, `R/efpt.R`, `R/ppass.R` | summary routines in `msm_ctmc.f90` |
| `R/deltamethod.R` | `src/msm_inference.f90` |

## Array conventions

The original C source emulates R's column-major arrays with `MI`, `MI3`, and
`MI4` macros.  Fortran is naturally column-major, so those index calculations
mostly disappear.  Public matrices use conventional Fortran `(row,column)`
indices.  State numbers are 1-based, matching the public R package API.

## Matrix exponential and derivatives

The old package contains several matrix-exponential implementations and many
closed-form special cases.  The Fortran translation uses one Pade
scaling-and-squaring implementation.  Transition-matrix derivatives are
computed as Frechet derivatives using

`exp([[A,E],[0,A]])(1:n,n+1:2n)`.

This replaces the original eigensystem / repeated-eigenvalue power-series
branch while preserving the derivative mathematically.

## Likelihood scaling

The hidden and censored forward recursions normalize at every observation and
accumulate the logarithm of the scaling constant.  `src/lik.c` normalizes by
the mean of the current vector rather than its sum.  Both give the identical
final log likelihood; sum normalization is used here because it makes the
filtered vectors probability distributions directly.

## HMM parameterization

R's hidden-model objects contain model IDs, parameter offsets, covariate
transformations, and derivative maps constructed from formulas.  Fortran
represents each outcome/state distribution explicitly as an `emission_model`
with a `kind` and parameter vector.  This removes the R metaprogramming layer
while preserving the distribution kernels.

For categorical models, `pars` contains the category probabilities directly,
without the legacy metadata slots used inside the C interface.

## Optimization

`msm` delegates numerical optimization to R's `optim` and optionally `minqa`.
These are external algorithms and are not copied.  The Fortran library returns
objective values and analytic Q gradients that can be passed to any FPM
optimization package.

## Regression philosophy

Where the original R runtime was unavailable, tests use independent closed
forms, explicit enumeration of hidden paths, and centered finite differences.
This is stronger for the translated kernels than comparing two implementations
which share the same internal approximation.

## First-passage systems

Hitting-probability equations are solved only on states with a directed path to the target set.  This avoids singular subgenerator systems when a model contains unrelated absorbing or recurrent classes.  `expected_first_passage` reports `HUGE()` when eventual passage is not certain, matching the mathematical fact that the unconditional first-passage time is infinite in that case.
