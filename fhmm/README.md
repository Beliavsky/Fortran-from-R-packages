# fhmm-fortran

A self-contained modern Fortran/FPM translation of the numerical algorithms in
[fHMM 1.4.3](original/fHMM-1.4.3), an R package for fitting ordinary and
hierarchical hidden Markov models to financial time series.

The project preserves the upstream GPL-3 license and attribution. It does not
require R, Rcpp, Armadillo, BLAS, LAPACK, or another statistics package.

## Main capabilities

- HMM and hierarchical-HMM log likelihoods in stable log form.
- Scaled forward filtering, prediction, backward smoothing, and state
  probabilities.
- Global Viterbi decoding for ordinary and hierarchical models.
- Normal, lognormal, Student-t, gamma, and Poisson state-dependent
  distributions.
- Native PDFs, CDFs, quantiles, and random generation for all five families.
- Constrained/unconstrained transition, mean, scale, and degrees-of-freedom
  parameter transformations.
- Ordinary and hierarchical simulation with deterministic optional seeds.
- Multi-start Nelder-Mead maximum-likelihood fitting.
- Numerical gradients, Hessians, covariance estimates, standard errors, AIC,
  and BIC.
- Quantile-based initial values.
- Pseudo-residuals, state reordering, model comparison, and forecasts.
- Fixed or randomly generated hierarchical chunk lengths.

## Basic example

```fortran
use fhmm
implicit none

type(hmm_parameters) :: par
type(inference_result) :: inf
real(dp) :: y(8)

par%distribution = dist_normal
allocate(par%gamma(2,2), par%mu(2), par%sigma(2), par%df(2))
par%gamma = reshape([0.95_dp, 0.10_dp, 0.05_dp, 0.90_dp], [2,2])
par%mu = [-0.5_dp, 0.8_dp]
par%sigma = [0.7_dp, 1.2_dp]
par%df = 10.0_dp

y = [-0.8_dp, -0.4_dp, -0.6_dp, 0.2_dp, &
      1.1_dp,  0.7_dp,  1.4_dp, 0.4_dp]

inf = forward_backward(y, par)
if (.not. inf%ok) error stop trim(inf%message)

print *, inf%log_likelihood
print *, inf%filtered(:, size(y))
```

## Building with FPM

```text
fpm build
fpm test
fpm run fhmm_demo
fpm run --example basic_hmm
fpm run --example hierarchical_hmm
```

The direct validation scripts are useful when FPM is unavailable:

```text
scripts/validate.sh
scripts/validate.bat
```

## Array conventions

- Transition matrices are row stochastic: `gamma(i,j)` is the probability of
  moving from state `i` to state `j`.
- Observation vectors have length `T`.
- Probability arrays are shaped `(states,T)`.
- In hierarchical data, `fine_observations(t,j)` stores fine observation `j`
  in coarse period `t`, and `chunk_lengths(t)` gives the active row length.
- Hierarchical fine models are stored in `parameters%fine(coarse_state)`.

## Forecast quantiles

The upstream R package forecasts a mixture by averaging component quantiles.
That is not generally the quantile of the mixture. `forecast_hmm` therefore
uses the actual mixture CDF by default. Set `upstream_quantiles=.true.` to
reproduce the upstream weighted-component calculation.

## Thread safety

Filtering, simulation, and distribution routines are reentrant. The fitting
routines hold objective data in module-level context while Nelder-Mead and the
finite-difference routines execute. Concurrent fitting calls should therefore
be externally serialized.

See [COVERAGE.md](COVERAGE.md), [PORTING_NOTES.md](PORTING_NOTES.md), and
[VALIDATION.md](VALIDATION.md) for details.
