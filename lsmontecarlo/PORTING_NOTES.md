# Porting notes

## Algorithmic correspondence

The original R code stores a matrix of simulated paths, works backward through
exercise dates, regresses discounted future cash flow on polynomial state
variables, and retains the first exercised payoff. The Fortran translation
implements the same Longstaff-Schwartz method directly with one cash-flow value
and one exercise index per path. This avoids large temporary matrices while
making discount timing explicit.

The regression bases preserve the original package choices:

- plain American put: `1`, `S`, `S^2`
- Asian American put: `1`, `S`, `S^2`, with exercise payoff based on the
  running arithmetic average
- quanto American put: `1`, `S`, `S^2`, `G`, `G^2`, and `S*G`

## Numerical robustness changes

1. Regressions use scaled column-pivoted QR rather than R `lm` calls wrapped in
   `try`. Dependent columns are removed by rank detection. This gives the
   correct reduced regression for deterministic or nearly deterministic state
   variables.
2. The control variate is applied to each path before computing the estimate
   and standard error. Its expectation is the same as the original aggregate
   correction.
3. Antithetic standard errors are calculated from antithetic pair averages,
   not by treating the two members of each pair as independent observations.
4. Black-Scholes routines explicitly support zero maturity and zero
   volatility.
5. All stochastic public routines accept an optional deterministic seed.
6. Input dimensions and parameter domains are checked and produce clear
   `error stop` messages.

## Original-source issues clarified

- `AmerPutLSM_AV` and `QuantoAmerPutLSM_AV` place `n` twice in their returned R
  lists. `option_result` instead records `n_paths` as the number of original
  paths and `effective_paths` as the total number simulated.
- `EuCallBS` calculates both call and put values but returns only the call. The
  Fortran call and put functions compute only the requested result.
- The original `firstValueRow` relies on cumulative sums and therefore assumes
  nonnegative payoff matrices. The translated routine explicitly retains the
  first strictly positive value in each row, matching its use in the package.
- The original surface plot is presentation code. The Fortran surface routines
  return data and leave plotting to the caller.

## Quanto convention

The translation preserves the package's convention: the first asset follows
`rate - dividend`, the multiplier asset follows `rate2 - dividend2`, and option
cash flows are discounted at `rate`. This is a direct computational port, not a
change to a separately derived quanto pricing measure.

## Memory use

Path matrices require approximately `8*n_paths*n_steps` bytes per asset in
double precision. Antithetic and quanto calculations multiply that amount by
two as expected. For very large simulations, users should choose path and step
counts with available memory in mind.
