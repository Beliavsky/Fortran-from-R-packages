# Porting notes

## Design

The translation uses small focused modules and one umbrella module,
`derivmkts`. Public calculations use `real(dp)` with
`dp = kind(1.0d0)`, explicit interfaces, `implicit none`, allocatable arrays,
and typed results.

R list and data-frame results are represented by derived types. R function
objects used by `greeks` are represented by a `pricing_function` procedure
interface. A wrapper can capture extra parameters, such as a barrier level,
before being passed to `numerical_greeks`.

## Numerical support

The project includes native implementations of:

- normal PDF and CDF
- bivariate normal CDF through adaptive Simpson integration
- bisection root finding in implied-value and yield routines
- Cholesky decomposition
- deterministic uniform, normal, Poisson, and binomial random generation
- means, sample variances, covariance, and interpolation

No external statistics or linear-algebra library is required.

## Intentional changes and corrections

1. Implied volatility is constrained to be positive, and the upper bracket
   expands automatically when required. The upstream default upper bound of
   `1e6` is accepted but is not evaluated unless necessary.
2. Black-Scholes time theta is returned with the same convention as the
   upstream numerical `greeks` routine: minus the derivative with respect to
   time to maturity, divided by 365.
3. Binomial exercise flags mark strictly optimal positive intrinsic exercise.
   The upstream equality test can label zero-value terminal nodes as exercise
   nodes because both continuation and intrinsic value are zero.
4. Monte Carlo routines use explicit deterministic seeds instead of changing
   and restoring R's global random state.
5. Multi-asset simulation takes an explicit covariance matrix and parameter
   vectors, avoiding R recycling rules.
6. Compound-option bivariate normal probabilities are calculated internally
   rather than through the R package `mnormt`.

## Compatibility details

The following upstream conventions are retained:

- Merton jump pricing uses `lambda * exp(alphaj)` as the Poisson-mixture
  intensity and the original adjusted-rate formula.
- Asian Monte Carlo `sd` outputs are payoff standard deviations, not Monte
  Carlo standard errors.
- The arithmetic Asian control variate uses 250 pilot observations and then
  applies the estimated coefficient to the complete simulated sample.
- Binomial theta follows the package's local two-step lattice approximation
  and is reported per day.

## Plotting

`binomplot` and the animated `quincunx` display are not compiled. Equivalent
numerical inputs for external plotting are returned in `binomial_result` and
`quincunx_result`.
