# Porting notes

## Marginal estimation

The R implementation delegates to `rugarch::ugarchroll`, which performs joint
model estimation and returns rich S4 objects. The supplied Fortran `rugarch`
port has a numerical GARCH fitter and ARMA utilities but no exact counterpart
to every `ugarchroll` class path. The translation therefore:

1. estimates the ARMA mean equation by iterated conditional least squares;
2. fits the chosen GARCH recursion to the resulting innovations;
3. combines the fitted mean and volatility parameters;
4. filters the last vine-training segment and the following forecast segment
   with fixed parameters;
5. applies the fitted innovation CDF to standardized residuals.

This preserves the rolling information flow and the main model equations. It
will not produce coefficients identical to R's hybrid solver.

## Vine estimation

The native vine dependency implements parametric C-vines and D-vines. The
D-vine path closely follows the package and supports conditional simulation.
The package's unrestricted regular-vine selection is represented by the C-vine
path because arbitrary R-vine matrices and Dissmann tree selection are not in
the current dependency port.

The D-vine ordering transforms each series to normal scores and greedily adds
the candidate with the largest sum of absolute partial correlations across the
requested depth. One or two conditioning variables are fixed at the beginning
of the order.

## Conditional simulation

For one conditioning variable, its raw copula value is fixed as the first
Rosenblatt coordinate. For two variables:

- quantile mode fixes the first unconditional quantile and the second
  conditional quantile;
- raw mode converts the second raw copula value to its conditional CDF before
  inverse Rosenblatt simulation.

This is algebraically equivalent to the auxiliary-matrix logic in the Rcpp
implementation for a D-vine ordered with conditioning variables first.

## Risk measures

`empirical_quantile()` uses R's default type-7 interpolation. Mean and median ES
use all observations less than or equal to the corresponding VaR. Monte Carlo
ES averages empirical quantiles at probabilities sampled uniformly between
zero and alpha.

When conditional estimation is requested, `overall` risk uses the combined
sample from all requested stress levels plus the observed/prior residual case,
matching the R implementation.

## Numerical and language differences

- Arrays are variables by observations.
- Procedures are serial and deterministic except for simulation.
- Call `rugarch::seed_rng` or the intrinsic `random_seed` before simulations for
  reproducible results.
- No missing-value semantics or R vector recycling are provided.
- Status integers and `result%message` replace R exceptions and warnings.
