# Validation

The automated test suite covers:

- Parameter-name generation, packing, unpacking, and default domains
- Independent reference values for `A(T)` and the two-factor covariance matrix
- Normal quantiles and analytical futures forecasts
- Correlated futures simulation
- Kalman filtering with a deliberately missing contract observation
- Theoretical and empirical futures volatility term structures
- European and Longstaff-Schwartz American option valuation
- Contract-number and maturity-matched stitching
- End-to-end one-factor simulated-data maximum-likelihood estimation

The checked build uses strict Fortran 2018, explicit-interface errors, array and
bounds checking, backtraces, and floating-point traps. An optimized build runs
the same tests.

## Release build

The release was compiled with GNU Fortran 14.2.0 in both checked and optimized
configurations. The checked configuration includes
`-Werror=implicit-interface`, `-fcheck=all`, and floating-point traps for
invalid operations, division by zero, and overflow. The example program is
compiled and run in both configurations.
