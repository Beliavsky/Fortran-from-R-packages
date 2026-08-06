# Validation

The test suite is deterministic and covers:

1. AR/PACF round trips and stationary-polynomial checks.
2. Exact AR(1), MA(1), FDWN, FGN, PLA, and seasonal-AR autocovariances.
3. Equality of Durbin-Levinson and direct covariance-matrix likelihoods.
4. Exact covariance-recursion simulation with supplied innovations.
5. Ordinary and seasonal differencing/integration round trips.
6. Exact AR(1) and ARIMA(1,1,0) forecasts and forecast variances.
7. Numerical Fisher information and identifiability.
8. Transfer-function and psi-weight formulas.
9. Recovery of simulated AR(1) and fractional-d parameters.
10. Recovery of static-regression and transfer-function coefficients.

Both checked and optimized builds use all warnings as errors. The checked build
also enables bounds, allocation, floating-point exception, backtrace, and
uninitialized-real diagnostics.
