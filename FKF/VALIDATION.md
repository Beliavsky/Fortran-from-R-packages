# Validation

The test suite covers:

1. Hand-calculated scalar local-level filter states, covariances, and likelihood.
2. Partial and fully missing multivariate observations, including NaN output layout.
3. Deterministic time-varying intercepts, transition/loadings, and covariances.
4. Durbin-Koopman smoothed states/covariances versus an independently coded RTS smoother.
5. Source-compatible versus corrected missing-data likelihood constants and singular-F status handling.

Run `make checked` and `make optimized` to reproduce the validation.
