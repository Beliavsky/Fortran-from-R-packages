# Porting notes

## Numerical design

The upstream local PCA diagonalizes a `T x T` matrix `X_r X_r'`. The Fortran implementation diagonalizes the smaller of `X_r X_r'` and `X_r' X_r`, then reconstructs the corresponding singular vectors. This is algebraically equivalent for nonzero singular values and is substantially faster when `p < T`.

The package is self-contained and does not require BLAS or LAPACK. Symmetric eigenproblems use Jacobi rotations; linear systems use pivoted Gaussian elimination; positive-definite repair floors eigenvalues exactly as the R code does.

## Expected-return forecasting

The R package fits `(0,0,0)`, `(1,0,0)`, `(0,0,1)`, and `(1,0,1)` models through `stats::arima`, chooses the smallest AIC, and averages the horizon forecasts. The Fortran port keeps those four candidates and AIC selection but estimates them by bounded conditional sum of squares with a deterministic pattern search. Results will not be bit-for-bit identical to R's exact likelihood optimizer.

## Upstream compatibility modes

### Boundary normalization

The R call `boundary_kernel(r, t, ...)` passes the local-PCA target as the first argument, but the function's boundary branch is based on its second argument. The default Fortran mode reproduces this. Set `source_compatible_boundary=.false.` to base the normalization on the target index.

### Expanding-window sample refresh

In the upstream R6 method, `est_data` and `bandwidth` are assigned inside the yearly factor-update branch. They therefore remain stale between factor updates. The default Fortran mode reproduces that behavior. Set `source_compatible_expanding=.false.` to use all observations available at every rebalance and recompute the bandwidth each time.

## R features omitted

R6 mutable objects are replaced by immutable/allocatable derived-type results. Plotting, progress reporting, `cli`, `prettyunits`, `dplyr`, `tidyr`, `ggplot2`, and R object-size helpers are not computational and are omitted.
