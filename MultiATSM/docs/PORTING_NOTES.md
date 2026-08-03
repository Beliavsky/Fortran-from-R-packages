# Porting and numerical notes

## Orientation and indexing

R code frequently switches between time-by-variable and variable-by-time
representations. The Fortran API standardizes on variables-by-time. Country
blocks are contiguous, and all indexing is one-based.

## Linear algebra

Dense least squares, eigensystems, SVD, Cholesky factors, and linear solves use
BLAS/LAPACK. Small ridge or jitter terms are added only after a failed solve or
factorization. The project does not vendor a separate linear-algebra library.

## Covariance normalization

`fit_var` reports residual covariance divided by the number of transition
observations, matching the Gaussian likelihood convention used by MultiATSM.
Callers needing the unbiased sample covariance can rescale by `T/(T-k)` or use
a custom estimator.

## PCA signs

Eigenvector signs are not mathematically identified. The same level, slope,
and curvature sign conventions used by the upstream helper are applied to the
first three components.

## Affine recursions

For each period `h`, the port uses

```text
B(h) = d + B(h-1) K1Q
A(h) = r0 + A(h-1) - 0.5 B(h-1) Sigma B(h-1)'
```

and divides by maturity to obtain yield loadings. Inputs and outputs are in the
caller's chosen period and rate units; no implicit annualization is applied.

## Likelihood

The likelihood combines Gaussian measurement portfolios and Gaussian state
innovations, and includes the PCA Jacobian term. Measurement variances are
concentrated country by country. A positive ridge/jitter is used for nearly
singular covariance matrices.

## Optimization

The optimizers minimize scalar callback functions. To maximize a likelihood,
return its negative value. BFGS uses central numerical gradients and an Armijo
backtracking line search. Nelder-Mead is derivative-free. Neither routine
performs hidden scaling or model-specific list reconstruction.

## Stationarity

`stabilize_transition` scales a matrix uniformly only when its spectral radius
exceeds the requested bound. It preserves eigenvectors and relative
eigenvalues, but it is not byte-for-byte equivalent to the upstream smooth
Jordan auxiliary transformation.

## Bootstrap

Random resampling uses Fortran's intrinsic generator after deterministic seed
initialization. IID, Rademacher wild, and overlapping moving-block schemes are
available. `bootstrap_var` refits an unrestricted VAR(1); restricted GVAR/JLL
bootstrap loops can be assembled from the lower-level public routines.

## JLL restriction difference

The PI orthogonalization and feedback zero restrictions are directly ported.
For the innovation covariance, the current implementation masks a Cholesky
factor and rebuilds the covariance. The R package may optimize free restricted
entries conditionally. Users requiring exact replication of that step should
wrap the public likelihood and `nelder_mead_minimize` routines.

## Error reporting

R exceptions are replaced with integer `info` outputs. Zero denotes success;
negative values indicate invalid dimensions or arguments; positive values may
indicate a usable fallback or a boundary condition. Examples stop on nonzero
status for clarity.
