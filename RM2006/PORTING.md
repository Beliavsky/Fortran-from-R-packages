# Porting notes

## Source mapping

The original package contains one exported R function:

```text
RM2006(data, tau0=1560, tau1=4, kmax=14, rho=sqrt(2))
```

It maps to the Fortran generic `rm2006` and its implementation procedure
`rm2006_covariance`.

The translation preserves the following algorithm exactly for normal valid
inputs:

1. Convert every return row into an outer-product matrix.
2. Construct geometrically spaced time scales.
3. Construct and normalize the RiskMetrics 2006 scale weights.
4. Build a finite exponentially weighted backcast for each scale.
5. Run one EWMA covariance recursion per scale.
6. Combine those recursions using the multiscale weights.

The Fortran implementation computes each outer product as needed rather than
materializing the complete `K x K x T` temporary array used by the R code.
This reduces peak memory while leaving the mathematical result unchanged.

## Small-sample edge case

The R expression for the backcast endpoint is:

```text
max(min(floor(log(0.01)/log(mu)), T), k)
```

When `k > T`, this can produce an endpoint greater than the number of
observations and lead to out-of-range indexing. The Fortran version caps the
minimum scale-dependent endpoint at `T`:

```text
max(min(raw_endpoint, T), min(k, T))
```

This agrees with the R result whenever `T >= kmax`, including the normal
intended use, while allowing short samples to be processed safely.

## Error handling

R generally reports errors through indexing or arithmetic failures. The
Fortran API validates dimensions, parameters, finite input data, and scale
weight normalization. It returns status codes instead of stopping the calling
program.

## Omitted R infrastructure

There are no plotting, class, dataset, or compiled-extension components in the
original package. The `.Rd` help pages, package metadata, and original R source
are retained under `original/` but are not part of the compiled library.

## Numerical ordering

The output keeps the original ordering `(K,K,T+1)`. This differs from the input
ordering `(T,K)`, so callers should use the last array index for time.
