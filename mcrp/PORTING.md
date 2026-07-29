# Porting notes

## Scope

The R package exports 22 numerical routines. Every one is represented in this
Fortran translation. No computational routine was omitted.

The FRAPO dependency is used by the R package only to construct an S4 `PortSol`
result. It is not required for the numerical algorithm. The Fortran port uses
`type(mcrp_result)` instead.

## Exact source behavior retained

### Different divisors by moment order

The source uses:

- `T-1` for the second centered moment (`M2`), because it calls the sample
  covariance convention.
- `T` for third and fourth centered moments (`M3`, `M4`).

This asymmetry is preserved.

### Tensor ordering

The columns of `M3` and `M4` follow the exact nested R `kronecker` ordering.
This matters when users contract the matrices manually.

### Standardized skewness and kurtosis decomposition formulas

The source routines named `PortSkewDeriv` and `PortKurtDeriv` are not the
ordinary calculus gradients of standardized skewness and kurtosis. Their
coefficients are chosen so that Euler-style weighted contributions sum to the
portfolio statistic:

- `sum(w * PortSkewDeriv) = PortSkew`
- `sum(w * PortKurtDeriv) = PortKurt`

The same formulas are used inside the optimization objective. They are
preserved exactly rather than silently replaced by conventional gradients.

### Variance criterion scaling

Inside `mcrp`, the variance contribution vector is divided by `2*variance`
before its sample variance is calculated. Consequently, the variance part of
the objective is scaled by one quarter relative to using ordinary percentage
variance contributions. This affects the interpretation of `lambda` and is
preserved.

### Raw-parameter bounds and final normalization

The source optimizes raw parameters and only afterward computes:

```text
w = x / sum(abs(x))
```

Bounds therefore constrain raw parameters, not final normalized weights. The
Fortran port preserves this behavior. With nonnegative raw parameters the final
weights sum to one; unrestricted parameters produce an L1-normalized long-short
portfolio.

## Numerical substitutions

### `stats::nlminb`

The R package delegates optimization to `nlminb`, whose implementation is not
part of mcrp. The Fortran library uses a self-contained bounded Nelder-Mead
method. It supports lower and upper bounds and preserves the exact objective,
but it need not follow the same iterations or select the same local minimum.

The scale-invariant radial direction can make any optimizer poorly identified.
The returned normalized weights are generally more meaningful than the raw
parameter scale.

### Efficient objective contractions

The public `M3` and `M4` routines materialize the same flattened tensors as R.
During optimization, the Fortran implementation evaluates equivalent moment
contractions directly from centered observations. This avoids allocating
`N x N^2` and `N x N^3` arrays at every fit and is algebraically identical.

## R-only infrastructure omitted

- S4 `PortSol` construction and methods
- R call capture
- column and row names
- `timeSeries` metadata
- bundled FRAPO demonstration data
- formatted printing and plotting

The original package tree is retained under `original/mcrp-master`.
