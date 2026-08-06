# Porting notes

## Matrix conventions

Matrices use ordinary Fortran column-major storage. Coefficients and formulas
follow the upstream R implementation. Covariance matrices must be finite,
symmetric, and have positive diagonal elements for minimum-torsion methods.

## Exact minimum torsion

The iteration follows the upstream sequence:

1. Form the correlation matrix and its symmetric square root.
2. Initialize the diagonal scaling to one.
3. Compute the polar factor from the scaled correlation root.
4. Update the diagonal from `diag(Q C^(1/2))`.
5. Stop on the upstream normalized relative objective change.

The result decorrelates the covariance matrix. Unlike the approximate method,
the exact method need not preserve each original marginal variance after the
torsion transform.

## Effective bets

The diversification probabilities are evaluated as

```text
solve(transpose(T), b) * (T * Sigma * b) / (b' * Sigma * b)
```

where the multiplication between the two vectors is elementwise. The entropy
threshold `p > 1e-5` is preserved from the R source.

## Optimization

The R function calls `NlcOptim::solnl`. Statically including a translation of
that GPL solver would change the licensing of this MIT package. The Fortran
port therefore uses an independently written spectral projected-gradient
method on the probability simplex. It uses the analytical ENB gradient,
Armijo backtracking, and Barzilai-Borwein step estimates. Numerical Hessian and
KKT multiplier estimates are returned for compatibility diagnostics.

Results should agree on the optimum but iteration counts and multipliers are
not expected to match `NlcOptim` exactly.
