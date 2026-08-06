# Validation

The test suite checks:

1. Symmetric matrix square-root reconstruction.
2. PCA orthogonality and covariance diagonalization.
3. Approximate and exact minimum-torsion decorrelation.
4. Diversification probabilities summing to one and valid ENB bounds.
5. Analytical ENB gradients against central finite differences.
6. Long-only fully invested ENB maximization to the theoretical dimension.
7. Invalid shapes, indefinite matrices, and dimension mismatches.

The checked build uses warnings as errors, bounds/array checking, backtraces,
and floating-point traps for invalid, divide-by-zero, and overflow operations.
The optimized build uses the same warning policy at `-O3`.
