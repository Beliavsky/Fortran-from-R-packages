# Testing

The permanent regression suite covers:

1. Weighted means, unbiased variances, and standardization.
2. Compact SVD, pseudoinverse identities, matrix powers, rank/condition,
   positive-definiteness checks, and Higham repair.
3. Symmetric packing/indexing and covariance/precision decomposition.
4. Independently calculated shrinkage intensities and covariance values.
5. Partial-correlation round trips and partial-variance estimates.
6. A `p > n` covariance/precision problem and low-rank matrix-power product.

`build_checked` enables bounds, allocation, and floating-point runtime checks.
`build_optimized` compiles the same source and executes the same tests at `-O3`.
