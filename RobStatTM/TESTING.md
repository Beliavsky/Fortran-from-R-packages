# Testing

The test suite uses deterministic contaminated datasets and checks allocation, dimensions, finite outputs, probability ranges, positive covariance diagonals, robust resistance to outliers, PCA orthonormality, DCML mixing bounds, RFPE decomposition, and nested-model degrees of freedom.

Test programs:

- `test_psi_location.f90`
- `test_regression.f90`
- `test_logistic.f90`
- `test_multivariate.f90`
- `test_pca.f90`

Run with:

```text
./scripts/build_checked.sh
./scripts/build_optimized.sh
```

The checked script enables `-std=f2018`, `-Wall`, `-Wextra`, `-Werror`, `-Wimplicit-interface`, `-fimplicit-none`, `-fcheck=all`, and backtraces.
