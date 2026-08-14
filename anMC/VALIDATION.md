# Validation

Compiler used for release validation:

`GNU Fortran 14.2.0`

Strict flags:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

## Permanent tests

- `test/test_core.f90`
  - multivariate-normal mean/covariance smoke checks;
  - covariance-input and upper-Cholesky-input sampling modes;
  - truncated standard-normal moment/bound checks;
  - deterministic active-dimension selection;
  - full-active-set `ProbaMax` and `ProbaMin` on the equicorrelation-0.5 family;
  - ordinary MC and ANMC remainder estimates on an independent conditional
    problem with exact probability 0.5;
  - conservative-estimate edge case where all marginal probabilities exceed
    `alpha`.
- `test/test_active.f90`
  - methods 1 through 5 return unique in-range active dimensions;
  - `select_q_dims` respects its q limits and returns conforming submatrices.
- `test/test_orthant.f90`
  - full bias-corrected MC and ANMC on a 20-dimensional zero-mean Gaussian with
    equicorrelation 0.5.  The exact value is `P(max X > 0)=20/21`.
- `test/test_conservative.f90`
  - independent Gaussian marginals with analytically known joint excursion
    probabilities; verifies selected set, level, and probability.

All permanent tests pass with runtime bounds/allocation checking enabled.

## Independent Gaussian-probability validation

Development differential tests generated random positive-definite correlation
matrices and compared the native Fortran rectangle-probability routine with
SciPy's multivariate-normal CDF at high integration accuracy.

- 40 one-sided rectangles, dimensions 2 through 6:
  maximum absolute difference `2.1021e-5`.
- 30 finite lower/upper boxes, dimensions 2 through 5:
  maximum absolute difference `8.0207e-6`.

There were zero failures under a tolerance of `max(0.004, 5*reported_error)`.

The separately supplied `mvtnorm-fortran` package was also run on the same
40 one-sided fixtures as an independent cross-check.  Its maximum discrepancy
from the SciPy reference was `4.4769e-6`; the maximum direct difference between
this port's native integrator and `mvtnorm-fortran` on those 40 cases was
`1.5996e-5` (mean `1.0985e-6`).  It is not linked into this release.

## End-to-end known-probability check

For `d=20`, zero means, unit variances and common correlation 0.5, symmetry gives

`P(max X > 0) = d/(d+1) = 0.95238095238...`.

With five equally spaced active dimensions in the strict optimized validation
run:

- ordinary MC result: about `0.948656`;
- ANMC result: about `0.950951`;
- exact value: `0.952381`.

Both are within the permanent test tolerance of 0.04.  Monte Carlo values can
change slightly with compiler/runtime timing because the upstream algorithm
allocates work according to a wall-clock budget.
