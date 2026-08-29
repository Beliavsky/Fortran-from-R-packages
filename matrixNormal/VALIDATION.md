# Validation

Validation was performed with GNU Fortran 14.2.0 using Fortran 2018, runtime checking, and implicit external interfaces promoted to errors.

Compiler flags used for the strict validation build included:

```
-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0
```

## Retained tests

* `test_utils` -- trace, `vec`, corrected `vech`, identity/ones matrices, Kronecker product, symmetry, PSD and PD checks.
* `test_density` -- reduction to a one-dimensional Normal and an independent SciPy 4-D multivariate-Normal log-density reference.
* `test_probability` -- independent diagonal-covariance rectangle probability, the documented `V kron U` ordering, legacy `U kron V` ordering, upper-only CDF, and full-support probability.
* `test_random` -- seeded one/many-draw reproducibility plus empirical means, variances, row covariance, and column covariance.

## Independent reference values

For

```
A = [[ 0.4, -0.7],
     [ 1.2,  0.3]]
M = [[ 0.1, -0.2],
     [ 0.5,  0.4]]
U = [[2.0, 0.3],
     [0.3, 1.5]]
V = [[1.2, 0.2],
     [0.2, 0.8]]
```

SciPy's multivariate Normal evaluated on column-major `vec(A)` with covariance `V kron U` gives

```
log density = -4.908844744462413
```

and the Fortran test reproduces this to approximately `2e-12` or better.

For the diagonal probability regression case, direct products of one-dimensional Normal probabilities give

```
V kron U : 0.006074900306399315
U kron V : 0.0062069338267269954
```

Both the corrected default and explicit legacy branch are tested.
