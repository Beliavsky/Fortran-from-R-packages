# Validation

The source tree was compiled with GNU Fortran 14.2.0 using:

```text
gfortran -std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all
```

The following test programs are supplied:

- `test_core`: normal quantile identity, finite/infinite Romberg integration,
  TOMS614 finite integration, two-dimensional integration, Gauss-Hermite
  moments, RK4, matrix-exponential ODE solution, `gettvc`, and contrasts.
- `test_distributions`: independent numerical reference values for inverse
  Gaussian, Laplace, Pareto, generalized gamma, power exponential,
  beta-binomial and gamma-count; normalization checks for overdispersed count
  families; GIG/simplex/quantile identities; PVF-CDF endpoint behavior.
- `test_roundtrip`: CDF/quantile consistency across every continuous family and
  smoke tests for all discrete quantiles.
- `test_pkpd`: representative one- and two-response PK/PD mean functions.

A separate TOMS614 improper-integral check evaluated
`integral_0^infinity exp(-x) dx` as approximately
`0.999999999987776` at requested tolerance `1e-9`.

Reference values used by `test_distributions` were independently computed with
SciPy/Python where corresponding standard distributions/special functions were
available.

All `.f90` source/test/example lines are checked to remain within the standard
132-column free-form limit.
