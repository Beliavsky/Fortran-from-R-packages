# Testing

The project ships four test programs:

- `test_core`: array helpers, matrix algebra, polynomial operations, and special
  functions.
- `test_calculus_optimization`: differentiation, quadrature, roots,
  unconstrained/constrained optimization, least squares, LP, and QP.
- `test_interpolation_ode_signal`: spline/PCHIP/Akima interpolation, ODE/PDE
  stepping, convolution/Fourier transforms, smoothing, and peak detection.
- `test_geometry_combinatorics`: planar/spherical geometry, circle fitting,
  number theory, combinatorics, random/structured matrices, and QP regression
  cases.

Two compiler configurations are validated:

```sh
./scripts/test_gfortran.sh
./scripts/test_optimized.sh
```

The strict configuration uses Fortran 2018 conformance, all common warnings as
errors, array bounds checking, and floating-point traps.  The optimized
configuration uses `-O3` with the same warning policy.  Each script builds and
runs all tests, all examples, and the demo.

The final release archive is extracted into a fresh directory and both scripts
are rerun from those exact contents.
