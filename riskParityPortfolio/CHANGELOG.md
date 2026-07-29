# Changelog

## 0.2.2-fortran.1 - 2026-07-27

- Initial modern Fortran 2018/FPM port.
- Added Spinu, Roncalli, Choi, Newton-Nesterov, diagonal, and active-risk solvers.
- Added all eight risk-concentration objectives, Jacobians, and gradients.
- Added SCA with expected-return and variance terms.
- Added box, equality, inequality, projection, and active-set QP support.
- Added examples, API documentation, provenance files, and strict tests.
- Corrected the nonuniform-budget gradient for the standard-deviation-scaled
  formulation by using `2*transpose(A)*g`.
