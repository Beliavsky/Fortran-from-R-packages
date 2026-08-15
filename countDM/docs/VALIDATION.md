# Validation

The release is compiled with GNU Fortran 14.2 using strict Fortran 2018,
implicit-interface errors, trampoline errors, and runtime checking.

Permanent tests:

- `test_core`: Touchard/Bell identities, Bell--Touchard PMF/CDF/QF/RNG,
  zero-inflated identities, and bundled datasets.
- `test_mle`: Bell closed-form MLE, Poisson/Borel analytic checks, and every
  numerical MLE model including parameter-domain and standard-error checks.

Independent SciPy/Nelder--Mead likelihood checks on the upstream stillbirth
sample give representative optima:

- Bell--Touchard: `(lambda, theta) = (1.19249644, 0.11077972)`
- ZIP: `(alpha, theta) = (0.72419095, 1.57835065)`
- ZIBell: `(alpha, lambda) = (0.61790144, 0.61558797)`
- ZIBell--Touchard: `(lambda, theta, pi) = (0.76270385, 0.61294635, 0.56569430)`
- ZOIP: `(alpha, beta, theta) = (0.76855200, 0.08891939, 2.43041718)`
- ZOIBell: `(alpha, beta, theta) = (0.78109453, 0.11940298, 0.76271954)`

The Fortran regression tests check these values with tolerances chosen well
above optimizer rounding noise but well below statistical differences.
