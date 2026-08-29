# Validation

The retained test suite exercises the translated algorithms rather than only
checking that the project compiles.

- `test_core`: hidden-layer squared-error network.  The analytic gradient is
  compared with central finite differences of the objective and the translated
  analytic Hessian with finite differences of the gradient.
- `test_modes`: exact binary-logistic prediction identity, softmax normalization,
  censored-softmax likelihood identity, and finite-difference Hessian checks for
  entropy, ordinary softmax, and censored softmax modes.
- `test_fit`: deterministic zero-hidden linear-output fit recovers intercept
  1.25 and slope 2.5 to numerical precision.
- `test_multinom`: three-class label fit, prediction, Fisher information, and
  generalized-inverse covariance checks.
- `test_multinom_counts`: response-count fitting, offsets, censored fitting, and
  an independent SciPy baseline-category multinomial BFGS reference.  The
  reference minimum is 61.44892130753712 with class-2 coefficients
  `(0.53347279, 0.40661046)` and class-3 coefficients
  `(0.41497049, 1.08252913)`.
- `test_utils`: one-hot encoding and duplicate-design-row consolidation.

The exact supplied `r_mod.f90` and the build copy are also compared after
removing whitespace and Fortran continuation markers; they are identical under
that formatting-only normalization.
