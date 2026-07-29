# Independent reference generation

Fixed numerical references in the tests were generated independently from the
Fortran implementation using small Python/NumPy translations of the published
and source formulas.

The reference calculations cover:

- the fractional binomial-weight convolution;
- the Haslett-Raftery recursions for `amk`, `ak`, `vk`, `phi`, and the tail
  approximation after `M` terms;
- the exact `fdsim` fractional-noise and inverse-ARMA recursion;
- autocovariance periodograms, OLS slopes, and standard-error expressions for
  GPH and Sperio.

The ARMA Jacobian is not validated against copied constants: each column is
compared at runtime with a central finite-difference derivative of the residual
recursion. Model-fit tests use newly generated series and verify statistical
recovery and matrix identities rather than a single optimizer path.
