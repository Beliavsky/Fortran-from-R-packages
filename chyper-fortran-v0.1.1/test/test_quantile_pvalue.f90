! SPDX-License-Identifier: MIT
program test_quantile_pvalue
  use chyper, only : dp, qchyper, pvalchyper, pchyper
  implicit none
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]

  if (qchyper(0.0_dp, 10, n, m) /= 0) error stop 1
  if (qchyper(0.90_dp, 10, n, m) /= 1) error stop 2
  if (qchyper(0.99_dp, 10, n, m) /= 2) error stop 3
  if (qchyper(1.0_dp, 10, n, m) /= 7) error stop 4
  ! Keep this point away from a CDF jump.  The previous v0.1.0 test used
  ! 0.9478252815000214, only about 1.2e-16 below F(1), which made the
  ! expected result compiler-rounding-dependent under the upstream strict
  ! CDF > p quantile convention.
  if (qchyper(0.95_dp, 10, n, m) /= 2) error stop 5
  if (abs(pvalchyper(2, 10, n, m) - &
      (1.0_dp - pchyper(1, 10, n, m))) > 1.0e-15_dp) error stop 6
  if (abs(pvalchyper(2, 10, n, m, .false.) - &
      pchyper(2, 10, n, m)) > 1.0e-15_dp) error stop 7
  print *, 'test_quantile_pvalue: PASS'
end program test_quantile_pvalue
