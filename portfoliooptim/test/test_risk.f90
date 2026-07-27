! SPDX-License-Identifier: GPL-3.0-only
program test_risk
  use portfoliooptim, only : dp, risk_result, risk_post, risk_code, &
    risk_cvar, risk_dcvar, risk_lsad, risk_mad
  implicit none
  real(dp) :: losses(4), probabilities(4)
  type(risk_result) :: result

  losses = [2.0_dp, 3.0_dp, 1.0_dp, 1.0_dp]
  probabilities = 0.25_dp
  result = risk_post(losses, probabilities, 0.90_dp)
  call assert_true(result%ok)
  call assert_close(result%var, 3.0_dp, 1.0e-13_dp)
  call assert_close(result%cvar, 3.0_dp, 1.0e-13_dp)
  call assert_close(result%mean, 1.75_dp, 1.0e-13_dp)
  call assert_close(result%mad, 0.75_dp, 1.0e-13_dp)
  call assert_true(risk_code('cvar') == risk_cvar)
  call assert_true(risk_code('DCVAR') == risk_dcvar)
  call assert_true(risk_code('lsad') == risk_lsad)
  call assert_true(risk_code('MAD') == risk_mad)
  call assert_true(risk_code('unknown') == 0)
  print '(a)', 'test_risk: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_risk
