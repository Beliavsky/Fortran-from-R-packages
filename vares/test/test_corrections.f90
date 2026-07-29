! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program test_corrections
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use vares
  implicit none
  integer :: failures
  real(dp) :: x, h, lower, upper

  failures = 0

  call check_close('CLG lower inversion', &
    pclg(varclg(0.20_dp, 1.4_dp, 1.2_dp, 0.3_dp), &
      1.4_dp, 1.2_dp, 0.3_dp), 0.20_dp, 3.0e-12_dp, failures)

  x = -0.7_dp
  lower = plaplace(x, 0.2_dp, 1.3_dp)
  upper = plaplace(x, 0.2_dp, 1.3_dp, lower_tail=.false.)
  call check_close('Laplace tail complement', lower+upper, 1.0_dp, &
    3.0e-14_dp, failures)
  call check_close('Laplace lower log', &
    plaplace(x, 0.2_dp, 1.3_dp, log_p=.true.), log(lower), &
    3.0e-14_dp, failures)
  call check_close('Laplace upper log', &
    plaplace(x, 0.2_dp, 1.3_dp, log_p=.true., lower_tail=.false.), &
    log(upper), 3.0e-14_dp, failures)

  call check_close('log-Laplace upper inversion', &
    ploglaplace(varloglaplace(0.83_dp, 1.4_dp, 2.2_dp, 1.7_dp), &
      1.4_dp, 2.2_dp, 1.7_dp), 0.83_dp, 4.0e-12_dp, failures)

  call check_close('asymmetric-power left inversion', &
    pasypower(varasypower(0.20_dp, 0.4_dp, 1.3_dp, 1.2_dp), &
      0.4_dp, 1.3_dp, 1.2_dp), 0.20_dp, 2.0e-8_dp, failures)
  call check_close('asymmetric-power right inversion', &
    pasypower(varasypower(0.80_dp, 0.4_dp, 1.3_dp, 1.2_dp), &
      0.4_dp, 1.3_dp, 1.2_dp), 0.80_dp, 2.0e-8_dp, failures)

  x = 0.6_dp
  call check_flags('generalized logistic III flags', &
    pgenlogis3(x, 1.6_dp, 0.2_dp, 1.3_dp), &
    pgenlogis3(x, 1.6_dp, 0.2_dp, 1.3_dp, lower_tail=.false.), &
    pgenlogis3(x, 1.6_dp, 0.2_dp, 1.3_dp, log_p=.true.), &
    failures)
  call check_flags('generalized logistic IV flags', &
    pgenlogis4(x, 1.2_dp, 1.6_dp, 0.2_dp, 1.3_dp), &
    pgenlogis4(x, 1.2_dp, 1.6_dp, 0.2_dp, 1.3_dp, &
      lower_tail=.false.), &
    pgenlogis4(x, 1.2_dp, 1.6_dp, 0.2_dp, 1.3_dp, &
      log_p=.true.), failures)
  call check_flags('Stacy gamma flags', &
    pstacygamma(1.2_dp, 1.4_dp, 1.3_dp, 1.1_dp), &
    pstacygamma(1.2_dp, 1.4_dp, 1.3_dp, 1.1_dp, lower_tail=.false.), &
    pstacygamma(1.2_dp, 1.4_dp, 1.3_dp, 1.1_dp, log_p=.true.), &
    failures)
  call check_flags('Nakagami flags', &
    pnakagami(1.2_dp, 1.7_dp, 1.4_dp), &
    pnakagami(1.2_dp, 1.7_dp, 1.4_dp, lower_tail=.false.), &
    pnakagami(1.2_dp, 1.7_dp, 1.4_dp, log_p=.true.), failures)

  call check_close('Kumaraswamy log-logistic inversion', &
    pkumloglogis(varkumloglogis(0.37_dp, 1.2_dp, 1.5_dp, &
      1.1_dp, 1.4_dp), 1.2_dp, 1.5_dp, 1.1_dp, 1.4_dp), &
    0.37_dp, 4.0e-12_dp, failures)

  lower = pburr7(1.3_dp, 1.4_dp, 1.2_dp)
  call check_close('Burr XII lower log', &
    pburr7(1.3_dp, 1.4_dp, 1.2_dp, log_p=.true.), &
    log(lower), 3.0e-14_dp, failures)

  call check_close('Dagum inversion', &
    pdagum(vardagum(0.37_dp, 1.3_dp, 1.1_dp, 1.4_dp), &
      1.3_dp, 1.1_dp, 1.4_dp), 0.37_dp, 4.0e-12_dp, failures)

  call check_close('Kumaraswamy Weibull inversion', &
    pkumweibull(varkumweibull(0.37_dp, 1.2_dp, 1.5_dp, &
      1.4_dp, 1.1_dp), 1.2_dp, 1.5_dp, 1.4_dp, 1.1_dp), &
    0.37_dp, 5.0e-12_dp, failures)
  x = 0.9_dp
  h = 1.0e-6_dp
  call check_close('Kumaraswamy Weibull density identity', &
    (pkumweibull(x+h, 1.2_dp, 1.5_dp, 1.4_dp, 1.1_dp) - &
     pkumweibull(x-h, 1.2_dp, 1.5_dp, 1.4_dp, 1.1_dp))/(2.0_dp*h), &
    dkumweibull(x, 1.2_dp, 1.5_dp, 1.4_dp, 1.1_dp), &
    2.0e-7_dp, failures)

  call check_close('exponentiated Weibull inversion', &
    pexpweibull(varexpweibull(0.37_dp, 1.3_dp, 1.4_dp, 1.1_dp), &
      1.3_dp, 1.4_dp, 1.1_dp), 0.37_dp, 5.0e-12_dp, failures)
  call check_close('exponentiated Weibull density identity', &
    (pexpweibull(x+h, 1.3_dp, 1.4_dp, 1.1_dp) - &
     pexpweibull(x-h, 1.3_dp, 1.4_dp, 1.1_dp))/(2.0_dp*h), &
    dexpweibull(x, 1.3_dp, 1.4_dp, 1.1_dp), &
    2.0e-7_dp, failures)

  lower = pinvexpexp(1.4_dp, 1.3_dp, 1.2_dp)
  call check_close('inverse exponentiated exponential lower log', &
    pinvexpexp(1.4_dp, 1.3_dp, 1.2_dp, log_p=.true.), &
    log(lower), 3.0e-14_dp, failures)

  if (failures /= 0) error stop 'VaRES correction tests failed'
  print '(a)', 'All VaRES correction tests passed.'
contains
  subroutine check_flags(label, low, up, loglow, failures)
    character(*), intent(in) :: label
    real(dp), intent(in) :: low, up, loglow
    integer, intent(inout) :: failures
    call check_close(trim(label)//' complement', low+up, 1.0_dp, &
      3.0e-12_dp, failures)
    call check_close(trim(label)//' log', loglow, log(low), &
      3.0e-12_dp, failures)
  end subroutine check_flags

  subroutine check_close(label, actual, expected, tol, failures)
    character(*), intent(in) :: label
    real(dp), intent(in) :: actual, expected, tol
    integer, intent(inout) :: failures
    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > tol*(1.0_dp+abs(expected))) then
      failures = failures + 1
      print '(a,2es24.14)', trim(label)//' failed: ', actual, expected
    end if
  end subroutine check_close
end program test_corrections
