! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program test_all
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use vares
  implicit none
  integer :: failures
  real(dp) :: p(5), q(5), back(5)

  failures = 0
  call check_close('normal cdf zero', pnormal(0.0_dp), 0.5_dp, 3.0e-14_dp, failures)
  call check_close('normal quantile median', varnormal(0.5_dp), 0.0_dp, 3.0e-14_dp, failures)
  call check_close('normal ES', esnormal(0.05_dp), &
    -dnormal(varnormal(0.05_dp))/0.05_dp, 2.0e-8_dp, failures)
  call check_close('exponential inversion', &
    pexponential(varexponential(0.25_dp, 2.0_dp), 2.0_dp), &
    0.25_dp, 2.0e-13_dp, failures)
  call check_close('uniform ES', esuniform(0.4_dp, 2.0_dp, 6.0_dp), &
    2.8_dp, 3.0e-13_dp, failures)
  call check_close('beta inversion', &
    pbetadist(varbetadist(0.37_dp, 1.4_dp, 2.3_dp), 1.4_dp, 2.3_dp), &
    0.37_dp, 2.0e-9_dp, failures)
  call check_close('gamma inversion', &
    pgamma(vargamma(0.37_dp, 2.2_dp, 1.4_dp), 2.2_dp, 1.4_dp), &
    0.37_dp, 2.0e-9_dp, failures)
  call check_close('student t inversion', &
    pT(varT(0.37_dp, 6.0_dp), 6.0_dp), 0.37_dp, 3.0e-9_dp, failures)
  call check_close('F inversion', &
    pF(varF(0.37_dp, 5.0_dp, 8.0_dp), 5.0_dp, 8.0_dp), &
    0.37_dp, 3.0e-9_dp, failures)
  call check_close('logistic upper log', &
    plogistic(1.2_dp, log_p=.true., lower_tail=.false.), &
    log(1.0_dp-plogistic(1.2_dp)), 2.0e-13_dp, failures)

  p = [0.01_dp, 0.05_dp, 0.25_dp, 0.50_dp, 0.95_dp]
  q = varnormal(p, mu=1.0_dp, sigma=2.0_dp)
  back = pnormal(q, mu=1.0_dp, sigma=2.0_dp)
  if (any(abs(back-p) > 3.0e-9_dp)) then
    failures = failures + 1
    print '(a)', 'elemental array inversion failed'
  end if

  if (failures /= 0) error stop 'VaRES core tests failed'
  print '(a)', 'All VaRES core tests passed.'
contains
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
end program test_all
