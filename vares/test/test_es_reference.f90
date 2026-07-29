! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program test_es_reference
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use vares
  implicit none
  integer :: failures
  real(dp), parameter :: p = 0.05_dp

  failures = 0
  call check_close('normal ES', esnormal(p, 0.2_dp, 1.3_dp), &
    -2.4815266497596533_dp, 3.0e-8_dp, failures)
  call check_close('Student t ES', est(p, 6.0_dp), &
    -2.7107385614799786_dp, 8.0e-8_dp, failures)
  call check_close('F ES', esf(p, 5.0_dp, 8.0_dp), &
    0.14045402504995963_dp, 8.0e-8_dp, failures)
  call check_close('gamma ES', esgamma(p, 2.2_dp, 1.4_dp), &
    0.20807113419669598_dp, 8.0e-8_dp, failures)
  call check_close('beta ES', esbetadist(p, 1.4_dp, 2.3_dp), &
    0.03297941532125819_dp, 8.0e-8_dp, failures)
  call check_close('lognormal ES', eslognorm(p, 0.2_dp, 1.1_dp), &
    0.13540489198223665_dp, 8.0e-8_dp, failures)
  call check_close('logistic ES', eslogistic(p, 0.2_dp, 1.3_dp), &
    -4.961396326992687_dp, 3.0e-8_dp, failures)
  call check_close('uniform ES', esuniform(p, 2.0_dp, 6.0_dp), &
    2.1_dp, 3.0e-13_dp, failures)
  call check_close('exponential ES', esexponential(p, 2.0_dp), &
    0.012713703318269935_dp, 3.0e-10_dp, failures)
  call check_close('Weibull ES', esweibull(p, 1.5_dp, 2.0_dp), &
    0.1645985741579519_dp, 3.0e-8_dp, failures)
  call check_close('Pareto ES', espareto(p, 1.0_dp, 3.0_dp), &
    1.0085241065536212_dp, 3.0e-10_dp, failures)
  call check_close('compound Laplace gamma ES', &
    esclg(p, 2.4_dp, 1.2_dp, 0.3_dp), &
    -2.595462689070291_dp, 5.0e-8_dp, failures)
  call check_close('Laplace ES', eslaplace(p, 0.2_dp, 1.3_dp), &
    -4.09336062089226_dp, 3.0e-8_dp, failures)
  call check_close('log-Laplace ES', &
    esloglaplace(p, 1.4_dp, 2.2_dp, 1.7_dp), &
    0.46003120489707555_dp, 5.0e-8_dp, failures)
  call check_close('asymmetric power ES', &
    esasypower(p, 0.4_dp, 1.3_dp, 1.2_dp), &
    -0.7214275319345309_dp, 8.0e-8_dp, failures)
  call check_close('generalized logistic III ES', &
    esgenlogis3(p, 1.6_dp, 0.2_dp, 1.3_dp), &
    -3.459833509292168_dp, 8.0e-8_dp, failures)
  call check_close('generalized logistic IV ES', &
    esgenlogis4(p, 1.2_dp, 1.6_dp, 0.2_dp, 1.3_dp), &
    -3.124962374988795_dp, 8.0e-8_dp, failures)
  call check_close('Stacy gamma ES', &
    esstacygamma(p, 1.4_dp, 1.3_dp, 1.1_dp), &
    0.15912610298033733_dp, 8.0e-8_dp, failures)
  call check_close('Kumaraswamy log-logistic ES', &
    eskumloglogis(p, 1.2_dp, 1.5_dp, 1.1_dp, 1.4_dp), &
    0.09397374126420492_dp, 8.0e-8_dp, failures)
  call check_close('Burr XII ES', esburr7(p, 1.4_dp, 1.2_dp), &
    0.03475589811928008_dp, 8.0e-8_dp, failures)
  call check_close('Dagum ES', esdagum(p, 1.3_dp, 1.1_dp, 1.4_dp), &
    0.14610467921992876_dp, 8.0e-8_dp, failures)
  call check_close('Kumaraswamy Weibull ES', &
    eskumweibull(p, 1.2_dp, 1.5_dp, 1.4_dp, 1.1_dp), &
    0.0926422877691206_dp, 8.0e-8_dp, failures)
  call check_close('exponentiated Weibull ES', &
    esexpweibull(p, 1.3_dp, 1.4_dp, 1.1_dp), &
    0.14030259825101493_dp, 8.0e-8_dp, failures)
  call check_close('inverse exponentiated exponential ES', &
    esinvexpexp(p, 1.3_dp, 1.2_dp), &
    0.32521026750487336_dp, 8.0e-8_dp, failures)
  call check_close('Nakagami ES', esnakagami(p, 1.7_dp, 1.4_dp), &
    0.34021741493037894_dp, 8.0e-8_dp, failures)

  if (failures /= 0) error stop 'VaRES ES reference tests failed'
  print '(a)', 'All VaRES ES reference tests passed.'
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
end program test_es_reference
