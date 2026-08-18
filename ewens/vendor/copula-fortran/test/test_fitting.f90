! SPDX-License-Identifier: GPL-3.0-or-later
program test_fitting
  use copula
  implicit none
  type(copula_model) :: model
  type(fit_result) :: fit
  real(dp), allocatable :: sample(:,:)
  real(dp) :: correlation(2,2)
  logical :: ok

  model = fgm_copula(0.7_dp)
  call rCopula(1200,model,sample,ok,20260727_i8)
  call assert_true(ok)
  fit = fit_copula(sample,family_fgm,'itau')
  call assert_true(fit%ok .and. fit%converged)
  call assert_close(fit%model%theta,0.7_dp,0.14_dp)
  fit = fit_copula(sample,family_fgm,'mpl')
  call assert_true(fit%ok .and. fit%converged)
  call assert_close(fit%model%theta,0.7_dp,0.14_dp)
  call assert_true(fit%standard_error > 0.0_dp)

  correlation = reshape([1.0_dp,0.55_dp,0.55_dp,1.0_dp],[2,2])
  model = normal_copula(correlation)
  call rCopula(1500,model,sample,ok,10101_i8)
  fit = fit_copula(sample,family_gaussian,'itau')
  call assert_true(fit%ok)
  call assert_close(fit%model%correlation(1,2),0.55_dp,0.06_dp)
  call assert_close(iTau(family_gaussian,tau(model)),0.55_dp,1.0e-13_dp)
  call assert_close(iRho(family_gaussian,rho(model)),0.55_dp,1.0e-13_dp)

  print '(a)', 'test_fitting: PASS'
contains
  subroutine assert_close(actual, reference, tolerance)
    real(dp), intent(in) :: actual, reference, tolerance
    if (abs(actual-reference) > tolerance*(1.0_dp+abs(reference))) then
      print '(a,3es25.16)', 'mismatch: ', actual, reference, abs(actual-reference)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true
end program test_fitting
