! SPDX-License-Identifier: GPL-3.0-only
program test_gaussian_ar1
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
  use imputefin
  use imputefin_rng, only : rng_state, rng_seed, rng_normal
  implicit none
  integer, parameter :: n=400
  real(dp) :: y(n),ym(n),nanv
  type(rng_state) :: rng
  type(ar1_fit_result) :: fit
  type(imputation_result) :: imp
  type(ar1_options) :: opt
  type(imputation_options) :: io
  integer :: i
  nanv=ieee_value(0.0_dp,ieee_quiet_nan)
  call rng_seed(rng,12345_8)
  y(1)=0.1_dp
  do i=2,n
    y(i)=0.2_dp+0.72_dp*y(i-1)+0.3_dp*rng_normal(rng)
  end do
  ym=y;ym(80:90)=nanv;ym(170)=nanv;ym(250:260)=nanv
  opt%tol=1.0e-7_dp;opt%maxiter=300
  call fit_ar1_gaussian(ym,fit,opt,return_iterates=.true.,return_conditional=.true.)
  call check(fit%status==impute_ok.or.fit%status==impute_not_converged,'Gaussian fit status')
  call check(abs(fit%phi0-0.2_dp)<0.08_dp,'Gaussian phi0')
  call check(abs(fit%phi1-0.72_dp)<0.10_dp,'Gaussian phi1')
  call check(abs(fit%sigma2-0.09_dp)<0.035_dp,'Gaussian sigma2')
  call check(allocated(fit%cond_mean).and.allocated(fit%cond_cov),'conditional moments')
  io%n_samples=2;io%seed=777_8
  call impute_ar1_gaussian(ym,imp,opt,io)
  call check(imp%status==impute_ok,'Gaussian imputation status')
  call check(.not.any(ieee_is_nan(imp%values(1:n,1,1))),'Gaussian imputation complete')
  do i=1,n
    if(.not.ieee_is_nan(ym(i)))call check(abs(imp%values(i,1,1)-ym(i))<1.0e-14_dp,'observed unchanged')
  end do
  print '(a)', 'test_gaussian_ar1: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine check
end program test_gaussian_ar1
