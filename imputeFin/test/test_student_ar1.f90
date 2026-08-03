! SPDX-License-Identifier: GPL-3.0-only
program test_student_ar1
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
  use imputefin
  use imputefin_rng, only : rng_state, rng_seed, rng_student_t
  implicit none
  integer,parameter::n=500
  real(dp)::y(n),ym(n),yo(n),nanv
  type(rng_state)::rng
  type(ar1_fit_result)::fit,fit_saem
  type(imputation_result)::imp
  type(ar1_options)::opt
  type(imputation_options)::io
  integer::i
  nanv=ieee_value(0.0_dp,ieee_quiet_nan);call rng_seed(rng,4567_8)
  y(1)=0.0_dp
  do i=2,n;y(i)=0.05_dp+0.60_dp*y(i-1)+0.22_dp*rng_student_t(rng,5.0_dp);end do
  ym=y;ym(100:108)=nanv;ym(220)=nanv;ym(330:336)=nanv
  opt%maxiter=250;opt%tol=1.0e-6_dp;opt%fast_and_heuristic=.true.
  call fit_ar1_t(ym,fit,opt,return_iterates=.true.)
  call check(fit%status==impute_ok.or.fit%status==impute_not_converged,'Student fit status')
  call check(abs(fit%phi1-0.60_dp)<0.15_dp,'Student phi1')
  call check(fit%sigma2>0.0_dp.and.fit%nu>0.2_dp,'Student scale and nu')
  io%n_samples=2;io%n_burn=30;io%n_thin=2;io%seed=991_8
  call impute_ar1_t(ym,imp,opt,io)
  call check(imp%status==impute_ok,'Student imputation status')
  call check(.not.any(ieee_is_nan(imp%values(:,1,1))),'Student imputation complete')
  yo=ym;yo(180)=25.0_dp;opt%remove_outliers=.true.;opt%outlier_prob_th=1.0e-4_dp
  call fit_ar1_t(yo,fit,opt)
  call check(allocated(fit%index_outliers),'outlier index allocated')
  call check(any(fit%index_outliers==180),'outlier detected')
  opt%remove_outliers=.false.;opt%fast_and_heuristic=.false.;opt%n_chain=3;opt%n_thin=1
  opt%saem_burn=8;opt%maxiter=20;opt%tol=1.0e-3_dp
  call fit_ar1_t(ym,fit_saem,opt)
  call check(fit_saem%status==impute_ok.or.fit_saem%status==impute_not_converged,'SAEM path')
  call check(fit_saem%sigma2>0.0_dp,'SAEM positive scale')
  print '(a)', 'test_student_ar1: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine check
end program test_student_ar1
