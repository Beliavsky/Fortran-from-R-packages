! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_fitting
  use fmultivar, only : dp, i8, skew_fit_result, fit_multivariate_normal, fit_skew_normal, &
    fit_skew_t, fit_skew_cauchy, mv_fit, mvsnorm_rng, mvst_rng, is_positive_definite
  implicit none
  real(dp),allocatable::x(:,:)
  type(skew_fit_result)::fit
  real(dp)::mu(1),omega(1,1),alpha(1),mu2(2),omega2(2,2),alpha2(2)
  logical::ok,pd_ok
  integer::failures
  failures=0;mu=[-0.3_dp];omega=reshape([1.4_dp],[1,1]);alpha=[3.0_dp]
  call mvsnorm_rng(550,mu,omega,alpha,x,112233_i8,ok)
  fit=fit_multivariate_normal(x)
  pd_ok=is_positive_definite(fit%omega)
  call check_true(fit%converged.and.pd_ok,'normal fit',failures)
  fit=fit_skew_normal(x,1200,2.0e-6_dp)
  pd_ok=is_positive_definite(fit%omega)
  call check_true(fit%loglik>-huge(1.0_dp)/2.0_dp.and.pd_ok, &
    'skew normal finite fit',failures)
  call check_true(fit%alpha(1)>0.2_dp,'skew normal direction',failures)
  call check_true(allocated(fit%hessian).and.allocated(fit%covariance),'skew normal inference',failures)

  call mvst_rng(650,mu,omega,[2.0_dp],6.0_dp,x,223344_i8,ok)
  fit=fit_skew_t(x,fixed_nu=6.0_dp,max_iter=1400,tol=2.0e-6_dp)
  call check_true(abs(fit%nu-6.0_dp)<1.0e-12_dp.and.fit%alpha(1)>0.1_dp,'fixed skew t fit',failures)
  fit=fit_skew_t(x,max_iter=1600,tol=5.0e-6_dp)
  call check_true(fit%nu>0.2_dp.and.fit%nu<200.0_dp.and.fit%alpha(1)>0.05_dp, &
    'free skew t fit',failures)
  call mvst_rng(500,mu,omega,[1.5_dp],1.0_dp,x,334455_i8,ok)
  fit=fit_skew_cauchy(x,1100,8.0e-6_dp)
  pd_ok=is_positive_definite(fit%omega)
  call check_true(abs(fit%nu-1.0_dp)<1.0e-12_dp.and.pd_ok, &
    'skew cauchy fit',failures)
  fit=mv_fit(x,'normal')
  call check_true(fit%converged,'mv_fit normal dispatch',failures)
  fit=mv_fit(x,'st',fixed_nu=6.0_dp,max_iter=1000,tol=5.0e-6_dp)
  call check_true(abs(fit%nu-6.0_dp)<1.0e-12_dp,'mv_fit skew t dispatch',failures)

  mu2=[-0.2_dp,0.3_dp]
  omega2=reshape([1.0_dp,0.35_dp,0.35_dp,1.4_dp],[2,2])
  alpha2=[2.0_dp,-1.0_dp]
  call mvsnorm_rng(200,mu2,omega2,alpha2,x,445566_i8,ok)
  fit=fit_skew_normal(x,1800,1.0e-5_dp);pd_ok=is_positive_definite(fit%omega)
  call check_true(fit%converged.and.pd_ok.and.fit%alpha(1)>0.1_dp, &
    'two-dimensional skew normal fit',failures)
  call mvst_rng(220,mu2,omega2,alpha2,6.0_dp,x,556677_i8,ok)
  fit=fit_skew_t(x,6.0_dp,1800,1.0e-5_dp);pd_ok=is_positive_definite(fit%omega)
  call check_true(fit%converged.and.pd_ok.and.abs(fit%nu-6.0_dp)<1.0e-12_dp, &
    'two-dimensional skew t fit',failures)
  call mvst_rng(200,mu2,omega2,alpha2,1.0_dp,x,667788_i8,ok)
  fit=fit_skew_cauchy(x,1800,2.0e-5_dp);pd_ok=is_positive_definite(fit%omega)
  call check_true(fit%converged.and.pd_ok.and.abs(fit%nu-1.0_dp)<1.0e-12_dp, &
    'two-dimensional skew cauchy fit',failures)

  if(failures>0)then;write(*,'(a,i0)')'Fitting test failures: ',failures;error stop 1;end if
  write(*,'(a)')'Multivariate Normal, skew-Normal, skew-t, and skew-Cauchy fitting tests passed.'
contains
  subroutine check_true(cond,name,nfail)
    logical,intent(in)::cond
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.cond)then;write(*,'(a)')trim(name)//' failed';nfail=nfail+1;end if
  end subroutine check_true
end program test_fitting
