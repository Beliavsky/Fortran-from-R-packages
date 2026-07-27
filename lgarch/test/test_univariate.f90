! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program test_univariate
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use lgarch_kinds, only : dp,pi
  use lgarch_rng, only : seed_rng
  use lgarch_univariate, only : lgarch_simulate,lgarch_arma_recursion,lgarch_objective,fit_lgarch,lgarch_fit_result, &
    LGARCH_LS,LGARCH_ML,LGARCH_CEX2,lgarch_is_stable
  implicit none
  integer,parameter::n=240
  real(dp)::z5(5),y5(5),ls5(5),xr5(5),expected(5),prevls(2),prevlz(2),phivec(2),innov
  real(dp)::yr(3),u(3),ly(3),pars(4),manual,eln,fit
  real(dp),allocatable::y(:),z(:),xreg(:,:),sd(:)
  type(lgarch_fit_result)::fls,fml,fce,fcustom
  integer::t

  z5=[0.5_dp,-1.0_dp,2.0_dp,-0.7_dp,1.2_dp]
  xr5=[0.01_dp,-0.02_dp,0.03_dp,0.0_dp,-0.01_dp]
  call lgarch_simulate(5,-0.1_dp,[0.2_dp,-0.05_dp],[0.6_dp],y5,xreg=xr5,innovations=z5, &
    backcast_lnsigma2=[-0.4_dp,-0.5_dp],backcast_lnz2=[-0.3_dp,-0.2_dp],log_sigma2=ls5)
  prevls=[-0.5_dp,-0.4_dp]; prevlz=[-0.2_dp,-0.3_dp]; phivec=[0.8_dp,-0.05_dp]
  do t=1,5
    innov=-0.1_dp+xr5(t)+0.2_dp*prevlz(1)-0.05_dp*prevlz(2)
    expected(t)=innov+phivec(1)*prevls(1)+phivec(2)*prevls(2)
    prevls=[expected(t),prevls(1)]; prevlz=[log(z5(t)*z5(t)),prevlz(1)]
  end do
  call assert_close(ls5,expected,1.0e-12_dp,"arbitrary-order simulation recursion")
  call assert_close(y5,exp(0.5_dp*expected)*z5,1.0e-12_dp,"simulation returns")
  if(.not.lgarch_is_stable([0.2_dp,-0.05_dp],[0.6_dp])) error stop "stable model rejected"
  if(lgarch_is_stable([0.7_dp,0.1_dp],[0.5_dp])) error stop "unstable model accepted"

  yr=[1.0_dp,0.0_dp,2.0_dp]
  call lgarch_arma_recursion(yr,0.1_dp,0.4_dp,-0.2_dp,u,ly)
  eln=log(2.0_dp)
  manual=0.0_dp-(0.1_dp+0.4_dp*eln)
  call assert_scalar(u(1),manual,1.0e-12_dp,"recursion residual 1")
  fit=0.1_dp+0.4_dp*0.0_dp-0.2_dp*manual
  call assert_scalar(ly(2),fit,1.0e-12_dp,"zero imputation")
  call assert_scalar(u(2),0.0_dp,1.0e-14_dp,"zero residual")
  manual=log(4.0_dp)-(0.1_dp+0.4_dp*fit)
  call assert_scalar(u(3),manual,1.0e-12_dp,"recursion residual 3")

  pars=[0.1_dp,0.4_dp,-0.2_dp,2.0_dp]
  manual=-0.5_dp*2.0_dp*(log(2.0_dp)+log(2.0_dp*pi))-0.5_dp*(u(1)*u(1)+u(3)*u(3))/2.0_dp
  call assert_scalar(lgarch_objective(yr,pars,1,1,LGARCH_ML),manual,1.0e-12_dp,"Gaussian objective")

  allocate(y(n),z(n),xreg(n,1),sd(n)); call seed_rng(9182)
  do t=1,n; xreg(t,1)=sin(0.03_dp*real(t,dp)); end do
  call lgarch_simulate(n,-0.15_dp,[0.12_dp],[0.76_dp],y,xreg=0.04_dp*xreg(:,1),z=z,sigma=sd)
  y(37)=0.0_dp; y(119)=0.0_dp
  call fit_lgarch(y,1,1,LGARCH_LS,fls,xreg=xreg,compute_vcov=.true.,max_iter=2500,tol=2.0e-7_dp)
  call check_fit(fls,"LS")
  call fit_lgarch(y,1,1,LGARCH_ML,fml,xreg=xreg,compute_vcov=.true.,max_iter=3000,tol=2.0e-7_dp)
  call check_fit(fml,"ML")
  call fit_lgarch(y,1,1,LGARCH_CEX2,fce,xreg=xreg,compute_vcov=.true.,max_iter=3000,tol=2.0e-7_dp)
  call check_fit(fce,"CEX2")
  if(.not.fce%mean_correction) error stop "CEX2 did not force mean correction"
  call fit_lgarch(y(:80),0,0,LGARCH_LS,fcustom,compute_vcov=.false.,initial_values=[-5.0_dp], &
    lower_bounds=[-20.0_dp],upper_bounds=[20.0_dp],objective_penalty=9.0e29_dp,max_iter=500)
  if(fcustom%arma_par(1)<-20.0_dp .or. fcustom%arma_par(1)>20.0_dp) error stop "custom bounds failed"
  print '(a)', 'Univariate simulation, recursion, objective, fitting, and covariance tests passed.'
contains
  subroutine check_fit(r,label)
    type(lgarch_fit_result),intent(in)::r
    character(len=*),intent(in)::label
    if(size(r%lgarch_par)/=5) error stop label//" lgarch parameter count"
    if(r%method==LGARCH_LS) then
      if(size(r%arma_par)/=4) error stop label//" arma parameter count"
    else
      if(size(r%arma_par)/=5) error stop label//" arma parameter count"
    end if
    if(.not.all(ieee_is_finite(r%arma_par))) error stop label//" nonfinite parameters"
    if(.not.all(ieee_is_finite(r%fitted_sd)) .or. any(r%fitted_sd<=0.0_dp)) error stop label//" invalid fitted sd"
    if(.not.ieee_is_finite(r%loglik_model) .or. .not.ieee_is_finite(r%rss)) error stop label//" invalid statistics"
    if(any(shape(r%vcov_arma)/=[size(r%arma_par),size(r%arma_par)])) error stop label//" vcov size"
  end subroutine check_fit
  subroutine assert_close(actual,expected,tol,label)
    real(dp),intent(in)::actual(:),expected(:),tol
    character(len=*),intent(in)::label
    if(maxval(abs(actual-expected))>tol) then; print *,trim(label),actual,expected; error stop 1; end if
  end subroutine assert_close
  subroutine assert_scalar(actual,expected,tol,label)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::label
    if(abs(actual-expected)>tol) then; print *,trim(label),actual,expected; error stop 1; end if
  end subroutine assert_scalar
end program test_univariate
