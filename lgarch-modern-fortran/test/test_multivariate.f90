! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program test_multivariate
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use lgarch_kinds, only : dp,pi
  use lgarch_rng, only : seed_rng
  use lgarch_linalg, only : covariance_matrix,inverse_matrix,logdet_spd
  use lgarch_multivariate, only : rmnorm,mlgarch_simulate,mlgarch_varma_recursion,mlgarch_objective,fit_mlgarch, &
    mlgarch_fit_result,mlgarch_is_stable
  implicit none
  integer,parameter::n=180,m=2
  real(dp)::mu(2),cov(2,2),draws(20000,2),cm(2,2)
  real(dp)::constant(2),arch(2,2),garch(2,2),zz(4,2),yy(4,2),ls(4,2),expect(4,2),prev(2),prevlz(2)
  real(dp)::yr(3,2),u(3,2),ly(3,2),phi(2,2),theta(2,2),pars(13),manual,sinv(2,2),ld,q
  real(dp),allocatable::y(:,:),z(:,:)
  type(mlgarch_fit_result)::fitres,fitcustom
  integer::info,t,idx
  logical::st,st2

  call seed_rng(131)
  mu=[0.2_dp,-0.1_dp]; cov=reshape([1.0_dp,0.35_dp,0.35_dp,2.0_dp],[2,2])
  call rmnorm(20000,mu,cov,draws,info); if(info/=0) error stop "rmnorm failed"
  call assert_close(sum(draws,dim=1)/20000.0_dp,mu,0.035_dp,"rmnorm mean")
  cm=covariance_matrix(draws); call assert_close(reshape(cm,[4]),reshape(cov,[4]),0.07_dp,"rmnorm covariance")

  constant=[-0.1_dp,-0.2_dp]; arch=reshape([0.1_dp,0.02_dp,0.01_dp,0.08_dp],[2,2]); garch=reshape([0.65_dp,0.01_dp,0.02_dp,0.7_dp],[2,2])
  zz=reshape([0.5_dp,-0.8_dp,1.2_dp,0.7_dp,-1.0_dp,0.6_dp,0.9_dp,-1.3_dp],[4,2])
  call mlgarch_simulate(4,constant,arch,garch,yy,innovations=zz,backcast_lnsigma2=[-0.4_dp,-0.3_dp], &
    backcast_lnz2=[-0.2_dp,-0.1_dp],log_sigma2=ls,stable=st)
  st2=mlgarch_is_stable(arch,garch)
  if(.not.st .or. .not.st2) error stop "multivariate stability"
  prev=[-0.4_dp,-0.3_dp]; prevlz=[-0.2_dp,-0.1_dp]
  do t=1,4
    expect(t,:)=constant+matmul(arch+garch,prev)+matmul(arch,prevlz)
    prev=expect(t,:); prevlz=log(zz(t,:)*zz(t,:))
  end do
  call assert_close(reshape(ls,[8]),reshape(expect,[8]),1.0e-12_dp,"multivariate simulation")
  call assert_close(reshape(yy,[8]),reshape(exp(0.5_dp*expect)*zz,[8]),1.0e-12_dp,"multivariate returns")

  yr=reshape([1.0_dp,0.0_dp,2.0_dp,0.8_dp,1.1_dp,0.0_dp],[3,2]); phi=0.0_dp; theta=0.0_dp
  phi(1,1)=0.3_dp; phi(2,2)=0.25_dp; theta(1,1)=-0.1_dp; theta(2,2)=-0.15_dp
  call mlgarch_varma_recursion(yr,[0.1_dp,-0.05_dp],phi,theta,u,ly)
  if(abs(u(2,1))>1.0e-14_dp .or. abs(u(3,2))>1.0e-14_dp) error stop "multivariate zero residual"
  if(.not.all(ieee_is_finite(u)) .or. .not.all(ieee_is_finite(ly))) error stop "multivariate recursion finite"

  idx=1; pars(idx:idx+1)=[0.1_dp,-0.05_dp]; idx=idx+2; pars(idx:idx+3)=reshape(phi,[4]); idx=idx+4
  pars(idx:idx+3)=reshape(theta,[4]); idx=idx+4; pars(idx:idx+1)=[1.8_dp,2.2_dp]; idx=idx+2; pars(idx)=0.25_dp
  manual=mlgarch_objective(yr,pars,1,1)
  if(.not.ieee_is_finite(manual)) error stop "multivariate objective nonfinite"
  cov=reshape([1.8_dp,0.25_dp,0.25_dp,2.2_dp],[2,2]); call inverse_matrix(cov,sinv,info); call logdet_spd(cov,ld,info)
  q=dot_product(u(1,:),matmul(sinv,u(1,:)))
  call assert_scalar(manual,-0.5_dp*2.0_dp*log(2.0_dp*pi)-0.5_dp*ld-0.5_dp*q,1.0e-11_dp,"multivariate objective reference")

  allocate(y(n,m),z(n,m)); call seed_rng(551)
  cov=reshape([1.0_dp,0.3_dp,0.3_dp,1.0_dp],[2,2])
  call mlgarch_simulate(n,constant,arch,garch,y,innovations_vcov=cov,z=z)
  y(41,1)=0.0_dp; y(127,2)=0.0_dp
  call fit_mlgarch(y,1,1,fitres,compute_vcov=.true.,max_iter=1600,tol=5.0e-6_dp)
  if(size(fitres%varma_par)/=13) error stop "multivariate fit parameter count"
  if(.not.all(ieee_is_finite(fitres%varma_par))) error stop "multivariate fit nonfinite parameters"
  if(.not.all(ieee_is_finite(fitres%fitted_sd)) .or. any(fitres%fitted_sd<=0.0_dp)) error stop "multivariate fitted sd"
  if(.not.ieee_is_finite(fitres%objective_varma) .or. .not.ieee_is_finite(fitres%loglik_model)) error stop "multivariate fit likelihood"
  if(any(shape(fitres%vcov_varma)/=[13,13])) error stop "multivariate vcov size"
  call fit_mlgarch(y(:60,:),0,0,fitcustom,compute_vcov=.false.,initial_values=[-5.0_dp,-5.0_dp,4.0_dp,4.0_dp,0.0_dp], &
    lower_bounds=[-20.0_dp,-20.0_dp,1.0e-8_dp,1.0e-8_dp,-10.0_dp], &
    upper_bounds=[20.0_dp,20.0_dp,100.0_dp,100.0_dp,10.0_dp],objective_penalty=8.0e29_dp,max_iter=600)
  if(size(fitcustom%varma_par)/=5 .or. .not.all(ieee_is_finite(fitcustom%varma_par))) error stop "custom multivariate fit"
  print '(a)', 'Multivariate simulation, recursion, objective, fitting, and covariance tests passed.'
contains
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
end program test_multivariate
