! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_elliptical
  use fcopulae, only : dp, i8, pi, copula_fit_result, elliptical_marginal_pdf, elliptical_marginal_cdf, &
    elliptical_marginal_quantile, elliptical_copula_density, elliptical_copula_cdf, elliptical_rng, &
    elliptical_tau, elliptical_rho, elliptical_tail_coeff, elliptical_fit
  use test_support
  implicit none
  character(len=16),parameter :: names(7)=[character(len=16)::'norm','cauchy','t','logistic','laplace','kotz','epower']
  integer :: i
  real(dp) :: rho,p,q,c,d,expected,lt,ut,a1,a2
  real(dp),allocatable :: x(:,:),u(:),v(:)
  type(copula_fit_result) :: fit

  rho=0.45_dp;p=0.73_dp
  expected=0.25_dp+asin(rho)/(2.0_dp*pi)
  call assert_close(elliptical_copula_cdf(0.5_dp,0.5_dp,rho,'norm'),expected,2.0e-7_dp, &
    'bivariate normal quadrant probability')
  call assert_close(elliptical_tau(rho),2.0_dp*asin(rho)/pi,1.0e-14_dp,'elliptical Kendall tau')
  call assert_close(elliptical_rho(rho,'norm'),6.0_dp*asin(rho/2.0_dp)/pi,1.0e-14_dp,'normal Spearman rho')

  do i=1,size(names)
    a1=1.0_dp;a2=1.0_dp
    if(trim(names(i))=='t')a1=5.0_dp
    if(trim(names(i))=='kotz')a1=1.4_dp
    if(trim(names(i))=='epower')then;a1=1.2_dp;a2=0.9_dp;end if
    q=elliptical_marginal_quantile(p,names(i),a1,a2)
    call assert_close(elliptical_marginal_cdf(q,names(i),a1,a2),p,3.0e-3_dp,'elliptical marginal quantile inversion')
    call assert_true(elliptical_marginal_pdf(q,names(i),a1,a2)>0.0_dp,'elliptical marginal density positive')
    c=elliptical_copula_cdf(0.41_dp,0.67_dp,rho,names(i),a1,a2)
    d=elliptical_copula_density(0.41_dp,0.67_dp,rho,names(i),a1,a2)
    call assert_true(c>=0.0_dp.and.c<=0.41_dp+1.0e-8_dp,'elliptical CDF bounds')
    call assert_true(d>=0.0_dp,'elliptical copula density nonnegative')
    call assert_finite(d,'elliptical copula density finite')
    call elliptical_rng(80,rho,names(i),x,int(12000+i,i8),a1,a2)
    call assert_all_finite(x,'elliptical simulation finite')
    call assert_true(minval(x)>0.0_dp.and.maxval(x)<1.0_dp,'elliptical simulation in open unit square')
    call assert_true(abs(sum(x(:,1))/real(size(x,1),dp)-0.5_dp)<0.13_dp,'elliptical first margin approximately uniform')
  end do

  call elliptical_tail_coeff(rho,'norm',lt,ut)
  call assert_close(lt,0.0_dp,1.0e-14_dp,'normal lower tail dependence')
  call assert_close(ut,0.0_dp,1.0e-14_dp,'normal upper tail dependence')
  call elliptical_tail_coeff(rho,'t',lt,ut,5.0_dp)
  call assert_true(lt>0.0_dp.and.abs(lt-ut)<1.0e-12_dp,'Student t symmetric tail dependence')

  allocate(u(5),v(5))
  u=[0.10_dp,0.27_dp,0.43_dp,0.69_dp,0.88_dp]
  v=[0.14_dp,0.31_dp,0.51_dp,0.73_dp,0.92_dp]
  do i=1,size(names)
    fit=elliptical_fit(u,v,names(i),max_iter=1)
    call assert_true(size(fit%param)>=1,'elliptical fit parameter count')
    call assert_finite(fit%loglik,'elliptical fit log likelihood')
    call assert_true(abs(fit%param(1))<1.0_dp,'elliptical fit correlation bound')
  end do

  print '(a)','Elliptical copula tests passed.'
end program test_elliptical
