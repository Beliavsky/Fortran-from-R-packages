! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_extreme_empirical
  use fcopulae, only : dp, i8, copula_fit_result, copula_grid, ev_default_param, ev_check, ev_dependence, &
    ev_cdf, ev_density, ev_rng, ev_tau, ev_rho, ev_tail_coeff, ev_fit, archm_cdf, empirical_copula_cdf, &
    empirical_copula_grid, empirical_density_grid, frechet_copula_cdf, marshall_olkin_cdf, debye_function
  use test_support
  implicit none
  character(len=16),parameter :: names(5)=[character(len=16)::'gumbel','galambos','husler.reiss','tawn','bb5']
  integer :: i,j
  real(dp),allocatable :: param(:),sample(:,:),u(:),v(:)
  real(dp) :: a,c,d,tau,rho,lt,ut,mass
  type(copula_fit_result) :: fit
  type(copula_grid) :: cg,dg

  do i=1,size(names)
    call ev_default_param(names(i),param)
    call assert_true(ev_check(param,names(i)),'extreme-value default parameter valid')
    do j=1,9
      a=ev_dependence(real(j,dp)/10.0_dp,param,names(i))
      call assert_true(a>=max(real(j,dp)/10.0_dp,1.0_dp-real(j,dp)/10.0_dp)-1.0e-10_dp.and.a<=1.0_dp+1.0e-10_dp, &
        'Pickands dependence bounds')
    end do
    c=ev_cdf(0.38_dp,0.72_dp,param,names(i));d=ev_density(0.38_dp,0.72_dp,param,names(i))
    call assert_true(c>=0.10_dp.and.c<=0.38_dp+1.0e-8_dp,'extreme-value CDF bounds')
    call assert_true(d>=0.0_dp,'extreme-value density nonnegative')
    call assert_finite(d,'extreme-value density finite')
    call ev_rng(50,param,names(i),sample,int(15000+i,i8))
    call assert_all_finite(sample,'extreme-value simulation finite')
    call assert_true(minval(sample)>0.0_dp.and.maxval(sample)<1.0_dp,'extreme-value simulation unit square')
    tau=ev_tau(param,names(i));rho=ev_rho(param,names(i));call ev_tail_coeff(param,names(i),lt,ut)
    call assert_true(tau>=-1.0_dp.and.tau<=1.0_dp,'extreme-value Kendall tau bounds')
    call assert_true(rho>=-1.0_dp.and.rho<=1.0_dp,'extreme-value Spearman rho bounds')
    call assert_true(lt>=0.0_dp.and.ut>=0.0_dp.and.ut<=1.0_dp,'extreme-value tail bounds')
  end do

  param=[1.2_dp]
  call assert_true(ev_dependence(0.3_dp,param,'gumbelII')>=0.7_dp,'supplemental Gumbel-II Pickands bound')
  call assert_close(ev_dependence(0.3_dp,param,'pi'),1.0_dp,1.0e-14_dp,'independence Pickands function')
  call assert_close(ev_dependence(0.3_dp,param,'m'),0.7_dp,1.0e-14_dp,'comonotonic Pickands function')

  call ev_default_param('gumbel',param);param=2.0_dp
  call assert_close(ev_cdf(0.38_dp,0.72_dp,param,'gumbel'),archm_cdf(0.38_dp,0.72_dp,2.0_dp,4),2.0e-12_dp, &
    'Gumbel extreme-value and Archimedean CDF equivalence')
  call assert_close(ev_tau(param,'gumbel'),0.5_dp,3.0e-4_dp,'Gumbel extreme-value Kendall tau')
  call ev_tail_coeff(param,'gumbel',lt,ut)
  call assert_close(ut,2.0_dp-sqrt(2.0_dp),2.0e-10_dp,'Gumbel upper-tail coefficient')

  allocate(u(8),v(8))
  u=[0.07_dp,0.16_dp,0.29_dp,0.41_dp,0.55_dp,0.67_dp,0.82_dp,0.94_dp]
  v=[0.09_dp,0.21_dp,0.34_dp,0.46_dp,0.61_dp,0.75_dp,0.86_dp,0.97_dp]
  do i=1,size(names)
    fit=ev_fit(u,v,names(i),max_iter=1)
    call assert_true(size(fit%param)>=1,'extreme-value fit parameter count')
    call assert_finite(fit%loglik,'extreme-value fit log likelihood')
    call assert_true(ev_check(fit%param,names(i)),'extreme-value fit bounds')
  end do

  call assert_close(empirical_copula_cdf(0.5_dp,0.5_dp,u,v),0.5_dp,1.0e-14_dp,'empirical copula point CDF')
  cg=empirical_copula_grid(u,v,8)
  call assert_close(cg%z(1,1),0.0_dp,1.0e-14_dp,'empirical copula lower corner')
  call assert_close(cg%z(9,9),1.0_dp,1.0e-14_dp,'empirical copula upper corner')
  dg=empirical_density_grid(u,v,8)
  mass=sum(dg%z)/64.0_dp
  call assert_close(mass,1.0_dp,1.0e-12_dp,'empirical copula density mass')
  call assert_close(frechet_copula_cdf(0.4_dp,0.7_dp,'m'),0.4_dp,1.0e-14_dp,'Frechet upper copula')
  call assert_close(frechet_copula_cdf(0.4_dp,0.7_dp,'pi'),0.28_dp,1.0e-14_dp,'independence copula')
  call assert_close(frechet_copula_cdf(0.4_dp,0.7_dp,'w'),0.1_dp,1.0e-14_dp,'Frechet lower copula')
  call assert_close(marshall_olkin_cdf(0.4_dp,0.7_dp,0.0_dp,0.0_dp),0.28_dp,1.0e-14_dp, &
    'Marshall-Olkin independence limit')
  call assert_close(debye_function(1.0_dp,1),0.777504634112248_dp,2.0e-7_dp,'Debye D1 reference value')

  print '(a)','Extreme-value and empirical copula tests passed.'
end program test_extreme_empirical
