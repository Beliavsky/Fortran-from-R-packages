! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_archimedean
  use fcopulae, only : dp, i8, copula_fit_result, archm_default_alpha, archm_range, archm_check, &
    archm_phi, archm_inv_phi, archm_cdf, archm_density, archm_rng, archm_tau, archm_rho, &
    archm_tail_coeff, archm_fit
  use test_support
  implicit none
  integer :: type_id
  real(dp) :: alpha,lower,upper,u,v,c,d,rho,lt,ut
  real(dp),allocatable :: sample(:,:),fit_u(:),fit_v(:)
  type(copula_fit_result) :: fit

  u=0.37_dp;v=0.68_dp
  do type_id=1,22
    alpha=archm_default_alpha(type_id)
    call archm_range(type_id,lower,upper)
    call assert_true(archm_check(alpha,type_id),'default Archimedean parameter is valid')
    call assert_close(archm_inv_phi(archm_phi(u,alpha,type_id),alpha,type_id),u,2.0e-8_dp, &
      'Archimedean generator inverse identity')
    c=archm_cdf(u,v,alpha,type_id)
    d=archm_density(u,v,alpha,type_id)
    call assert_true(c>=max(0.0_dp,u+v-1.0_dp)-1.0e-10_dp.and.c<=min(u,v)+1.0e-10_dp, &
      'Archimedean CDF obeys Frechet bounds')
    call assert_true(d>=0.0_dp,'Archimedean density is nonnegative')
    call assert_finite(d,'Archimedean density is finite')
    call archm_rng(30,alpha,type_id,sample,int(9000+type_id,i8))
    call assert_all_finite(sample,'Archimedean simulation is finite')
    call assert_true(minval(sample)>=0.0_dp.and.maxval(sample)<=1.0_dp,'Archimedean simulation lies in unit square')
  end do

  call assert_close(archm_cdf(u,v,0.0_dp,1),u*v,2.0e-12_dp,'type 1 independence limit')
  alpha=2.0_dp
  c=(u**(-alpha)+v**(-alpha)-1.0_dp)**(-1.0_dp/alpha)
  call assert_close(archm_cdf(u,v,alpha,1),c,2.0e-12_dp,'Clayton CDF formula')
  call assert_close(archm_tau(alpha,1),0.5_dp,3.0e-4_dp,'Clayton Kendall tau')
  call assert_close(archm_tau(2.0_dp,4),0.5_dp,3.0e-4_dp,'Gumbel Kendall tau')
  rho=archm_rho(2.0_dp,1)
  call assert_true(rho>0.6_dp.and.rho<0.8_dp,'Clayton Spearman rho range')
  call archm_tail_coeff(2.0_dp,1,lt,ut)
  call assert_close(lt,2.0_dp**(-0.5_dp),3.0e-4_dp,'Clayton lower-tail coefficient')
  call assert_true(ut<2.0e-3_dp,'Clayton upper-tail coefficient is zero')

  allocate(fit_u(6),fit_v(6))
  fit_u=[0.08_dp,0.19_dp,0.33_dp,0.52_dp,0.71_dp,0.89_dp]
  fit_v=[0.12_dp,0.25_dp,0.40_dp,0.58_dp,0.77_dp,0.93_dp]
  do type_id=1,22
    fit=archm_fit(fit_u,fit_v,type_id,max_iter=1)
    call assert_true(size(fit%param)==1,'Archimedean fit parameter count')
    call assert_finite(fit%loglik,'Archimedean fit log likelihood')
    call assert_true(archm_check(fit%param(1),type_id),'Archimedean fit respects bounds')
  end do

  print '(a)','Archimedean copula tests passed.'
end program test_archimedean
