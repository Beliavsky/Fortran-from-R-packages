! SPDX-License-Identifier: GPL-3.0-only
program test_mvt
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use fitheavytail
   use test_support, only: check, make_data, sample_cov_mle
   implicit none
   integer :: i
   real(dp) :: x(240,3), xna(240,3), target(3,3), nan_value
   type(heavy_tail_fit) :: gaussian, robust, factor_fit, missing_fit

   call make_data(x)
   target=sample_cov_mle(x)
   call fit_mvt(x,gaussian,fixed_nu=1.0e12_dp, &
      scale_covmat=.false.,ptol=1.0e-7_dp,max_iter=50)
   call check(gaussian%status==ht_success.or. &
      gaussian%status==ht_no_convergence,'Gaussian status')
   call check(maxval(abs(gaussian%mu-sum(x,dim=1)/240.0_dp))<1.0e-6_dp, &
      'Gaussian mean')
   call check(maxval(abs(gaussian%covariance-target))<1.0e-5_dp, &
      'Gaussian covariance')

   call fit_mvt(x,robust,nu_method='iterative', &
      nu_iterative_method='POP',max_iter=150,ptol=5.0e-4_dp)
   call check(robust%nu>=2.5_dp.and.robust%nu<=100.0_dp, &
      'iterative nu bounds')
   call check(allocated(robust%latent_weights),'latent weights')

   call fit_mvt(x,factor_fit,fixed_nu=7.0_dp,factors=2, &
      max_iter=80,ptol=1.0e-4_dp)
   call check(allocated(factor_fit%loadings),'factor loadings')
   call check(all(shape(factor_fit%loadings)==[3,2]), &
      'factor loading shape')
   call check(allocated(factor_fit%psi),'factor psi')
   call check(all(factor_fit%psi>=0.0_dp),'factor psi nonnegative')

   xna=x
   nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
   xna(1,1)=nan_value
   xna(2,2)=nan_value
   xna(3,3)=nan_value
   call fit_mvt(xna,missing_fit,fixed_nu=7.0_dp,na_rm=.false., &
      max_iter=80,ptol=1.0e-4_dp)
   call check(allocated(missing_fit%covariance),'missing covariance')
   call check(all([(missing_fit%covariance(i,i)>0.0_dp,i=1,3)]), &
      'missing-data positive diagonal')
   write(*,'(a)') 'test_mvt: PASS'
end program test_mvt
