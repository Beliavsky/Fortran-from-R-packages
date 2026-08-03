! SPDX-License-Identifier: GPL-3.0-only
program test_mvst
   use fitheavytail
   use test_support, only: check, make_data
   implicit none
   integer :: i
   real(dp) :: x(120,2), zero_gamma(2), gamma_fixed(2)
   type(heavy_tail_fit) :: symmetric_fit, skew_fit, mvt_fit

   call make_data(x)
   zero_gamma=0.0_dp
   call fit_mvt(x,mvt_fit,fixed_nu=8.0_dp,scale_covmat=.false., &
      max_iter=100,ptol=1.0e-5_dp)
   call fit_mvst(x,symmetric_fit,fixed_nu=8.0_dp, &
      fixed_gamma=zero_gamma,max_iter=100,ptol=1.0e-5_dp)
   call check(maxval(abs(mvt_fit%mu-symmetric_fit%mu))<2.0e-4_dp, &
      'symmetric skew-t mean')
   call check(maxval(abs(mvt_fit%covariance- &
      symmetric_fit%covariance))<2.0e-3_dp, &
      'symmetric skew-t covariance')

   gamma_fixed=[0.15_dp,-0.08_dp]
   call fit_mvst(x,skew_fit,fixed_nu=8.0_dp, &
      initial_gamma=gamma_fixed,max_iter=80,ptol=1.0e-4_dp)
   call check(allocated(skew_fit%gamma),'skewness vector')
   call check(allocated(skew_fit%covariance),'skew covariance')
   call check(all([(skew_fit%scatter(i,i)>0.0_dp,i=1,2)]), &
      'skew scatter diagonal')
   write(*,'(a)') 'test_mvst: PASS'
end program test_mvst
