! SPDX-License-Identifier: GPL-3.0-only
program test_estimator_automatic
   use, intrinsic :: iso_fortran_env, only : int64
   use ghyp, only : ghyp_model_type, student_t_uv, rghyp, normal_quantile
   use ghyp_kinds, only : i8
   use sharpe_rratio, only : dp, snr_result, estimate_snr, gaussian_nu
   implicit none
   integer, parameter :: n = 100
   real(dp) :: normal_sample(n)
   real(dp), allocatable :: heavy_sample(:,:)
   type(ghyp_model_type) :: model
   type(snr_result) :: normal_result, heavy_result
   logical :: ok
   integer :: i

   do i = 1, n
      normal_sample(i) = normal_quantile((real(i,dp)-0.5_dp)/real(n,dp))+0.05_dp
   end do
   normal_result = estimate_snr(normal_sample,num_perm=12,seed=123_int64)
   call assert_true(normal_result%ok,'automatic Gaussian estimate')
   call assert_true(normal_result%gaussian_selected,'Gaussian branch selected')
   call assert_close(normal_result%nu,gaussian_nu,0.0_dp,'Gaussian nu sentinel')

   model = student_t_uv(3.0_dp,mu=0.05_dp)
   call rghyp(80,model,heavy_sample,ok,12345_i8)
   call assert_true(ok,'heavy-tail simulation')
   heavy_result = estimate_snr(heavy_sample(:,1),num_perm=12,seed=321_int64, &
      max_fit_iterations=300)
   call assert_true(heavy_result%ok,'automatic heavy-tail estimate')
   call assert_true(heavy_result%nu_estimated,'nu estimated flag')
   call assert_true(.not. heavy_result%gaussian_selected,'heavy-tail branch selected')
   call assert_true(heavy_result%nu > 2.0_dp .and. heavy_result%nu < 100.0_dp, &
      'finite fitted tail exponent')

   print '(a)', 'test_estimator_automatic: PASS'

contains

   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop message
   end subroutine assert_true

   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance) error stop message
   end subroutine assert_close

end program test_estimator_automatic
