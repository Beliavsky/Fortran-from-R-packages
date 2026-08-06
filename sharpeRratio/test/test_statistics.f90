! SPDX-License-Identifier: GPL-3.0-only
program test_statistics
   use ghyp, only : ghyp_model_type, student_t_uv, rghyp, normal_quantile
   use ghyp_kinds, only : i8
   use sharpe_rratio, only : dp, test_n, quantile_type7, sample_variance
   implicit none
   integer, parameter :: n = 100
   real(dp) :: normal_sample(n), statistic
   real(dp), allocatable :: sample(:,:)
   type(ghyp_model_type) :: model
   logical :: ok
   integer :: i

   do i = 1, n
      normal_sample(i) = normal_quantile((real(i,dp)-0.5_dp)/real(n,dp))
   end do
   statistic = test_n(normal_sample)
   call assert_true(abs(statistic) < 3.0_dp,'normal sample classification')
   call assert_close(quantile_type7([1.0_dp,2.0_dp,3.0_dp,4.0_dp],0.25_dp), &
      1.75_dp,1.0e-14_dp,'R type-7 quantile')
   call assert_close(sample_variance([1.0_dp,2.0_dp,3.0_dp]),1.0_dp, &
      1.0e-14_dp,'sample variance')

   model = student_t_uv(3.0_dp)
   call rghyp(n,model,sample,ok,24680_i8)
   call assert_true(ok,'Student simulation')
   statistic = test_n(sample(:,1))
   call assert_true(abs(statistic) > 3.0_dp,'heavy-tail sample classification')

   print '(a)', 'test_statistics: PASS'

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

end program test_statistics
