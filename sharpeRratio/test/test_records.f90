! SPDX-License-Identifier: GPL-3.0-only
program test_records
   use, intrinsic :: iso_fortran_env, only : int64
   use sharpe_rratio, only : dp, r0_result, num_records_up, num_records_down, compute_r0bar
   implicit none
   real(dp) :: increasing(4), positive(12), negative(12), mixed(8)
   type(r0_result) :: r1, r2, r3

   increasing = [1.0_dp,2.0_dp,2.0_dp,3.0_dp]
   call assert_true(num_records_up(increasing) == 3,'upper-record count')
   call assert_true(num_records_down(increasing) == 1,'lower-record count')
   call assert_true(num_records_up([3.0_dp,2.0_dp,2.0_dp,1.0_dp]) == 1, &
      'upper records decreasing')
   call assert_true(num_records_down([3.0_dp,2.0_dp,2.0_dp,1.0_dp]) == 3, &
      'lower records decreasing')

   positive = 0.1_dp
   negative = -0.1_dp
   r1 = compute_r0bar(positive,num_perm=7,seed=123_int64)
   r2 = compute_r0bar(negative,num_perm=7,seed=123_int64)
   call assert_true(r1%ok .and. r2%ok,'constant-sign computations')
   call assert_close(r1%mean,11.0_dp,1.0e-14_dp,'positive R0')
   call assert_close(r2%mean,-11.0_dp,1.0e-14_dp,'negative R0')
   call assert_close(r1%q1,r1%q2,0.0_dp,'positive deterministic quantiles')

   mixed = [0.4_dp,-0.2_dp,0.1_dp,-0.3_dp,0.5_dp,-0.1_dp,0.2_dp,-0.4_dp]
   r1 = compute_r0bar(mixed,num_perm=25,seed=9876_int64)
   r2 = compute_r0bar(mixed,num_perm=25,seed=9876_int64)
   r3 = compute_r0bar(mixed,num_perm=25,q1=0.1_dp,q2=0.9_dp, &
      seed=9876_int64,source_compatible=.false.)
   call assert_true(r1%ok .and. r2%ok .and. r3%ok,'mixed computations')
   call assert_close(r1%mean,r2%mean,0.0_dp,'seed reproducibility')
   call assert_close(r1%q1,r2%q1,0.0_dp,'seed quantile reproducibility')
   call assert_true(r3%q1 <= r3%q2,'corrected quantile ordering')

   print '(a)', 'test_records: PASS'

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

end program test_records
