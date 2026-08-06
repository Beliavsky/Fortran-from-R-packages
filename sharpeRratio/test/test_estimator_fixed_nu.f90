! SPDX-License-Identifier: GPL-3.0-only
program test_estimator_fixed_nu
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use sharpe_rratio, only : dp, snr_result, estimate_snr
   implicit none
   real(dp) :: positive(40), negative(40), with_nan(41)
   type(snr_result) :: pos, neg, filtered

   positive = 0.01_dp
   negative = -0.01_dp
   pos = estimate_snr(positive,num_perm=10,nu=5.0_dp,seed=12_int64)
   neg = estimate_snr(negative,num_perm=10,nu=5.0_dp,seed=12_int64)
   call assert_true(pos%ok .and. neg%ok,'fixed-nu estimates')
   call assert_close(pos%r0bar,39.0_dp,0.0_dp,'positive R0bar')
   call assert_close(neg%r0bar,-39.0_dp,0.0_dp,'negative R0bar')
   call assert_true(pos%snr > 0.0_dp .and. neg%snr < 0.0_dp,'SNR signs')
   call assert_close(pos%snr,-neg%snr,2.0e-14_dp,'SNR sign symmetry')
   call assert_close(pos%ci_lower,pos%ci_upper,0.0_dp,'deterministic interval')

   with_nan(1:40) = positive
   with_nan(41) = ieee_value(0.0_dp,ieee_quiet_nan)
   filtered = estimate_snr(with_nan,num_perm=5,nu=5.0_dp,seed=99_int64)
   call assert_true(filtered%ok .and. filtered%n == 40,'non-finite filtering')

   print '(a)', 'test_estimator_fixed_nu: PASS'

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

end program test_estimator_fixed_nu
