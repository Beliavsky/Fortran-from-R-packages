! SPDX-License-Identifier: GPL-3.0-only
program test_calibration
   use sharpe_rratio, only : dp, calibration_a, calibration_a_medium, &
      calibration_f, a_full, f_full, correction_b, theta_snr
   implicit none
   real(dp) :: positive, negative

   call assert_close(calibration_a(0.0_dp),0.0_dp,0.0_dp,'a at zero')
   call assert_close(calibration_a(0.41_dp),0.36245061149957286_dp, &
      2.0e-15_dp,'a interior knot')
   call assert_close(calibration_a(0.82_dp),1.0520620689598335_dp, &
      2.0e-15_dp,'a final knot')
   call assert_close(calibration_a_medium(0.6186242400000002_dp), &
      0.629940788348712_dp,2.0e-15_dp,'a-medium first knot')
   call assert_close(calibration_a_medium(0.9682733699999996_dp), &
      1.889822365046136_dp,2.0e-15_dp,'a-medium interior knot')
   call assert_close(calibration_f(0.246_dp),1.7623132950234297_dp, &
      2.0e-15_dp,'f first knot')
   call assert_close(calibration_f(3.115_dp),0.519974154175771_dp, &
      2.0e-15_dp,'f interior knot')
   call assert_close(f_full(4.0_dp),f_full(3.0_dp),0.0_dp,'f upper clamp')
   call assert_true(a_full(0.75_dp) > a_full(0.70_dp),'a-full monotonic segment')
   call assert_true(correction_b(0.2_dp,100) > 0.0_dp,'positive correction b')

   positive = theta_snr(0.25_dp,100,5.0_dp,.true.)
   negative = theta_snr(-0.25_dp,100,5.0_dp,.true.)
   call assert_close(positive,-negative,2.0e-15_dp,'theta sign symmetry')
   call assert_true(positive > 0.0_dp,'positive theta')

   print '(a)', 'test_calibration: PASS'

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

end program test_calibration
