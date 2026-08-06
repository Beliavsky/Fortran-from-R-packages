program test_distribution
   use adgoftest, only : dp, ad_distribution_cdf, ad_p_value
   implicit none

   real(dp) :: value

   value = ad_distribution_cdf(0.5_dp, 10)
   call assert_close(value, 0.2573659942339062_dp, 2.0e-14_dp, 'cdf 0.5')

   value = ad_p_value(2.0_dp, 10)
   call assert_close(value, 0.0930646507877626_dp, 2.0e-14_dp, 'p 2.0')

   value = ad_distribution_cdf(0.1_dp, 10)
   call assert_true(value < 0.0_dp, 'source-compatible cdf should preserve approximation undershoot')

   value = ad_distribution_cdf(0.1_dp, 10, clamp_probability=.true.)
   call assert_close(value, 0.0_dp, 0.0_dp, 'clamped cdf')

   value = ad_p_value(0.1_dp, 10, clamp_probability=.true.)
   call assert_close(value, 1.0_dp, 0.0_dp, 'clamped p')

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es24.16))') trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_distribution
