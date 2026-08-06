program test_uniform_api
   use adgoftest, only : dp, ad_test_result, ad_test_uniform, ad_success
   implicit none

   real(dp), parameter :: x(5) = [0.05_dp, 0.2_dp, 0.4_dp, 0.7_dp, 0.9_dp]
   type(ad_test_result) :: result

   call ad_test_uniform(x, result)
   call assert_true(result%status == ad_success, 'status')
   call assert_true(result%n == 5, 'sample size')
   call assert_close(result%statistic, 0.2685490104415731_dp, 2.0e-14_dp, 'statistic')
   call assert_close(result%p_value, 0.9603953677778795_dp, 2.0e-14_dp, 'p-value')

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

end program test_uniform_api
