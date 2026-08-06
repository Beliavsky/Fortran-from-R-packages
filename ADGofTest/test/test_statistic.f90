program test_statistic
   use adgoftest, only : dp, ad_statistic, ad_success
   implicit none

   real(dp), parameter :: x1(5) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
   real(dp), parameter :: x2(5) = [0.9_dp, 0.05_dp, 0.7_dp, 0.2_dp, 0.4_dp]
   real(dp) :: value
   integer :: status

   value = ad_statistic(x1, status=status)
   call assert_true(status == ad_success, 'x1 status')
   call assert_close(value, 1.4644741743512215_dp, 2.0e-14_dp, 'x1 statistic')

   value = ad_statistic(x2, status=status)
   call assert_true(status == ad_success, 'x2 status')
   call assert_close(value, 0.2685490104415731_dp, 2.0e-14_dp, 'x2 statistic')

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

end program test_statistic
