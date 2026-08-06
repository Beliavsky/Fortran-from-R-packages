module test_cdf_functions
   use adgoftest, only : dp
   implicit none
contains
   pure function standard_normal_cdf(value) result(probability)
      real(dp), intent(in) :: value
      real(dp) :: probability
      probability = 0.5_dp * erfc(-value / sqrt(2.0_dp))
   end function standard_normal_cdf
end module test_cdf_functions

program test_cdf_api
   use adgoftest, only : dp, ad_test, ad_test_result, ad_success
   use test_cdf_functions, only : standard_normal_cdf
   implicit none

   real(dp), parameter :: x(5) = [-1.2815515655446004_dp, -0.8416212335729143_dp, &
      -0.5244005127080409_dp, -0.2533471031357997_dp, 0.0_dp]
   type(ad_test_result) :: result

   call ad_test(x, standard_normal_cdf, result)
   call assert_true(result%status == ad_success, 'status')
   call assert_close(result%statistic, 1.4644741743512215_dp, 5.0e-12_dp, 'callback statistic')
   call assert_close(result%p_value, 0.1855700680847705_dp, 5.0e-12_dp, 'callback p-value')

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

end program test_cdf_api
