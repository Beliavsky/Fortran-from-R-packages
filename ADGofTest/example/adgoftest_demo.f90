module demo_cdf_functions
   use adgoftest, only : dp
   implicit none
contains
   pure function standard_normal_cdf(x) result(probability)
      real(dp), intent(in) :: x
      real(dp) :: probability
      probability = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function standard_normal_cdf
end module demo_cdf_functions

program adgoftest_demo
   use adgoftest, only : dp, ad_test, ad_test_result, ad_success
   use demo_cdf_functions, only : standard_normal_cdf
   implicit none

   real(dp), parameter :: sample(10) = [-1.3_dp, -0.8_dp, -0.4_dp, -0.1_dp, 0.0_dp, &
      0.2_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.4_dp]
   type(ad_test_result) :: result

   call ad_test(sample, standard_normal_cdf, result, clamp_p_value=.true.)
   if (result%status /= ad_success) then
      write(*, '(a)') trim(result%message)
      error stop 1
   end if

   write(*, '(a,i0)') 'n         = ', result%n
   write(*, '(a,f12.6)') 'AD        = ', result%statistic
   write(*, '(a,f12.6)') 'p-value   = ', result%p_value
end program adgoftest_demo
