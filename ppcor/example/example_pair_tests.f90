program example_pair_tests
   use ppcor, only : dp, ppcor_test_result, pcor_test, spcor_test, ppcor_spearman
   implicit none
   real(dp) :: x(20), y(20), z(20,2)
   type(ppcor_test_result) :: result
   integer :: i
   real(dp) :: t

   do i = 1, size(x)
      t = real(i,dp)
      z(i,1) = sin(0.27_dp*t)
      z(i,2) = cos(0.16_dp*t)
      x(i) = 0.6_dp*z(i,1) + sin(0.93_dp*t)
      y(i) = 0.5_dp*x(i) + 0.4_dp*z(i,2) + cos(0.61_dp*t)
   end do

   call pcor_test(x, y, z, result, ppcor_spearman)
   print '(a,f10.6,a,es11.4)', 'Spearman partial r = ', result%estimate, &
         ', p = ', result%p_value
   call spcor_test(x, y, z, result, ppcor_spearman)
   print '(a,f10.6,a,es11.4)', 'Spearman semi-partial r = ', result%estimate, &
         ', p = ', result%p_value
end program example_pair_tests
