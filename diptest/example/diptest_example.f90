! SPDX-License-Identifier: GPL-2.0-or-later
program diptest_example
   use diptest_api, only : dp, dip_test, dip_test_result, dip_ok
   implicit none
   type(dip_test_result) :: result
   real(dp) :: x(8), uniforms(8, 40)
   integer :: i, j
   x = [-3.0_dp, -2.9_dp, -2.8_dp, -2.7_dp, 2.7_dp, 2.8_dp, 2.9_dp, 3.0_dp]
   call dip_test(x, result)
   if (result%status /= dip_ok) error stop 'dip test failed'
   print '(a,f12.8,a,f12.8)', 'Dip = ', result%statistic, ', table p = ', result%p_value
   ! Fixed caller-supplied samples illustrate the reproducible API, not a good RNG.
   do j = 1, size(uniforms, 2)
      do i = 1, size(uniforms, 1)
         uniforms(i, j) = real(modulo(97 * i + 73 * j + 11 * i * j, 509), dp) / 509.0_dp
      end do
   end do
   call dip_test(x, result, uniform_samples=uniforms)
   if (result%status /= dip_ok) error stop 'simulated dip test failed'
   print '(a,f12.8,a,i0)', 'Fixed-sample p = ', result%p_value, ', replicates = ', result%replicates
end program diptest_example
