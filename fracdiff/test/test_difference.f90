! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program test_difference
   use fracdiff, only : dp, fd_ok, diffseries, diffseries_direct, fractional_weights, &
                        polynomial_multiply, minimum_ar_root_modulus
   use test_support
   implicit none

   real(dp) :: x(8), fast(8), direct(8), expected(8), weights(8)
   real(dp) :: linear(20), dlinear(20), expected_linear(20)
   real(dp) :: a(3), b(2), product(4), root_modulus
   integer :: status, i

   x = [3.0_dp,1.0_dp,4.0_dp,1.0_dp,5.0_dp,9.0_dp,2.0_dp,6.0_dp]
   expected = [-0.875_dp,-2.55125_dp,1.29073125_dp,-2.5307589375_dp, &
               2.39267401421875_dp,5.182070143666563_dp, &
               -3.618196411177042_dp,2.349594546622669_dp]
   call diffseries(x,0.37_dp,fast,status)
   call assert_true(status==fd_ok,"FFT differencing status")
   call diffseries_direct(x,0.37_dp,direct,status)
   call assert_true(status==fd_ok,"direct differencing status")
   call assert_vector_close(fast,expected,2.0e-12_dp,"FFT differencing reference")
   call assert_vector_close(direct,expected,2.0e-14_dp,"direct differencing reference")
   call assert_vector_close(fast,direct,2.0e-12_dp,"FFT/direct agreement")

   linear=[(real(i,dp),i=1,20)]
   expected_linear(1)=-9.5_dp
   expected_linear(2:)=1.0_dp
   call diffseries(linear,1.0_dp,dlinear,status)
   call assert_vector_close(dlinear,expected_linear,2.0e-12_dp,"integer differencing identity")

   call fractional_weights(8,0.37_dp,weights,status)
   call assert_close(weights(1),1.0_dp,0.0_dp,"weight zero")
   call assert_close(weights(2),-0.37_dp,1.0e-15_dp,"weight one")

   a=[1.0_dp,2.0_dp,3.0_dp]; b=[-1.0_dp,4.0_dp]
   call polynomial_multiply(a,b,product)
   call assert_vector_close(product,[-1.0_dp,2.0_dp,5.0_dp,12.0_dp],1.0e-15_dp,"polynomial multiply")

   root_modulus=minimum_ar_root_modulus([0.5_dp],status)
   call assert_true(status==0,"AR root status")
   call assert_close(root_modulus,2.0_dp,1.0e-12_dp,"AR root modulus")

   write(*,'(a)') "test_difference: PASS"

end program test_difference
