! SPDX-License-Identifier: GPL-2.0-or-later
program moments_example
   use moments, only : dp, all_moments, skewness, kurtosis, geary
   implicit none

   real(dp) :: x(8)
   real(dp), allocatable :: mu(:)
   integer :: k

   x = [1.2_dp, 0.8_dp, 1.5_dp, 2.1_dp, 0.7_dp, 1.0_dp, 3.2_dp, 1.4_dp]
   mu = all_moments(x, order_max=4, central=.true.)

   do k = 0, 4
      write(*, '(a,i0,a,f12.6)') 'central moment ', k, ': ', mu(k + 1)
   end do
   write(*, '(a,f12.6)') 'skewness: ', skewness(x)
   write(*, '(a,f12.6)') 'Pearson kurtosis: ', kurtosis(x)
   write(*, '(a,f12.6)') 'Geary kurtosis: ', geary(x)
end program moments_example
