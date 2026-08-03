! SPDX-License-Identifier: GPL-2.0-or-later
program density_cdf_example
   use gnorm
   implicit none
   real(dp) :: x
   integer :: i

   print '(a)', '     x             density             cdf'
   do i = -4, 4
      x = 0.5_dp * real(i, dp)
      print '(f7.2,2es20.10)', x, dgnorm(x, alpha=sqrt(2.0_dp), beta=2.0_dp), &
         pgnorm(x, alpha=sqrt(2.0_dp), beta=2.0_dp)
   end do
end program density_cdf_example
