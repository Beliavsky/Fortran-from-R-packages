! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
program basic_roll
   use roll, only : dp, roll_mean, roll_sd
   implicit none
   real(dp) :: x(8)
   real(dp), allocatable :: mean_x(:), sd_x(:)

   x = [1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp]
   call roll_mean(x, 4, mean_x, min_obs=1)
   call roll_sd(x, 4, sd_x, min_obs=2)
   write (*, '(a,8f8.3)') 'rolling mean:', mean_x
   write (*, '(a,8f8.3)') 'rolling sd:  ', sd_x
end program basic_roll
