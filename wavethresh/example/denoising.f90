! SPDX-License-Identifier: GPL-2.0-or-later
program denoising
   use wavethresh
   implicit none
   real(dp) :: x(32)
   real(dp), allocatable :: y(:)
   type(wd_t) :: object
   type(wd_t) :: shrunk
   integer :: i
   x = [(sin(0.18_dp * real(i, dp)) + 0.15_dp * cos(1.2_dp * real(i, dp)), i = 1, 32)]
   object = wd(x, filter_number=1.0_dp, family="DaubExPhase")
   shrunk = threshold_wd(object, threshold_type="soft", policy="universal")
   y = wr_wd(shrunk)
   print *, dof(shrunk), y(1:4)
end program denoising
