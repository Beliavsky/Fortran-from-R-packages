! SPDX-License-Identifier: GPL-2.0-or-later
program wavethresh_demo
   use wavethresh
   implicit none
   real(dp) :: x(8)
   real(dp), allocatable :: y(:)
   type(wd_t) :: object
   x = [1.0_dp, 2.0_dp, 1.5_dp, -0.5_dp, 3.0_dp, 2.5_dp, 0.0_dp, -1.0_dp]
   object = wd(x, filter_number=2.0_dp, family="DaubExPhase")
   y = wr_wd(object)
   print '(a,es12.4)', 'maximum reconstruction error = ', maxval(abs(y - x))
end program wavethresh_demo
