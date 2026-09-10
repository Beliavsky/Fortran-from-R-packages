! SPDX-License-Identifier: GPL-2.0-or-later
program basic_dwt
   use wavethresh
   implicit none
   real(dp) :: x(16)
   real(dp), allocatable :: y(:)
   type(wd_t) :: object
   integer :: i
   x = [(sin(0.3_dp * real(i, dp)), i = 1, 16)]
   object = wd(x, filter_number=4.0_dp, family="DaubLeAsymm")
   y = wr_wd(object)
   print *, object%nlevels, maxval(abs(y - x))
end program basic_dwt
