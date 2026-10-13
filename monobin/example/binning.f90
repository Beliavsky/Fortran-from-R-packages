! SPDX-License-Identifier: GPL-3.0-or-later
program binning_example
   use monobin
   implicit none
   real(dp) :: x(300), y(300)
   type(bin_result) :: fit
   integer, allocatable :: group(:)
   real(dp), allocatable :: risk(:)
   integer :: i
   do i = 1, size(x)
      x(i) = real((i - 1) / 10, dp)
      y(i) = merge(1.0_dp, 0.0_dp, modulo(i * 17, 100) < 10 + int(2.0_dp * x(i)))
   end do
   call iso_bin(x, y, fit)
   if (fit%status /= bin_ok) error stop 'Binning failed'
   call apply_bins(fit, [2.0_dp, 15.0_dp, 50.0_dp], group, risk)
   print '(a,i0)', 'Complete bins: ', size(fit%bins) - fit%n_special
   print '(a,3i5)', 'Predicted groups: ', group
   print '(a,3f10.5)', 'Predicted event rates: ', risk
end program binning_example
