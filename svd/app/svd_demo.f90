! SPDX-License-Identifier: GPL-2.0-or-later
program svd_demo
   use svd
   implicit none

   real(dp) :: a(3, 3)
   type(trlan_eigen_result) :: fit

   a = 0.0_dp
   a(1, 1) = 4.0_dp
   a(2, 2) = 2.0_dp
   a(3, 3) = -1.0_dp
   fit = trlan_eigen(a, 2)
   if (fit%info /= 0) error stop "svd_demo: eigensolver failed"
   print '(a,*(1x,f8.4))', 'leading eigenvalues:', fit%d
end program svd_demo
