! SPDX-License-Identifier: GPL-2.0-or-later
program dense_svd
   use svd
   implicit none

   real(dp) :: a(4, 3)
   type(propack_svd_result) :: fit

   a = reshape([3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [4, 3])
   fit = propack_svd(a, 2)
   if (fit%info /= 0) error stop "dense_svd: solver failed"
   print '(a,*(1x,f8.4))', 'leading singular values:', fit%d
end program dense_svd
