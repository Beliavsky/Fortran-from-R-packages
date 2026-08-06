! SPDX-License-Identifier: GPL-3.0-only
program nearest_correlation
   use matrix, only : dp, near_positive_definite, symmetric_eigen, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), nearest(:,:), values(:), vectors(:,:)
   integer :: info, i

   a = reshape([1.0_dp, 0.9_dp, 1.1_dp, 0.9_dp, 1.0_dp, -0.8_dp, &
                1.1_dp, -0.8_dp, 1.0_dp], [3, 3])
   call near_positive_definite(a, nearest, info, corr=.true.)
   if (info /= matrix_success) error stop 'nearPD failed'
   call symmetric_eigen(nearest, values, vectors, info)
   print '(a)', 'nearest correlation matrix:'
   do i = 1, 3
      print '(3(1x,f10.6))', nearest(i, :)
   end do
   print '(a,*(1x,es12.4))', 'eigenvalues:', values
end program nearest_correlation
