! SPDX-License-Identifier: GPL-3.0-only
program dense_solve
   use matrix, only : dp, solve_linear, determinant, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), x(:,:)
   real(dp) :: det
   integer :: info

   a = reshape([3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
   b = reshape([9.0_dp, 8.0_dp], [2, 1])
   call solve_linear(a, b, x, info)
   if (info /= matrix_success) error stop 'solve failed'
   call determinant(a, det, info)
   print '(a,2(1x,f10.5))', 'x =', x(:, 1)
   print '(a,1x,f10.5)', 'det(A) =', det
end program dense_solve
