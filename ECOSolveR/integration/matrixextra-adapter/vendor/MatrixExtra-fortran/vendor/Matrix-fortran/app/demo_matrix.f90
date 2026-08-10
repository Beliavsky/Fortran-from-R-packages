! SPDX-License-Identifier: GPL-3.0-only
program demo_matrix
   use matrix, only : dp, csr_matrix, csr_from_dense, csr_matvec, solve_linear, &
      near_positive_definite, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), x(:,:), adjusted(:,:), y(:)
   type(csr_matrix) :: sparse_a
   integer :: info

   a = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
   b = reshape([1.0_dp, 2.0_dp], [2, 1])
   call solve_linear(a, b, x, info)
   if (info /= matrix_success) error stop 'dense solve failed'
   print '(a,2(1x,f10.6))', 'solution:', x(:, 1)

   call csr_from_dense(a, sparse_a)
   call csr_matvec(sparse_a, [1.0_dp, 2.0_dp], y, info)
   print '(a,2(1x,f10.6))', 'sparse product:', y

   a = reshape([1.0_dp, 1.2_dp, 1.2_dp, 1.0_dp], [2, 2])
   call near_positive_definite(a, adjusted, info, corr=.true.)
   if (info /= matrix_success) error stop 'near_positive_definite failed'
   print '(a)', 'nearest positive-definite correlation matrix:'
   print '(2(1x,f10.6))', adjusted(1, :)
   print '(2(1x,f10.6))', adjusted(2, :)
end program demo_matrix
