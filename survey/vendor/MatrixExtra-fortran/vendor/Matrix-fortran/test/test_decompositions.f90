! SPDX-License-Identifier: GPL-3.0-only
program test_decompositions
   use matrix, only : dp, solve_linear, inverse_matrix, determinant, log_determinant, &
      cholesky_factor, cholesky_solve, qr_factor, least_squares, symmetric_eigen, &
      singular_value_decomposition, rank_matrix, pseudoinverse, ldlt_factor, ldlt_solve, &
      eye, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), x(:,:), ainv(:,:), l(:,:), q(:,:), r(:,:)
   real(dp), allocatable :: values(:), vectors(:,:), u(:,:), s(:), vt(:,:), ap(:,:), d(:)
   real(dp) :: det, logdet, sign_det
   integer :: info, rnk

   a = reshape([4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], [2, 2])
   b = reshape([1.0_dp, 2.0_dp], [2, 1])
   call solve_linear(a, b, x, info)
   call check(info == matrix_success, 'solve status')
   call check_close(matmul(a, x), b, 1.0e-12_dp, 'solve residual')
   call inverse_matrix(a, ainv, info)
   call check_close(matmul(a, ainv), eye(2), 1.0e-12_dp, 'inverse')
   call determinant(a, det, info)
   call check(abs(det - 10.0_dp) < 1.0e-12_dp, 'determinant')
   call log_determinant(a, logdet, sign_det, info)
   call check(abs(logdet - log(10.0_dp)) < 1.0e-12_dp .and. sign_det > 0.0_dp, 'logdet')

   a = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
   call cholesky_factor(a, l, info)
   call check(info == matrix_success, 'cholesky')
   call check_close(matmul(l, transpose(l)), a, 1.0e-12_dp, 'cholesky reconstruction')
   call cholesky_solve(l, b, x, info)
   call check_close(matmul(a, x), b, 1.0e-12_dp, 'cholesky solve')

   a = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [3, 2])
   call qr_factor(a, q, r, info)
   call check(info == matrix_success, 'qr')
   call check_close(matmul(q, r), a, 1.0e-11_dp, 'qr reconstruction')
   call check_close(matmul(transpose(q), q), eye(2), 1.0e-11_dp, 'q orthogonal')
   b = reshape([1.0_dp, 2.0_dp, 2.0_dp], [3, 1])
   call least_squares(a, b, x, info)
   call check_close(matmul(transpose(a), matmul(a, x) - b), reshape([0.0_dp, 0.0_dp], [2, 1]), &
      1.0e-10_dp, 'least squares normal equations')

   a = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
   call symmetric_eigen(a, values, vectors, info)
   call check(info == matrix_success, 'eigen')
   call check_close(matmul(a, vectors), matmul(vectors, diagonal(values)), 1.0e-10_dp, 'eigen equation')

   a = reshape([3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [3, 2])
   call singular_value_decomposition(a, u, s, vt, info)
   call check(info == matrix_success, 'svd')
   call check_close(matmul(u, matmul(diagonal(s), vt)), a, 1.0e-10_dp, 'svd reconstruction')
   rnk = rank_matrix(a, info=info)
   call check(rnk == 2, 'rank')
   call pseudoinverse(a, ap, info)
   call check_close(matmul(a, matmul(ap, a)), a, 1.0e-10_dp, 'pseudoinverse')

   a = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
   call ldlt_factor(a, l, d, info)
   call check(info == matrix_success, 'ldlt')
   call check_close(matmul(l, matmul(diagonal(d), transpose(l))), a, 1.0e-12_dp, 'ldlt reconstruction')
   call ldlt_solve(l, d, b(:2, :), x, info)
   call check_close(matmul(a, x), b(:2, :), 1.0e-11_dp, 'ldlt solve')

   print '(a)', 'test_decompositions: PASS'
contains
   function diagonal(v) result(mat)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: mat(:,:)
      integer :: i
      allocate(mat(size(v), size(v)), source=0.0_dp)
      do i = 1, size(v)
         mat(i, i) = v(i)
      end do
   end function diagonal

   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not. condition) then
         print '(a)', 'FAIL: ' // name
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tol, name)
      real(dp), intent(in) :: actual(:,:), expected(:,:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: name
      call check(all(shape(actual) == shape(expected)) .and. maxval(abs(actual - expected)) <= tol, name)
   end subroutine check_close
end program test_decompositions
