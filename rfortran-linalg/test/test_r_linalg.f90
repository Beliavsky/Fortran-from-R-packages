! SPDX-License-Identifier: MIT
program test_r_linalg
   use iso_fortran_env, only : real64
   use r_linalg, only : cholesky_factor, r_linalg_invalid_shape
   use r_linalg, only : solve_system, spd_inverse_logdet
   use r_linalg, only : symmetric_eigen, symmetric_eigenvalues, symmetrize
   implicit none

   call test_intrinsic_norm2
   call test_symmetrize
   call test_solve
   call test_eigen
   call test_cholesky_and_inverse
   call test_shapes
   print *, 'test_r_linalg: PASS'

contains

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop message
   end subroutine check

   subroutine test_intrinsic_norm2
      real(real64) :: x(2)
      intrinsic :: norm2

      x = [3.0_real64, 4.0_real64]
      call check(abs(norm2(x) - 5.0_real64) < 1.0e-14_real64, 'norm2')
   end subroutine test_intrinsic_norm2

   subroutine test_symmetrize
      real(real64) :: a(2, 2), expected(2, 2)

      a = reshape([1.0_real64, 2.0_real64, 4.0_real64, 3.0_real64], [2, 2])
      expected = reshape([1.0_real64, 3.0_real64, 3.0_real64, 3.0_real64], [2, 2])
      call check(maxval(abs(symmetrize(a) - expected)) < 1.0e-14_real64, 'symmetrize')
   end subroutine test_symmetrize

   subroutine test_solve
      real(real64) :: a(2, 2), b(2), b_matrix(2, 2), x(2), x_matrix(2, 2)
      integer :: info

      a = reshape([4.0_real64, 1.0_real64, 1.0_real64, 3.0_real64], [2, 2])
      b = [1.0_real64, 2.0_real64]
      call solve_system(a, b, x, info)
      call check(info == 0, 'solve info')
      call check(maxval(abs(matmul(a, x) - b)) < 1.0e-13_real64, 'solve residual')
      b_matrix = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], [2, 2])
      call solve_system(a, b_matrix, x_matrix, info)
      call check(info == 0, 'matrix solve info')
      call check(maxval(abs(matmul(a, x_matrix) - b_matrix)) < 1.0e-13_real64, 'matrix solve residual')
   end subroutine test_solve

   subroutine test_eigen
      real(real64), allocatable :: values(:), values_only(:), vectors(:, :)
      real(real64) :: a(2, 2), reconstructed(2, 2)
      integer :: info

      a = reshape([2.0_real64, 1.0_real64, 1.0_real64, 2.0_real64], [2, 2])
      call symmetric_eigen(a, values, vectors, info, descending=.true.)
      call check(info == 0, 'eigen info')
      call check(maxval(abs(values - [3.0_real64, 1.0_real64])) < 1.0e-13_real64, 'eigenvalues')
      reconstructed = matmul(vectors, matmul(diagonal(values), transpose(vectors)))
      call check(maxval(abs(reconstructed - a)) < 1.0e-13_real64, 'eigenvectors')
      call symmetric_eigenvalues(a, values_only, info)
      call check(info == 0, 'eigenvalues-only info')
      call check(maxval(abs(values_only - [1.0_real64, 3.0_real64])) < 1.0e-13_real64, 'eigenvalues only')
   end subroutine test_eigen

   subroutine test_cholesky_and_inverse
      real(real64), allocatable :: factor(:, :), inverse(:, :)
      real(real64) :: a(2, 2), identity(2, 2), logdet
      integer :: info

      a = reshape([4.0_real64, 1.0_real64, 1.0_real64, 3.0_real64], [2, 2])
      identity = 0.0_real64
      identity(1, 1) = 1.0_real64
      identity(2, 2) = 1.0_real64
      call cholesky_factor(a, factor, info)
      call check(info == 0, 'cholesky info')
      call check(maxval(abs(matmul(factor, transpose(factor)) - a)) < 1.0e-13_real64, 'cholesky')
      call spd_inverse_logdet(a, inverse, logdet, info)
      call check(info == 0, 'inverse info')
      call check(maxval(abs(matmul(a, inverse) - identity)) < 1.0e-13_real64, 'inverse')
      call check(abs(logdet - log(11.0_real64)) < 1.0e-13_real64, 'log determinant')
   end subroutine test_cholesky_and_inverse

   subroutine test_shapes
      real(real64) :: a(2, 3), b(2), x(2)
      integer :: info

      a = 0.0_real64
      b = 0.0_real64
      call solve_system(a, b, x, info)
      call check(info == r_linalg_invalid_shape, 'invalid solve shape')
   end subroutine test_shapes

   pure function diagonal(x) result(a)
      real(real64), intent(in) :: x(:)
      real(real64) :: a(size(x), size(x))
      integer :: i

      a = 0.0_real64
      do i = 1, size(x)
         a(i, i) = x(i)
      end do
   end function diagonal

end program test_r_linalg
