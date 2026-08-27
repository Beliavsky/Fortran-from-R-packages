! SPDX-License-Identifier: MIT
program test_r_linalg
   use iso_fortran_env, only : real64
   use r_linalg, only : cholesky_factor, r_linalg_invalid_shape
   use r_linalg, only : balance_matrix, complex_schur
   use r_linalg, only : complex_thin_svd, thin_qr
   use r_linalg, only : full_svd, inverse_matrix, numerical_rank, singular_values
   use r_linalg, only : general_complex_eigen, general_real_eigen, general_real_eigenvalues
   use r_linalg, only : least_squares, least_squares_svd
   use r_linalg, only : solve_cholesky, solve_spd, solve_system
   use r_linalg, only : spd_inverse_logdet
   use r_linalg, only : spectral_radius
   use r_linalg, only : symmetric_eigen, symmetric_eigenvalues, symmetrize
   use r_linalg, only : rank_revealing_qr, thin_svd
   use r_linalg, only : real_schur
   implicit none

   call test_intrinsic_norm2
   call test_symmetrize
   call test_thin_qr
   call test_solve
   call test_complex_solve
   call test_general_inverse
   call test_least_squares
   call test_least_squares_svd
   call test_eigen
   call test_cholesky_and_inverse
   call test_spd_solves
   call test_rank_revealing_qr
   call test_thin_svd
   call test_complex_thin_svd
   call test_full_svd
   call test_singular_values_and_rank
   call test_spectral_radius
   call test_general_eigen
   call test_real_schur
   call test_complex_schur
   call test_balance
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

   subroutine test_thin_qr
      real(real64), allocatable :: q(:, :)
      real(real64) :: a(3, 2), identity2(2, 2)
      integer :: info

      a = reshape([1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 1.0_real64], [3, 2])
      identity2 = diagonal([1.0_real64, 1.0_real64])
      call thin_qr(a, q, info)
      call check(info == 0, 'thin QR info')
      call check(all(shape(q) == [3, 2]), 'thin QR shape')
      call check(maxval(abs(matmul(transpose(q), q) - identity2)) < 1.0e-13_real64, 'thin QR orthogonality')
      call check(maxval(abs(matmul(q, matmul(transpose(q), a)) - a)) < 1.0e-13_real64, 'thin QR column space')
   end subroutine test_thin_qr

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

   subroutine test_complex_solve
      complex(real64) :: a(2, 2), b(2), bm(2, 2), x(2), xm(2, 2)
      integer :: info

      a(1, 1) = cmplx(2.0_real64, 1.0_real64, real64)
      a(2, 1) = cmplx(1.0_real64, 0.0_real64, real64)
      a(1, 2) = cmplx(0.0_real64, -1.0_real64, real64)
      a(2, 2) = cmplx(3.0_real64, 0.5_real64, real64)
      b = [cmplx(1.0_real64, 2.0_real64, real64), cmplx(-1.0_real64, 0.5_real64, real64)]
      call solve_system(a, b, x, info)
      call check(info == 0, 'complex solve info')
      call check(maxval(abs(matmul(a, x) - b)) < 1.0e-13_real64, 'complex solve residual')
      bm(:, 1) = b
      bm(:, 2) = 2.0_real64 * b
      call solve_system(a, bm, xm, info)
      call check(info == 0, 'complex matrix solve info')
      call check(maxval(abs(matmul(a, xm) - bm)) < 1.0e-13_real64, 'complex matrix solve residual')
   end subroutine test_complex_solve

   subroutine test_general_inverse
      real(real64), allocatable :: inverse(:, :)
      real(real64) :: a(2, 2), identity(2, 2)
      integer :: info

      a = reshape([2.0_real64, 1.0_real64, 3.0_real64, 2.0_real64], [2, 2])
      identity = 0.0_real64
      identity(1, 1) = 1.0_real64
      identity(2, 2) = 1.0_real64
      call inverse_matrix(a, inverse, info)
      call check(info == 0, 'general inverse info')
      call check(maxval(abs(matmul(a, inverse) - identity)) < 1.0e-13_real64, 'general inverse')
   end subroutine test_general_inverse

   subroutine test_least_squares
      real(real64) :: a(3, 2), b(3), b_matrix(3, 2), x(2), x_matrix(2, 2)
      integer :: info

      a = reshape([1.0_real64, 1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, 2.0_real64], [3, 2])
      b = [1.0_real64, 3.0_real64, 5.0_real64]
      call least_squares(a, b, x, info)
      call check(info == 0, 'QR least-squares vector info')
      call check(maxval(abs(x - [1.0_real64, 2.0_real64])) < 1.0e-13_real64, 'QR least-squares vector')
      b_matrix(:, 1) = b
      b_matrix(:, 2) = 2.0_real64 * b
      call least_squares(a, b_matrix, x_matrix, info)
      call check(info == 0, 'QR least-squares matrix info')
      call check(maxval(abs(x_matrix(:, 1) - x)) < 1.0e-13_real64, 'QR least-squares first RHS')
      call check(maxval(abs(x_matrix(:, 2) - 2.0_real64 * x)) < 1.0e-13_real64, 'QR least-squares second RHS')
   end subroutine test_least_squares

   subroutine test_least_squares_svd
      real(real64) :: a(3, 2), b(3), x(2)
      integer :: info, rank

      a(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
      a(:, 2) = 2.0_real64 * a(:, 1)
      b = 5.0_real64 * a(:, 1)
      call least_squares_svd(a, b, x, rank, info)
      call check(info == 0, 'SVD least-squares info')
      call check(rank == 1, 'SVD least-squares rank')
      call check(maxval(abs(matmul(a, x) - b)) < 1.0e-13_real64, 'SVD least-squares residual')
      call check(maxval(abs(x - [1.0_real64, 2.0_real64])) < 1.0e-13_real64, 'SVD minimum-norm solution')
   end subroutine test_least_squares_svd

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

   subroutine test_spd_solves
      real(real64), allocatable :: lower(:, :), upper(:, :)
      real(real64) :: a(2, 2), b(2), b_matrix(2, 2), x(2), x_matrix(2, 2)
      integer :: info

      a = reshape([4.0_real64, 1.0_real64, 1.0_real64, 3.0_real64], [2, 2])
      b = [1.0_real64, 2.0_real64]
      b_matrix = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], [2, 2])
      call solve_spd(a, b, x, info)
      call check(info == 0, 'SPD vector solve info')
      call check(maxval(abs(matmul(a, x) - b)) < 1.0e-13_real64, 'SPD vector solve')
      call solve_spd(a, b_matrix, x_matrix, info, upper=.true.)
      call check(info == 0, 'SPD matrix solve info')
      call check(maxval(abs(matmul(a, x_matrix) - b_matrix)) < 1.0e-13_real64, 'SPD matrix solve')

      call cholesky_factor(a, lower, info)
      call solve_cholesky(lower, b, x, info)
      call check(info == 0, 'lower Cholesky solve info')
      call check(maxval(abs(matmul(a, x) - b)) < 1.0e-13_real64, 'lower Cholesky solve')
      call cholesky_factor(a, upper, info, upper=.true.)
      call solve_cholesky(upper, b_matrix, x_matrix, info, upper=.true.)
      call check(info == 0, 'upper Cholesky solve info')
      call check(maxval(abs(matmul(a, x_matrix) - b_matrix)) < 1.0e-13_real64, 'upper Cholesky solve')
   end subroutine test_spd_solves

   subroutine test_shapes
      real(real64) :: a(2, 3), b(2), x(2)
      real(real64), allocatable :: inverse(:, :)
      integer :: info

      a = 0.0_real64
      b = 0.0_real64
      call solve_system(a, b, x, info)
      call check(info == r_linalg_invalid_shape, 'invalid solve shape')
      call inverse_matrix(a, inverse, info)
      call check(info == r_linalg_invalid_shape, 'invalid inverse shape')
      call check(all(shape(inverse) == [0, 0]), 'invalid inverse result shape')
   end subroutine test_shapes

   subroutine test_rank_revealing_qr
      real(real64) :: a(3, 3)
      integer, allocatable :: pivots(:)
      integer :: info, rank

      a(:, 1) = [1.0_real64, 0.0_real64, 1.0_real64]
      a(:, 2) = [0.0_real64, 1.0_real64, 1.0_real64]
      a(:, 3) = a(:, 1) + a(:, 2)
      call rank_revealing_qr(a, pivots, rank, info)
      call check(info == 0, 'pivoted QR info')
      call check(rank == 2, 'pivoted QR rank')
      call check(size(pivots) == 3, 'pivoted QR permutation size')
      call check(all(sort_integers(pivots) == [1, 2, 3]), 'pivoted QR permutation')
   end subroutine test_rank_revealing_qr

   subroutine test_thin_svd
      real(real64), allocatable :: s(:), u(:, :), vt(:, :)
      real(real64) :: a(3, 2), reconstructed(3, 2)
      integer :: info

      a = reshape([3.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, 0.0_real64], [3, 2])
      call thin_svd(a, u, s, vt, info)
      call check(info == 0, 'thin SVD info')
      call check(all(shape(u) == [3, 2]), 'thin SVD U shape')
      call check(all(shape(vt) == [2, 2]), 'thin SVD VT shape')
      call check(maxval(abs(s - [3.0_real64, 2.0_real64])) < 1.0e-13_real64, 'thin SVD values')
      reconstructed = matmul(u, matmul(diagonal(s), vt))
      call check(maxval(abs(reconstructed - a)) < 1.0e-13_real64, 'thin SVD reconstruction')
   end subroutine test_thin_svd

   subroutine test_complex_thin_svd
      complex(real64), allocatable :: u(:, :), vt(:, :)
      complex(real64) :: a(2, 2), reconstructed(2, 2), sigma(2, 2)
      real(real64), allocatable :: s(:)
      integer :: info

      a = cmplx(0.0_real64, 0.0_real64, real64)
      a(1, 1) = cmplx(1.0_real64, 1.0_real64, real64)
      a(2, 2) = cmplx(2.0_real64, -1.0_real64, real64)
      call complex_thin_svd(a, u, s, vt, info)
      call check(info == 0, 'complex thin SVD info')
      sigma = cmplx(diagonal(s), 0.0_real64, real64)
      reconstructed = matmul(u, matmul(sigma, vt))
      call check(maxval(abs(reconstructed - a)) < 1.0e-13_real64, 'complex thin SVD reconstruction')
   end subroutine test_complex_thin_svd

   subroutine test_full_svd
      real(real64), allocatable :: s(:), u(:, :), vt(:, :)
      real(real64) :: a(3, 2), identity3(3, 3), reconstructed(3, 2)
      integer :: info

      a = reshape([3.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, 0.0_real64], [3, 2])
      call full_svd(a, u, s, vt, info)
      call check(info == 0, 'full SVD info')
      call check(all(shape(u) == [3, 3]), 'full SVD U shape')
      call check(all(shape(vt) == [2, 2]), 'full SVD VT shape')
      reconstructed = matmul(u(:, 1:2), matmul(diagonal(s), vt))
      call check(maxval(abs(reconstructed - a)) < 1.0e-13_real64, 'full SVD reconstruction')
      identity3 = diagonal([1.0_real64, 1.0_real64, 1.0_real64])
      call check(maxval(abs(matmul(transpose(u), u) - identity3)) < 1.0e-13_real64, 'full SVD orthogonal U')
   end subroutine test_full_svd

   subroutine test_singular_values_and_rank
      real(real64), allocatable :: values(:)
      real(real64) :: a(3, 3)
      integer :: info, rank

      a = 0.0_real64
      a(1, 1) = 3.0_real64
      a(2, 2) = 1.0e-8_real64
      call singular_values(a, values, info)
      call check(info == 0, 'singular values info')
      call check(maxval(abs(values - [3.0_real64, 1.0e-8_real64, 0.0_real64])) < 1.0e-14_real64, 'singular values')
      call numerical_rank(a, rank, info)
      call check(info == 0, 'numerical rank info')
      call check(rank == 2, 'default numerical rank')
      call numerical_rank(a, rank, info, tolerance=1.0e-7_real64)
      call check(info == 0, 'toleranced numerical rank info')
      call check(rank == 1, 'toleranced numerical rank')
   end subroutine test_singular_values_and_rank

   subroutine test_spectral_radius
      real(real64) :: a(2, 2), radius
      integer :: info

      a = reshape([0.0_real64, 2.0_real64, -2.0_real64, 0.0_real64], [2, 2])
      call spectral_radius(a, radius, info)
      call check(info == 0, 'spectral radius info')
      call check(abs(radius - 2.0_real64) < 1.0e-13_real64, 'spectral radius')
   end subroutine test_spectral_radius

   subroutine test_general_eigen
      real(real64), allocatable :: wr(:), wi(:), vr(:, :)
      complex(real64), allocatable :: values(:), complex_vectors(:, :)
      real(real64) :: a(2, 2)
      complex(real64) :: complex_a(2, 2), vector(2)
      integer :: info, j

      a = reshape([0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64], [2, 2])
      call general_real_eigenvalues(a, wr, wi, info)
      call check(info == 0, 'general real eigenvalues info')
      call check(maxval(abs(wr)) < 1.0e-13_real64, 'general real eigenvalue real parts')
      call check(maxval(abs(abs(wi) - 1.0_real64)) < 1.0e-13_real64, 'general real eigenvalue imaginary parts')
      call general_real_eigen(a, wr, wi, vr, info)
      call check(info == 0, 'general real eigen info')
      do j = 1, 2
         if (wi(j) <= 0.0_real64) cycle
         vector = cmplx(vr(:, j), vr(:, j + 1), real64)
         call check(maxval(abs(matmul(a, vector) - cmplx(wr(j), wi(j), real64) * vector)) < 1.0e-13_real64, 'real eigenvector')
      end do

      complex_a = cmplx(0.0_real64, 0.0_real64, real64)
      complex_a(1, 1) = cmplx(1.0_real64, 2.0_real64, real64)
      complex_a(1, 2) = cmplx(3.0_real64, -1.0_real64, real64)
      complex_a(2, 2) = cmplx(-2.0_real64, 0.5_real64, real64)
      call general_complex_eigen(complex_a, values, complex_vectors, info)
      call check(info == 0, 'general complex eigen info')
      do j = 1, 2
         vector = complex_vectors(:, j)
         call check(maxval(abs(matmul(complex_a, vector) - values(j) * vector)) < 1.0e-12_real64, 'complex eigenvector')
      end do
   end subroutine test_general_eigen

   subroutine test_real_schur
      real(real64), allocatable :: t(:, :), rp(:), ip(:), q(:, :)
      real(real64) :: a(2, 2), identity2(2, 2), reconstructed(2, 2)
      integer :: info

      a = reshape([0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64], [2, 2])
      identity2 = diagonal([1.0_real64, 1.0_real64])
      call real_schur(a, t, rp, ip, q, info)
      call check(info == 0, 'real Schur info')
      reconstructed = matmul(q, matmul(t, transpose(q)))
      call check(maxval(abs(reconstructed - a)) < 1.0e-13_real64, 'real Schur reconstruction')
      call check(maxval(abs(matmul(transpose(q), q) - identity2)) < 1.0e-13_real64, 'real Schur orthogonality')
      call check(maxval(abs(abs(ip) - 1.0_real64)) < 1.0e-13_real64, 'real Schur complex pair')
      call check(maxval(abs(rp)) < 1.0e-13_real64, 'real Schur real parts')
   end subroutine test_real_schur

   subroutine test_complex_schur
      complex(real64), allocatable :: t(:, :), values(:), q(:, :)
      complex(real64) :: a(2, 2), identity2(2, 2), reconstructed(2, 2)
      integer :: info

      a = cmplx(0.0_real64, 0.0_real64, real64)
      a(1, 1) = cmplx(1.0_real64, 2.0_real64, real64)
      a(1, 2) = cmplx(3.0_real64, -1.0_real64, real64)
      a(2, 1) = cmplx(-0.5_real64, 0.25_real64, real64)
      a(2, 2) = cmplx(-2.0_real64, 0.5_real64, real64)
      identity2 = cmplx(diagonal([1.0_real64, 1.0_real64]), 0.0_real64, real64)
      call complex_schur(a, t, values, q, info)
      call check(info == 0, 'complex Schur info')
      reconstructed = matmul(q, matmul(t, transpose(conjg(q))))
      call check(maxval(abs(reconstructed - a)) < 1.0e-12_real64, 'complex Schur reconstruction')
      call check(maxval(abs(matmul(transpose(conjg(q)), q) - identity2)) < 1.0e-13_real64, 'complex Schur unitary vectors')
      call check(maxval(abs(values - [t(1, 1), t(2, 2)])) < 1.0e-13_real64, 'complex Schur eigenvalues')
   end subroutine test_complex_schur

   subroutine test_balance
      real(real64), allocatable :: rb(:, :), scale(:)
      complex(real64), allocatable :: cb(:, :)
      real(real64) :: ra(2, 2), real_det
      complex(real64) :: ca(2, 2), complex_det
      integer :: ihi, ilo, info

      ra = reshape([1.0e-8_real64, 2.0_real64, 1.0_real64, 1.0e8_real64], [2, 2])
      call balance_matrix(ra, rb, scale, ilo, ihi, info)
      call check(info == 0, 'real balance info')
      call check(all(shape(rb) == [2, 2]) .and. size(scale) == 2, 'real balance shapes')
      real_det = ra(1, 1) * ra(2, 2) - ra(1, 2) * ra(2, 1)
      call check(abs(rb(1, 1) + rb(2, 2) - ra(1, 1) - ra(2, 2)) < 1.0e-7_real64, 'real balance trace')
      call check(abs(rb(1, 1) * rb(2, 2) - rb(1, 2) * rb(2, 1) - real_det) < 1.0e-12_real64, 'real balance determinant')

      ca = cmplx(ra, 0.0_real64, real64)
      ca(1, 2) = cmplx(1.0_real64, 2.0_real64, real64)
      call balance_matrix(ca, cb, scale, ilo, ihi, info)
      call check(info == 0, 'complex balance info')
      call check(all(shape(cb) == [2, 2]) .and. size(scale) == 2, 'complex balance shapes')
      complex_det = ca(1, 1) * ca(2, 2) - ca(1, 2) * ca(2, 1)
      call check(abs(cb(1, 1) + cb(2, 2) - ca(1, 1) - ca(2, 2)) < 1.0e-7_real64, 'complex balance trace')
      call check(abs(cb(1, 1) * cb(2, 2) - cb(1, 2) * cb(2, 1) - complex_det) < 1.0e-12_real64, 'complex balance determinant')
   end subroutine test_balance

   pure function diagonal(x) result(a)
      real(real64), intent(in) :: x(:)
      real(real64) :: a(size(x), size(x))
      integer :: i

      a = 0.0_real64
      do i = 1, size(x)
         a(i, i) = x(i)
      end do
   end function diagonal

   pure function sort_integers(x) result(sorted)
      integer, intent(in) :: x(:)
      integer :: sorted(size(x))
      integer :: i, j, value

      sorted = x
      do i = 2, size(sorted)
         value = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= value) exit
            sorted(j + 1) = sorted(j)
            j = j - 1
         end do
         sorted(j + 1) = value
      end do
   end function sort_integers

end program test_r_linalg
