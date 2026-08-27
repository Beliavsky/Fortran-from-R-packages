! SPDX-License-Identifier: MIT
module r_linalg
   use iso_fortran_env, only : real64
   use la_lapack, only : dgeev => geev, dgels => gels, dgelss => gelss
   use la_lapack, only : dgees => gees
   use la_lapack, only : dgebal => gebal
   use la_lapack, only : dgeqrf => geqrf, dorgqr => orgqr
   use la_lapack, only : dgesdd => gesdd
   use la_lapack_d, only : dgeqp3 => la_dgeqp3
   use la_lapack, only : dgesv => gesv
   use la_lapack, only : dposv => posv
   use la_lapack, only : dpotrf => potrf, dpotri => potri, dpotrs => potrs
   use la_lapack, only : dsyev => syev
   implicit none
   private

   integer, parameter, public :: r_linalg_dp = real64
   integer, parameter, public :: r_linalg_invalid_shape = -1001

   public :: cholesky_factor
   public :: balance_matrix
   public :: complex_schur
   public :: complex_thin_svd
   public :: full_svd
   public :: general_complex_eigen
   public :: general_real_eigen
   public :: general_real_eigenvalues
   public :: inverse_matrix
   public :: least_squares
   public :: least_squares_svd
   public :: numerical_rank
   public :: rank_revealing_qr
   public :: real_schur
   public :: singular_values
   public :: solve_cholesky
   public :: solve_spd
   public :: solve_system
   public :: spectral_radius
   public :: spd_inverse_logdet
   public :: symmetric_eigen
   public :: symmetric_eigenvalues
   public :: symmetrize
   public :: thin_svd
   public :: thin_qr

   interface solve_system
      module procedure solve_system_vector
      module procedure solve_system_matrix
      module procedure solve_system_complex_vector
      module procedure solve_system_complex_matrix
   end interface solve_system

   interface balance_matrix
      module procedure balance_matrix_real
      module procedure balance_matrix_complex
   end interface balance_matrix

   interface least_squares
      module procedure least_squares_vector
      module procedure least_squares_matrix
   end interface least_squares

   interface least_squares_svd
      module procedure least_squares_svd_vector
      module procedure least_squares_svd_matrix
   end interface least_squares_svd

   interface solve_spd
      module procedure solve_spd_vector
      module procedure solve_spd_matrix
   end interface solve_spd

   interface solve_cholesky
      module procedure solve_cholesky_vector
      module procedure solve_cholesky_matrix
   end interface solve_cholesky

contains

   subroutine thin_qr(a, q, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: q(:, :)
      integer, intent(out) :: info
      real(real64), allocatable :: matrix(:, :), tau(:), work(:)
      real(real64) :: work_query(1)
      integer :: k, lwork, m, n

      m = size(a, 1)
      n = size(a, 2)
      k = min(m, n)
      allocate(q(m, k))
      if (k == 0) then
         info = 0
         return
      end if

      allocate(matrix(m, n), tau(k))
      matrix = a
      call dgeqrf(m, n, matrix, m, tau, work_query, -1, info)
      if (info /= 0) then
         q = 0.0_real64
         return
      end if
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      call dgeqrf(m, n, matrix, m, tau, work, lwork, info)
      if (info /= 0) then
         q = 0.0_real64
         return
      end if

      deallocate(work)
      call dorgqr(m, k, k, matrix, m, tau, work_query, -1, info)
      if (info /= 0) then
         q = 0.0_real64
         return
      end if
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      call dorgqr(m, k, k, matrix, m, tau, work, lwork, info)
      if (info == 0) then
         q = matrix(:, 1:k)
      else
         q = 0.0_real64
      end if
   end subroutine thin_qr

   pure function symmetrize(a) result(s)
      real(real64), intent(in) :: a(:, :)
      real(real64) :: s(size(a, 1), size(a, 2))

      s = 0.5_real64 * (a + transpose(a))
   end function symmetrize

   subroutine solve_system_vector(a, b, x, info)
      real(real64), intent(in) :: a(:, :), b(:)
      real(real64), intent(out) :: x(:)
      integer, intent(out) :: info
      real(real64), allocatable :: rhs(:, :)
      integer, allocatable :: pivots(:)
      real(real64), allocatable :: work(:, :)
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         info = 0
         return
      end if

      allocate(work(n, n), rhs(n, 1), pivots(n))
      work = a
      rhs(:, 1) = b
      call dgesv(n, 1, work, n, pivots, rhs, n, info)
      if (info == 0) x = rhs(:, 1)
   end subroutine solve_system_vector

   subroutine solve_system_matrix(a, b, x, info)
      real(real64), intent(in) :: a(:, :), b(:, :)
      real(real64), intent(out) :: x(:, :)
      integer, intent(out) :: info
      integer, allocatable :: pivots(:)
      real(real64), allocatable :: work(:, :), rhs(:, :)
      integer :: n, nrhs

      n = size(a, 1)
      nrhs = size(b, 2)
      if (size(a, 2) /= n .or. size(b, 1) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0 .or. nrhs == 0) then
         info = 0
         return
      end if

      allocate(work(n, n), rhs(n, nrhs), pivots(n))
      work = a
      rhs = b
      call dgesv(n, nrhs, work, n, pivots, rhs, n, info)
      if (info == 0) x = rhs
   end subroutine solve_system_matrix

   subroutine solve_system_complex_vector(a, b, x, info)
      complex(real64), intent(in) :: a(:, :), b(:)
      complex(real64), intent(out) :: x(:)
      integer, intent(out) :: info
      complex(real64), allocatable :: rhs(:, :), work(:, :)
      integer, allocatable :: pivots(:)
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         info = 0
         return
      end if

      allocate(work(n, n), rhs(n, 1), pivots(n))
      work = a
      rhs(:, 1) = b
      call dgesv(n, 1, work, n, pivots, rhs, n, info)
      if (info == 0) x = rhs(:, 1)
   end subroutine solve_system_complex_vector

   subroutine solve_system_complex_matrix(a, b, x, info)
      complex(real64), intent(in) :: a(:, :), b(:, :)
      complex(real64), intent(out) :: x(:, :)
      integer, intent(out) :: info
      complex(real64), allocatable :: rhs(:, :), work(:, :)
      integer, allocatable :: pivots(:)
      integer :: n, nrhs

      n = size(a, 1)
      nrhs = size(b, 2)
      if (size(a, 2) /= n .or. size(b, 1) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0 .or. nrhs == 0) then
         info = 0
         return
      end if

      allocate(work(n, n), rhs(n, nrhs), pivots(n))
      work = a
      rhs = b
      call dgesv(n, nrhs, work, n, pivots, rhs, n, info)
      if (info == 0) x = rhs
   end subroutine solve_system_complex_matrix

   subroutine inverse_matrix(a, inverse, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: inverse(:, :)
      integer, intent(out) :: info
      real(real64), allocatable :: identity(:, :)
      integer :: i, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(inverse(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(inverse(n, n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(identity(n, n))
      identity = 0.0_real64
      do i = 1, n
         identity(i, i) = 1.0_real64
      end do
      inverse = 0.0_real64
      call solve_system(a, identity, inverse, info)
   end subroutine inverse_matrix

   subroutine least_squares_vector(a, b, x, info)
      real(real64), intent(in) :: a(:, :), b(:)
      real(real64), intent(out) :: x(:)
      integer, intent(out) :: info
      real(real64), allocatable :: b_matrix(:, :), x_matrix(:, :)
      integer :: m, n

      m = size(a, 1)
      n = size(a, 2)
      if (size(b) /= m .or. size(x) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      allocate(b_matrix(m, 1), x_matrix(n, 1))
      b_matrix(:, 1) = b
      call least_squares_matrix(a, b_matrix, x_matrix, info)
      if (info == 0) x = x_matrix(:, 1)
   end subroutine least_squares_vector

   subroutine least_squares_matrix(a, b, x, info)
      real(real64), intent(in) :: a(:, :), b(:, :)
      real(real64), intent(out) :: x(:, :)
      integer, intent(out) :: info
      real(real64), allocatable :: matrix(:, :), rhs(:, :), work(:)
      real(real64) :: work_query(1)
      integer :: lda, ldb, lwork, m, n, nrhs

      m = size(a, 1)
      n = size(a, 2)
      nrhs = size(b, 2)
      if (size(b, 1) /= m) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (m == 0 .or. n == 0 .or. nrhs == 0) then
         x = 0.0_real64
         info = 0
         return
      end if

      lda = max(1, m)
      ldb = max(m, n)
      allocate(matrix(lda, n), rhs(ldb, nrhs))
      matrix = a
      rhs = 0.0_real64
      rhs(1:m, :) = b
      call dgels('N', m, n, nrhs, matrix, lda, rhs, ldb, work_query, -1, info)
      if (info /= 0) return
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      matrix = a
      rhs = 0.0_real64
      rhs(1:m, :) = b
      call dgels('N', m, n, nrhs, matrix, lda, rhs, ldb, work, lwork, info)
      if (info == 0) x = rhs(1:n, :)
   end subroutine least_squares_matrix

   subroutine least_squares_svd_vector(a, b, x, rank, info, rcond)
      real(real64), intent(in) :: a(:, :), b(:)
      real(real64), intent(out) :: x(:)
      integer, intent(out) :: rank, info
      real(real64), intent(in), optional :: rcond
      real(real64), allocatable :: b_matrix(:, :), x_matrix(:, :)
      integer :: m, n

      m = size(a, 1)
      n = size(a, 2)
      if (size(b) /= m .or. size(x) /= n) then
         rank = 0
         info = r_linalg_invalid_shape
         return
      end if
      allocate(b_matrix(m, 1), x_matrix(n, 1))
      b_matrix(:, 1) = b
      if (present(rcond)) then
         call least_squares_svd_matrix(a, b_matrix, x_matrix, rank, info, rcond)
      else
         call least_squares_svd_matrix(a, b_matrix, x_matrix, rank, info)
      end if
      if (info == 0) x = x_matrix(:, 1)
   end subroutine least_squares_svd_vector

   subroutine least_squares_svd_matrix(a, b, x, rank, info, rcond)
      real(real64), intent(in) :: a(:, :), b(:, :)
      real(real64), intent(out) :: x(:, :)
      integer, intent(out) :: rank, info
      real(real64), intent(in), optional :: rcond
      real(real64), allocatable :: matrix(:, :), rhs(:, :), values(:), work(:)
      real(real64) :: rank_tolerance, work_query(1)
      integer :: lda, ldb, lwork, m, n, nrhs

      m = size(a, 1)
      n = size(a, 2)
      nrhs = size(b, 2)
      rank = 0
      if (size(b, 1) /= m) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (m == 0 .or. n == 0 .or. nrhs == 0) then
         x = 0.0_real64
         info = 0
         return
      end if

      rank_tolerance = -1.0_real64
      if (present(rcond)) rank_tolerance = rcond
      lda = max(1, m)
      ldb = max(m, n)
      allocate(matrix(lda, n), rhs(ldb, nrhs), values(min(m, n)))
      matrix = a
      rhs = 0.0_real64
      rhs(1:m, :) = b
      call dgelss(m, n, nrhs, matrix, lda, rhs, ldb, values, rank_tolerance, rank, work_query, -1, info)
      if (info /= 0) return
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      matrix = a
      rhs = 0.0_real64
      rhs(1:m, :) = b
      call dgelss(m, n, nrhs, matrix, lda, rhs, ldb, values, rank_tolerance, rank, work, lwork, info)
      if (info == 0) x = rhs(1:n, :)
   end subroutine least_squares_svd_matrix

   subroutine solve_spd_vector(a, b, x, info, upper)
      real(real64), intent(in) :: a(:, :), b(:)
      real(real64), intent(out) :: x(:)
      integer, intent(out) :: info
      logical, intent(in), optional :: upper
      real(real64), allocatable :: factor(:, :), rhs(:, :)
      character(len=1) :: triangle
      logical :: use_upper
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         info = 0
         return
      end if

      use_upper = .false.
      if (present(upper)) use_upper = upper
      triangle = merge('U', 'L', use_upper)
      allocate(factor(n, n), rhs(n, 1))
      factor = symmetrize(a)
      rhs(:, 1) = b
      call dposv(triangle, n, 1, factor, n, rhs, n, info)
      if (info == 0) x = rhs(:, 1)
   end subroutine solve_spd_vector

   subroutine solve_spd_matrix(a, b, x, info, upper)
      real(real64), intent(in) :: a(:, :), b(:, :)
      real(real64), intent(out) :: x(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: upper
      real(real64), allocatable :: factor(:, :), rhs(:, :)
      character(len=1) :: triangle
      logical :: use_upper
      integer :: n, nrhs

      n = size(a, 1)
      nrhs = size(b, 2)
      if (size(a, 2) /= n .or. size(b, 1) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0 .or. nrhs == 0) then
         info = 0
         return
      end if

      use_upper = .false.
      if (present(upper)) use_upper = upper
      triangle = merge('U', 'L', use_upper)
      allocate(factor(n, n), rhs(n, nrhs))
      factor = symmetrize(a)
      rhs = b
      call dposv(triangle, n, nrhs, factor, n, rhs, n, info)
      if (info == 0) x = rhs
   end subroutine solve_spd_matrix

   subroutine solve_cholesky_vector(factor, b, x, info, upper)
      real(real64), intent(in) :: factor(:, :), b(:)
      real(real64), intent(out) :: x(:)
      integer, intent(out) :: info
      logical, intent(in), optional :: upper
      real(real64), allocatable :: rhs(:, :)
      character(len=1) :: triangle
      logical :: use_upper
      integer :: n

      n = size(factor, 1)
      if (size(factor, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         info = 0
         return
      end if

      use_upper = .false.
      if (present(upper)) use_upper = upper
      triangle = merge('U', 'L', use_upper)
      allocate(rhs(n, 1))
      rhs(:, 1) = b
      call dpotrs(triangle, n, 1, factor, n, rhs, n, info)
      if (info == 0) x = rhs(:, 1)
   end subroutine solve_cholesky_vector

   subroutine solve_cholesky_matrix(factor, b, x, info, upper)
      real(real64), intent(in) :: factor(:, :), b(:, :)
      real(real64), intent(out) :: x(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: upper
      real(real64), allocatable :: rhs(:, :)
      character(len=1) :: triangle
      logical :: use_upper
      integer :: n, nrhs

      n = size(factor, 1)
      nrhs = size(b, 2)
      if (size(factor, 2) /= n .or. size(b, 1) /= n) then
         info = r_linalg_invalid_shape
         return
      end if
      if (size(x, 1) /= n .or. size(x, 2) /= nrhs) then
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0 .or. nrhs == 0) then
         info = 0
         return
      end if

      use_upper = .false.
      if (present(upper)) use_upper = upper
      triangle = merge('U', 'L', use_upper)
      allocate(rhs(n, nrhs))
      rhs = b
      call dpotrs(triangle, n, nrhs, factor, n, rhs, n, info)
      if (info == 0) x = rhs
   end subroutine solve_cholesky_matrix

   subroutine rank_revealing_qr(a, pivots, rank, info, tolerance)
      real(real64), intent(in) :: a(:, :)
      integer, allocatable, intent(out) :: pivots(:)
      integer, intent(out) :: rank, info
      real(real64), intent(in), optional :: tolerance
      real(real64), allocatable :: matrix(:, :), tau(:), work(:)
      real(real64) :: threshold, work_query(1)
      integer :: i, lwork, m, n

      m = size(a, 1)
      n = size(a, 2)
      allocate(pivots(n))
      pivots = 0
      rank = 0
      if (m == 0 .or. n == 0) then
         do i = 1, n
            pivots(i) = i
         end do
         info = 0
         return
      end if

      allocate(matrix(m, n), tau(min(m, n)))
      matrix = a
      call dgeqp3(m, n, matrix, m, pivots, tau, work_query, -1, info)
      if (info /= 0) return
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      matrix = a
      pivots = 0
      call dgeqp3(m, n, matrix, m, pivots, tau, work, lwork, info)
      if (info /= 0) return

      threshold = real(max(m, n), real64) * epsilon(1.0_real64) * abs(matrix(1, 1))
      if (present(tolerance)) threshold = max(0.0_real64, tolerance)
      do i = 1, min(m, n)
         if (abs(matrix(i, i)) > threshold) rank = rank + 1
      end do
   end subroutine rank_revealing_qr

   subroutine thin_svd(a, u, singular_values, vt, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: u(:, :), singular_values(:), vt(:, :)
      integer, intent(out) :: info
      integer, allocatable :: integer_work(:)
      real(real64), allocatable :: matrix(:, :), work(:)
      real(real64) :: work_query(1)
      integer :: lwork, m, n, r

      m = size(a, 1)
      n = size(a, 2)
      r = min(m, n)
      allocate(u(m, r), singular_values(r), vt(r, n))
      if (r == 0) then
         info = 0
         return
      end if

      allocate(matrix(m, n), integer_work(8 * r))
      matrix = a
      call dgesdd('S', m, n, matrix, m, singular_values, u, m, vt, r, work_query, -1, integer_work, info)
      if (info /= 0) return
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      matrix = a
      call dgesdd('S', m, n, matrix, m, singular_values, u, m, vt, r, work, lwork, integer_work, info)
   end subroutine thin_svd

   subroutine complex_thin_svd(a, u, singular_values, vt, info)
      complex(real64), intent(in) :: a(:, :)
      complex(real64), allocatable, intent(out) :: u(:, :), vt(:, :)
      real(real64), allocatable, intent(out) :: singular_values(:)
      integer, intent(out) :: info
      complex(real64), allocatable :: ac(:, :), work(:)
      complex(real64) :: query(1)
      integer, allocatable :: iw(:)
      real(real64), allocatable :: rw(:)
      integer :: lwork, m, n, r, rwsize

      m = size(a, 1)
      n = size(a, 2)
      r = min(m, n)
      allocate(u(m, r), singular_values(r), vt(r, n))
      if (r == 0) then
         info = 0
         return
      end if

      rwsize = max(1, 5 * r * r + 7 * r)
      allocate(ac(m, n), rw(rwsize), iw(8 * r))
      ac = a
      call dgesdd('S', m, n, ac, m, singular_values, u, m, vt, r, query, -1, rw, iw, info)
      if (info /= 0) return
      lwork = max(1, int(real(query(1), real64)))
      allocate(work(lwork))
      ac = a
      call dgesdd('S', m, n, ac, m, singular_values, u, m, vt, r, work, lwork, rw, iw, info)
   end subroutine complex_thin_svd

   subroutine full_svd(a, u, singular_values, vt, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: u(:, :), singular_values(:), vt(:, :)
      integer, intent(out) :: info
      integer, allocatable :: integer_work(:)
      real(real64), allocatable :: matrix(:, :), work(:)
      real(real64) :: work_query(1)
      integer :: lwork, m, n, r

      m = size(a, 1)
      n = size(a, 2)
      r = min(m, n)
      allocate(u(m, m), singular_values(r), vt(n, n))
      if (r == 0) then
         u = 0.0_real64
         vt = 0.0_real64
         info = 0
         return
      end if

      allocate(matrix(m, n), integer_work(8 * r))
      matrix = a
      call dgesdd('A', m, n, matrix, m, singular_values, u, m, vt, n, work_query, -1, integer_work, info)
      if (info /= 0) return
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      matrix = a
      call dgesdd('A', m, n, matrix, m, singular_values, u, m, vt, n, work, lwork, integer_work, info)
   end subroutine full_svd

   subroutine singular_values(a, values, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: values(:)
      integer, intent(out) :: info
      integer, allocatable :: integer_work(:)
      real(real64), allocatable :: matrix(:, :), work(:)
      real(real64) :: unused_u(1, 1), unused_vt(1, 1), work_query(1)
      integer :: lwork, m, n, r

      m = size(a, 1)
      n = size(a, 2)
      r = min(m, n)
      allocate(values(r))
      if (r == 0) then
         info = 0
         return
      end if

      allocate(matrix(m, n), integer_work(8 * r))
      matrix = a
      call dgesdd('N', m, n, matrix, m, values, unused_u, 1, unused_vt, 1, work_query, -1, integer_work, info)
      if (info /= 0) return
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      matrix = a
      call dgesdd('N', m, n, matrix, m, values, unused_u, 1, unused_vt, 1, work, lwork, integer_work, info)
   end subroutine singular_values

   subroutine numerical_rank(a, rank, info, tolerance)
      real(real64), intent(in) :: a(:, :)
      integer, intent(out) :: rank, info
      real(real64), intent(in), optional :: tolerance
      real(real64), allocatable :: values(:)
      real(real64) :: threshold
      integer :: m, n

      m = size(a, 1)
      n = size(a, 2)
      call singular_values(a, values, info)
      if (info /= 0) then
         rank = 0
         return
      end if
      if (size(values) == 0) then
         rank = 0
         return
      end if

      threshold = real(max(m, n), real64) * epsilon(1.0_real64) * max(1.0_real64, values(1))
      if (present(tolerance)) threshold = tolerance
      rank = count(values > threshold)
   end subroutine numerical_rank

   subroutine balance_matrix_real(a, balanced, scale, ilo, ihi, info, job)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: balanced(:, :), scale(:)
      integer, intent(out) :: ilo, ihi, info
      character(len=1), intent(in), optional :: job
      character(len=1) :: operation
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(balanced(0, 0), scale(0))
         ilo = 1
         ihi = 0
         info = r_linalg_invalid_shape
         return
      end if
      allocate(balanced(n, n), scale(n))
      if (n == 0) then
         ilo = 1
         ihi = 0
         info = 0
         return
      end if

      operation = 'B'
      if (present(job)) operation = job
      balanced = a
      call dgebal(operation, n, balanced, n, ilo, ihi, scale, info)
   end subroutine balance_matrix_real

   subroutine balance_matrix_complex(a, balanced, scale, ilo, ihi, info, job)
      complex(real64), intent(in) :: a(:, :)
      complex(real64), allocatable, intent(out) :: balanced(:, :)
      real(real64), allocatable, intent(out) :: scale(:)
      integer, intent(out) :: ilo, ihi, info
      character(len=1), intent(in), optional :: job
      character(len=1) :: operation
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(balanced(0, 0), scale(0))
         ilo = 1
         ihi = 0
         info = r_linalg_invalid_shape
         return
      end if
      allocate(balanced(n, n), scale(n))
      if (n == 0) then
         ilo = 1
         ihi = 0
         info = 0
         return
      end if

      operation = 'B'
      if (present(job)) operation = job
      balanced = a
      call dgebal(operation, n, balanced, n, ilo, ihi, scale, info)
   end subroutine balance_matrix_complex

   subroutine general_real_eigenvalues(a, wr, wi, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: wr(:), wi(:)
      integer, intent(out) :: info
      real(real64), allocatable :: matrix(:, :), work(:)
      real(real64) :: vl(1, 1), vr(1, 1), work_query(1)
      integer :: lwork, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(wr(0), wi(0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(wr(n), wi(n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(matrix(n, n))
      matrix = a
      call dgeev('N', 'N', n, matrix, n, wr, wi, vl, 1, vr, 1, work_query, -1, info)
      if (info /= 0) return
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      matrix = a
      call dgeev('N', 'N', n, matrix, n, wr, wi, vl, 1, vr, 1, work, lwork, info)
   end subroutine general_real_eigenvalues

   subroutine general_real_eigen(a, wr, wi, right_vectors, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: wr(:), wi(:), right_vectors(:, :)
      integer, intent(out) :: info
      real(real64), allocatable :: matrix(:, :), work(:)
      real(real64) :: left_vectors(1, 1), work_query(1)
      integer :: lwork, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(wr(0), wi(0), right_vectors(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(wr(n), wi(n), right_vectors(n, n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(matrix(n, n))
      matrix = a
      call dgeev('N', 'V', n, matrix, n, wr, wi, left_vectors, 1, right_vectors, n, work_query, -1, info)
      if (info /= 0) return
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      matrix = a
      call dgeev('N', 'V', n, matrix, n, wr, wi, left_vectors, 1, right_vectors, n, work, lwork, info)
   end subroutine general_real_eigen

   subroutine general_complex_eigen(a, values, right_vectors, info)
      complex(real64), intent(in) :: a(:, :)
      complex(real64), allocatable, intent(out) :: values(:), right_vectors(:, :)
      integer, intent(out) :: info
      complex(real64), allocatable :: matrix(:, :), work(:)
      complex(real64) :: left_vectors(1, 1), work_query(1)
      real(real64), allocatable :: real_work(:)
      integer :: lwork, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(values(0), right_vectors(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(values(n), right_vectors(n, n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(matrix(n, n), real_work(max(1, 2 * n)))
      matrix = a
      call dgeev('N', 'V', n, matrix, n, values, left_vectors, 1, right_vectors, n, work_query, -1, real_work, info)
      if (info /= 0) return
      lwork = max(1, ceiling(real(work_query(1), real64)))
      allocate(work(lwork))
      matrix = a
      call dgeev('N', 'V', n, matrix, n, values, left_vectors, 1, right_vectors, n, work, lwork, real_work, info)
   end subroutine general_complex_eigen

   subroutine complex_schur(a, t, values, q, info)
      complex(real64), intent(in) :: a(:, :)
      complex(real64), allocatable, intent(out) :: t(:, :), values(:), q(:, :)
      integer, intent(out) :: info
      complex(real64), allocatable :: work(:)
      complex(real64) :: query(1)
      real(real64), allocatable :: rwork(:)
      logical, allocatable :: bwork(:)
      integer :: lwork, n, sdim

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(t(0, 0), values(0), q(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(t(n, n), values(n), q(n, n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(rwork(max(1, n)), bwork(max(1, n)))
      t = a
      call dgees('V', 'N', zselect_none, n, t, n, sdim, values, q, n, query, -1, rwork, bwork, info)
      if (info /= 0) return
      lwork = max(1, ceiling(real(query(1), real64)))
      allocate(work(lwork))
      t = a
      call dgees('V', 'N', zselect_none, n, t, n, sdim, values, q, n, work, lwork, rwork, bwork, info)
   end subroutine complex_schur

   pure logical function zselect_none(value)
      complex(real64), intent(in) :: value

      zselect_none = real(value, real64) > huge(1.0_real64)
   end function zselect_none

   subroutine spectral_radius(a, radius, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), intent(out) :: radius
      integer, intent(out) :: info
      real(real64), allocatable :: rp(:), ip(:)
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         radius = huge(1.0_real64)
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         radius = 0.0_real64
         info = 0
         return
      end if

      call general_real_eigenvalues(a, rp, ip, info)
      if (info == 0) then
         radius = maxval(hypot(rp, ip))
      else
         radius = huge(1.0_real64)
      end if
   end subroutine spectral_radius

   subroutine real_schur(a, t, wr, wi, q, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: t(:, :), wr(:), wi(:), q(:, :)
      integer, intent(out) :: info
      real(real64), allocatable :: work(:)
      real(real64) :: work_query(1)
      logical, allocatable :: selected_work(:)
      integer :: lwork, n, selected_count

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(t(0, 0), wr(0), wi(0), q(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      allocate(t(n, n), wr(n), wi(n), q(n, n))
      if (n == 0) then
         info = 0
         return
      end if

      allocate(selected_work(max(1, n)))
      t = a
      call dgees('V', 'N', schur_select_none, n, t, n, selected_count, wr, wi, q, n, work_query, -1, selected_work, info)
      if (info /= 0) return
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      t = a
      call dgees('V', 'N', schur_select_none, n, t, n, selected_count, wr, wi, q, n, work, lwork, selected_work, info)
   end subroutine real_schur

   pure logical function schur_select_none(real_part, imaginary_part)
      real(real64), intent(in) :: real_part, imaginary_part

      schur_select_none = real_part + imaginary_part < -huge(1.0_real64)
   end function schur_select_none

   subroutine symmetric_eigen(a, values, vectors, info, descending)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: values(:), vectors(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: descending
      real(real64), allocatable :: ascending(:), matrix(:, :), work(:)
      real(real64) :: work_query(1)
      logical :: reverse_order
      integer :: j, lwork, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(values(0), vectors(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         allocate(values(0), vectors(0, 0))
         info = 0
         return
      end if

      allocate(ascending(n), matrix(n, n))
      matrix = symmetrize(a)
      call dsyev('V', 'U', n, matrix, n, ascending, work_query, -1, info)
      if (info /= 0) then
         allocate(values(0), vectors(0, 0))
         return
      end if
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      call dsyev('V', 'U', n, matrix, n, ascending, work, lwork, info)
      if (info /= 0) then
         allocate(values(0), vectors(0, 0))
         return
      end if

      reverse_order = .false.
      if (present(descending)) reverse_order = descending
      allocate(values(n), vectors(n, n))
      if (reverse_order) then
         do j = 1, n
            values(j) = ascending(n - j + 1)
            vectors(:, j) = matrix(:, n - j + 1)
         end do
      else
         values = ascending
         vectors = matrix
      end if
   end subroutine symmetric_eigen

   subroutine symmetric_eigenvalues(a, values, info, descending)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: values(:)
      integer, intent(out) :: info
      logical, intent(in), optional :: descending
      real(real64), allocatable :: vectors(:, :)

      call symmetric_eigen(a, values, vectors, info, descending)
   end subroutine symmetric_eigenvalues

   subroutine cholesky_factor(a, factor, info, upper)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: factor(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: upper
      character(len=1) :: triangle
      logical :: use_upper
      integer :: i, j, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(factor(0, 0))
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         allocate(factor(0, 0))
         info = 0
         return
      end if

      use_upper = .false.
      if (present(upper)) use_upper = upper
      triangle = merge('U', 'L', use_upper)
      allocate(factor(n, n))
      factor = symmetrize(a)
      call dpotrf(triangle, n, factor, n, info)
      if (info /= 0) return
      if (use_upper) then
         do j = 1, n
            do i = j + 1, n
               factor(i, j) = 0.0_real64
            end do
         end do
      else
         do j = 1, n
            do i = 1, j - 1
               factor(i, j) = 0.0_real64
            end do
         end do
      end if
   end subroutine cholesky_factor

   subroutine spd_inverse_logdet(a, inverse, logdet, info)
      real(real64), intent(in) :: a(:, :)
      real(real64), allocatable, intent(out) :: inverse(:, :)
      real(real64), intent(out) :: logdet
      integer, intent(out) :: info
      integer :: i, j, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(inverse(0, 0))
         logdet = huge(1.0_real64)
         info = r_linalg_invalid_shape
         return
      end if
      if (n == 0) then
         allocate(inverse(0, 0))
         logdet = 0.0_real64
         info = 0
         return
      end if

      allocate(inverse(n, n))
      inverse = symmetrize(a)
      call dpotrf('L', n, inverse, n, info)
      if (info /= 0) then
         logdet = huge(1.0_real64)
         return
      end if
      logdet = 0.0_real64
      do i = 1, n
         logdet = logdet + 2.0_real64 * log(inverse(i, i))
      end do
      call dpotri('L', n, inverse, n, info)
      if (info /= 0) return
      do j = 1, n
         do i = 1, j - 1
            inverse(i, j) = inverse(j, i)
         end do
      end do
   end subroutine spd_inverse_logdet

end module r_linalg
