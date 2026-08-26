! SPDX-License-Identifier: MIT
module r_linalg
   use iso_fortran_env, only : real64
   use la_lapack, only : dgesv => gesv
   use la_lapack, only : dpotrf => potrf, dpotri => potri
   use la_lapack, only : dsyev => syev
   implicit none
   private

   integer, parameter, public :: r_linalg_dp = real64
   integer, parameter, public :: r_linalg_invalid_shape = -1001

   public :: cholesky_factor
   public :: solve_system
   public :: spd_inverse_logdet
   public :: symmetric_eigen
   public :: symmetric_eigenvalues
   public :: symmetrize

   interface solve_system
      module procedure solve_system_vector
      module procedure solve_system_matrix
   end interface solve_system

contains

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
