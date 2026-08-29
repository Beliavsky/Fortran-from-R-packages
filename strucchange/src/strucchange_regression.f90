! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_regression
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix, least_squares_svd, symmetric_eigen
   implicit none
   private
   public :: inverse_crossprod
   public :: ols_fit
   public :: root_matrix
contains
   subroutine ols_fit(x, y, beta, residuals, rss, rank, info)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: beta(:), residuals(:)
      real(dp), intent(out) :: rss
      integer, intent(out) :: rank, info
      integer :: n, k

      n = size(x, 1)
      k = size(x, 2)
      if (size(y) /= n) then
         allocate(beta(0), residuals(0))
         rss = 0.0_dp
         rank = 0
         info = -1
         return
      end if
      allocate(beta(k), residuals(n))
      if (k == 0) then
         residuals = y
         rss = sum(residuals ** 2)
         rank = 0
         info = 0
         return
      end if
      call least_squares_svd(x, y, beta, rank, info)
      if (info /= 0) then
         beta = 0.0_dp
         residuals = y
         rss = huge(1.0_dp)
         return
      end if
      residuals = y - matmul(x, beta)
      rss = sum(residuals ** 2)
   end subroutine ols_fit

   subroutine inverse_crossprod(x, inverse, info)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: inverse(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: xtx(:, :)

      allocate(xtx(size(x, 2), size(x, 2)))
      xtx = matmul(transpose(x), x)
      call inverse_matrix(xtx, inverse, info)
   end subroutine inverse_crossprod

   subroutine root_matrix(a, root, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
      real(dp) :: tolerance
      integer :: j, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(root(0, 0))
         info = -1
         return
      end if
      call symmetric_eigen(a, values, vectors, info)
      if (info /= 0) then
         allocate(root(0, 0))
         return
      end if
      tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))
      if (any(values < -tolerance)) then
         allocate(root(0, 0))
         info = -2
         return
      end if
      where (values < 0.0_dp)
         values = 0.0_dp
      end where
      allocate(scaled(n, n), root(n, n))
      scaled = vectors
      do j = 1, n
         scaled(:, j) = scaled(:, j) * sqrt(values(j))
      end do
      root = matmul(scaled, transpose(vectors))
      root = 0.5_dp * (root + transpose(root))
   end subroutine root_matrix
end module strucchange_regression
