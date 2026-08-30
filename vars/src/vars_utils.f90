! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_utils
   use r_kinds, only : dp
   use r_linalg, only : full_svd, inverse_matrix
   implicit none
   private

   public :: identity_matrix, kronecker_product, commutation_matrix
   public :: determinant_logabs, covariance_matrix, residual_standardize
   public :: duplication_matrix, vech_lower, null_space, matrix_rank_svd
   public :: sample_mean, sample_sd, quantile_linear, center_columns
   public :: safe_inverse, trace_matrix, outer_product

contains

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x), size(y))
      integer :: i, j

      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_product

   function kronecker_product(a, b) result(c)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable :: c(:, :)
      integer :: i, j, mb, nb

      mb = size(b, 1)
      nb = size(b, 2)
      allocate(c(size(a, 1) * mb, size(a, 2) * nb))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            c((i - 1) * mb + 1:i * mb, (j - 1) * nb + 1:j * nb) = a(i, j) * b
         end do
      end do
   end function kronecker_product

   pure function commutation_matrix(n) result(kmat)
      integer, intent(in) :: n
      real(dp) :: kmat(n * n, n * n)
      integer :: i, j, source, target

      kmat = 0.0_dp
      do j = 1, n
         do i = 1, n
            source = i + (j - 1) * n
            target = j + (i - 1) * n
            kmat(target, source) = 1.0_dp
         end do
      end do
   end function commutation_matrix

   subroutine determinant_logabs(a, logabsdet, sign_det, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: logabsdet
      integer, intent(out) :: sign_det, info
      real(dp), allocatable :: work(:, :), row_tmp(:)
      real(dp) :: pivot_abs, factor
      integer :: i, j, k, n, pivot

      n = size(a, 1)
      if (size(a, 2) /= n) then
         info = -1
         logabsdet = -huge(1.0_dp)
         sign_det = 0
         return
      end if
      allocate(work(n, n), row_tmp(n))
      work = a
      sign_det = 1
      logabsdet = 0.0_dp
      info = 0
      do k = 1, n
         pivot = k
         pivot_abs = abs(work(k, k))
         do i = k + 1, n
            if (abs(work(i, k)) > pivot_abs) then
               pivot_abs = abs(work(i, k))
               pivot = i
            end if
         end do
         if (pivot_abs <= tiny(1.0_dp)) then
            info = k
            sign_det = 0
            logabsdet = -huge(1.0_dp)
            return
         end if
         if (pivot /= k) then
            row_tmp = work(k, :)
            work(k, :) = work(pivot, :)
            work(pivot, :) = row_tmp
            sign_det = -sign_det
         end if
         if (work(k, k) < 0.0_dp) sign_det = -sign_det
         logabsdet = logabsdet + log(abs(work(k, k)))
         do i = k + 1, n
            factor = work(i, k) / work(k, k)
            work(i, k) = 0.0_dp
            do j = k + 1, n
               work(i, j) = work(i, j) - factor * work(k, j)
            end do
         end do
      end do
   end subroutine determinant_logabs

   subroutine safe_inverse(a, ainv, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: info

      call inverse_matrix(a, ainv, info)
   end subroutine safe_inverse

   pure function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: value
      integer :: i, n

      n = min(size(a, 1), size(a, 2))
      value = 0.0_dp
      do i = 1, n
         value = value + a(i, i)
      end do
   end function trace_matrix

   subroutine covariance_matrix(x, cov, denominator)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: cov(:, :)
      real(dp), intent(in), optional :: denominator
      real(dp), allocatable :: centered(:, :)
      real(dp) :: denom
      integer :: n

      n = size(x, 1)
      allocate(centered(size(x, 1), size(x, 2)))
      call center_columns(x, centered)
      denom = real(max(1, n - 1), dp)
      if (present(denominator)) denom = denominator
      cov = matmul(transpose(centered), centered) / denom
   end subroutine covariance_matrix

   subroutine center_columns(x, centered)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: centered(:, :)
      integer :: j

      do j = 1, size(x, 2)
         centered(:, j) = x(:, j) - sample_mean(x(:, j))
      end do
   end subroutine center_columns

   subroutine residual_standardize(x, z)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: z(:, :)
      real(dp) :: mu, sd
      integer :: j

      do j = 1, size(x, 2)
         mu = sample_mean(x(:, j))
         sd = sample_sd(x(:, j))
         if (sd > 0.0_dp) then
            z(:, j) = (x(:, j) - mu) / sd
         else
            z(:, j) = 0.0_dp
         end if
      end do
   end subroutine residual_standardize

   pure function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function sample_mean

   pure function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, mu

      if (size(x) <= 1) then
         value = 0.0_dp
      else
         mu = sample_mean(x)
         value = sqrt(sum((x - mu) ** 2) / real(size(x) - 1, dp))
      end if
   end function sample_sd

   function quantile_linear(x, probability) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: probability
      real(dp) :: value, h, frac, tmp
      real(dp), allocatable :: y(:)
      integer :: i, j, lo, hi, n

      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(y(n))
      y = x
      do i = 2, n
         tmp = y(i)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= tmp) exit
            y(j + 1) = y(j)
            j = j - 1
         end do
         y(j + 1) = tmp
      end do
      if (probability <= 0.0_dp) then
         value = y(1)
         return
      end if
      if (probability >= 1.0_dp) then
         value = y(n)
         return
      end if
      h = 1.0_dp + real(n - 1, dp) * probability
      lo = floor(h)
      hi = ceiling(h)
      frac = h - real(lo, dp)
      value = (1.0_dp - frac) * y(lo) + frac * y(hi)
   end function quantile_linear

   pure subroutine vech_lower(a, v)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: v(:)
      integer :: i, j, idx, n

      n = size(a, 1)
      idx = 0
      do j = 1, n
         do i = j, n
            idx = idx + 1
            v(idx) = a(i, j)
         end do
      end do
   end subroutine vech_lower

   subroutine duplication_matrix(n, d)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: d(:, :)
      integer :: i, j, idx, col1, col2

      allocate(d(n * n, n * (n + 1) / 2))
      d = 0.0_dp
      idx = 0
      do j = 1, n
         do i = j, n
            idx = idx + 1
            col1 = i + (j - 1) * n
            d(col1, idx) = 1.0_dp
            if (i /= j) then
               col2 = j + (i - 1) * n
               d(col2, idx) = 1.0_dp
            end if
         end do
      end do
   end subroutine duplication_matrix

   subroutine null_space(a, basis, info, tolerance)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: basis(:, :)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: u(:, :), s(:), vt(:, :)
      real(dp) :: tol
      integer :: i, n, rank

      n = size(a, 2)
      call full_svd(a, u, s, vt, info)
      if (info /= 0) then
         allocate(basis(0, 0))
         return
      end if
      if (size(s) == 0) then
         allocate(basis(n, n))
         basis = identity_matrix(n)
         return
      end if
      tol = real(max(size(a, 1), size(a, 2)), dp) * epsilon(1.0_dp) * max(1.0_dp, s(1))
      if (present(tolerance)) tol = tolerance
      rank = count(s > tol)
      allocate(basis(n, n - rank))
      if (n - rank == 0) return
      do i = 1, n - rank
         basis(:, i) = vt(rank + i, :)
      end do
   end subroutine null_space

   subroutine matrix_rank_svd(a, rank, info, tolerance)
      real(dp), intent(in) :: a(:, :)
      integer, intent(out) :: rank, info
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: u(:, :), s(:), vt(:, :)
      real(dp) :: tol

      call full_svd(a, u, s, vt, info)
      if (info /= 0) then
         rank = 0
         return
      end if
      if (size(s) == 0) then
         rank = 0
         return
      end if
      tol = real(max(size(a, 1), size(a, 2)), dp) * epsilon(1.0_dp) * max(1.0_dp, s(1))
      if (present(tolerance)) tol = tolerance
      rank = count(s > tol)
   end subroutine matrix_rank_svd

end module vars_utils
