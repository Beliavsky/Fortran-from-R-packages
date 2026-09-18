! SPDX-License-Identifier: GPL-3.0-only
! Matrix/statistical helpers for the dccmidas translation.
module dccmidas_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_linalg, only : determinant, inverse_matrix, signed_log_determinant
   use dccmidas_types, only : DCCMIDAS_SUCCESS, DCCMIDAS_INVALID_INPUT, DCCMIDAS_LINALG_ERROR
   implicit none
   private

   public :: identity_matrix, diagonal_matrix, outer_product
   public :: sample_covariance_rows, sample_covariance_series
   public :: inv, det, logdet_quadratic, covariance_to_correlation
   public :: lower_triangle_squared_distance

contains

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Matrix order; must be nonnegative.
      real(dp) :: a(max(0, n), max(0, n))
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function diagonal_matrix(x) result(a)
      real(dp), intent(in) :: x(:) !! Diagonal entries of the returned square matrix.
      real(dp) :: a(size(x), size(x))
      integer :: i

      a = 0.0_dp
      do i = 1, size(x)
         a(i, i) = x(i)
      end do
   end function diagonal_matrix

   pure function outer_product(x) result(a)
      real(dp), intent(in) :: x(:) !! Vector whose outer product is required.
      real(dp) :: a(size(x), size(x))
      integer :: i, j

      do j = 1, size(x)
         do i = 1, size(x)
            a(i, j) = x(i) * x(j)
         end do
      end do
   end function outer_product

   pure subroutine sample_covariance_rows(x, cov, status)
      real(dp), intent(in) :: x(:, :) !! Observations by variables matrix with shape `(n, k)`.
      real(dp), intent(out) :: cov(:, :) !! Sample covariance matrix with shape `(k, k)` and `n-1` denominator.
      integer, intent(out) :: status !! Zero on success; invalid-input code for incompatible dimensions or too few rows.
      real(dp) :: mean(size(x, 2)), centered(size(x, 2))
      integer :: i, n, k

      n = size(x, 1)
      k = size(x, 2)
      cov = 0.0_dp
      if (size(cov, 1) /= k .or. size(cov, 2) /= k .or. n < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      mean = sum(x, dim=1) / real(n, dp)
      do i = 1, n
         centered = x(i, :) - mean
         cov = cov + outer_product(centered)
      end do
      cov = cov / real(n - 1, dp)
      status = DCCMIDAS_SUCCESS
   end subroutine sample_covariance_rows

   pure subroutine sample_covariance_series(x, cov, status)
      real(dp), intent(in) :: x(:, :) !! Variables by time matrix with shape `(k, n)`.
      real(dp), intent(out) :: cov(:, :) !! Sample covariance across time with shape `(k, k)`.
      integer, intent(out) :: status !! Zero on success; invalid-input code for incompatible dimensions or too few times.
      real(dp) :: mean(size(x, 1)), centered(size(x, 1))
      integer :: t, n, k

      k = size(x, 1)
      n = size(x, 2)
      cov = 0.0_dp
      if (size(cov, 1) /= k .or. size(cov, 2) /= k .or. n < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      mean = sum(x, dim=2) / real(n, dp)
      do t = 1, n
         centered = x(:, t) - mean
         cov = cov + outer_product(centered)
      end do
      cov = cov / real(n - 1, dp)
      status = DCCMIDAS_SUCCESS
   end subroutine sample_covariance_series

   pure subroutine inv(x, inverse, status)
      real(dp), intent(in) :: x(:, :) !! Square matrix to invert.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Allocated inverse matrix, or an empty matrix on invalid shape.
      integer, intent(out) :: status !! Zero on success or a dccmidas linear-algebra/input error code.
      integer :: info

      if (size(x, 1) /= size(x, 2)) then
         allocate(inverse(0, 0))
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call inverse_matrix(x, inverse, info)
      if (info == 0) then
         status = DCCMIDAS_SUCCESS
      else
         status = DCCMIDAS_LINALG_ERROR
      end if
   end subroutine inv

   pure subroutine det(x, value, status)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose determinant is requested.
      real(dp), intent(out) :: value !! Determinant of `x`, or zero on failure.
      integer, intent(out) :: status !! Zero on success or a dccmidas linear-algebra/input error code.
      integer :: info

      if (size(x, 1) /= size(x, 2)) then
         value = 0.0_dp
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call determinant(x, value, info)
      if (info == 0) then
         status = DCCMIDAS_SUCCESS
      else
         status = DCCMIDAS_LINALG_ERROR
      end if
   end subroutine det

   pure subroutine logdet_quadratic(a, x, logdet, quadratic, status)
      real(dp), intent(in) :: a(:, :) !! Positive-definite covariance/correlation matrix.
      real(dp), intent(in) :: x(:) !! Observation vector conformable with `a`.
      real(dp), intent(out) :: logdet !! Natural logarithm of the determinant of `a`.
      real(dp), intent(out) :: quadratic !! Quadratic form `x' inv(a) x`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: sign
      integer :: info

      logdet = ieee_value(0.0_dp, ieee_quiet_nan)
      quadratic = huge(1.0_dp)
      if (size(a, 1) /= size(a, 2) .or. size(a, 1) /= size(x)) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call signed_log_determinant(a, sign, logdet, info)
      if (info /= 0 .or. sign <= 0.0_dp .or. .not. ieee_is_finite(logdet)) then
         status = DCCMIDAS_LINALG_ERROR
         return
      end if
      call inverse_matrix(a, inverse, info)
      if (info /= 0) then
         status = DCCMIDAS_LINALG_ERROR
         return
      end if
      quadratic = dot_product(x, matmul(inverse, x))
      if (.not. ieee_is_finite(quadratic)) then
         status = DCCMIDAS_LINALG_ERROR
         return
      end if
      status = DCCMIDAS_SUCCESS
   end subroutine logdet_quadratic

   pure subroutine covariance_to_correlation(q, r, status)
      real(dp), intent(in) :: q(:, :) !! Symmetric covariance-like matrix with positive diagonal.
      real(dp), intent(out) :: r(:, :) !! Correlation-normalized matrix with unit diagonal.
      integer, intent(out) :: status !! Zero on success or invalid-input code when normalization is impossible.
      real(dp) :: scale(size(q, 1))
      integer :: i, j, n

      n = size(q, 1)
      r = 0.0_dp
      if (size(q, 2) /= n .or. size(r, 1) /= n .or. size(r, 2) /= n) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      do i = 1, n
         if (q(i, i) <= 0.0_dp .or. .not. ieee_is_finite(q(i, i))) then
            status = DCCMIDAS_LINALG_ERROR
            return
         end if
         scale(i) = sqrt(q(i, i))
      end do
      do j = 1, n
         do i = 1, n
            r(i, j) = q(i, j) / (scale(i) * scale(j))
         end do
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine covariance_to_correlation

   pure real(dp) function lower_triangle_squared_distance(a, b) result(value)
      real(dp), intent(in) :: a(:, :) !! First square matrix.
      real(dp), intent(in) :: b(:, :) !! Second square matrix of the same shape.
      integer :: i, j, n

      value = 0.0_dp
      n = min(size(a, 1), size(a, 2), size(b, 1), size(b, 2))
      do j = 1, n
         do i = j, n
            value = value + (a(i, j) - b(i, j)) ** 2
         end do
      end do
   end function lower_triangle_squared_distance

end module dccmidas_matrix
