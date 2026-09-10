! SPDX-License-Identifier: GPL-2.0-only
module marss_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   implicit none
   private
   public :: zscore
   public :: ldiag
   public :: identity_matrix
   public :: is_symmetric
   public :: count_observed
   public :: normal_quantile

   interface zscore
      module procedure zscore_vector
      module procedure zscore_rows
   end interface zscore

   interface ldiag
      module procedure ldiag_scalar
      module procedure ldiag_vector
   end interface ldiag

contains

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Order of the square identity matrix, must be nonnegative.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure logical function is_symmetric(a, tol)
      real(dp), intent(in) :: a(:, :) !! Matrix to test for symmetry.
      real(dp), intent(in), optional :: tol !! Absolute comparison tolerance, defaults to 1e-10.
      real(dp) :: eps

      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      if (size(a, 1) /= size(a, 2)) then
         is_symmetric = .false.
      else if (size(a) == 0) then
         is_symmetric = .true.
      else
         is_symmetric = maxval(abs(a - transpose(a))) <= eps
      end if
   end function is_symmetric

   pure integer function count_observed(y)
      real(dp), intent(in) :: y(:, :) !! Observation matrix whose NaNs represent missing values.
      integer :: i
      integer :: t

      count_observed = 0
      do t = 1, size(y, 2)
         do i = 1, size(y, 1)
            if (.not. ieee_is_nan(y(i, t))) count_observed = count_observed + 1
         end do
      end do
   end function count_observed

   pure function zscore_vector(x, mean_only) result(z)
      real(dp), intent(in) :: x(:) !! Input vector, NaNs are ignored when estimating mean and variance and preserved in output.
      logical, intent(in), optional :: mean_only !! If true, subtract only the sample mean without variance scaling.
      real(dp) :: z(size(x))
      real(dp) :: mean_x
      real(dp) :: sd_x
      real(dp) :: ss
      integer :: i
      integer :: n
      logical :: center_only

      mean_x = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            mean_x = mean_x + x(i)
            n = n + 1
         end if
      end do
      if (n > 0) mean_x = mean_x / real(n, dp)
      ss = 0.0_dp
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) ss = ss + (x(i) - mean_x)**2
      end do
      if (n > 1) then
         sd_x = sqrt(ss / real(n - 1, dp))
      else
         sd_x = 0.0_dp
      end if
      center_only = .false.
      if (present(mean_only)) center_only = mean_only
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            z(i) = x(i)
         else if (center_only) then
            z(i) = x(i) - mean_x
         else if (sd_x > 0.0_dp) then
            z(i) = (x(i) - mean_x) / sd_x
         else
            z(i) = 0.0_dp
         end if
      end do
   end function zscore_vector

   pure function zscore_rows(x, mean_only) result(z)
      real(dp), intent(in) :: x(:, :) !! Matrix whose rows are standardized independently, matching MARSS zscore matrix behavior.
      logical, intent(in), optional :: mean_only !! If true, subtract each row mean without variance scaling.
      real(dp) :: z(size(x, 1), size(x, 2))
      integer :: i

      do i = 1, size(x, 1)
         if (present(mean_only)) then
            z(i, :) = zscore_vector(x(i, :), mean_only)
         else
            z(i, :) = zscore_vector(x(i, :))
         end if
      end do
   end function zscore_rows

   pure function ldiag_scalar(x, nrow) result(a)
      real(dp), intent(in) :: x !! Scalar value recycled through the full matrix, matching R matrix(x,nrow,nrow).
      integer, intent(in) :: nrow !! Number of rows and columns of the returned square matrix.
      real(dp) :: a(nrow, nrow)

      a = x
   end function ldiag_scalar

   pure function ldiag_vector(x, nrow) result(a)
      real(dp), intent(in) :: x(:) !! Values placed on the diagonal, values beyond the matrix order are ignored.
      integer, intent(in), optional :: nrow !! Matrix order, defaults to size(x).
      real(dp), allocatable :: a(:, :)
      integer :: i
      integer :: n

      n = size(x)
      if (present(nrow)) n = nrow
      allocate(a(n, n))
      a = 0.0_dp
      do i = 1, min(n, size(x))
         a(i, i) = x(i)
      end do
   end function ldiag_vector

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p !! Lower-tail standard-normal probability strictly between zero and one.
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 = 2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 = 1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 = 2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 = 1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 = 6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 = 4.374664141464968_dp
      real(dp), parameter :: c6 = 2.938163982698783_dp
      real(dp), parameter :: d1 = 7.784695709041462e-3_dp
      real(dp), parameter :: d2 = 3.224671290700398e-1_dp
      real(dp), parameter :: d3 = 2.445134137142996_dp
      real(dp), parameter :: d4 = 3.754408661907416_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q
      real(dp) :: r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
            ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
            (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
            ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      end if
   end function normal_quantile

end module marss_utils
