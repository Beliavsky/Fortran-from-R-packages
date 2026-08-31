! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Small numerical helpers used by the mitml computational translation.
module mitml_numeric
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: mean_real
   public :: sample_variance
   public :: sample_covariance
   public :: covariance_columns
   public :: mean_columns
   public :: mean_cube_slices
   public :: trace_matrix

contains

   pure real(dp) function mean_real(x) result(value)
      real(dp), intent(in) :: x(:) !! Finite vector whose arithmetic mean is requested.

      if (size(x) == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function mean_real

   pure real(dp) function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose sample variance uses denominator n - 1.
      real(dp) :: center

      if (size(x) < 2) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      center = mean_real(x)
      value = sum((x - center)**2) / real(size(x) - 1, dp)
   end function sample_variance

   pure real(dp) function sample_covariance(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First vector in the sample covariance calculation.
      real(dp), intent(in) :: y(:) !! Second vector; must have the same length as x.
      real(dp) :: mean_x
      real(dp) :: mean_y

      if (size(x) /= size(y) .or. size(x) < 2) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mean_x = mean_real(x)
      mean_y = mean_real(y)
      value = sum((x - mean_x) * (y - mean_y)) / real(size(x) - 1, dp)
   end function sample_covariance

   pure subroutine mean_columns(x, values)
      real(dp), intent(in) :: x(:, :) !! Matrix with observations in columns, shape p by m.
      real(dp), intent(out) :: values(:) !! Row means, length p.
      integer :: i

      if (size(values) /= size(x, 1)) error stop "mean_columns: output length mismatch"
      do i = 1, size(x, 1)
         values(i) = mean_real(x(i, :))
      end do
   end subroutine mean_columns

   pure subroutine covariance_columns(x, covariance)
      real(dp), intent(in) :: x(:, :) !! Matrix with variables in rows and observations in columns.
      real(dp), intent(out) :: covariance(:, :) !! Sample covariance matrix of rows of x.
      real(dp) :: means(size(x, 1))
      integer :: i
      integer :: j
      integer :: m

      if (size(covariance, 1) /= size(x, 1) .or. size(covariance, 2) /= size(x, 1)) &
         error stop "covariance_columns: output shape mismatch"
      m = size(x, 2)
      if (m < 2) then
         covariance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      call mean_columns(x, means)
      do j = 1, size(x, 1)
         do i = 1, size(x, 1)
            covariance(i, j) = sum((x(i, :) - means(i)) * (x(j, :) - means(j))) / real(m - 1, dp)
         end do
      end do
   end subroutine covariance_columns

   pure subroutine mean_cube_slices(x, mean_matrix)
      real(dp), intent(in) :: x(:, :, :) !! Three-dimensional array whose third dimension is averaged.
      real(dp), intent(out) :: mean_matrix(:, :) !! Mean of x over its third dimension.

      if (size(mean_matrix, 1) /= size(x, 1) .or. size(mean_matrix, 2) /= size(x, 2)) &
         error stop "mean_cube_slices: output shape mismatch"
      if (size(x, 3) == 0) then
         mean_matrix = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         mean_matrix = sum(x, dim=3) / real(size(x, 3), dp)
      end if
   end subroutine mean_cube_slices

   pure real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal sum is requested.
      integer :: i

      if (size(a, 1) /= size(a, 2)) error stop "trace_matrix: matrix must be square"
      value = 0.0_dp
      do i = 1, size(a, 1)
         value = value + a(i, i)
      end do
   end function trace_matrix

end module mitml_numeric
