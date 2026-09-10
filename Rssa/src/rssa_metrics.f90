! SPDX-License-Identifier: GPL-2.0-or-later
! Weighted norms and correlations translated from Rssa 1.1.
module rssa_metrics
   use rssa_kinds, only : dp
   use rssa_matrices, only : hankel_weights, hankel_weights_2d
   use rssa_reconstruction, only : elementary_series_ssa
   use rssa_types, only : cssa_result, mssa_result, ssa2d_result, ssa_result
   implicit none
   private
   public :: wnorm, wnorm_default, wnorm_complex, wnorm_ssa, wnorm_mssa, wnorm_2d
   public :: wcor_default, wcor_ssa, frobenius_cor
   interface wnorm
      module procedure wnorm_default
      module procedure wnorm_complex
      module procedure wnorm_ssa
      module procedure wnorm_mssa
      module procedure wnorm_2d
      module procedure wnorm_cssa
   end interface wnorm
contains
   pure real(dp) function wnorm_default(series, window) result(value)
      real(dp), intent(in) :: series(:) !! Real series whose squared values receive Hankel overlap weights.
      integer, intent(in), optional :: window !! Window length; defaults to (N+1)/2.
      real(dp), allocatable :: weights(:)
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      weights = hankel_weights(size(series), l)
      value = sqrt(sum(weights * series * series))
   end function wnorm_default
   pure real(dp) function wnorm_complex(series, window) result(value)
      complex(dp), intent(in) :: series(:) !! Complex series whose squared moduli receive Hankel weights.
      integer, intent(in), optional :: window !! Window length; defaults to (N+1)/2.
      real(dp), allocatable :: weights(:)
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      weights = hankel_weights(size(series), l)
      value = sqrt(sum(weights * abs(series) ** 2))
   end function wnorm_complex
   pure real(dp) function wnorm_ssa(object) result(value)
      type(ssa_result), intent(in) :: object !! Real SSA object containing source series and window.
      value = 0.0_dp
      if (allocated(object%series)) value = wnorm_default(object%series, object%window)
   end function wnorm_ssa
   pure real(dp) function wnorm_cssa(object) result(value)
      type(cssa_result), intent(in) :: object !! Complex SSA object containing source series and window.
      value = 0.0_dp
      if (allocated(object%series)) value = wnorm_complex(object%series, object%window)
   end function wnorm_cssa
   pure real(dp) function wnorm_mssa(object) result(value)
      type(mssa_result), intent(in) :: object !! MSSA object with active per-channel lengths.
      real(dp), allocatable :: weights(:)
      integer :: j
      value = 0.0_dp
      if (.not. allocated(object%series)) return
      do j = 1, size(object%series, 2)
         weights = hankel_weights(object%lengths(j), object%window)
         value = value + sum(weights * object%series(1:object%lengths(j), j) ** 2)
      end do
      value = sqrt(value)
   end function wnorm_mssa
   pure real(dp) function wnorm_2d(object) result(value)
      type(ssa2d_result), intent(in) :: object !! Rectangular 2-D SSA object.
      real(dp), allocatable :: weights(:, :)
      value = 0.0_dp
      if (.not. allocated(object%field)) return
      weights = hankel_weights_2d(shape(object%field), object%window)
      value = sqrt(sum(weights * object%field * object%field))
   end function wnorm_2d
   pure real(dp) function wcor_default(x, y, window) result(value)
      real(dp), intent(in) :: x(:) !! First real series.
      real(dp), intent(in) :: y(:) !! Second real series of the same length.
      integer, intent(in), optional :: window !! Hankel window; defaults to (N+1)/2.
      real(dp), allocatable :: weights(:)
      real(dp) :: nx, ny
      integer :: l
      value = 0.0_dp
      if (size(x) /= size(y)) return
      l = (size(x) + 1) / 2
      if (present(window)) l = window
      weights = hankel_weights(size(x), l)
      nx = sqrt(sum(weights * x * x))
      ny = sqrt(sum(weights * y * y))
      if (nx > 0.0_dp .and. ny > 0.0_dp) value = sum(weights * x * y) / (nx * ny)
   end function wcor_default
   function wcor_ssa(object, indices) result(matrix)
      type(ssa_result), intent(in) :: object !! SSA object whose elementary reconstructed components are correlated.
      integer, intent(in), optional :: indices(:) !! Selected components; defaults to all.
      real(dp), allocatable :: matrix(:, :), components(:, :)
      integer, allocatable :: idx(:)
      integer :: i, j, n
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      n = size(idx)
      allocate(components(size(object%series), n), matrix(n, n))
      do i = 1, n
         components(:, i) = elementary_series_ssa(object, idx(i))
      end do
      do j = 1, n
         do i = 1, n
            matrix(i, j) = wcor_default(components(:, i), components(:, j), object%window)
         end do
      end do
   end function wcor_ssa
   function frobenius_cor(components) result(matrix)
      real(dp), intent(in) :: components(:, :) !! Columns containing vectorized matrices or series components.
      real(dp), allocatable :: matrix(:, :)
      integer :: i, j, n
      real(dp) :: ni, nj
      n = size(components, 2)
      allocate(matrix(n, n))
      do j = 1, n
         nj = sqrt(sum(components(:, j) ** 2))
         do i = 1, n
            ni = sqrt(sum(components(:, i) ** 2))
            if (ni > 0.0_dp .and. nj > 0.0_dp) then
               matrix(i, j) = dot_product(components(:, i), components(:, j)) / (ni * nj)
            else
               matrix(i, j) = 0.0_dp
            end if
         end do
      end do
   end function frobenius_cor
end module rssa_metrics
