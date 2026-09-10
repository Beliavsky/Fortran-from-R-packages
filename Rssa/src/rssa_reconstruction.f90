! SPDX-License-Identifier: GPL-2.0-or-later
! Reconstruction and residual kernels translated from Rssa 1.1.
module rssa_reconstruction
   use rssa_kinds, only : dp
   use rssa_matrices, only : complex_hankelize_matrix, hankelize_2d, hankelize_matrix
   use rssa_types, only : cssa_result, mssa_result, ssa2d_result, ssa_result
   implicit none
   private
   public :: reconstruct_ssa, reconstruct_mssa, reconstruct_2d, reconstruct_complex
   public :: elementary_series_ssa, elementary_series_mssa, elementary_field_2d, elementary_series_complex
   public :: residuals_ssa, residuals_mssa, residuals_2d, residuals_complex
contains
   function elementary_series_ssa(object, index) result(series)
      type(ssa_result), intent(in) :: object !! Real SSA decomposition.
      integer, intent(in) :: index !! One-based eigentriple index.
      real(dp), allocatable :: series(:), matrix(:, :)
      if (index < 1 .or. index > size(object%sigma)) then
         allocate(series(0))
         return
      end if
      matrix = object%sigma(index) * spread(object%u(:, index), 2, size(object%v, 1)) * &
               spread(object%v(:, index), 1, size(object%u, 1))
      series = hankelize_matrix(matrix)
   end function elementary_series_ssa

   function reconstruct_ssa(object, indices) result(series)
      type(ssa_result), intent(in) :: object !! Real or Toeplitz SSA decomposition.
      integer, intent(in), optional :: indices(:) !! Components to reconstruct; defaults to all stored eigentriples.
      real(dp), allocatable :: series(:), one(:)
      integer, allocatable :: idx(:)
      integer :: i
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      allocate(series(size(object%series)))
      series = 0.0_dp
      do i = 1, size(idx)
         one = elementary_series_ssa(object, idx(i))
         if (size(one) == size(series)) series = series + one
      end do
   end function reconstruct_ssa

   function residuals_ssa(object, indices) result(residual)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      integer, intent(in), optional :: indices(:) !! Components removed from the original series.
      real(dp), allocatable :: residual(:)
      if (present(indices)) then
         residual = object%series - reconstruct_ssa(object, indices)
      else
         residual = object%series - reconstruct_ssa(object)
      end if
   end function residuals_ssa

   function elementary_series_mssa(object, index) result(series)
      type(mssa_result), intent(in) :: object !! MSSA decomposition.
      integer, intent(in) :: index !! One-based eigentriple index.
      real(dp), allocatable :: series(:, :), matrix(:, :), one(:)
      integer :: j, first, last, k
      allocate(series(size(object%series, 1), size(object%series, 2)))
      series = 0.0_dp
      if (index < 1 .or. index > size(object%sigma)) return
      first = 1
      do j = 1, size(object%series, 2)
         k = object%lengths(j) - object%window + 1
         last = first + k - 1
         matrix = object%sigma(index) * spread(object%u(:, index), 2, k) * &
                  spread(object%v(first:last, index), 1, object%window)
         one = hankelize_matrix(matrix)
         series(1:object%lengths(j), j) = one
         first = last + 1
      end do
   end function elementary_series_mssa

   function reconstruct_mssa(object, indices) result(series)
      type(mssa_result), intent(in) :: object !! MSSA decomposition.
      integer, intent(in), optional :: indices(:) !! Components to reconstruct; defaults to all.
      real(dp), allocatable :: series(:, :), one(:, :)
      integer, allocatable :: idx(:)
      integer :: i
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      allocate(series(size(object%series, 1), size(object%series, 2)))
      series = 0.0_dp
      do i = 1, size(idx)
         one = elementary_series_mssa(object, idx(i))
         series = series + one
      end do
   end function reconstruct_mssa

   function residuals_mssa(object, indices) result(residual)
      type(mssa_result), intent(in) :: object !! MSSA object.
      integer, intent(in), optional :: indices(:) !! Components removed from source channels.
      real(dp), allocatable :: residual(:, :)
      if (present(indices)) then
         residual = object%series - reconstruct_mssa(object, indices)
      else
         residual = object%series - reconstruct_mssa(object)
      end if
   end function residuals_mssa

   function elementary_field_2d(object, index) result(field)
      type(ssa2d_result), intent(in) :: object !! Two-dimensional SSA decomposition.
      integer, intent(in) :: index !! One-based eigentriple index.
      real(dp), allocatable :: field(:, :), matrix(:, :)
      if (index < 1 .or. index > size(object%sigma)) then
         allocate(field(0, 0))
         return
      end if
      matrix = object%sigma(index) * spread(object%u(:, index), 2, size(object%v, 1)) * &
               spread(object%v(:, index), 1, size(object%u, 1))
      field = hankelize_2d(matrix, shape(object%field), object%window)
   end function elementary_field_2d

   function reconstruct_2d(object, indices) result(field)
      type(ssa2d_result), intent(in) :: object !! Two-dimensional SSA decomposition.
      integer, intent(in), optional :: indices(:) !! Components to reconstruct; defaults to all.
      real(dp), allocatable :: field(:, :), one(:, :)
      integer, allocatable :: idx(:)
      integer :: i
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      allocate(field(size(object%field, 1), size(object%field, 2)))
      field = 0.0_dp
      do i = 1, size(idx)
         one = elementary_field_2d(object, idx(i))
         if (all(shape(one) == shape(field))) field = field + one
      end do
   end function reconstruct_2d

   function residuals_2d(object, indices) result(residual)
      type(ssa2d_result), intent(in) :: object !! Two-dimensional SSA object.
      integer, intent(in), optional :: indices(:) !! Components removed from source field.
      real(dp), allocatable :: residual(:, :)
      if (present(indices)) then
         residual = object%field - reconstruct_2d(object, indices)
      else
         residual = object%field - reconstruct_2d(object)
      end if
   end function residuals_2d

   function elementary_series_complex(object, index) result(series)
      type(cssa_result), intent(in) :: object !! Complex SSA decomposition.
      integer, intent(in) :: index !! One-based eigentriple index.
      complex(dp), allocatable :: series(:), matrix(:, :)
      if (index < 1 .or. index > size(object%sigma)) then
         allocate(series(0))
         return
      end if
      matrix = cmplx(object%sigma(index), 0.0_dp, kind=dp) * &
               spread(object%u(:, index), 2, size(object%v, 1)) * &
               spread(conjg(object%v(:, index)), 1, size(object%u, 1))
      series = complex_hankelize_matrix(matrix)
   end function elementary_series_complex

   function reconstruct_complex(object, indices) result(series)
      type(cssa_result), intent(in) :: object !! Complex SSA decomposition.
      integer, intent(in), optional :: indices(:) !! Components to reconstruct; defaults to all.
      complex(dp), allocatable :: series(:), one(:)
      integer, allocatable :: idx(:)
      integer :: i
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      allocate(series(size(object%series)))
      series = (0.0_dp, 0.0_dp)
      do i = 1, size(idx)
         one = elementary_series_complex(object, idx(i))
         if (size(one) == size(series)) series = series + one
      end do
   end function reconstruct_complex

   function residuals_complex(object, indices) result(residual)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      integer, intent(in), optional :: indices(:) !! Components removed from the source series.
      complex(dp), allocatable :: residual(:)
      if (present(indices)) then
         residual = object%series - reconstruct_complex(object, indices)
      else
         residual = object%series - reconstruct_complex(object)
      end if
   end function residuals_complex
end module rssa_reconstruction
