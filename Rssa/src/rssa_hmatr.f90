! SPDX-License-Identifier: GPL-2.0-or-later
! Heterogeneity matrix calculation translated from Rssa 1.1.
module rssa_hmatr
   use rssa_kinds, only : dp
   use rssa_reconstruction, only : elementary_series_ssa
   use rssa_types, only : ssa_result
   implicit none
   private
   public :: hmatr
contains
   function hmatr(object, indices) result(matrix)
      type(ssa_result), intent(in) :: object !! Decomposed SSA object supplying elementary reconstructed series.
      integer, intent(in), optional :: indices(:) !! Components to compare; defaults to all stored components.
      real(dp), allocatable :: matrix(:, :), components(:, :), one(:)
      integer, allocatable :: idx(:)
      real(dp) :: denom
      integer :: i, j, n
      if (.not. allocated(object%sigma) .or. .not. allocated(object%series)) then
         allocate(matrix(0, 0))
         return
      end if
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i=1, size(idx))]
      end if
      n = size(idx)
      allocate(components(size(object%series), n), matrix(n, n))
      do i = 1, n
         one = elementary_series_ssa(object, idx(i))
         components(:, i) = one
      end do
      do j = 1, n
         do i = 1, n
            denom = sqrt(sum(components(:, i)**2) * sum(components(:, j)**2))
            if (denom > 0.0_dp) then
               matrix(i, j) = 1.0_dp - abs(dot_product(components(:, i), components(:, j))) / denom
            else
               matrix(i, j) = 0.0_dp
            end if
         end do
      end do
   end function hmatr
end module rssa_hmatr
