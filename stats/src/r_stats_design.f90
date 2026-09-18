! SPDX-License-Identifier: MIT
! SPDX-FileComment: Shared dense-design utilities for translated R stats models.
module r_stats_design
   use r_kinds, only: dp
   implicit none
   private

   public :: independent_columns

contains

   pure subroutine independent_columns(matrix, relative_tolerance, indices)
      !! Selects the first numerically independent design columns by modified Gram-Schmidt.
      real(dp), intent(in) :: matrix(:, :) !! Weighted design matrix with shape `(n,p)`.
      real(dp), intent(in) :: relative_tolerance !! Relative norm threshold.
      integer, allocatable, intent(out) :: indices(:) !! One-based independent column indices.
      real(dp), allocatable :: basis(:, :), candidate(:)
      real(dp) :: column_norm, scale, threshold
      integer :: j, rank

      allocate (basis(size(matrix, 1), size(matrix, 2)))
      allocate (indices(size(matrix, 2)))
      scale = maxval(sqrt(sum(matrix**2, dim=1)))
      threshold = relative_tolerance*max(1.0_dp, scale)
      rank = 0
      do j = 1, size(matrix, 2)
         candidate = matrix(:, j)
         if (rank > 0) then
            candidate = candidate - matmul(basis(:, :rank), &
                                            matmul(transpose(basis(:, :rank)), candidate))
         end if
         column_norm = norm2(candidate)
         if (column_norm > threshold) then
            rank = rank + 1
            indices(rank) = j
            basis(:, rank) = candidate/column_norm
         end if
      end do
      indices = indices(:rank)
   end subroutine independent_columns

end module r_stats_design
