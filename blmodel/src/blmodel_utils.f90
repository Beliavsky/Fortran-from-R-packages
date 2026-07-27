! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_utils
  use blmodel_kinds, only : dp
  implicit none
  private

  public :: diag_of, make_diag

contains

  pure function diag_of(matrix) result(diagonal)
    real(dp), intent(in) :: matrix(:,:)
    real(dp) :: diagonal(min(size(matrix, 1), size(matrix, 2)))
    integer :: i

    do i = 1, size(diagonal)
      diagonal(i) = matrix(i, i)
    end do
  end function diag_of

  pure function make_diag(diagonal) result(matrix)
    real(dp), intent(in) :: diagonal(:)
    real(dp) :: matrix(size(diagonal), size(diagonal))
    integer :: i

    matrix = 0.0_dp
    do i = 1, size(diagonal)
      matrix(i, i) = diagonal(i)
    end do
  end function make_diag

end module blmodel_utils
