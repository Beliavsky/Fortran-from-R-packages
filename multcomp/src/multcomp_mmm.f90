! SPDX-License-Identifier: GPL-2.0-only
module multcomp_mmm
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type
  use multcomp_parm, only : make_parm
  implicit none
  private

  public :: mmm_parm_from_iid
  public :: block_diagonal_matrix

contains

  subroutine mmm_parm_from_iid(blocks, iid, result)
    type(parm_type), intent(in) :: blocks(:) !! Marginal model parameter blocks supplying estimates and model-based variances.
    real(dp), intent(in) :: iid(:, :) !! Zero-filled IID influence rows for all concatenated coefficients by observation.
    type(parm_type), intent(out) :: result !! Combined asymptotic parameter object with cross-model sandwich correlations.

    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: marginal_sd(:)
    real(dp), allocatable :: raw_sd(:)
    integer :: first
    integer :: i
    integer :: j
    integer :: last
    integer :: nobs
    integer :: total

    result%ok = .false.
    result%message = ''
    total = 0
    do i = 1, size(blocks)
      if (.not. blocks(i)%ok) then
        result%message = 'every marginal parameter block must be valid'
        return
      end if
      total = total + size(blocks(i)%coef)
    end do
    nobs = size(iid, 2)
    if (total < 1 .or. size(iid, 1) /= total .or. nobs < 1) then
      result%message = 'iid dimensions must match all concatenated coefficients and observations'
      return
    end if

    allocate(coef(total), marginal_sd(total))
    first = 1
    do i = 1, size(blocks)
      last = first + size(blocks(i)%coef) - 1
      coef(first:last) = blocks(i)%coef
      do j = first, last
        marginal_sd(j) = sqrt(blocks(i)%vcov(j - first + 1, j - first + 1))
      end do
      first = last + 1
    end do

    covariance = matmul(iid, transpose(iid)) / real(nobs * nobs, dp)
    allocate(raw_sd(total), correlation(total, total))
    do i = 1, total
      if (covariance(i, i) <= 0.0_dp) then
        result%message = 'IID covariance has a nonpositive marginal variance'
        return
      end if
      raw_sd(i) = sqrt(covariance(i, i))
    end do
    do i = 1, total
      do j = 1, total
        correlation(i, j) = covariance(i, j) / (raw_sd(i) * raw_sd(j))
        covariance(i, j) = correlation(i, j) * marginal_sd(i) * marginal_sd(j)
      end do
    end do
    call make_parm(coef, covariance, result, df=0.0_dp)
  end subroutine mmm_parm_from_iid

  subroutine block_diagonal_matrix(blocks, block_rows, block_cols, output, ok)
    real(dp), intent(in) :: blocks(:, :, :) !! Padded matrix blocks stored as (max_rows,max_cols,nblocks).
    integer, intent(in) :: block_rows(:) !! Active row count for each padded block.
    integer, intent(in) :: block_cols(:) !! Active column count for each padded block.
    real(dp), allocatable, intent(out) :: output(:, :) !! Dense block-diagonal matrix containing active block regions.
    logical, intent(out) :: ok !! True when block dimensions are valid and output is constructed.

    integer :: col_first
    integer :: col_last
    integer :: i
    integer :: row_first
    integer :: row_last
    integer :: total_cols
    integer :: total_rows

    ok = .false.
    if (size(block_rows) /= size(blocks, 3) .or. size(block_cols) /= size(blocks, 3)) return
    if (any(block_rows < 0) .or. any(block_rows > size(blocks, 1))) return
    if (any(block_cols < 0) .or. any(block_cols > size(blocks, 2))) return

    total_rows = sum(block_rows)
    total_cols = sum(block_cols)
    allocate(output(total_rows, total_cols))
    output = 0.0_dp
    row_first = 1
    col_first = 1
    do i = 1, size(blocks, 3)
      row_last = row_first + block_rows(i) - 1
      col_last = col_first + block_cols(i) - 1
      if (block_rows(i) > 0 .and. block_cols(i) > 0) then
        output(row_first:row_last, col_first:col_last) = &
          blocks(1:block_rows(i), 1:block_cols(i), i)
      end if
      row_first = row_last + 1
      col_first = col_last + 1
    end do
    ok = .true.
  end subroutine block_diagonal_matrix

end module multcomp_mmm
