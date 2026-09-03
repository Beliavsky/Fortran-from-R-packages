module grbase_reductions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use r_kinds, only : dp
  implicit none
  private

  public :: row_sums
  public :: column_sums
  public :: columnwise_product
  public :: matrix_nonzero_indices

contains

  pure function row_sums(x) result(out)
    real(dp), intent(in) :: x(:, :) !! Numeric matrix whose rows are reduced across all columns.
    real(dp) :: out(size(x, 1))

    out = sum(x, dim=2)
  end function row_sums

  pure function column_sums(x) result(out)
    real(dp), intent(in) :: x(:, :) !! Numeric matrix whose columns are reduced across all rows.
    real(dp) :: out(size(x, 2))

    out = sum(x, dim=1)
  end function column_sums

  pure function columnwise_product(weights, x) result(out)
    real(dp), intent(in) :: weights(:) !! Nonempty column multipliers, recycled cyclically across matrix columns.
    real(dp), intent(in) :: x(:, :) !! Numeric matrix whose columns are multiplied by recycled weights.
    real(dp) :: out(size(x, 1), size(x, 2))
    integer :: j

    if (size(weights) == 0) then
      out = 0.0_dp
      return
    end if

    do j = 1, size(x, 2)
      out(:, j) = weights(mod(j - 1, size(weights)) + 1) * x(:, j)
    end do
  end function columnwise_product

  pure function matrix_nonzero_indices(x) result(indices)
    real(dp), intent(in) :: x(:, :) !! Numeric matrix scanned by row then column for exact nonzero entries.
    integer, allocatable :: indices(:, :)
    integer :: i
    integer :: j
    integer :: k
    integer :: n

    n = count(ieee_is_nan(x) .or. abs(x) > 0.0_dp)
    allocate(indices(n, 2))
    k = 0
    do i = 1, size(x, 1)
      do j = 1, size(x, 2)
        if (ieee_is_nan(x(i, j)) .or. abs(x(i, j)) > 0.0_dp) then
          k = k + 1
          indices(k, :) = [i, j]
        end if
      end do
    end do
  end function matrix_nonzero_indices

end module grbase_reductions
