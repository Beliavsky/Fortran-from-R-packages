! SPDX-License-Identifier: MIT
module mfgarch_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use mfgarch_kinds, only : dp
  use mfgarch_status, only : mfgarch_success, mfgarch_dimension_error, &
    mfgarch_singular_matrix, mfgarch_invalid_argument
  use r_rolling, only : r_roll_mean_right
  implicit none
  private

  real(dp), parameter :: pi_dp = acos(-1.0_dp)

  public :: sample_mean, sample_variance, variance_ignore_nan
  public :: invert_matrix, identity_matrix, vector_norm
  public :: normal_cdf, rolling_mean, finite_value, pi_dp

contains

  pure logical function finite_value(x)
    real(dp), intent(in) :: x
    finite_value = ieee_is_finite(x)
  end function finite_value

  pure function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function sample_mean

  pure function sample_variance(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp) :: mean_x

    if (size(x) <= 1) then
      value = 0.0_dp
      return
    end if
    mean_x = sample_mean(x)
    value = sum((x - mean_x)**2) / real(size(x) - 1, dp)
  end function sample_variance

  pure function variance_ignore_nan(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp) :: mean_x
    integer :: i, n

    mean_x = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (finite_value(x(i))) then
        n = n + 1
        mean_x = mean_x + (x(i) - mean_x) / real(n, dp)
      end if
    end do
    if (n <= 1) then
      value = 0.0_dp
      return
    end if
    value = 0.0_dp
    do i = 1, size(x)
      if (finite_value(x(i))) value = value + (x(i) - mean_x)**2
    end do
    value = value / real(n - 1, dp)
  end function variance_ignore_nan

  pure function vector_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i

    a = 0.0_dp
    do i = 1, n
      a(i,i) = 1.0_dp
    end do
  end function identity_matrix

  subroutine invert_matrix(a, inverse, status, tolerance)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: inverse(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: aug(:,:), row_tmp(:)
    real(dp) :: pivot, tol
    integer :: n, i, j, pivot_row

    status = mfgarch_success
    if (size(a,1) /= size(a,2)) then
      status = mfgarch_dimension_error
      allocate(inverse(0,0))
      return
    end if
    n = size(a,1)
    allocate(aug(n,2*n), row_tmp(2*n), inverse(n,n))
    aug(:,1:n) = a
    aug(:,n+1:2*n) = identity_matrix(n)
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance

    do i = 1, n
      pivot_row = i
      do j = i + 1, n
        if (abs(aug(j,i)) > abs(aug(pivot_row,i))) pivot_row = j
      end do
      if (abs(aug(pivot_row,i)) <= tol) then
        status = mfgarch_singular_matrix
        inverse = 0.0_dp
        return
      end if
      if (pivot_row /= i) then
        row_tmp = aug(i,:)
        aug(i,:) = aug(pivot_row,:)
        aug(pivot_row,:) = row_tmp
      end if
      pivot = aug(i,i)
      aug(i,:) = aug(i,:) / pivot
      do j = 1, n
        if (j /= i) aug(j,:) = aug(j,:) - aug(j,i) * aug(i,:)
      end do
    end do
    inverse = aug(:,n+1:2*n)
  end subroutine invert_matrix

  pure function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  subroutine rolling_mean(x, window, values, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: window
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(out) :: status
    if (window <= 0) then
      status = mfgarch_invalid_argument
      allocate(values(0))
      return
    end if
    call r_roll_mean_right(x, window, values)
    status = mfgarch_success
  end subroutine rolling_mean

end module mfgarch_math
