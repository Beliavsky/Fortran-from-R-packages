! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_hac
  use intrinsicfrp_kinds, only: dp, status_ok, status_invalid
  use intrinsicfrp_linalg, only: solve_linear, inverse_matrix, identity_matrix
  implicit none
  private
  public :: hac_covariance, hac_variance, hac_standard_errors, newey_west_lags

contains

  pure integer function newey_west_lags(n) result(n_lags)
    integer, intent(in) :: n
    if (n > 5) then
      n_lags = int(floor(4.0_dp * (0.01_dp * real(n, dp)) ** (2.0_dp / 9.0_dp)))
    else
      n_lags = 0
    end if
  end function newey_west_lags

  subroutine hac_covariance(series, covariance, status, prewhite)
    real(dp), intent(in) :: series(:, :)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status
    logical, intent(in), optional :: prewhite
    real(dp), allocatable :: work(:, :), coefficients(:, :), inv_temp(:, :)
    real(dp), allocatable :: xtx(:, :), rhs(:, :)
    real(dp) :: weight
    integer :: n, p, lag, st
    logical :: pw

    n = size(series, 1)
    p = size(series, 2)
    allocate(covariance(p, p))
    covariance = 0.0_dp
    status = status_invalid
    if (n < 2 .or. p < 1) return
    pw = .false.
    if (present(prewhite)) pw = prewhite
    work = series

    if (pw .and. n > 2) then
      xtx = matmul(transpose(work(1:n - 1, :)), work(1:n - 1, :))
      rhs = matmul(transpose(work(1:n - 1, :)), work(2:n, :))
      call solve_linear(xtx, rhs, coefficients, st)
      work(2:n, :) = work(2:n, :) - matmul(work(1:n - 1, :), coefficients)
    else
      allocate(coefficients(0, 0))
    end if

    covariance = matmul(transpose(work), work) / real(n, dp)
    do lag = 1, newey_west_lags(n)
      weight = 1.0_dp - real(lag, dp) / real(newey_west_lags(n) + 1, dp)
      covariance = covariance + weight * ( &
        matmul(transpose(work(lag + 1:n, :)), work(1:n - lag, :)) + &
        matmul(transpose(work(1:n - lag, :)), work(lag + 1:n, :))) / real(n, dp)
    end do

    if (pw .and. n > 2) then
      call inverse_matrix(identity_matrix(p) - transpose(coefficients), inv_temp, st)
      covariance = matmul(inv_temp, matmul(covariance, transpose(inv_temp)))
    end if
    covariance = 0.5_dp * (covariance + transpose(covariance))
    status = status_ok
  end subroutine hac_covariance

  subroutine hac_standard_errors(series, standard_errors, status, prewhite)
    real(dp), intent(in) :: series(:, :)
    real(dp), allocatable, intent(out) :: standard_errors(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: prewhite
    real(dp), allocatable :: cov(:, :)
    integer :: i
    call hac_covariance(series, cov, status, prewhite)
    allocate(standard_errors(size(series, 2)))
    do i = 1, size(standard_errors)
      standard_errors(i) = sqrt(max(0.0_dp, cov(i, i)))
    end do
  end subroutine hac_standard_errors

  subroutine hac_variance(series, variance, status, prewhite)
    real(dp), intent(in) :: series(:)
    real(dp), intent(out) :: variance
    integer, intent(out) :: status
    logical, intent(in), optional :: prewhite
    real(dp), allocatable :: x(:, :), cov(:, :)
    allocate(x(size(series), 1))
    x(:, 1) = series
    call hac_covariance(x, cov, status, prewhite)
    variance = cov(1, 1)
  end subroutine hac_variance

end module intrinsicfrp_hac
