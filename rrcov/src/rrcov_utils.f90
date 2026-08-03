! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_utils
  use rrcov_kinds, only : dp
  use rrcov_types, only : covariance_result, rrcov_success, rrcov_dimension_error
  use rrcov_stats, only : chi_square_quantile
  implicit none
  private
  public :: vecnorm, covariance_to_correlation, outlier_flags, ilr_transform
contains
  pure function vecnorm(x, p) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: p
    real(dp) :: value, exponent
    exponent = 2.0_dp
    if (present(p)) exponent = p
    if (exponent <= 0.0_dp) then
      value = 0.0_dp
    else
      value = sum(abs(x) ** exponent) ** (1.0_dp / exponent)
    end if
  end function vecnorm

  function covariance_to_correlation(covariance, status) result(correlation)
    real(dp), intent(in) :: covariance(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: correlation(:, :), scale(:)
    integer :: i, j, p
    p = size(covariance, 1)
    allocate(correlation(p, p), scale(p))
    correlation = 0.0_dp
    if (size(covariance, 2) /= p) then
      if (present(status)) status = rrcov_dimension_error
      return
    end if
    do i = 1, p
      scale(i) = sqrt(max(covariance(i, i), 0.0_dp))
    end do
    do j = 1, p
      do i = 1, p
        if (scale(i) > tiny(1.0_dp) .and. scale(j) > tiny(1.0_dp)) then
          correlation(i, j) = covariance(i, j) / (scale(i) * scale(j))
        end if
      end do
      correlation(j, j) = 1.0_dp
    end do
    correlation = max(-1.0_dp, min(1.0_dp, correlation))
    if (present(status)) status = rrcov_success
  end function covariance_to_correlation

  function outlier_flags(estimate, probability) result(flags)
    type(covariance_result), intent(in) :: estimate
    real(dp), intent(in), optional :: probability
    logical, allocatable :: flags(:)
    real(dp) :: cutoff, prob
    prob = 0.975_dp
    if (present(probability)) prob = min(1.0_dp - epsilon(1.0_dp), max(epsilon(1.0_dp), probability))
    cutoff = chi_square_quantile(prob, real(max(1, estimate%rank), dp))
    allocate(flags(size(estimate%distances)))
    flags = estimate%distances > cutoff
  end function outlier_flags

  subroutine ilr_transform(composition, transformed, status)
    real(dp), intent(in) :: composition(:, :)
    real(dp), allocatable, intent(out) :: transformed(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: logx(:, :)
    real(dp) :: left_mean, right_value, coefficient
    integer :: n, p, i, j
    n = size(composition, 1)
    p = size(composition, 2)
    if (p < 2 .or. any(composition <= 0.0_dp)) then
      allocate(transformed(n, max(0, p - 1)))
      transformed = 0.0_dp
      status = rrcov_dimension_error
      return
    end if
    allocate(logx(n, p), transformed(n, p - 1))
    logx = log(composition)
    do i = 1, n
      do j = 1, p - 1
        left_mean = sum(logx(i, 1:j)) / real(j, dp)
        right_value = logx(i, j + 1)
        coefficient = sqrt(real(j, dp) / real(j + 1, dp))
        transformed(i, j) = coefficient * (left_mean - right_value)
      end do
    end do
    status = rrcov_success
  end subroutine ilr_transform
end module rrcov_utils
