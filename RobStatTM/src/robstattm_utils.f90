! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_utils
  use robstattm_kinds, only : dp
  use rrcov_linalg, only : general_inverse, determinant, mahalanobis_squared, make_positive_definite
  use rrcov_stats, only : median
  implicit none
  private
  public :: add_intercept, covariance_to_correlation, weighted_center_covariance
  public :: safe_logistic, logistic_deviance_residuals, median_absolute
  public :: normalize_determinant, matrix_inverse, squared_mahalanobis
  public :: lower_string, mean_value, outer_product
contains
  subroutine add_intercept(x, with_intercept, design)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in) :: with_intercept
    real(dp), allocatable, intent(out) :: design(:, :)
    if (with_intercept) then
      allocate(design(size(x, 1), size(x, 2) + 1))
      design(:, 1) = 1.0_dp
      design(:, 2:) = x
    else
      design = x
    end if
  end subroutine add_intercept

  pure function lower_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code
    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lower_string

  pure function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    a = spread(x, 2, size(y)) * spread(y, 1, size(x))
  end function outer_product

  elemental function safe_logistic(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    if (x >= 0.0_dp) then
      value = 1.0_dp / (1.0_dp + exp(-min(x, 700.0_dp)))
    else
      value = exp(max(x, -700.0_dp)) / (1.0_dp + exp(max(x, -700.0_dp)))
    end if
  end function safe_logistic

  subroutine logistic_deviance_residuals(y, fitted, residuals)
    real(dp), intent(in) :: y(:), fitted(:)
    real(dp), allocatable, intent(out) :: residuals(:)
    real(dp), allocatable :: p(:), dev(:)
    allocate(residuals(size(y)), p(size(y)), dev(size(y)))
    p = max(1.0e-14_dp, min(1.0_dp - 1.0e-14_dp, fitted))
    dev = -2.0_dp * (y * log(p) + (1.0_dp - y) * log(1.0_dp - p))
    residuals = sign(sqrt(max(dev, 0.0_dp)), y - fitted)
  end subroutine logistic_deviance_residuals

  function median_absolute(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = median(abs(x))
  end function median_absolute

  subroutine covariance_to_correlation(covariance, correlation)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: correlation(:, :)
    real(dp), allocatable :: scale(:)
    integer :: i, j, p
    p = size(covariance, 1)
    allocate(correlation(p, p), scale(p))
    do i = 1, p
      scale(i) = sqrt(max(covariance(i, i), tiny(1.0_dp)))
    end do
    do j = 1, p
      do i = 1, p
        correlation(i, j) = covariance(i, j) / (scale(i) * scale(j))
      end do
    end do
    do i = 1, p
      correlation(i, i) = 1.0_dp
    end do
  end subroutine covariance_to_correlation

  subroutine weighted_center_covariance(x, weights, center, covariance, denominator)
    real(dp), intent(in) :: x(:, :), weights(:)
    real(dp), allocatable, intent(out) :: center(:), covariance(:, :)
    real(dp), intent(in), optional :: denominator
    real(dp), allocatable :: z(:, :)
    real(dp) :: sw, den
    integer :: j, p
    p = size(x, 2)
    allocate(center(p), covariance(p, p), z(size(x, 1), p))
    sw = sum(max(weights, 0.0_dp))
    if (sw <= tiny(1.0_dp)) then
      center = sum(x, dim=1) / real(size(x, 1), dp)
      covariance = 0.0_dp
      return
    end if
    do j = 1, p
      center(j) = sum(max(weights, 0.0_dp) * x(:, j)) / sw
      z(:, j) = x(:, j) - center(j)
    end do
    den = sw
    if (present(denominator)) den = max(denominator, tiny(1.0_dp))
    covariance = matmul(transpose(z), z * spread(max(weights, 0.0_dp), 2, p)) / den
    covariance = make_positive_definite(0.5_dp * (covariance + transpose(covariance)))
  end subroutine weighted_center_covariance

  subroutine normalize_determinant(covariance, status)
    real(dp), intent(inout) :: covariance(:, :)
    integer, intent(out) :: status
    real(dp) :: detv, factor
    integer :: p
    p = size(covariance, 1)
    detv = determinant(covariance, status)
    if (status /= 0 .or. detv <= tiny(1.0_dp)) then
      covariance = make_positive_definite(covariance)
      detv = determinant(covariance, status)
    end if
    if (status == 0 .and. detv > tiny(1.0_dp)) then
      factor = detv ** (1.0_dp / real(p, dp))
      covariance = covariance / factor
    end if
  end subroutine normalize_determinant

  function matrix_inverse(a, status) result(inverse)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: inverse(:, :)
    inverse = general_inverse(a, status)
  end function matrix_inverse

  subroutine squared_mahalanobis(x, center, covariance, distances, status)
    real(dp), intent(in) :: x(:, :), center(:), covariance(:, :)
    real(dp), allocatable, intent(out) :: distances(:)
    integer, intent(out) :: status
    call mahalanobis_squared(x, center, covariance, distances, status)
  end subroutine squared_mahalanobis
end module robstattm_utils
