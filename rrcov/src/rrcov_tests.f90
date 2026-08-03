! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_tests
  use rrcov_kinds, only : dp
  use rrcov_types, only : test_result, covariance_result, rrcov_success, &
    rrcov_invalid_argument, rrcov_dimension_error
  use rrcov_stats, only : mean_vector, covariance_matrix, f_cdf, chi_square_cdf
  use rrcov_sort, only : rank_values
  use rrcov_linalg, only : symmetric_inverse, make_positive_definite, determinant
  use rrcov_robust, only : robust_covariance
  implicit none
  private
  public :: hotelling_t2_one_sample, hotelling_t2_two_sample, wilks_test
contains
  subroutine hotelling_t2_one_sample(x, mu, result, method, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :), mu(:)
    type(test_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    type(covariance_result) :: estimate
    real(dp), allocatable :: inverse(:, :), difference(:)
    real(dp) :: f_value
    character(len=16) :: covariance_method
    integer :: n, p, status
    n = size(x, 1)
    p = size(x, 2)
    if (size(mu) /= p .or. n <= p) then
      result%status = rrcov_dimension_error
      result%method = "Hotelling one-sample T2"
      return
    end if
    covariance_method = "classic"
    if (present(method)) covariance_method = method
    call robust_covariance(x, covariance_method, estimate, alpha=alpha, nsamp=nsamp, seed=seed)
    inverse = symmetric_inverse(make_positive_definite(estimate%covariance, 1.0e-10_dp), status)
    difference = estimate%center - mu
    result%statistic = real(n, dp) * dot_product(difference, matmul(inverse, difference))
    result%df1 = real(p, dp)
    result%df2 = real(n - p, dp)
    f_value = result%df2 * result%statistic / (result%df1 * real(n - 1, dp))
    result%p_value = 1.0_dp - f_cdf(f_value, result%df1, result%df2)
    result%status = estimate%status
    result%method = "Hotelling one-sample T2 (" // trim(covariance_method) // ")"
  end subroutine hotelling_t2_one_sample

  subroutine hotelling_t2_two_sample(x, y, result, method, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :), y(:, :)
    type(test_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    type(covariance_result) :: x_estimate, y_estimate
    real(dp), allocatable :: pooled(:, :), inverse(:, :), difference(:)
    real(dp) :: f_value
    character(len=16) :: covariance_method
    integer :: n1, n2, p, status
    n1 = size(x, 1)
    n2 = size(y, 1)
    p = size(x, 2)
    if (size(y, 2) /= p .or. n1 + n2 <= p + 1 .or. n1 < 2 .or. n2 < 2) then
      result%status = rrcov_dimension_error
      result%method = "Hotelling two-sample T2"
      return
    end if
    covariance_method = "classic"
    if (present(method)) covariance_method = method
    call robust_covariance(x, covariance_method, x_estimate, alpha=alpha, nsamp=nsamp, seed=seed)
    call robust_covariance(y, covariance_method, y_estimate, alpha=alpha, nsamp=nsamp, seed=seed)
    pooled = (real(n1 - 1, dp) * x_estimate%covariance + &
      real(n2 - 1, dp) * y_estimate%covariance) / real(n1 + n2 - 2, dp)
    pooled = make_positive_definite(pooled, 1.0e-10_dp)
    inverse = symmetric_inverse(pooled, status)
    difference = x_estimate%center - y_estimate%center
    result%statistic = real(n1 * n2, dp) / real(n1 + n2, dp) * &
      dot_product(difference, matmul(inverse, difference))
    result%df1 = real(p, dp)
    result%df2 = real(n1 + n2 - p - 1, dp)
    f_value = result%df2 * result%statistic / &
      (result%df1 * real(n1 + n2 - 2, dp))
    result%p_value = 1.0_dp - f_cdf(f_value, result%df1, result%df2)
    result%status = max(x_estimate%status, y_estimate%status)
    result%method = "Hotelling two-sample T2 (" // trim(covariance_method) // ")"
  end subroutine hotelling_t2_two_sample

  subroutine wilks_test(x, grouping, result, method, approximation, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(test_result), intent(out) :: result
    character(len=*), intent(in), optional :: method, approximation
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    type(covariance_result) :: total_estimate, within_estimate, group_estimate
    real(dp), allocatable :: work(:, :), residuals(:, :), group_data(:, :), &
      group_center(:), ranks(:), w(:, :), t(:, :)
    integer, allocatable :: labels(:), rows(:)
    character(len=16) :: covariance_method, approx
    real(dp) :: log_lambda, chi_value, exponent, transformed, df1, df2, f_value
    integer :: n, p, ng, g, i, j, count_rows, status
    n = size(x, 1)
    p = size(x, 2)
    if (size(grouping) /= n .or. n < 2 .or. p < 1) then
      result%status = rrcov_dimension_error
      result%method = "Wilks MANOVA"
      return
    end if
    covariance_method = "classic"
    if (present(method)) covariance_method = method
    approx = "bartlett"
    if (present(approximation)) approx = lowercase(approximation)
    call unique_labels(grouping, labels)
    ng = size(labels)
    if (ng < 2) then
      result%status = rrcov_invalid_argument
      result%method = "Wilks MANOVA"
      return
    end if
    allocate(work(n, p))
    work = x
    if (lowercase(covariance_method) == "rank") then
      allocate(ranks(n))
      do j = 1, p
        call rank_values(work(:, j), ranks)
        work(:, j) = ranks
      end do
      covariance_method = "classic"
    end if
    allocate(residuals(n, p))
    do g = 1, ng
      count_rows = count(grouping == labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == labels(g))
      group_data = work(rows, :)
      call robust_covariance(group_data, covariance_method, group_estimate, &
        alpha=alpha, nsamp=nsamp, seed=seed)
      group_center = group_estimate%center
      do i = 1, count_rows
        residuals(rows(i), :) = group_data(i, :) - group_center
      end do
      deallocate(rows, group_data, group_center)
    end do
    call robust_covariance(residuals, covariance_method, within_estimate, &
      alpha=alpha, nsamp=nsamp, seed=seed)
    call robust_covariance(work, covariance_method, total_estimate, &
      alpha=alpha, nsamp=nsamp, seed=seed)
    w = make_positive_definite(within_estimate%covariance, 1.0e-10_dp)
    t = make_positive_definite(total_estimate%covariance, 1.0e-10_dp)
    result%lambda = max(tiny(1.0_dp), min(1.0_dp, determinant(w, status) / &
      max(determinant(t, status), tiny(1.0_dp))))
    log_lambda = log(result%lambda)
    select case (trim(approx))
    case ("rao", "f")
      exponent = sqrt(max(1.0_dp, &
        (real(p * p * (ng - 1) * (ng - 1) - 4, dp)) / &
        max(1.0_dp, real(p * p + (ng - 1) * (ng - 1) - 5, dp))))
      df1 = real(p * (ng - 1), dp)
      df2 = exponent * (real(n - 1, dp) - 0.5_dp * real(p + ng, dp)) - &
        0.5_dp * real(p * (ng - 1) - 2, dp)
      transformed = result%lambda ** (1.0_dp / exponent)
      f_value = df2 * (1.0_dp - transformed) / max(df1 * transformed, tiny(1.0_dp))
      result%statistic = f_value
      result%df1 = df1
      result%df2 = df2
      result%p_value = 1.0_dp - f_cdf(f_value, df1, df2)
      result%method = "One-way MANOVA (Rao F approximation)"
    case default
      chi_value = -(real(n - 1, dp) - 0.5_dp * real(p + ng, dp)) * log_lambda
      result%statistic = chi_value
      result%df1 = real(p * (ng - 1), dp)
      result%df2 = 0.0_dp
      result%p_value = 1.0_dp - chi_square_cdf(chi_value, result%df1)
      result%method = "One-way MANOVA (Bartlett chi-square)"
    end select
    result%status = max(within_estimate%status, total_estimate%status)
  end subroutine wilks_test

  subroutine unique_labels(grouping, labels)
    integer, intent(in) :: grouping(:)
    integer, allocatable, intent(out) :: labels(:)
    integer, allocatable :: work(:)
    integer :: i, count_labels
    allocate(work(size(grouping)))
    count_labels = 0
    do i = 1, size(grouping)
      if (count_labels == 0 .or. .not. any(work(1:count_labels) == grouping(i))) then
        count_labels = count_labels + 1
        work(count_labels) = grouping(i)
      end if
    end do
    allocate(labels(count_labels))
    labels = work(1:count_labels)
  end subroutine unique_labels

  pure function lowercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, code
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        value(i:i) = achar(code + 32)
      else
        value(i:i) = text(i:i)
      end if
    end do
  end function lowercase
end module rrcov_tests
