! SPDX-License-Identifier: GPL-3.0-only
module pwev_metrics
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use pwev_kinds, only : dp
  implicit none
  private
  public :: pwev_accuracy_table, pwev_metric_vector
contains

  subroutine pwev_metric_vector(actual, predicted, metrics)
    real(dp), intent(in) :: actual(:), predicted(:)
    real(dp), intent(out) :: metrics(9)
    real(dp), allocatable :: error(:), absolute_error(:), sorted(:)
    real(dp) :: mean_actual, denominator, corr, nan_value
    integer :: n

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    n = size(actual)
    metrics = nan_value
    if (n <= 0 .or. size(predicted) /= n) return
    allocate(error(n), absolute_error(n), sorted(n))
    error = actual - predicted
    absolute_error = abs(error)
    mean_actual = sum(actual) / real(n, dp)

    metrics(1) = sqrt(sum(error**2) / real(n, dp))
    if (any(abs(actual) <= tiny(1.0_dp))) then
      metrics(2) = nan_value
    else
      metrics(2) = sum(absolute_error / abs(actual)) / real(n, dp)
    end if
    metrics(3) = sum(absolute_error) / real(n, dp)
    denominator = sum((actual - mean_actual)**2)
    if (denominator > 0.0_dp) metrics(4) = sqrt(sum(error**2) / denominator)
    sorted = absolute_error
    call insertion_sort(sorted)
    if (mod(n, 2) == 1) then
      metrics(5) = sorted((n + 1) / 2)
    else
      metrics(5) = 0.5_dp * (sorted(n / 2) + sorted(n / 2 + 1))
    end if
    if (all(actual > -1.0_dp) .and. all(predicted > -1.0_dp)) then
      metrics(6) = sqrt(sum((log1p_safe(actual) - log1p_safe(predicted))**2) / real(n, dp))
    end if
    denominator = sum(abs(actual - mean_actual))
    if (denominator > 0.0_dp) metrics(7) = sum(absolute_error) / denominator
    if (all(abs(actual) + abs(predicted) > 0.0_dp)) then
      metrics(8) = sum(2.0_dp * absolute_error / (abs(actual) + abs(predicted))) / real(n, dp)
    end if
    corr = correlation(actual, predicted)
    if (ieee_is_finite(corr)) metrics(9) = corr * corr
  end subroutine pwev_metric_vector

  subroutine pwev_accuracy_table(train_actual, train_models, test_actual, test_models, accuracy, round_values)
    real(dp), intent(in) :: train_actual(:), train_models(:, :)
    real(dp), intent(in) :: test_actual(:), test_models(:, :)
    real(dp), allocatable, intent(out) :: accuracy(:, :)
    logical, intent(in), optional :: round_values
    real(dp) :: values(9)
    integer :: j
    logical :: do_round

    do_round = .true.
    if (present(round_values)) do_round = round_values
    allocate(accuracy(size(train_models, 2), 18))
    do j = 1, size(train_models, 2)
      call pwev_metric_vector(train_actual, train_models(:, j), values)
      accuracy(j, 1:9) = values
      call pwev_metric_vector(test_actual, test_models(:, j), values)
      accuracy(j, 10:18) = values
    end do
    if (do_round) then
      where (ieee_is_finite(accuracy)) accuracy = anint(10000.0_dp * accuracy) / 10000.0_dp
    end if
  end subroutine pwev_accuracy_table

  pure elemental real(dp) function log1p_safe(x) result(value)
    real(dp), intent(in) :: x
    value = log(1.0_dp + x)
  end function log1p_safe

  real(dp) function correlation(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: mx, my, sx, sy
    mx = sum(x) / real(size(x), dp)
    my = sum(y) / real(size(y), dp)
    sx = sum((x - mx)**2)
    sy = sum((y - my)**2)
    if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = sum((x - mx) * (y - my)) / sqrt(sx * sy)
    end if
  end function correlation

  subroutine insertion_sort(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: key
    integer :: i, j
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine insertion_sort

end module pwev_metrics
