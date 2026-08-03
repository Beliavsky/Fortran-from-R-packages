! SPDX-License-Identifier: BSD-3-Clause
module metrics_time_series
    use metrics_kinds, only : dp, metrics_success, metrics_invalid_size, metrics_invalid_argument
    use metrics_utils, only : quiet_nan, positive_infinity
    implicit none
    private

    public :: mase

contains

    real(dp) function mase(actual, predicted, step_size, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(in), optional :: step_size
        integer, intent(out), optional :: stat
        integer :: step, n, naive_end
        real(dp) :: sum_errors, naive_errors, denominator

        step = 1
        if (present(step_size)) step = step_size
        n = size(actual)
        if (n /= size(predicted) .or. n == 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        if (step < 1 .or. step >= n) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if

        naive_end = n - step
        sum_errors = sum(abs(actual - predicted))
        naive_errors = sum(abs(actual(step + 1:n) - actual(1:naive_end)))
        denominator = real(n, dp) * naive_errors / real(naive_end, dp)
        if (denominator == 0.0_dp) then
            if (sum_errors == 0.0_dp) then
                value = quiet_nan()
            else
                value = positive_infinity()
            end if
        else
            value = sum_errors / denominator
        end if
        if (present(stat)) stat = metrics_success
    end function mase

end module metrics_time_series
