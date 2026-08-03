! SPDX-License-Identifier: BSD-3-Clause
module metrics_regression
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use metrics_kinds, only : dp, metrics_success, metrics_invalid_size, metrics_invalid_argument
    use metrics_utils, only : mean_value, median_value, quiet_nan, positive_infinity, negative_infinity
    implicit none
    private

    public :: bias, percent_bias, se, sse, mse, rmse
    public :: ae, mae, mdae, ape, mape, smape
    public :: sle, msle, rmsle, rse, rrse, rae, explained_variation

contains

    elemental real(dp) function se(actual, predicted) result(value)
        real(dp), intent(in) :: actual, predicted
        value = (actual - predicted)**2
    end function se

    elemental real(dp) function ae(actual, predicted) result(value)
        real(dp), intent(in) :: actual, predicted
        value = abs(actual - predicted)
    end function ae

    elemental real(dp) function ape(actual, predicted) result(value)
        real(dp), intent(in) :: actual, predicted
        real(dp) :: numerator

        numerator = abs(actual - predicted)
        if (actual /= 0.0_dp) then
            value = numerator / abs(actual)
        else if (numerator == 0.0_dp) then
            value = quiet_nan()
        else
            value = positive_infinity()
        end if
    end function ape

    elemental real(dp) function sle(actual, predicted) result(value)
        real(dp), intent(in) :: actual, predicted
        value = (log(1.0_dp + actual) - log(1.0_dp + predicted))**2
    end function sle

    real(dp) function bias(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = mean_value(actual - predicted)
    end function bias

    real(dp) function percent_bias(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        real(dp), allocatable :: terms(:)
        integer :: i

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        allocate(terms(size(actual)))
        do i = 1, size(actual)
            terms(i) = signed_ratio(actual(i) - predicted(i), abs(actual(i)))
        end do
        value = mean_value(terms)
    end function percent_bias

    real(dp) function sse(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = sum(se(actual, predicted))
    end function sse

    real(dp) function mse(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = mean_value(se(actual, predicted))
    end function mse

    real(dp) function rmse(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = sqrt(mse(actual, predicted, stat))
    end function rmse

    real(dp) function mae(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = mean_value(ae(actual, predicted))
    end function mae

    real(dp) function mdae(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = median_value(ae(actual, predicted))
    end function mdae

    real(dp) function mape(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        value = mean_value(ape(actual, predicted))
    end function mape

    real(dp) function smape(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        real(dp), allocatable :: terms(:)
        real(dp) :: denominator
        integer :: i

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        allocate(terms(size(actual)))
        do i = 1, size(actual)
            denominator = abs(actual(i)) + abs(predicted(i))
            if (denominator == 0.0_dp) then
                terms(i) = quiet_nan()
            else
                terms(i) = 2.0_dp * abs(actual(i) - predicted(i)) / denominator
            end if
        end do
        value = mean_value(terms)
    end function smape

    real(dp) function msle(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        if (any(actual < -1.0_dp) .or. any(predicted < -1.0_dp)) then
            if (present(stat)) stat = metrics_invalid_argument
        end if
        value = mean_value(sle(actual, predicted))
    end function msle

    real(dp) function rmsle(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = sqrt(msle(actual, predicted, stat))
    end function rmsle

    real(dp) function rse(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        real(dp) :: center, denominator

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        center = mean_value(actual)
        denominator = sum((actual - center)**2)
        value = signed_ratio(sum((actual - predicted)**2), denominator)
    end function rse

    real(dp) function rrse(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = sqrt(rse(actual, predicted, stat))
    end function rrse

    real(dp) function rae(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        real(dp) :: center, denominator

        if (.not. valid_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        center = mean_value(actual)
        denominator = sum(abs(actual - center))
        value = signed_ratio(sum(abs(actual - predicted)), denominator)
    end function rae

    real(dp) function explained_variation(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = 1.0_dp - rse(actual, predicted, stat)
    end function explained_variation

    logical function valid_pair(actual, predicted, stat) result(ok)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        ok = size(actual) == size(predicted) .and. size(actual) > 0
        if (present(stat)) then
            if (ok) then
                stat = metrics_success
            else
                stat = metrics_invalid_size
            end if
        end if
    end function valid_pair

    elemental real(dp) function signed_ratio(numerator, denominator) result(value)
        real(dp), intent(in) :: numerator, denominator

        if (denominator /= 0.0_dp) then
            value = numerator / denominator
        else if (numerator > 0.0_dp) then
            value = positive_infinity()
        else if (numerator < 0.0_dp) then
            value = negative_infinity()
        else
            value = quiet_nan()
        end if
    end function signed_ratio

end module metrics_regression
