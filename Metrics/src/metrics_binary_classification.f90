! SPDX-License-Identifier: BSD-3-Clause
module metrics_binary_classification
    use metrics_kinds, only : dp, metrics_success, metrics_invalid_size, metrics_invalid_argument
    use metrics_utils, only : average_ranks, quiet_nan, positive_infinity
    implicit none
    private

    public :: auc, ll, logloss, precision, recall, fbeta_score

contains

    real(dp) function auc(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:)
        real(dp), intent(in) :: predicted(:)
        integer, intent(out), optional :: stat
        real(dp), allocatable :: ranks(:)
        integer :: n_pos, n_neg

        if (size(actual) /= size(predicted) .or. size(actual) == 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        if (any(actual /= 0 .and. actual /= 1)) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if

        allocate(ranks(size(predicted)))
        call average_ranks(predicted, ranks)
        n_pos = count(actual == 1)
        n_neg = size(actual) - n_pos
        if (n_pos == 0 .or. n_neg == 0) then
            value = quiet_nan()
        else
            value = (sum(ranks, mask=actual == 1) - real(n_pos * (n_pos + 1), dp) / 2.0_dp) / &
                    real(n_pos * n_neg, dp)
        end if
        if (present(stat)) stat = metrics_success
    end function auc

    elemental real(dp) function ll(actual, predicted) result(value)
        integer, intent(in) :: actual
        real(dp), intent(in) :: predicted

        if ((actual /= 0 .and. actual /= 1) .or. predicted < 0.0_dp .or. predicted > 1.0_dp) then
            value = positive_infinity()
        else if (actual == 1) then
            if (predicted == 1.0_dp) then
                value = 0.0_dp
            else if (predicted == 0.0_dp) then
                value = positive_infinity()
            else
                value = -log(predicted)
            end if
        else
            if (predicted == 0.0_dp) then
                value = 0.0_dp
            else if (predicted == 1.0_dp) then
                value = positive_infinity()
            else
                value = -log(1.0_dp - predicted)
            end if
        end if
    end function ll

    real(dp) function logloss(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:)
        real(dp), intent(in) :: predicted(:)
        integer, intent(out), optional :: stat

        if (size(actual) /= size(predicted) .or. size(actual) == 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        value = sum(ll(actual, predicted)) / real(size(actual), dp)
        if (present(stat)) stat = metrics_success
    end function logloss

    real(dp) function precision(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: denominator

        if (.not. valid_binary_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        denominator = count(predicted == 1)
        if (denominator == 0) then
            value = quiet_nan()
        else
            value = real(count(actual == 1 .and. predicted == 1), dp) / real(denominator, dp)
        end if
    end function precision

    real(dp) function recall(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        integer :: denominator

        if (.not. valid_binary_pair(actual, predicted, stat)) then
            value = quiet_nan()
            return
        end if
        denominator = count(actual == 1)
        if (denominator == 0) then
            value = quiet_nan()
        else
            value = real(count(actual == 1 .and. predicted == 1), dp) / real(denominator, dp)
        end if
    end function recall

    real(dp) function fbeta_score(actual, predicted, beta, stat) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        real(dp), intent(in), optional :: beta
        integer, intent(out), optional :: stat
        real(dp) :: beta_value, prec, rec, denominator

        beta_value = 1.0_dp
        if (present(beta)) beta_value = beta
        if (beta_value < 0.0_dp) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if
        prec = precision(actual, predicted, stat)
        rec = recall(actual, predicted)
        denominator = beta_value**2 * prec + rec
        if (denominator == 0.0_dp) then
            value = quiet_nan()
        else
            value = (1.0_dp + beta_value**2) * prec * rec / denominator
        end if
    end function fbeta_score

    logical function valid_binary_pair(actual, predicted, stat) result(ok)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat

        ok = size(actual) == size(predicted) .and. size(actual) > 0
        if (ok) ok = all((actual == 0 .or. actual == 1) .and. (predicted == 0 .or. predicted == 1))
        if (present(stat)) then
            if (size(actual) /= size(predicted) .or. size(actual) == 0) then
                stat = metrics_invalid_size
            else if (.not. ok) then
                stat = metrics_invalid_argument
            else
                stat = metrics_success
            end if
        end if
    end function valid_binary_pair

end module metrics_binary_classification
