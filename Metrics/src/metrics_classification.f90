! SPDX-License-Identifier: BSD-3-Clause
module metrics_classification
    use metrics_kinds, only : dp, metrics_success, metrics_invalid_size, metrics_invalid_argument
    use metrics_utils, only : quiet_nan
    implicit none
    private

    public :: ce, accuracy
    public :: score_quadratic_weighted_kappa, mean_quadratic_weighted_kappa
    public :: scorequadraticweightedkappa, meanquadraticweightedkappa

    interface ce
        module procedure ce_integer
        module procedure ce_real
        module procedure ce_character
    end interface ce

    interface accuracy
        module procedure accuracy_integer
        module procedure accuracy_real
        module procedure accuracy_character
    end interface accuracy

    interface scorequadraticweightedkappa
        module procedure score_quadratic_weighted_kappa
    end interface scorequadraticweightedkappa

    interface meanquadraticweightedkappa
        module procedure mean_quadratic_weighted_kappa
    end interface meanquadraticweightedkappa

contains

    real(dp) function ce_integer(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        if (.not. same_nonempty_size(size(actual), size(predicted), stat)) then
            value = quiet_nan()
        else
            value = real(count(actual /= predicted), dp) / real(size(actual), dp)
        end if
    end function ce_integer

    real(dp) function ce_real(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        if (.not. same_nonempty_size(size(actual), size(predicted), stat)) then
            value = quiet_nan()
        else
            value = real(count(actual /= predicted), dp) / real(size(actual), dp)
        end if
    end function ce_real

    real(dp) function ce_character(actual, predicted, stat) result(value)
        character(len=*), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        if (.not. same_nonempty_size(size(actual), size(predicted), stat)) then
            value = quiet_nan()
        else
            value = real(count(actual /= predicted), dp) / real(size(actual), dp)
        end if
    end function ce_character

    real(dp) function accuracy_integer(actual, predicted, stat) result(value)
        integer, intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = 1.0_dp - ce_integer(actual, predicted, stat)
    end function accuracy_integer

    real(dp) function accuracy_real(actual, predicted, stat) result(value)
        real(dp), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = 1.0_dp - ce_real(actual, predicted, stat)
    end function accuracy_real

    real(dp) function accuracy_character(actual, predicted, stat) result(value)
        character(len=*), intent(in) :: actual(:), predicted(:)
        integer, intent(out), optional :: stat
        value = 1.0_dp - ce_character(actual, predicted, stat)
    end function accuracy_character

    real(dp) function score_quadratic_weighted_kappa(rater_a, rater_b, min_rating, max_rating, stat) result(value)
        integer, intent(in) :: rater_a(:), rater_b(:)
        integer, intent(in), optional :: min_rating, max_rating
        integer, intent(out), optional :: stat
        integer :: lower, upper, ncat, i, ia, ib
        real(dp), allocatable :: observed(:,:), expected(:,:), hist_a(:), hist_b(:), weights(:,:)
        real(dp) :: numerator, denominator

        if (size(rater_a) /= size(rater_b) .or. size(rater_a) == 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        lower = min(minval(rater_a), minval(rater_b))
        upper = max(maxval(rater_a), maxval(rater_b))
        if (present(min_rating)) lower = min_rating
        if (present(max_rating)) upper = max_rating
        if (upper < lower .or. any(rater_a < lower) .or. any(rater_a > upper) .or. &
            any(rater_b < lower) .or. any(rater_b > upper)) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if

        ncat = upper - lower + 1
        allocate(observed(ncat, ncat), expected(ncat, ncat), hist_a(ncat), hist_b(ncat), weights(ncat, ncat))
        observed = 0.0_dp
        hist_a = 0.0_dp
        hist_b = 0.0_dp
        do i = 1, size(rater_a)
            ia = rater_a(i) - lower + 1
            ib = rater_b(i) - lower + 1
            observed(ia, ib) = observed(ia, ib) + 1.0_dp
            hist_a(ia) = hist_a(ia) + 1.0_dp
            hist_b(ib) = hist_b(ib) + 1.0_dp
        end do
        observed = observed / real(size(rater_a), dp)
        hist_a = hist_a / real(size(rater_a), dp)
        hist_b = hist_b / real(size(rater_b), dp)
        do ia = 1, ncat
            do ib = 1, ncat
                expected(ia, ib) = hist_a(ia) * hist_b(ib)
                weights(ia, ib) = real((ia - ib)**2, dp)
            end do
        end do
        numerator = sum(weights * observed)
        denominator = sum(weights * expected)
        if (denominator == 0.0_dp) then
            value = quiet_nan()
        else
            value = 1.0_dp - numerator / denominator
        end if
        if (present(stat)) stat = metrics_success
    end function score_quadratic_weighted_kappa

    real(dp) function mean_quadratic_weighted_kappa(kappas, weights, stat) result(value)
        real(dp), intent(in) :: kappas(:)
        real(dp), intent(in), optional :: weights(:)
        integer, intent(out), optional :: stat
        real(dp), allocatable :: normalized_weights(:), clipped(:), z(:)
        real(dp) :: mean_weight, zmean
        integer :: i

        if (size(kappas) == 0) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_size
            return
        end if
        allocate(normalized_weights(size(kappas)), clipped(size(kappas)), z(size(kappas)))
        if (present(weights)) then
            if (size(weights) /= size(kappas)) then
                value = quiet_nan()
                if (present(stat)) stat = metrics_invalid_size
                return
            end if
            normalized_weights = weights
        else
            normalized_weights = 1.0_dp
        end if
        mean_weight = sum(normalized_weights) / real(size(normalized_weights), dp)
        if (mean_weight == 0.0_dp) then
            value = quiet_nan()
            if (present(stat)) stat = metrics_invalid_argument
            return
        end if
        normalized_weights = normalized_weights / mean_weight

        do i = 1, size(kappas)
            if (kappas(i) == 0.0_dp) then
                clipped(i) = 0.0_dp
            else
                clipped(i) = sign(min(0.999_dp, abs(kappas(i))), kappas(i))
                clipped(i) = sign(max(0.001_dp, abs(clipped(i))), clipped(i))
            end if
            z(i) = 0.5_dp * log((1.0_dp + clipped(i)) / (1.0_dp - clipped(i)))
        end do
        zmean = sum(z * normalized_weights) / real(size(z), dp)
        value = tanh(zmean)
        if (present(stat)) stat = metrics_success
    end function mean_quadratic_weighted_kappa

    logical function same_nonempty_size(n1, n2, stat) result(ok)
        integer, intent(in) :: n1, n2
        integer, intent(out), optional :: stat
        ok = n1 == n2 .and. n1 > 0
        if (present(stat)) then
            if (ok) then
                stat = metrics_success
            else
                stat = metrics_invalid_size
            end if
        end if
    end function same_nonempty_size

end module metrics_classification
