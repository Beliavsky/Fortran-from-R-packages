! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 multivariate amputation helpers.
module mice_ampute
    use r_kinds, only : dp
    use mice_rng, only : mice_rng_state, rng_normal, rng_uniform
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape
    implicit none
    private

    integer, parameter, public :: ampute_right = 1
    integer, parameter, public :: ampute_left = 2
    integer, parameter, public :: ampute_mid = 3
    integer, parameter, public :: ampute_tail = 4

    public :: ampute_mcar
    public :: ampute_continuous
    public :: ampute_discrete

contains

    subroutine ampute_mcar(assignments, patterns, prop, rng, missing, info)
        integer, intent(in) :: assignments(:) !! One-based candidate-pattern assignment for each case.
        integer, intent(in) :: patterns(:, :) !! Pattern-by-variable matrix with zero marking cells to amputate and one to retain.
        real(dp), intent(in), value :: prop !! Probability that a candidate case receives its assigned missingness pattern.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for Bernoulli amputation draws.
        logical, allocatable, intent(out) :: missing(:, :) !! Generated case-by-variable missingness mask.
        integer, intent(out) :: info !! `mice_ok` on success or an argument/shape status code.
        integer :: i, j, pat

        if (prop < 0.0_dp .or. prop > 1.0_dp .or. size(patterns, 1) < 1 .or. size(patterns, 2) < 1) then
            info = mice_invalid_argument
            return
        end if
        if (any(assignments < 1) .or. any(assignments > size(patterns, 1))) then
            info = mice_invalid_shape
            return
        end if
        allocate(missing(size(assignments), size(patterns, 2)))
        missing = .false.
        do i = 1, size(assignments)
            if (rng_uniform(rng) > prop) cycle
            pat = assignments(i)
            do j = 1, size(patterns, 2)
                missing(i, j) = patterns(pat, j) == 0
            end do
        end do
        info = mice_ok
    end subroutine ampute_mcar

    subroutine ampute_continuous(assignments, patterns, scores, prop, types, rng, missing, info)
        integer, intent(in) :: assignments(:) !! One-based candidate-pattern assignment for each case.
        integer, intent(in) :: patterns(:, :) !! Pattern-by-variable matrix with zero marking cells to amputate.
        real(dp), intent(in) :: scores(:, :) !! Weighted sum scores; column `j` is used for candidates assigned to pattern `j`.
        real(dp), intent(in), value :: prop !! Target probability of amputation used to calibrate each probability curve.
        integer, intent(in) :: types(:) !! Curve type per pattern: RIGHT=1, LEFT=2, MID=3, TAIL=4.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for normal calibration and Bernoulli amputation draws.
        logical, allocatable, intent(out) :: missing(:, :) !! Generated case-by-variable missingness mask.
        integer, intent(out) :: info !! `mice_ok` on success or an argument/shape status code.

        real(dp), allocatable :: testset(:), candidate(:)
        real(dp) :: mean_test, sd_test, mean_score, shift, probability
        integer :: i, j, k, n_candidate, pat

        if (size(scores, 1) /= size(assignments) .or. size(scores, 2) /= size(patterns, 1) .or. &
            size(types) /= size(patterns, 1)) then
            info = mice_invalid_shape
            return
        end if
        if (prop < 0.0_dp .or. prop > 1.0_dp .or. any(types < ampute_right) .or. any(types > ampute_tail)) then
            info = mice_invalid_argument
            return
        end if
        if (any(assignments < 1) .or. any(assignments > size(patterns, 1))) then
            info = mice_invalid_shape
            return
        end if
        allocate(missing(size(assignments), size(patterns, 2)), testset(10000))
        missing = .false.
        do i = 1, size(testset)
            testset(i) = rng_normal(rng)
        end do
        mean_test = sum(testset) / real(size(testset), dp)
        sd_test = sqrt(sum((testset - mean_test)**2) / real(size(testset) - 1, dp))
        if (sd_test > 0.0_dp) testset = (testset - mean_test) / sd_test

        do pat = 1, size(patterns, 1)
            n_candidate = count(assignments == pat)
            if (n_candidate < 1) cycle
            allocate(candidate(n_candidate))
            k = 0
            do i = 1, size(assignments)
                if (assignments(i) /= pat) cycle
                k = k + 1
                candidate(k) = scores(i, pat)
            end do
            mean_score = sum(candidate) / real(n_candidate, dp)
            if (n_candidate == 1 .or. nearly_constant(candidate)) then
                shift = 0.0_dp
            else
                shift = calibrate_shift(testset, prop, types(pat))
            end if
            do i = 1, size(assignments)
                if (assignments(i) /= pat) cycle
                if (n_candidate == 1 .or. nearly_constant(candidate)) then
                    probability = prop
                else
                    probability = curve_probability(scores(i, pat), mean_score, shift, types(pat))
                end if
                if (rng_uniform(rng) <= probability) then
                    do j = 1, size(patterns, 2)
                        missing(i, j) = patterns(pat, j) == 0
                    end do
                end if
            end do
            deallocate(candidate)
        end do
        info = mice_ok
    end subroutine ampute_continuous

    subroutine ampute_discrete(assignments, patterns, scores, prop, odds, rng, missing, info)
        integer, intent(in) :: assignments(:) !! One-based candidate-pattern assignment for each case.
        integer, intent(in) :: patterns(:, :) !! Pattern-by-variable matrix with zero marking cells to amputate.
        real(dp), intent(in) :: scores(:, :) !! Weighted sum scores by case and pattern.
        real(dp), intent(in), value :: prop !! Target overall amputation proportion.
        real(dp), intent(in) :: odds(:, :) !! Positive relative missingness odds for equal-frequency score groups by pattern.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for Bernoulli amputation draws.
        logical, allocatable, intent(out) :: missing(:, :) !! Generated case-by-variable missingness mask.
        integer, intent(out) :: info !! `mice_ok` on success or an argument/shape status code.

        real(dp), allocatable :: candidate(:), sorted(:)
        real(dp) :: probability, sum_odds
        integer :: group, i, j, k, n_candidate, ng, pat, position

        if (size(scores, 1) /= size(assignments) .or. size(scores, 2) /= size(patterns, 1) .or. &
            size(odds, 1) /= size(patterns, 1)) then
            info = mice_invalid_shape
            return
        end if
        if (prop < 0.0_dp .or. prop > 1.0_dp .or. any(odds < 0.0_dp)) then
            info = mice_invalid_argument
            return
        end if
        allocate(missing(size(assignments), size(patterns, 2)))
        missing = .false.
        do pat = 1, size(patterns, 1)
            n_candidate = count(assignments == pat)
            if (n_candidate < 1) cycle
            ng = count(odds(pat, :) > 0.0_dp)
            if (ng < 1) cycle
            sum_odds = sum(odds(pat, 1:ng))
            if (sum_odds <= 0.0_dp) cycle
            allocate(candidate(n_candidate), sorted(n_candidate))
            k = 0
            do i = 1, size(assignments)
                if (assignments(i) /= pat) cycle
                k = k + 1
                candidate(k) = scores(i, pat)
            end do
            sorted = candidate
            call sort_real(sorted)
            k = 0
            do i = 1, size(assignments)
                if (assignments(i) /= pat) cycle
                k = k + 1
                position = rank_right(sorted, candidate(k))
                group = min(ng, max(1, 1 + int(real(position - 1, dp) * real(ng, dp) / real(n_candidate, dp))))
                probability = real(ng, dp) * prop * odds(pat, group) / sum_odds
                probability = min(1.0_dp, max(0.0_dp, probability))
                if (rng_uniform(rng) <= probability) then
                    do j = 1, size(patterns, 2)
                        missing(i, j) = patterns(pat, j) == 0
                    end do
                end if
            end do
            deallocate(candidate, sorted)
        end do
        info = mice_ok
    end subroutine ampute_discrete

    pure real(dp) function curve_probability(x, mean_x, shift, curve_type) result(probability)
        real(dp), intent(in), value :: x !! Weighted score at which the missingness probability is evaluated.
        real(dp), intent(in), value :: mean_x !! Mean score in the population used by the chosen probability curve.
        real(dp), intent(in), value :: shift !! Additive logit shift calibrated to the target proportion.
        integer, intent(in), value :: curve_type !! RIGHT, LEFT, MID, or TAIL integer curve code.
        real(dp) :: eta

        select case (curve_type)
        case (ampute_left)
            eta = mean_x - x + shift
        case (ampute_mid)
            eta = -abs(x - mean_x) + 0.75_dp + shift
        case (ampute_tail)
            eta = abs(x - mean_x) - 0.75_dp + shift
        case default
            eta = -mean_x + x + shift
        end select
        if (eta >= 0.0_dp) then
            probability = 1.0_dp / (1.0_dp + exp(-eta))
        else
            probability = exp(eta) / (1.0_dp + exp(eta))
        end if
    end function curve_probability

    pure real(dp) function calibrate_shift(testset, target, curve_type) result(shift)
        real(dp), intent(in) :: testset(:) !! Standardized normal calibration sample analogous to upstream `testset`.
        real(dp), intent(in), value :: target !! Desired mean missingness probability.
        integer, intent(in), value :: curve_type !! RIGHT, LEFT, MID, or TAIL integer curve code.
        real(dp) :: center, hi, lo, mean_test, value
        integer :: iter

        lo = -8.0_dp
        hi = 8.0_dp
        mean_test = sum(testset) / real(size(testset), dp)
        do iter = 1, 100
            center = 0.5_dp * (lo + hi)
            value = mean_curve_probability(testset, mean_test, center, curve_type)
            if (value < target) then
                lo = center
            else
                hi = center
            end if
            if (hi - lo <= 1.0e-3_dp) exit
        end do
        shift = 0.5_dp * (lo + hi)
    end function calibrate_shift

    pure real(dp) function mean_curve_probability(x, mean_x, shift, curve_type) result(value)
        real(dp), intent(in) :: x(:) !! Score vector used to average a probability curve.
        real(dp), intent(in), value :: mean_x !! Mean score used by LEFT/RIGHT/MID/TAIL transformations.
        real(dp), intent(in), value :: shift !! Candidate additive logit shift.
        integer, intent(in), value :: curve_type !! RIGHT, LEFT, MID, or TAIL integer curve code.
        integer :: i

        value = 0.0_dp
        do i = 1, size(x)
            value = value + curve_probability(x(i), mean_x, shift, curve_type)
        end do
        value = value / real(size(x), dp)
    end function mean_curve_probability

    pure logical function nearly_constant(values) result(constant)
        real(dp), intent(in) :: values(:) !! Values tested for numerical equality over their full range.
        real(dp) :: scale

        if (size(values) < 2) then
            constant = .true.
            return
        end if
        scale = max(1.0_dp, maxval(abs(values)))
        constant = maxval(values) - minval(values) <= 16.0_dp * epsilon(1.0_dp) * scale
    end function nearly_constant

    pure subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! Real vector sorted in ascending order in place.
        real(dp) :: key
        integer :: i, j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real

    pure integer function rank_right(sorted, value) result(position)
        real(dp), intent(in) :: sorted(:) !! Ascending score vector.
        real(dp), intent(in), value :: value !! Score whose rightmost rank is requested.
        integer :: i

        position = 1
        do i = 1, size(sorted)
            if (sorted(i) <= value) position = i
        end do
    end function rank_right

end module mice_ampute
