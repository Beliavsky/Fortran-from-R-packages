! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_metrics
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use classint_kinds, only: dp
    use classint_types, only: class_intervals, jenks_indices
    use classint_utils, only: mean_dp, sample_sd, unique_sorted
    implicit none
    private

    public :: find_cols, classify_intervals, class_counts, jenks_tests
    public :: gvf, tai, oai, classint_loglik, classint_aic, classint_n_partitions

contains

    pure subroutine find_cols(fit, cols)
        !! Assigns each retained observation to a fitted interval.
        type(class_intervals), intent(in) :: fit !! Fit defining the breaks and interval closure.
        integer, allocatable, intent(out) :: cols(:) !! One-based classes; zero marks an original non-finite value.
        integer :: i
        integer :: j
        integer :: k
        integer :: count_gt
        real(dp) :: value

        if (.not. allocated(fit%breaks)) error stop "find_cols: fit has no breaks"
        k = size(fit%breaks) - 1
        if (k < 1) error stop "find_cols: fewer than two breaks"
        allocate(cols(size(fit%values)), source=0)
        do i = 1, size(fit%values)
            value = fit%values(i)
            if (.not. ieee_is_finite(value)) cycle
            if (trim(fit%interval_closure) == "right") then
                count_gt = 0
                do j = 1, size(fit%breaks)
                    if (value > fit%breaks(j)) count_gt = count_gt + 1
                end do
                cols(i) = max(1, count_gt)
            else
                cols(i) = 1
                do j = 2, k
                    if (value >= fit%breaks(j)) cols(i) = j
                end do
                cols(i) = max(1, min(k, cols(i)))
            end if
        end do
    end subroutine find_cols

    pure subroutine classify_intervals(fit, cols)
        !! Returns the fitted class label for each retained observation.
        type(class_intervals), intent(in) :: fit !! Intervals applied to the retained observations.
        integer, allocatable, intent(out) :: cols(:) !! One-based labels; zero marks a non-finite observation.

        call find_cols(fit, cols)
    end subroutine classify_intervals

    pure function class_counts(fit) result(counts)
        !! Counts retained observations in each fitted interval.
        type(class_intervals), intent(in) :: fit !! Interval fit whose retained observations are counted by class.
        integer, allocatable :: counts(:)
        integer, allocatable :: cols(:)
        integer :: i
        integer :: k

        k = size(fit%breaks) - 1
        allocate(counts(k), source=0)
        call find_cols(fit, cols)
        do i = 1, size(cols)
            if (cols(i) >= 1 .and. cols(i) <= k) counts(cols(i)) = counts(cols(i)) + 1
        end do
    end function class_counts

    pure function gvf(var, cols) result(value)
        !! Computes Jenks' goodness of variance fit for a classification.
        real(dp), intent(in) :: var(:) !! Finite classified values used by the goodness-of-variance-fit calculation.
        integer, intent(in) :: cols(:) !! Positive one-based label for each value.
        real(dp) :: value
        real(dp) :: total_ss
        real(dp) :: within_ss
        real(dp), allocatable :: group(:)
        integer :: k
        integer :: c

        if (size(var) /= size(cols) .or. size(var) < 1) error stop "gvf: invalid shapes"
        k = maxval(cols)
        total_ss = sum((var - mean_dp(var))**2)
        within_ss = 0.0_dp
        do c = 1, k
            group = pack(var, cols == c)
            if (size(group) > 0) within_ss = within_ss + sum((group - mean_dp(group))**2)
        end do
        if (total_ss <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value = 1.0_dp - within_ss / total_ss
        end if
    end function gvf

    pure function tai(var, cols) result(value)
        !! Computes the tabular accuracy index for a classification.
        real(dp), intent(in) :: var(:) !! Finite classified numeric values used by the tabular accuracy index.
        integer, intent(in) :: cols(:) !! One-based class label for every value; same length as var.
        real(dp) :: value
        real(dp) :: total_abs
        real(dp) :: within_abs
        real(dp), allocatable :: group(:)
        integer :: k
        integer :: c

        if (size(var) /= size(cols) .or. size(var) < 1) error stop "tai: invalid shapes"
        k = maxval(cols)
        total_abs = sum(abs(var - mean_dp(var)))
        within_abs = 0.0_dp
        do c = 1, k
            group = pack(var, cols == c)
            if (size(group) > 0) within_abs = within_abs + sum(abs(group - mean_dp(group)))
        end do
        if (total_abs <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value = 1.0_dp - within_abs / total_abs
        end if
    end function tai

    pure function oai(var, cols, area) result(value)
        !! Computes the overview accuracy index using observation areas as weights.
        real(dp), intent(in) :: var(:) !! Finite classified numeric values used by the overview accuracy index.
        integer, intent(in) :: cols(:) !! One-based class label for each observation.
        real(dp), intent(in) :: area(:) !! Observation areas or weights; must have the same length as var.
        real(dp) :: value
        real(dp) :: total_abs
        real(dp) :: within_abs
        real(dp) :: mu
        integer :: k
        integer :: c
        integer :: i
        integer :: nclass

        if (size(var) /= size(cols) .or. size(var) /= size(area)) error stop "oai: invalid shapes"
        mu = mean_dp(var)
        total_abs = sum(abs(var - mu) * area)
        within_abs = 0.0_dp
        k = maxval(cols)
        do c = 1, k
            nclass = count(cols == c)
            if (nclass > 0) then
                mu = sum(var, mask=cols == c) / real(nclass, dp)
                do i = 1, size(var)
                    if (cols(i) == c) within_abs = within_abs + abs(var(i) - mu) * area(i)
                end do
            end if
        end do
        if (total_abs <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value = 1.0_dp - within_abs / total_abs
        end if
    end function oai

    pure function jenks_tests(fit, area) result(indices)
        !! Computes the Jenks and Armstrong accuracy indices for a fitted classification.
        type(class_intervals), intent(in) :: fit !! Classification assessed with Jenks/Armstrong indices.
        real(dp), intent(in), optional :: area(:) !! Polygon areas matching `fit%values`, when supplied.
        type(jenks_indices) :: indices
        integer, allocatable :: cols_all(:)
        integer, allocatable :: cols(:)
        real(dp), allocatable :: var(:)
        real(dp), allocatable :: weights(:)
        logical, allocatable :: keep(:)
        integer :: i

        call find_cols(fit, cols_all)
        allocate(keep(size(fit%values)))
        do i = 1, size(fit%values)
            keep(i) = ieee_is_finite(fit%values(i))
        end do
        var = pack(fit%values, keep)
        cols = pack(cols_all, keep)
        indices%n_classes = size(fit%breaks) - 1
        indices%goodness_of_fit = gvf(var, cols)
        indices%tabular_accuracy = tai(var, cols)
        if (present(area)) then
            if (size(area) /= size(fit%values)) error stop "jenks_tests: area length mismatch"
            weights = pack(area, keep)
            indices%overview_accuracy = oai(var, cols, weights)
            indices%has_overview = .true.
        end if
    end function jenks_tests

    pure elemental function classint_loglik(fit) result(loglik)
        !! Computes a Gaussian within-class log-likelihood for a fitted classification.
        type(class_intervals), intent(in) :: fit !! Fit evaluated with within-class Gaussian densities.
        real(dp) :: loglik
        real(dp), allocatable :: current(:)
        logical, allocatable :: mask(:)
        real(dp) :: lo
        real(dp) :: hi
        real(dp) :: mu
        real(dp) :: sd
        real(dp), allocatable :: unique_x(:)
        integer :: c
        integer :: k

        k = size(fit%breaks) - 1
        loglik = 0.0_dp
        allocate(mask(size(fit%values)))
        do c = 1, k
            lo = fit%breaks(c)
            hi = fit%breaks(c + 1)
            if ((c == 1 .and. trim(fit%interval_closure) == "right") .or. &
                (c == k .and. trim(fit%interval_closure) == "left")) then
                mask = ieee_is_finite(fit%values) .and. fit%values >= lo .and. fit%values <= hi
            else if (trim(fit%interval_closure) == "right") then
                mask = ieee_is_finite(fit%values) .and. fit%values > lo .and. fit%values <= hi
            else
                mask = ieee_is_finite(fit%values) .and. fit%values >= lo .and. fit%values < hi
            end if
            current = pack(fit%values, mask)
            if (size(current) == 0) cycle
            unique_x = unique_sorted(current)
            if (size(unique_x) == 1) cycle
            mu = mean_dp(current)
            sd = sample_sd(current)
            loglik = loglik - real(size(current), dp) * log(sd * sqrt(2.0_dp * acos(-1.0_dp))) - &
                     sum((current - mu)**2) / (2.0_dp * sd * sd)
        end do
    end function classint_loglik

    pure elemental function classint_aic(fit, penalty) result(aic)
        !! Computes an AIC-like score from the within-class Gaussian log-likelihood.
        type(class_intervals), intent(in) :: fit !! Fit providing the likelihood and class count.
        real(dp), intent(in), optional :: penalty !! Per-degree-of-freedom multiplier; defaults to two.
        real(dp) :: aic
        real(dp) :: kpen
        integer :: df

        kpen = 2.0_dp
        if (present(penalty)) kpen = penalty
        df = size(fit%breaks) - 1
        aic = -2.0_dp * classint_loglik(fit) + kpen * real(df, dp)
    end function classint_aic

    pure elemental function classint_n_partitions(fit) result(value)
        !! Computes the number of ways the fitted observations can be divided into the fitted classes.
        type(class_intervals), intent(in) :: fit !! Fit providing finite observation and class counts.
        real(dp) :: value
        integer :: n
        integer :: k
        integer :: r
        integer :: i

        n = fit%nobs
        k = size(fit%breaks) - 1
        if (n > 170) then
            value = ieee_value(0.0_dp, ieee_positive_inf)
            return
        end if
        if (k < 1 .or. k > n) then
            value = 0.0_dp
            return
        end if
        r = min(k - 1, n - k)
        value = 1.0_dp
        do i = 1, r
            value = value * real(n - i, dp) / real(i, dp)
        end do
    end function classint_n_partitions
end module classint_metrics
