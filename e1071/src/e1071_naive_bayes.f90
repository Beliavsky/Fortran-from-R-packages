module e1071_naive_bayes
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use e1071_kinds, only: dp
    use e1071_constants, only: e1071_pi
    implicit none
    private

    type, public :: naive_bayes_model
        integer, allocatable :: class_labels(:)
        real(dp), allocatable :: apriori(:)
        real(dp), allocatable :: numeric_mean(:, :)
        real(dp), allocatable :: numeric_sd(:, :)
        integer, allocatable :: categorical_levels(:)
        real(dp), allocatable :: categorical_prob(:, :, :)
        real(dp) :: laplace = 0.0_dp
    end type naive_bayes_model

    public :: naive_bayes_fit, naive_bayes_predict

contains

    subroutine naive_bayes_fit(y, model, x_numeric, x_categorical, categorical_levels, laplace)
        integer, intent(in) :: y(:) !! Integer class labels for training observations; arbitrary distinct integer values are
        !! accepted.
        type(naive_bayes_model), intent(out) :: model !! Fitted class priors and conditional numeric/categorical distributions.
        real(dp), intent(in), optional :: x_numeric(:, :) !! Optional numeric predictors; rows must correspond to y and NaN
        !! values are omitted.
        integer, intent(in), optional :: x_categorical(:, :) !! Optional categorical predictors encoded 1..nlevels; nonpositive
        !! values are missing.
        integer, intent(in), optional :: categorical_levels(:) !! Number of levels for each categorical predictor; required with
        !! x_categorical.
        real(dp), intent(in), optional :: laplace !! Additive smoothing applied to categorical counts; defaults to zero and must
        !! be nonnegative.
        integer :: nclass
        integer :: nnum
        integer :: ncat
        integer :: maxlevel
        integer :: c
        integer :: j
        integer :: i
        integer :: level
        integer :: count_class
        integer :: ncomplete
        real(dp) :: smoothing
        real(dp) :: mu
        real(dp) :: ss

        if (size(y) < 1) error stop "naive_bayes_fit: y must not be empty"
        if (present(x_numeric)) then
            if (size(x_numeric, 1) /= size(y)) error stop "naive_bayes_fit: x_numeric row mismatch"
            nnum = size(x_numeric, 2)
        else
            nnum = 0
        end if
        if (present(x_categorical)) then
            if (size(x_categorical, 1) /= size(y)) error stop "naive_bayes_fit: x_categorical row mismatch"
            if (.not. present(categorical_levels)) error stop "naive_bayes_fit: categorical_levels is required"
            if (size(categorical_levels) /= size(x_categorical, 2)) error stop "naive_bayes_fit: categorical level mismatch"
            if (any(categorical_levels < 1)) error stop "naive_bayes_fit: categorical levels must be positive"
            ncat = size(x_categorical, 2)
            maxlevel = maxval(categorical_levels)
        else
            ncat = 0
            maxlevel = 0
        end if
        smoothing = 0.0_dp
        if (present(laplace)) smoothing = laplace
        if (smoothing < 0.0_dp) error stop "naive_bayes_fit: laplace must be nonnegative"

        call unique_sorted_labels(y, model%class_labels)
        nclass = size(model%class_labels)
        allocate(model%apriori(nclass))
        model%apriori = 0.0_dp
        do c = 1, nclass
            model%apriori(c) = real(count(y == model%class_labels(c)), dp)
        end do
        model%laplace = smoothing

        allocate(model%numeric_mean(nclass, nnum), model%numeric_sd(nclass, nnum))
        model%numeric_mean = 0.0_dp
        model%numeric_sd = 0.0_dp
        if (present(x_numeric)) then
            do c = 1, nclass
                do j = 1, nnum
                    mu = 0.0_dp
                    ncomplete = 0
                    do i = 1, size(y)
                        if (y(i) /= model%class_labels(c)) cycle
                        if (ieee_is_nan(x_numeric(i, j))) cycle
                        mu = mu + x_numeric(i, j)
                        ncomplete = ncomplete + 1
                    end do
                    if (ncomplete > 0) mu = mu / real(ncomplete, dp)
                    model%numeric_mean(c, j) = mu
                    ss = 0.0_dp
                    do i = 1, size(y)
                        if (y(i) /= model%class_labels(c)) cycle
                        if (ieee_is_nan(x_numeric(i, j))) cycle
                        ss = ss + (x_numeric(i, j) - mu)**2
                    end do
                    if (ncomplete > 1) then
                        model%numeric_sd(c, j) = sqrt(ss / real(ncomplete - 1, dp))
                    else
                        model%numeric_sd(c, j) = 0.0_dp
                    end if
                end do
            end do
        end if

        allocate(model%categorical_levels(ncat), model%categorical_prob(nclass, ncat, maxlevel))
        model%categorical_prob = 0.0_dp
        if (ncat > 0) model%categorical_levels = categorical_levels
        if (present(x_categorical)) then
            do c = 1, nclass
                count_class = count(y == model%class_labels(c))
                do j = 1, ncat
                    do level = 1, model%categorical_levels(j)
                        model%categorical_prob(c, j, level) = &
                            (real(count((y == model%class_labels(c)) .and. (x_categorical(:, j) == level)), dp) + smoothing) &
                            / (real(count_class, dp) + smoothing * real(model%categorical_levels(j), dp))
                    end do
                end do
            end do
        end if
    end subroutine naive_bayes_fit

    subroutine naive_bayes_predict(model, class, probability, x_numeric, x_categorical, threshold, eps)
        type(naive_bayes_model), intent(in) :: model !! Fitted naive Bayes model containing class priors and conditional
        !! distributions.
        integer, allocatable, intent(out) :: class(:) !! Predicted original integer class label for each supplied observation.
        real(dp), allocatable, intent(out), optional :: probability(:, :) !! Optional normalized posterior probabilities by row
        !! and model class order.
        real(dp), intent(in), optional :: x_numeric(:, :) !! Numeric predictors with the same column count used in fitting; NaNs
        !! contribute no factor.
        integer, intent(in), optional :: x_categorical(:, :) !! Categorical predictors encoded 1..nlevels;
        !! nonpositive/out-of-range values are missing.
        real(dp), intent(in), optional :: threshold !! Replacement probability for tiny/degenerate conditionals; defaults to 0.001.
        real(dp), intent(in), optional :: eps !! Conditional probabilities or standard deviations at or below this trigger
        !! replacement; default zero.
        real(dp), allocatable :: logp(:, :)
        real(dp) :: use_threshold
        real(dp) :: use_eps
        real(dp) :: sd
        real(dp) :: p
        real(dp) :: maxlog
        real(dp) :: total
        integer :: n
        integer :: nclass
        integer :: i
        integer :: c
        integer :: j
        integer :: level

        call prediction_shape(model, x_numeric, x_categorical, n)
        use_threshold = 0.001_dp
        if (present(threshold)) use_threshold = threshold
        use_eps = 0.0_dp
        if (present(eps)) use_eps = eps
        if (use_threshold <= 0.0_dp) error stop "naive_bayes_predict: threshold must be positive"
        nclass = size(model%class_labels)
        allocate(logp(n, nclass), class(n))
        do i = 1, n
            do c = 1, nclass
                logp(i, c) = log(max(model%apriori(c), tiny(1.0_dp)))
                if (present(x_numeric)) then
                    do j = 1, size(x_numeric, 2)
                        if (ieee_is_nan(x_numeric(i, j))) cycle
                        sd = model%numeric_sd(c, j)
                        if (sd <= use_eps) sd = use_threshold
                        p = normal_density(x_numeric(i, j), model%numeric_mean(c, j), sd)
                        if (p <= use_eps) p = use_threshold
                        logp(i, c) = logp(i, c) + log(p)
                    end do
                end if
                if (present(x_categorical)) then
                    do j = 1, size(x_categorical, 2)
                        level = x_categorical(i, j)
                        if (level < 1 .or. level > model%categorical_levels(j)) cycle
                        p = model%categorical_prob(c, j, level)
                        if (p <= use_eps) p = use_threshold
                        logp(i, c) = logp(i, c) + log(p)
                    end do
                end if
            end do
            c = maxloc(logp(i, :), dim=1)
            class(i) = model%class_labels(c)
        end do

        if (present(probability)) then
            allocate(probability(n, nclass))
            do i = 1, n
                maxlog = maxval(logp(i, :))
                probability(i, :) = exp(logp(i, :) - maxlog)
                total = sum(probability(i, :))
                probability(i, :) = probability(i, :) / total
            end do
        end if
    end subroutine naive_bayes_predict

    subroutine prediction_shape(model, x_numeric, x_categorical, n)
        type(naive_bayes_model), intent(in) :: model !! Fitted model defining expected numeric and categorical feature counts.
        real(dp), intent(in), optional :: x_numeric(:, :) !! Optional numeric predictor matrix whose dimensions are validated.
        integer, intent(in), optional :: x_categorical(:, :) !! Optional categorical predictor matrix whose dimensions are
        !! validated.
        integer, intent(out) :: n !! Number of prediction rows inferred from the supplied matrices.

        n = -1
        if (present(x_numeric)) then
            if (size(x_numeric, 2) /= size(model%numeric_mean, 2)) error stop "naive_bayes_predict: numeric column mismatch"
            n = size(x_numeric, 1)
        else if (size(model%numeric_mean, 2) > 0) then
            error stop "naive_bayes_predict: numeric predictors are required"
        end if
        if (present(x_categorical)) then
            if (size(x_categorical, 2) /= size(model%categorical_levels)) then
                error stop "naive_bayes_predict: categorical column mismatch"
            end if
            if (n >= 0 .and. size(x_categorical, 1) /= n) error stop "naive_bayes_predict: row mismatch"
            n = size(x_categorical, 1)
        else if (size(model%categorical_levels) > 0) then
            error stop "naive_bayes_predict: categorical predictors are required"
        end if
        if (n < 0) error stop "naive_bayes_predict: no predictor matrix supplied"
    end subroutine prediction_shape

    pure function normal_density(x, mean_value, sd) result(value)
        real(dp), intent(in) :: x !! Scalar observation at which the normal density is evaluated.
        real(dp), intent(in) :: mean_value !! Mean of the normal conditional distribution.
        real(dp), intent(in) :: sd !! Positive standard deviation of the normal conditional distribution.
        real(dp) :: value
        real(dp) :: z

        z = (x - mean_value) / sd
        value = exp(-0.5_dp * z * z) / (sqrt(2.0_dp * e1071_pi) * sd)
    end function normal_density

    subroutine unique_sorted_labels(y, labels)
        integer, intent(in) :: y(:) !! Integer class labels from which sorted distinct values are extracted.
        integer, allocatable, intent(out) :: labels(:) !! Sorted distinct labels occurring in y.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: key

        work = y
        do i = 2, size(work)
            key = work(i)
            j = i - 1
            do while (j >= 1)
                if (work(j) <= key) exit
                work(j + 1) = work(j)
                j = j - 1
            end do
            work(j + 1) = key
        end do
        n = 1
        do i = 2, size(work)
            if (work(i) /= work(i - 1)) n = n + 1
        end do
        allocate(labels(n))
        labels(1) = work(1)
        j = 1
        do i = 2, size(work)
            if (work(i) /= work(i - 1)) then
                j = j + 1
                labels(j) = work(i)
            end if
        end do
    end subroutine unique_sorted_labels

end module e1071_naive_bayes
