module e1071_lca
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_uniform
    use e1071_special, only: normal_cdf, chi_square_sf
    use e1071_utils, only: bincombinations
    implicit none
    private

    type, public :: lca_model
        real(dp), allocatable :: class_probability(:)
        real(dp), allocatable :: item_probability(:, :)
        integer, allocatable :: pattern_class(:)
        real(dp) :: log_likelihood = 0.0_dp
        real(dp) :: saturated_log_likelihood = 0.0_dp
        real(dp) :: chi_square = 0.0_dp
        real(dp) :: likelihood_ratio = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: saturated_bic = 0.0_dp
        integer :: n = 0
        integer :: n_parameters = 0
    end type lca_model

    type, public :: lca_bootstrap_result
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: saturated_log_likelihood(:)
        real(dp), allocatable :: likelihood_ratio(:)
        real(dp), allocatable :: chi_square(:)
        real(dp) :: likelihood_ratio_mean = 0.0_dp
        real(dp) :: likelihood_ratio_sd = 0.0_dp
        real(dp) :: chi_square_mean = 0.0_dp
        real(dp) :: chi_square_sd = 0.0_dp
        real(dp) :: z_likelihood_ratio = 0.0_dp
        real(dp) :: z_chi_square = 0.0_dp
        real(dp) :: p_z_likelihood_ratio = 0.0_dp
        real(dp) :: p_z_chi_square = 0.0_dp
        real(dp) :: p_likelihood_ratio = 0.0_dp
        real(dp) :: p_chi_square = 0.0_dp
    end type lca_bootstrap_result

    public :: countpattern, lca_fit_counts, lca_fit_data, lca_predict
    public :: lca_bootstrap, lca_goodness_of_fit

contains

    subroutine countpattern(x, counts, matching)
        integer, intent(in) :: x(:, :) !! Binary observation-by-item matrix; every entry must be zero or one.
        integer, allocatable, intent(out) :: counts(:) !! Counts for all 2**nitems binary patterns in lexicographic binary order.
        integer, allocatable, intent(out), optional :: matching(:) !! Optional one-based pattern index for each input observation.
        integer :: npattern
        integer :: i
        integer :: j
        integer :: index

        if (any(x < 0) .or. any(x > 1)) error stop "countpattern: x must be binary"
        npattern = 2**size(x, 2)
        allocate(counts(npattern))
        counts = 0
        if (present(matching)) allocate(matching(size(x, 1)))
        do i = 1, size(x, 1)
            index = x(i, 1)
            do j = 2, size(x, 2)
                index = 2 * index + x(i, j)
            end do
            index = index + 1
            counts(index) = counts(index) + 1
            if (present(matching)) matching(i) = index
        end do
    end subroutine countpattern

    subroutine lca_fit_counts(counts, k, niter, rng, model)
        integer, intent(in) :: counts(:) !! Frequency vector over every binary response pattern; length must be a power of two.
        integer, intent(in) :: k !! Number of latent classes; must be positive.
        integer, intent(in) :: niter !! Number of EM updates performed exactly, matching the e1071 interface.
        type(rng_state), intent(inout) :: rng !! Mutable random generator used for uniform initialization of probabilities.
        type(lca_model), intent(out) :: model !! Fitted latent-class probabilities, item probabilities, assignments, and fit
        !! statistics.
        integer, allocatable :: patterns(:, :)
        real(dp), allocatable :: posterior(:, :)
        real(dp), allocatable :: joint(:, :)
        real(dp), allocatable :: weighted(:, :)
        real(dp), allocatable :: pattern_probability(:)
        real(dp), allocatable :: class_size(:)
        real(dp) :: total_class
        real(dp) :: likelihood
        real(dp) :: saturated
        real(dp) :: expected
        integer :: nvar
        integer :: npattern
        integer :: n
        integer :: i
        integer :: j
        integer :: v
        integer :: iter

        if (k < 1 .or. niter < 1) error stop "lca_fit_counts: k and niter must be positive"
        if (any(counts < 0)) error stop "lca_fit_counts: counts must be nonnegative"
        npattern = size(counts)
        nvar = exact_log2(npattern)
        if (nvar < 0) error stop "lca_fit_counts: count length must be a power of two"
        n = sum(counts)
        if (n <= 0) error stop "lca_fit_counts: total count must be positive"
        patterns = bincombinations(nvar)
        allocate(model%class_probability(k), model%item_probability(k, nvar))
        allocate(posterior(k, npattern), joint(k, npattern), weighted(k, npattern))
        allocate(pattern_probability(npattern), class_size(k), model%pattern_class(npattern))

        do j = 1, k
            model%class_probability(j) = rng_uniform(rng)
            do v = 1, nvar
                model%item_probability(j, v) = rng_uniform(rng)
            end do
        end do
        model%class_probability = model%class_probability / sum(model%class_probability)

        do iter = 1, niter
            call lca_joint(patterns, model%class_probability, model%item_probability, joint)
            pattern_probability = sum(joint, dim=1)
            do i = 1, npattern
                if (pattern_probability(i) <= 0.0_dp) error stop "lca_fit_counts: zero pattern probability"
                posterior(:, i) = joint(:, i) / pattern_probability(i)
            end do
            do j = 1, k
                weighted(j, :) = posterior(j, :) * real(counts, dp)
                class_size(j) = sum(weighted(j, :))
                model%class_probability(j) = class_size(j) / real(n, dp)
                do v = 1, nvar
                    total_class = sum(weighted(j, :) * real(patterns(:, v), dp))
                    if (class_size(j) > 0.0_dp) model%item_probability(j, v) = total_class / class_size(j)
                end do
            end do
        end do

        call lca_joint(patterns, model%class_probability, model%item_probability, joint)
        pattern_probability = sum(joint, dim=1)
        likelihood = 0.0_dp
        saturated = 0.0_dp
        model%chi_square = 0.0_dp
        do i = 1, npattern
            if (counts(i) > 0) then
                likelihood = likelihood + real(counts(i), dp) * log(max(pattern_probability(i), tiny(1.0_dp)))
                saturated = saturated + real(counts(i), dp) * log(real(counts(i), dp) / real(n, dp))
            end if
            expected = real(n, dp) * pattern_probability(i)
            if (expected > 0.0_dp) model%chi_square = model%chi_square + (real(counts(i), dp) - expected)**2 / expected
            model%pattern_class(i) = maxloc(joint(:, i), dim=1)
        end do
        model%log_likelihood = likelihood
        model%saturated_log_likelihood = saturated
        model%likelihood_ratio = 2.0_dp * (saturated - likelihood)
        model%n = n
        model%n_parameters = k * (nvar + 1) - 1
        model%bic = -2.0_dp * likelihood + log(real(n, dp)) * real(model%n_parameters, dp)
        model%saturated_bic = -2.0_dp * saturated + log(real(n, dp)) * real(npattern - 1, dp)
    end subroutine lca_fit_counts

    subroutine lca_fit_data(x, k, niter, rng, model, matching)
        integer, intent(in) :: x(:, :) !! Binary observation-by-item data matrix.
        integer, intent(in) :: k !! Number of latent classes.
        integer, intent(in) :: niter !! Number of EM iterations.
        type(rng_state), intent(inout) :: rng !! Mutable random generator used for model initialization.
        type(lca_model), intent(out) :: model !! Fitted latent-class model.
        integer, allocatable, intent(out), optional :: matching(:) !! Optional fitted latent-class assignment for each original
        !! observation.
        integer, allocatable :: counts(:)
        integer, allocatable :: pattern_index(:)
        integer :: i

        call countpattern(x, counts, pattern_index)
        call lca_fit_counts(counts, k, niter, rng, model)
        if (present(matching)) then
            allocate(matching(size(x, 1)))
            do i = 1, size(x, 1)
                matching(i) = model%pattern_class(pattern_index(i))
            end do
        end if
    end subroutine lca_fit_data

    subroutine lca_predict(model, x, class)
        type(lca_model), intent(in) :: model !! Fitted model whose pattern_class array supplies latent-class assignments.
        integer, intent(in) :: x(:, :) !! Binary observation-by-item matrix with the fitted model's item count.
        integer, allocatable, intent(out) :: class(:) !! One-based fitted latent-class assignment for every observation.
        integer, allocatable :: counts(:)
        integer, allocatable :: matching(:)
        integer :: i

        if (2**size(x, 2) /= size(model%pattern_class)) error stop "lca_predict: item count mismatch"
        call countpattern(x, counts, matching)
        allocate(class(size(x, 1)))
        do i = 1, size(x, 1)
            class(i) = model%pattern_class(matching(i))
        end do
    end subroutine lca_predict

    subroutine lca_goodness_of_fit(model, df, p_likelihood_ratio, p_chi_square)
        type(lca_model), intent(in) :: model !! Fitted latent-class model carrying likelihood-ratio and Pearson statistics.
        integer, intent(out) :: df !! Saturated-minus-fitted parameter degrees of freedom.
        real(dp), intent(out) :: p_likelihood_ratio !! Upper-tail chi-square approximation for the likelihood-ratio statistic.
        real(dp), intent(out) :: p_chi_square !! Upper-tail chi-square approximation for the Pearson statistic.
        integer :: nvar

        nvar = exact_log2(size(model%pattern_class))
        df = 2**nvar - 1 - model%n_parameters
        if (df <= 0) then
            p_likelihood_ratio = 0.0_dp
            p_chi_square = 0.0_dp
        else
            p_likelihood_ratio = chi_square_sf(model%likelihood_ratio, real(df, dp))
            p_chi_square = chi_square_sf(model%chi_square, real(df, dp))
        end if
    end subroutine lca_goodness_of_fit

    subroutine lca_bootstrap(model, nsamples, lca_iter, rng, result)
        type(lca_model), intent(in) :: model !! Fitted latent-class model defining the parametric bootstrap data-generating
        !! distribution.
        integer, intent(in) :: nsamples !! Number of parametric bootstrap samples; must be positive.
        integer, intent(in) :: lca_iter !! EM iterations used to refit each bootstrap sample; must be positive.
        type(rng_state), intent(inout) :: rng !! Mutable random generator used for class, item, and refit initialization draws.
        type(lca_bootstrap_result), intent(out) :: result !! Bootstrap fit statistics and e1071-style empirical/normal-tail
        !! summaries.
        integer, allocatable :: x(:, :)
        integer, allocatable :: counts(:)
        type(lca_model) :: fit
        integer :: sample_index
        integer :: i
        integer :: v
        integer :: c
        integer :: nclass
        integer :: nvar
        real(dp) :: u

        if (nsamples < 1 .or. lca_iter < 1) error stop "lca_bootstrap: invalid iteration count"
        nclass = size(model%class_probability)
        nvar = size(model%item_probability, 2)
        allocate(result%log_likelihood(nsamples), result%saturated_log_likelihood(nsamples))
        allocate(result%likelihood_ratio(nsamples), result%chi_square(nsamples))
        allocate(x(model%n, nvar))
        do sample_index = 1, nsamples
            do i = 1, model%n
                u = rng_uniform(rng)
                c = sample_class(model%class_probability, u)
                do v = 1, nvar
                    x(i, v) = merge(1, 0, rng_uniform(rng) < model%item_probability(c, v))
                end do
            end do
            call countpattern(x, counts)
            call lca_fit_counts(counts, nclass, lca_iter, rng, fit)
            result%log_likelihood(sample_index) = fit%log_likelihood
            result%saturated_log_likelihood(sample_index) = fit%saturated_log_likelihood
            result%likelihood_ratio(sample_index) = fit%likelihood_ratio
            result%chi_square(sample_index) = fit%chi_square
        end do
        call mean_sd(result%likelihood_ratio, result%likelihood_ratio_mean, result%likelihood_ratio_sd)
        call mean_sd(result%chi_square, result%chi_square_mean, result%chi_square_sd)
        result%z_likelihood_ratio = (model%likelihood_ratio - result%likelihood_ratio_mean) / result%likelihood_ratio_sd
        result%z_chi_square = (model%chi_square - result%chi_square_mean) / result%chi_square_sd
        result%p_z_likelihood_ratio = 1.0_dp - normal_cdf(result%z_likelihood_ratio)
        result%p_z_chi_square = 1.0_dp - normal_cdf(result%z_chi_square)
        result%p_likelihood_ratio = real(count(model%likelihood_ratio < result%likelihood_ratio), dp) / real(nsamples, dp)
        result%p_chi_square = real(count(model%chi_square < result%chi_square), dp) / real(nsamples, dp)
    end subroutine lca_bootstrap

    subroutine lca_joint(patterns, class_probability, item_probability, joint)
        integer, intent(in) :: patterns(:, :) !! Binary pattern matrix with one pattern per row.
        real(dp), intent(in) :: class_probability(:) !! Latent-class prior probabilities summing to one.
        real(dp), intent(in) :: item_probability(:, :) !! Bernoulli success probabilities by class and item.
        real(dp), intent(out) :: joint(:, :) !! Joint P(pattern,class) matrix by class and pattern.
        real(dp) :: logp
        real(dp) :: p
        integer :: c
        integer :: i
        integer :: v

        do c = 1, size(class_probability)
            do i = 1, size(patterns, 1)
                logp = log(max(class_probability(c), tiny(1.0_dp)))
                do v = 1, size(patterns, 2)
                    if (patterns(i, v) == 1) then
                        p = item_probability(c, v)
                    else
                        p = 1.0_dp - item_probability(c, v)
                    end if
                    logp = logp + log(max(p, tiny(1.0_dp)))
                end do
                joint(c, i) = exp(logp)
            end do
        end do
    end subroutine lca_joint

    pure function exact_log2(n) result(power)
        integer, intent(in) :: n !! Positive integer tested for exact power-of-two length.
        integer :: power
        integer :: value

        if (n < 1) then
            power = -1
            return
        end if
        value = n
        power = 0
        do while (value > 1 .and. modulo(value, 2) == 0)
            value = value / 2
            power = power + 1
        end do
        if (value /= 1) power = -1
    end function exact_log2

    pure function sample_class(probability, u) result(index)
        real(dp), intent(in) :: probability(:) !! Class probabilities summing to one.
        real(dp), intent(in) :: u !! Uniform variate in (0,1) used for inversion sampling.
        integer :: index
        real(dp) :: cumulative

        cumulative = 0.0_dp
        do index = 1, size(probability)
            cumulative = cumulative + probability(index)
            if (u <= cumulative) return
        end do
        index = size(probability)
    end function sample_class

    subroutine mean_sd(x, mean_value, sd_value)
        real(dp), intent(in) :: x(:) !! Numeric vector whose arithmetic mean and sample standard deviation are requested.
        real(dp), intent(out) :: mean_value !! Arithmetic mean of x.
        real(dp), intent(out) :: sd_value !! Sample standard deviation using denominator n-1, or zero for a singleton.

        mean_value = sum(x) / real(size(x), dp)
        if (size(x) > 1) then
            sd_value = sqrt(sum((x - mean_value)**2) / real(size(x) - 1, dp))
        else
            sd_value = 0.0_dp
        end if
    end subroutine mean_sd

end module e1071_lca
