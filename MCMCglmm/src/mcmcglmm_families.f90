! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_families
    use r_kinds, only : dp
    use mcmcglmm_math, only : logistic, log1pexp, normal_cdf, normal_logpdf, poisson_logpmf
    use mcmcglmm_distributions, only : pkk_probability
    implicit none
    private

    real(dp), parameter :: log_two_pi = 1.8378770664093454835606594728112353_dp

    public :: binomial_logpmf
    public :: poisson_cdf
    public :: student_t_logpdf
    public :: noncentral_t_logpdf
    public :: ordinal_probit_loglik
    public :: scalar_family_loglik
    public :: two_part_family_loglik
    public :: multinomial_log_kernel
    public :: ztmb_log_kernel
    public :: ztmultinomial_log_kernel

contains

    pure elemental real(dp) function log1p_safe(x) result(value)
        real(dp), intent(in) :: x !! Real argument greater than minus one for stable evaluation of log(1+x).
        real(dp) :: x2

        if (abs(x) < 1.0e-5_dp) then
            x2 = x * x
            value = x - 0.5_dp * x2 + x2 * x / 3.0_dp - 0.25_dp * x2 * x2 + 0.2_dp * x2 * x2 * x
        else
            value = log(1.0_dp + x)
        end if
    end function log1p_safe

    pure elemental real(dp) function expm1_safe(x) result(value)
        real(dp), intent(in) :: x !! Real argument for stable evaluation of exp(x)-1 near zero.
        real(dp) :: x2

        if (abs(x) < 1.0e-5_dp) then
            x2 = x * x
            value = x + 0.5_dp * x2 + x2 * x / 6.0_dp + x2 * x2 / 24.0_dp + x2 * x2 * x / 120.0_dp
        else
            value = exp(x) - 1.0_dp
        end if
    end function expm1_safe

    pure elemental real(dp) function binomial_logpmf(successes, trials, probability) result(value)
        integer, intent(in) :: successes !! Number of observed successes, required to lie between zero and trials.
        integer, intent(in) :: trials !! Number of Bernoulli trials, required to be nonnegative.
        real(dp), intent(in) :: probability !! Success probability in the closed interval from zero to one.

        if (trials < 0 .or. successes < 0 .or. successes > trials .or. &
            probability < 0.0_dp .or. probability > 1.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        if (probability <= 0.0_dp) then
            if (successes == 0) then
                value = 0.0_dp
            else
                value = -huge(1.0_dp)
            end if
            return
        end if
        if (probability >= 1.0_dp) then
            if (successes == trials) then
                value = 0.0_dp
            else
                value = -huge(1.0_dp)
            end if
            return
        end if
        value = log_gamma(real(trials + 1, dp)) - log_gamma(real(successes + 1, dp)) - &
            log_gamma(real(trials - successes + 1, dp)) + real(successes, dp) * log(probability) + &
            real(trials - successes, dp) * log1p_safe(-probability)
    end function binomial_logpmf

    pure real(dp) function regularized_gamma_q(shape, x) result(value)
        real(dp), intent(in) :: shape !! Positive gamma shape parameter.
        real(dp), intent(in) :: x !! Nonnegative gamma argument.
        integer, parameter :: max_iterations = 400
        real(dp), parameter :: tolerance = 1.0e-14_dp
        real(dp) :: a0
        real(dp) :: an
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: delta
        real(dp) :: factor
        real(dp) :: sum_value
        integer :: iteration

        if (shape <= 0.0_dp .or. x < 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (x <= 0.0_dp) then
            value = 1.0_dp
            return
        end if

        factor = exp(-x + shape * log(x) - log_gamma(shape))
        if (x < shape + 1.0_dp) then
            a0 = shape
            sum_value = 1.0_dp / shape
            an = sum_value
            do iteration = 1, max_iterations
                a0 = a0 + 1.0_dp
                an = an * x / a0
                sum_value = sum_value + an
                if (abs(an) <= abs(sum_value) * tolerance) exit
            end do
            value = max(0.0_dp, min(1.0_dp, 1.0_dp - factor * sum_value))
            return
        end if

        b = x + 1.0_dp - shape
        c = 1.0_dp / tiny(1.0_dp)
        d = 1.0_dp / b
        value = d
        do iteration = 1, max_iterations
            an = -real(iteration, dp) * (real(iteration, dp) - shape)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + an / c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp / d
            delta = d * c
            value = value * delta
            if (abs(delta - 1.0_dp) <= tolerance) exit
        end do
        value = max(0.0_dp, min(1.0_dp, factor * value))
    end function regularized_gamma_q

    pure elemental real(dp) function poisson_cdf(count, lambda) result(value)
        integer, intent(in) :: count !! Integer Poisson upper endpoint; negative values have cumulative probability zero.
        real(dp), intent(in) :: lambda !! Positive Poisson mean.

        if (count < 0) then
            value = 0.0_dp
        else if (lambda < 0.0_dp) then
            value = 0.0_dp
        else if (lambda <= tiny(1.0_dp)) then
            value = 1.0_dp
        else
            value = regularized_gamma_q(real(count + 1, dp), lambda)
        end if
    end function poisson_cdf

    pure elemental real(dp) function student_t_logpdf(x, degrees_freedom) result(value)
        real(dp), intent(in) :: x !! Standardized Student-t variate.
        real(dp), intent(in) :: degrees_freedom !! Positive Student-t degrees of freedom.

        if (degrees_freedom <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        value = log_gamma(0.5_dp * (degrees_freedom + 1.0_dp)) - &
            log_gamma(0.5_dp * degrees_freedom) - 0.5_dp * log(degrees_freedom) - &
            0.5_dp * log_two_pi + 0.5_dp * log(2.0_dp) - &
            0.5_dp * (degrees_freedom + 1.0_dp) * log1p_safe(x * x / degrees_freedom)
    end function student_t_logpdf


    pure real(dp) function noncentral_t_logpdf(x, degrees_freedom, noncentrality) result(value)
        real(dp), intent(in) :: x !! Noncentral Student-t variate.
        real(dp), intent(in) :: degrees_freedom !! Positive degrees of freedom.
        real(dp), intent(in) :: noncentrality !! Noncentrality parameter on the standard-normal numerator scale.
        integer, parameter :: n_steps = 512
        real(dp) :: a
        real(dp) :: constant_log
        real(dp) :: h
        real(dp) :: integrand
        real(dp) :: log_integrand
        real(dp) :: maximum_log
        real(dp) :: mode
        real(dp) :: s
        real(dp) :: smax
        real(dp) :: sum_value
        integer :: i
        integer :: weight

        if (degrees_freedom <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        a = degrees_freedom + x * x
        mode = (x * noncentrality + &
            sqrt(x * x * noncentrality * noncentrality + 4.0_dp * a * degrees_freedom)) / (2.0_dp * a)
        smax = max(1.0_dp, mode + 12.0_dp / sqrt(a) + 4.0_dp / sqrt(degrees_freedom))
        constant_log = log(2.0_dp * degrees_freedom) + &
            (0.5_dp * degrees_freedom - 1.0_dp) * log(degrees_freedom) - &
            0.5_dp * degrees_freedom * log(2.0_dp) - log_gamma(0.5_dp * degrees_freedom) - &
            0.5_dp * log_two_pi
        maximum_log = constant_log + degrees_freedom * log(max(mode, tiny(1.0_dp))) - &
            0.5_dp * (x * mode - noncentrality) ** 2 - 0.5_dp * degrees_freedom * mode * mode
        h = smax / real(n_steps, dp)
        sum_value = 0.0_dp
        do i = 0, n_steps
            if (i == 0) then
                integrand = 0.0_dp
            else
                s = real(i, dp) * h
                log_integrand = constant_log + degrees_freedom * log(s) - &
                    0.5_dp * (x * s - noncentrality) ** 2 - 0.5_dp * degrees_freedom * s * s
                integrand = exp(log_integrand - maximum_log)
            end if
            if (i == 0 .or. i == n_steps) then
                weight = 1
            else if (modulo(i, 2) == 0) then
                weight = 2
            else
                weight = 4
            end if
            sum_value = sum_value + real(weight, dp) * integrand
        end do
        if (sum_value <= 0.0_dp) then
            value = -huge(1.0_dp)
        else
            value = maximum_log + log(h * sum_value / 3.0_dp)
        end if
    end function noncentral_t_logpdf

    pure real(dp) function ordinal_probit_loglik(category, eta, cutpoints) result(value)
        integer, intent(in) :: category !! One-based observed ordinal category between one and size(cutpoints)-1.
        real(dp), intent(in) :: eta !! Latent ordered-probit link value before the unit-variance observation layer.
        real(dp), intent(in) :: cutpoints(:) !! Increasing category boundaries, normally with effectively infinite endpoints.
        real(dp) :: probability

        if (size(cutpoints) < 3 .or. category < 1 .or. category >= size(cutpoints)) then
            value = -huge(1.0_dp)
            return
        end if
        if (cutpoints(category + 1) <= cutpoints(category)) then
            value = -huge(1.0_dp)
            return
        end if
        if (cutpoints(category) < -1.0e32_dp) then
            probability = normal_cdf(cutpoints(category + 1) - eta)
        else if (cutpoints(category + 1) > 1.0e32_dp) then
            probability = normal_cdf(eta - cutpoints(category))
        else
            probability = normal_cdf(cutpoints(category + 1) - eta) - &
                normal_cdf(cutpoints(category) - eta)
        end if
        value = log(max(probability, tiny(1.0_dp)))
    end function ordinal_probit_loglik

    pure elemental real(dp) function log_positive_probability(log_complement) result(value)
        real(dp), intent(in) :: log_complement !! Log probability of the complementary zero event, required not to exceed zero.

        if (log_complement >= 0.0_dp) then
            value = -huge(1.0_dp)
        else
            value = log(-expm1_safe(log_complement))
        end if
    end function log_positive_probability

    pure elemental real(dp) function scalar_family_loglik(family, y, eta, additional, additional2) result(value)
        integer, intent(in) :: family !! MCMCglmm native family code for an active scalar-response likelihood.
        real(dp), intent(in) :: y !! Native response value after the same preprocessing used by MCMCglmm.
        real(dp), intent(in) :: eta !! Current latent/link value used by the native family likelihood.
        real(dp), intent(in) :: additional !! Family-specific auxiliary value such as an upper endpoint, trials, or scale.
        real(dp), intent(in) :: additional2 !! Second family-specific auxiliary value; used for Student-t degrees of freedom.
        real(dp) :: cdf_lower
        real(dp) :: cdf_upper
        real(dp) :: lambda
        real(dp) :: probability
        real(dp) :: scale
        integer :: trials

        select case (family)
        case (1)
            value = normal_logpdf(y, eta, sqrt(max(additional, tiny(1.0_dp))))
        case (2)
            value = poisson_logpmf(nint(y), exp(eta))
        case (3)
            trials = max(1, nint(additional))
            probability = logistic(eta)
            value = binomial_logpmf(nint(y), trials, probability)
        case (4, 5)
            if (y < 0.0_dp) then
                value = -huge(1.0_dp)
            else
                lambda = exp(eta)
                value = log(lambda) - lambda * y
            end if
        case (6)
            if (eta >= y .and. eta <= additional) then
                value = 0.0_dp
            else
                value = -huge(1.0_dp)
            end if
        case (7)
            lambda = exp(eta)
            cdf_lower = poisson_cdf(nint(y), lambda)
            cdf_upper = poisson_cdf(nint(additional), lambda)
            value = log(max(cdf_upper - cdf_lower, tiny(1.0_dp)))
        case (8, 9)
            if (additional < y .or. y < 0.0_dp) then
                value = -huge(1.0_dp)
            else
                lambda = exp(eta)
                cdf_lower = 1.0_dp - exp(-lambda * y)
                cdf_upper = 1.0_dp - exp(-lambda * additional)
                value = log(max(cdf_upper - cdf_lower, tiny(1.0_dp)))
            end if
        case (16)
            lambda = exp(eta)
            value = poisson_logpmf(nint(y), lambda) - log_positive_probability(-lambda)
        case (17)
            if (nint(y) < 0) then
                value = -huge(1.0_dp)
            else
                value = eta - log1pexp(eta) - real(nint(y), dp) * (eta + log1pexp(-eta))
            end if
        case (22)
            trials = nint(additional)
            probability = logistic(eta)
            if (trials < 1 .or. nint(y) < 0 .or. nint(y) > 1) then
                value = -huge(1.0_dp)
            else if (nint(y) == 0) then
                value = real(trials, dp) * log1p_safe(-probability)
            else
                value = log_positive_probability(real(trials, dp) * log1p_safe(-probability))
            end if
        case (23)
            scale = additional
            if (scale <= 0.0_dp .or. additional2 <= 0.0_dp) then
                value = -huge(1.0_dp)
            else
                value = noncentral_t_logpdf(y / scale, additional2, eta / scale) - log(scale)
            end if
        case (24)
            scale = additional
            if (scale <= 0.0_dp .or. additional2 <= 0.0_dp) then
                value = -huge(1.0_dp)
            else
                value = student_t_logpdf((y - eta) / scale, additional2) - log(scale)
            end if
        case default
            value = -huge(1.0_dp)
        end select
    end function scalar_family_loglik

    pure elemental real(dp) function two_part_family_loglik(family, y, eta_main, eta_zero, trials) result(value)
        integer, intent(in) :: family !! Native two-process family code: 11, 15, 18, 19, or 25.
        real(dp), intent(in) :: y !! Original count or success count before MCMCglmm expands the response into two traits.
        real(dp), intent(in) :: eta_main !! Main count/binomial latent link value.
        real(dp), intent(in) :: eta_zero !! Zero-process latent link value.
        integer, intent(in) :: trials !! Binomial trials for families 19 and 25; ignored for Poisson families.
        real(dp) :: log_main
        real(dp) :: log_positive
        real(dp) :: log_zero
        real(dp) :: lambda
        real(dp) :: p_main
        real(dp) :: p_zero
        real(dp) :: terms(2)

        select case (family)
        case (11)
            lambda = exp(eta_main)
            p_zero = logistic(eta_zero)
            log_main = poisson_logpmf(nint(y), lambda)
            if (nint(y) == 0) then
                terms = [log1p_safe(-p_zero) + log_main, log(p_zero)]
                value = maxval(terms) + log(sum(exp(terms - maxval(terms))))
            else
                value = log1p_safe(-p_zero) + log_main
            end if
        case (15)
            lambda = exp(eta_main)
            p_zero = logistic(eta_zero)
            if (nint(y) == 0) then
                value = log(p_zero)
            else
                log_main = poisson_logpmf(nint(y), lambda) - log_positive_probability(-lambda)
                value = log1p_safe(-p_zero) + log_main
            end if
        case (18)
            lambda = exp(eta_main)
            log_zero = -exp(eta_zero)
            if (nint(y) == 0) then
                value = log_zero
            else
                log_main = poisson_logpmf(nint(y), lambda) - log_positive_probability(-lambda)
                log_positive = log_positive_probability(log_zero)
                value = log_positive + log_main
            end if
        case (19)
            if (trials < 1) then
                value = -huge(1.0_dp)
                return
            end if
            p_main = logistic(eta_main)
            p_zero = logistic(eta_zero)
            log_main = binomial_logpmf(nint(y), trials, p_main)
            if (nint(y) == 0) then
                terms = [log1p_safe(-p_zero) + log_main, log(p_zero)]
                value = maxval(terms) + log(sum(exp(terms - maxval(terms))))
            else
                value = log1p_safe(-p_zero) + log_main
            end if
        case (25)
            if (trials < 1) then
                value = -huge(1.0_dp)
                return
            end if
            p_main = logistic(eta_main)
            p_zero = logistic(eta_zero)
            if (nint(y) == 0) then
                value = log(p_zero)
            else
                log_main = binomial_logpmf(nint(y), trials, p_main)
                log_main = log_main - log_positive_probability(real(trials, dp) * log1p_safe(-p_main))
                value = log1p_safe(-p_zero) + log_main
            end if
        case default
            value = -huge(1.0_dp)
        end select
    end function two_part_family_loglik

    pure real(dp) function multinomial_log_kernel(counts, eta) result(value)
        integer, intent(in) :: counts(:) !! Category counts including the reference category as the last element.
        real(dp), intent(in) :: eta(:) !! Logit values for all non-reference categories, of length size(counts)-1.
        real(dp) :: maximum_value
        real(dp) :: denominator
        integer :: category
        integer :: total

        if (size(counts) < 2 .or. size(eta) /= size(counts) - 1 .or. any(counts < 0)) then
            value = -huge(1.0_dp)
            return
        end if
        total = sum(counts)
        maximum_value = max(0.0_dp, maxval(eta))
        denominator = exp(-maximum_value) + sum(exp(eta - maximum_value))
        value = -real(total, dp) * (maximum_value + log(denominator))
        do category = 1, size(eta)
            value = value + real(counts(category), dp) * eta(category)
        end do
    end function multinomial_log_kernel

    pure real(dp) function ztmb_log_kernel(binary_response, eta) result(value)
        integer, intent(in) :: binary_response(:) !! Zero/one multiple-Bernoulli responses conditioned on at least one success.
        real(dp), intent(in) :: eta(:) !! Logits for the corresponding Bernoulli responses.
        real(dp) :: log_all_zero
        real(dp) :: probability
        integer :: category

        if (size(binary_response) < 2 .or. size(eta) /= size(binary_response) .or. &
            any(binary_response < 0) .or. any(binary_response > 1) .or. sum(binary_response) < 1) then
            value = -huge(1.0_dp)
            return
        end if
        value = 0.0_dp
        log_all_zero = 0.0_dp
        do category = 1, size(eta)
            probability = logistic(eta(category))
            if (binary_response(category) == 1) then
                value = value + log(probability)
            else
                value = value + log1p_safe(-probability)
            end if
            log_all_zero = log_all_zero + log1p_safe(-probability)
        end do
        value = value - log_positive_probability(log_all_zero)
    end function ztmb_log_kernel

    pure real(dp) function ztmultinomial_log_kernel(counts, eta) result(value)
        integer, intent(in) :: counts(:) !! Category counts including the reference category as the last element.
        real(dp), intent(in) :: eta(:) !! Logits for non-reference categories, matching the MCMCglmm reference-category convention.
        real(dp), allocatable :: probability(:)
        real(dp), allocatable :: weights(:)
        integer :: category
        integer :: n_present
        integer :: total
        real(dp) :: denominator

        if (size(counts) < 2 .or. size(eta) /= size(counts) - 1 .or. any(counts < 0)) then
            value = -huge(1.0_dp)
            return
        end if
        total = sum(counts)
        if (total < 1) then
            value = -huge(1.0_dp)
            return
        end if
        n_present = count(counts > 0)
        if (n_present < 1) then
            value = -huge(1.0_dp)
            return
        end if
        allocate(weights(n_present), probability(n_present))
        n_present = 0
        denominator = 0.0_dp
        value = 0.0_dp
        do category = 1, size(eta)
            if (counts(category) > 0) then
                n_present = n_present + 1
                weights(n_present) = exp(eta(category))
                denominator = denominator + weights(n_present)
                value = value + real(counts(category), dp) * eta(category)
            end if
        end do
        if (counts(size(counts)) > 0) then
            n_present = n_present + 1
            weights(n_present) = 1.0_dp
            denominator = denominator + 1.0_dp
        end if
        if (denominator <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        value = value - real(total, dp) * log(denominator)
        probability = weights / denominator
        value = value - log(max(pkk_probability(probability, real(total, dp)), tiny(1.0_dp)))
    end function ztmultinomial_log_kernel

end module mcmcglmm_families
