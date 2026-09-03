! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 simulate.MCMCglmm response transformations; see NOTICE.md and upstream/.
module mcmcglmm_simulation
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state, rng_uniform, rng_exponential, rng_chisq, rng_normal, rng_poisson
    use mcmcglmm_math, only : logistic, normal_cdf
    use mcmcglmm_matrix, only : sample_mvn_covariance, sample_mvn_precision
    implicit none
    private

    public :: simulate_scalar_response
    public :: simulate_ordinal_response
    public :: simulate_threshold_response
    public :: simulate_two_part_response
    public :: simulate_multinomial_response
    public :: simulate_multi_term_gaussian_latent

contains

    pure subroutine bernoulli_sample(state, probability, value)
        type(rng_state), intent(inout) :: state !! Generator state consumed by one uniform draw.
        real(dp), intent(in) :: probability !! Bernoulli success probability in the closed interval [0,1].
        integer, intent(out) :: value !! Sampled indicator, either zero or one.
        real(dp) :: u

        call rng_uniform(state, u)
        if (u < probability) then
            value = 1
        else
            value = 0
        end if
    end subroutine bernoulli_sample

    pure subroutine binomial_sample(state, trials, probability, value)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Bernoulli component draws.
        integer, intent(in) :: trials !! Nonnegative number of Bernoulli trials.
        real(dp), intent(in) :: probability !! Per-trial success probability in [0,1].
        integer, intent(out) :: value !! Number of sampled successes from zero through trials.
        integer :: draw
        integer :: i

        value = 0
        do i = 1, trials
            call bernoulli_sample(state, probability, draw)
            value = value + draw
        end do
    end subroutine binomial_sample

    pure subroutine zero_truncated_poisson_sample(state, lambda, value, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the conditional inverse-CDF draw.
        real(dp), intent(in) :: lambda !! Positive Poisson mean before conditioning on a nonzero count.
        integer, intent(out) :: value !! Positive Poisson count conditional on being greater than zero.
        integer, intent(out) :: info !! Zero on success; nonzero when lambda is nonpositive.
        real(dp) :: cumulative
        real(dp) :: denominator
        real(dp) :: probability
        real(dp) :: p0
        real(dp) :: u

        info = 0
        if (lambda <= 0.0_dp) then
            value = 0
            info = 1
            return
        end if
        p0 = exp(-lambda)
        if (lambda < 1.0e-5_dp) then
            denominator = lambda * (1.0_dp - 0.5_dp * lambda + lambda * lambda / 6.0_dp)
        else
            denominator = 1.0_dp - p0
        end if
        probability = p0 * lambda / max(denominator, tiny(1.0_dp))
        call rng_uniform(state, u)
        cumulative = probability
        value = 1
        do while (u > cumulative)
            value = value + 1
            probability = probability * lambda / real(value, dp)
            cumulative = cumulative + probability
            if (value > max(100, int(lambda + 20.0_dp * sqrt(lambda + 1.0_dp)))) exit
        end do
    end subroutine zero_truncated_poisson_sample

    pure subroutine zero_truncated_binomial_sample(state, trials, probability, value, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the conditional inverse-CDF draw.
        integer, intent(in) :: trials !! Positive number of Bernoulli trials.
        real(dp), intent(in) :: probability !! Per-trial success probability in [0,1].
        integer, intent(out) :: value !! Positive success count conditional on at least one success.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid trials or probability.
        real(dp) :: cumulative
        real(dp) :: denominator
        real(dp) :: mass
        real(dp) :: p0
        real(dp) :: ratio
        real(dp) :: u

        info = 0
        if (trials < 1 .or. probability < 0.0_dp .or. probability > 1.0_dp) then
            value = 0
            info = 1
            return
        end if
        if (probability <= epsilon(1.0_dp)) then
            value = 1
            return
        end if
        if (probability >= 1.0_dp - epsilon(1.0_dp)) then
            value = trials
            return
        end if
        p0 = (1.0_dp - probability) ** trials
        denominator = 1.0_dp - p0
        ratio = probability / (1.0_dp - probability)
        mass = p0 * real(trials, dp) * ratio / denominator
        call rng_uniform(state, u)
        cumulative = mass
        value = 1
        do while (u > cumulative .and. value < trials)
            mass = mass * real(trials - value, dp) / real(value + 1, dp) * ratio
            value = value + 1
            cumulative = cumulative + mass
        end do
    end subroutine zero_truncated_binomial_sample

    pure subroutine geometric_sample(state, probability, value, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by one uniform draw.
        real(dp), intent(in) :: probability !! Success probability of the geometric law on {0,1,...}.
        integer, intent(out) :: value !! Number of failures before the first success.
        integer, intent(out) :: info !! Zero on success; nonzero when probability is outside (0,1].
        real(dp) :: u

        info = 0
        if (probability <= 0.0_dp .or. probability > 1.0_dp) then
            value = 0
            info = 1
            return
        end if
        if (probability >= 1.0_dp - epsilon(1.0_dp)) then
            value = 0
            return
        end if
        call rng_uniform(state, u)
        value = int(log(max(1.0_dp - u, tiny(1.0_dp))) / log(1.0_dp - probability))
    end subroutine geometric_sample

    pure subroutine simulate_scalar_response(state, family, link_value, additional, additional2, response, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the requested response transformation.
        integer, intent(in) :: family !! Native scalar family code using MCMCglmm's response-simulation convention.
        real(dp), intent(in) :: link_value !! Latent/link value after Gaussian residual simulation where applicable.
        real(dp), intent(in) :: additional !! Family auxiliary value such as trials or Student-t scale.
        real(dp), intent(in) :: additional2 !! Second auxiliary value, used as Student-t degrees of freedom.
        real(dp), intent(out) :: response !! Simulated response on the observed-data scale.
        integer, intent(out) :: info !! Zero on success; nonzero for unsupported family or invalid auxiliaries.
        integer :: count
        integer :: indicator
        integer :: trials
        real(dp) :: chi_square
        real(dp) :: draw
        real(dp) :: rate
        real(dp) :: probability
        real(dp) :: z

        info = 0
        select case (family)
        case (1, 6)
            response = link_value
        case (2, 7)
            rate = exp(link_value)
            call rng_poisson(state, rate, count)
            response = real(count, dp)
        case (3)
            trials = nint(additional)
            if (trials < 1) then
                response = 0.0_dp
                info = 1
                return
            end if
            call binomial_sample(state, trials, logistic(link_value), count)
            response = real(count, dp)
        case (4, 5, 8, 9)
            rate = exp(link_value)
            call rng_exponential(state, 1.0_dp / rate, draw)
            response = draw
        case (16)
            call zero_truncated_poisson_sample(state, exp(link_value), count, info)
            response = real(count, dp)
        case (17)
            call geometric_sample(state, logistic(link_value), count, info)
            response = real(count, dp)
        case (22)
            trials = nint(additional)
            if (trials < 1) then
                response = 0.0_dp
                info = 1
                return
            end if
            probability = 1.0_dp - (1.0_dp - logistic(link_value)) ** trials
            call bernoulli_sample(state, probability, indicator)
            response = real(indicator, dp)
        case (23)
            if (additional <= 0.0_dp .or. additional2 <= 0.0_dp) then
                response = 0.0_dp
                info = 1
                return
            end if
            call rng_normal(state, z)
            call rng_chisq(state, additional2, chi_square)
            response = additional * (z + link_value / additional) / sqrt(chi_square / additional2)
        case (24)
            if (additional <= 0.0_dp .or. additional2 <= 0.0_dp) then
                response = 0.0_dp
                info = 1
                return
            end if
            call rng_normal(state, z)
            call rng_chisq(state, additional2, chi_square)
            response = link_value + additional * z / sqrt(chi_square / additional2)
        case default
            response = 0.0_dp
            info = 2
        end select
    end subroutine simulate_scalar_response

    pure subroutine simulate_ordinal_response(state, link_value, cutpoints, category, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by one categorical probability draw.
        real(dp), intent(in) :: link_value !! Family-14 latent ordered-probit link before the unit observation layer.
        real(dp), intent(in) :: cutpoints(:) !! Increasing category boundaries including effectively infinite endpoints.
        integer, intent(out) :: category !! One-based simulated ordinal category.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid cutpoints.
        real(dp) :: cumulative
        real(dp) :: probability
        real(dp) :: total
        real(dp) :: u
        integer :: j

        info = 0
        if (size(cutpoints) < 3) then
            category = 0
            info = 1
            return
        end if
        do j = 2, size(cutpoints)
            if (cutpoints(j) <= cutpoints(j - 1)) then
                category = 0
                info = 1
                return
            end if
        end do
        total = 0.0_dp
        do j = 1, size(cutpoints) - 1
            probability = normal_cdf(cutpoints(j + 1) - link_value) - &
                normal_cdf(cutpoints(j) - link_value)
            total = total + max(probability, 0.0_dp)
        end do
        call rng_uniform(state, u)
        u = u * max(total, tiny(1.0_dp))
        cumulative = 0.0_dp
        category = size(cutpoints) - 1
        do j = 1, size(cutpoints) - 1
            probability = normal_cdf(cutpoints(j + 1) - link_value) - &
                normal_cdf(cutpoints(j) - link_value)
            cumulative = cumulative + max(probability, 0.0_dp)
            if (u <= cumulative) then
                category = j
                exit
            end if
        end do
    end subroutine simulate_ordinal_response

    pure subroutine simulate_threshold_response(link_value, cutpoints, category, info)
        real(dp), intent(in) :: link_value !! Family-20 latent liability after residual simulation.
        real(dp), intent(in) :: cutpoints(:) !! Increasing category boundaries including effectively infinite endpoints.
        integer, intent(out) :: category !! One-based threshold category containing the liability.
        integer, intent(out) :: info !! Zero on success; nonzero when the cutpoint vector is invalid.
        integer :: j

        info = 0
        if (size(cutpoints) < 3) then
            category = 0
            info = 1
            return
        end if
        do j = 2, size(cutpoints)
            if (cutpoints(j) <= cutpoints(j - 1)) then
                category = 0
                info = 1
                return
            end if
        end do
        category = size(cutpoints) - 1
        do j = 1, size(cutpoints) - 1
            if (link_value >= cutpoints(j) .and. link_value < cutpoints(j + 1)) then
                category = j
                return
            end if
        end do
    end subroutine simulate_threshold_response

    pure subroutine simulate_two_part_response(state, family, eta_main, eta_zero, trials, response, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by zero-process and positive-response draws.
        integer, intent(in) :: family !! Two-process family code 11, 15, 18, 19, or 25.
        real(dp), intent(in) :: eta_main !! Main count/binomial link value.
        real(dp), intent(in) :: eta_zero !! Zero-process link value using MCMCglmm's native sign convention.
        integer, intent(in) :: trials !! Binomial trials for families 19/25; ignored for Poisson families.
        real(dp), intent(out) :: response !! Simulated nonnegative count or binomial response.
        integer, intent(out) :: info !! Zero on success; nonzero for unsupported family or invalid trials.
        integer :: active
        integer :: count
        real(dp) :: active_probability

        info = 0
        select case (family)
        case (11)
            active_probability = 1.0_dp - logistic(eta_zero)
            call bernoulli_sample(state, active_probability, active)
            if (active == 0) then
                response = 0.0_dp
            else
                call rng_poisson(state, exp(eta_main), count)
                response = real(count, dp)
            end if
        case (15)
            active_probability = 1.0_dp - logistic(eta_zero)
            call bernoulli_sample(state, active_probability, active)
            if (active == 0) then
                response = 0.0_dp
            else
                call zero_truncated_poisson_sample(state, exp(eta_main), count, info)
                response = real(count, dp)
            end if
        case (18)
            active_probability = 1.0_dp - exp(-exp(eta_zero))
            call bernoulli_sample(state, active_probability, active)
            if (active == 0) then
                response = 0.0_dp
            else
                call zero_truncated_poisson_sample(state, exp(eta_main), count, info)
                response = real(count, dp)
            end if
        case (19)
            if (trials < 1) then
                response = 0.0_dp
                info = 1
                return
            end if
            active_probability = 1.0_dp - logistic(eta_zero)
            call bernoulli_sample(state, active_probability, active)
            if (active == 0) then
                response = 0.0_dp
            else
                call binomial_sample(state, trials, logistic(eta_main), count)
                response = real(count, dp)
            end if
        case (25)
            if (trials < 1) then
                response = 0.0_dp
                info = 1
                return
            end if
            active_probability = 1.0_dp - logistic(eta_zero)
            call bernoulli_sample(state, active_probability, active)
            if (active == 0) then
                response = 0.0_dp
            else
                call zero_truncated_binomial_sample(state, trials, logistic(eta_main), count, info)
                response = real(count, dp)
            end if
        case default
            response = 0.0_dp
            info = 2
        end select
    end subroutine simulate_two_part_response

    pure subroutine simulate_multinomial_response(state, family, eta, trials, response, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by multinomial or Bernoulli category draws.
        integer, intent(in) :: family !! Grouped family code 3, 26, or 27.
        real(dp), intent(in) :: eta(:) !! Non-reference logits for 3/27, or Bernoulli logits for family 26.
        integer, intent(in) :: trials !! Multinomial size for families 3/27; ignored for family 26.
        integer, allocatable, intent(out) :: response(:) !! Allocated category counts or binary indicators.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid family or trial count.
        integer :: category
        integer :: draw
        integer :: j
        real(dp), allocatable :: probability(:)
        real(dp) :: cumulative
        real(dp) :: denominator
        real(dp) :: u

        info = 0
        select case (family)
        case (3, 27)
            if (trials < 1 .or. size(eta) < 1) then
                allocate(response(0))
                info = 1
                return
            end if
            allocate(response(size(eta) + 1), probability(size(eta) + 1))
            probability(1) = 1.0_dp
            probability(2:) = exp(eta - max(0.0_dp, maxval(eta)))
            probability(1) = exp(-max(0.0_dp, maxval(eta)))
            denominator = sum(probability)
            probability = probability / denominator
            response = 0
            do draw = 1, trials
                call rng_uniform(state, u)
                cumulative = 0.0_dp
                category = size(probability)
                do j = 1, size(probability)
                    cumulative = cumulative + probability(j)
                    if (u <= cumulative) then
                        category = j
                        exit
                    end if
                end do
                response(category) = response(category) + 1
            end do
        case (26)
            if (size(eta) < 1) then
                allocate(response(0))
                info = 1
                return
            end if
            allocate(response(size(eta)))
            do
                do j = 1, size(eta)
                    call bernoulli_sample(state, logistic(eta(j)), response(j))
                end do
                if (sum(response) > 0) exit
            end do
        case default
            allocate(response(0))
            info = 2
        end select
    end subroutine simulate_multinomial_response


    pure subroutine simulate_multi_term_gaussian_latent(state, x, z, random_term, a_inverse, beta, g_matrix, &
                                                        r_matrix, latent, random_effects, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by random-effect and residual MVN draws.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design matrix.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random-effect levels.
        real(dp), intent(in) :: beta(:, :) !! p by traits fixed-effect coefficients for the selected posterior draw.
        real(dp), intent(in) :: g_matrix(:, :, :) !! traits by traits by n_term random-effect covariance blocks.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance for the selected posterior draw.
        real(dp), allocatable, intent(out) :: latent(:, :) !! Allocated n by traits Gaussian latent/predictive response.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated q by traits newly sampled random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shapes or failed SPD sampling.
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: g_inverse(:, :)
        real(dp), allocatable :: mean_zero(:)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: residual_draw(:)
        real(dp), allocatable :: term_a_inverse(:, :)
        integer, allocatable :: term_indices(:)
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: n
        integer :: nterm
        integer :: p
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: term
        integer :: traits

        info = 0
        n = size(x, 1)
        p = size(x, 2)
        q = size(z, 2)
        traits = size(beta, 2)
        nterm = size(g_matrix, 3)
        if (n < 1 .or. p < 1 .or. q < 1 .or. traits < 1 .or. nterm < 1 .or. size(z, 1) /= n .or. &
            size(random_term) /= q .or. size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(beta, 1) /= p .or. size(g_matrix, 1) /= traits .or. size(g_matrix, 2) /= traits .or. &
            size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits .or. any(random_term < 1) .or. &
            any(random_term > nterm)) then
            allocate(latent(0, 0), random_effects(0, 0))
            info = 1
            return
        end if
        allocate(random_effects(q, traits))
        random_effects = 0.0_dp
        do term = 1, nterm
            qterm = count(random_term == term)
            if (qterm < 1) then
                allocate(latent(0, 0))
                info = 2
                return
            end if
            allocate(term_indices(qterm))
            i = 0
            do random_index = 1, q
                if (random_term(random_index) /= term) cycle
                i = i + 1
                term_indices(i) = random_index
            end do
            term_a_inverse = a_inverse(term_indices, term_indices)
            call inverse_matrix(g_matrix(:, :, term), g_inverse, info)
            if (info /= 0) then
                allocate(latent(0, 0))
                return
            end if
            allocate(precision(qterm * traits, qterm * traits))
            do b = 1, traits
                do a = 1, traits
                    do j = 1, qterm
                        do i = 1, qterm
                            precision((a - 1) * qterm + i, (b - 1) * qterm + j) = &
                                g_inverse(a, b) * term_a_inverse(i, j)
                        end do
                    end do
                end do
            end do
            allocate(mean_zero(qterm * traits))
            mean_zero = 0.0_dp
            call sample_mvn_precision(state, mean_zero, precision, packed, conditional_mean, info)
            if (info /= 0) then
                allocate(latent(0, 0))
                return
            end if
            do a = 1, traits
                random_effects(term_indices, a) = packed((a - 1) * qterm + 1:a * qterm)
            end do
            deallocate(term_indices, precision, mean_zero)
        end do

        latent = matmul(x, beta) + matmul(z, random_effects)
        allocate(mean_zero(traits))
        mean_zero = 0.0_dp
        do i = 1, n
            call sample_mvn_covariance(state, mean_zero, r_matrix, residual_draw, info)
            if (info /= 0) return
            latent(i, :) = latent(i, :) + residual_draw
        end do
    end subroutine simulate_multi_term_gaussian_latent

end module mcmcglmm_simulation
