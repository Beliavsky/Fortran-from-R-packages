! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 sampler update blocks; see NOTICE.md and upstream/.
module mcmcglmm_engine_features
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state, rng_exponential, rng_normal, rng_uniform
    use mcmcglmm_math, only : logistic, log1pexp, normal_cdf, normal_quantile
    use mcmcglmm_matrix, only : mvn_log_density
    use mcmcglmm_distributions, only : truncated_normal_sample
    implicit none
    private

    public :: optimal_acceptance_ratio
    public :: adaptive_mh_observe
    public :: adaptive_mh_decay
    public :: adaptive_mh_finalize
    public :: binary_slice_liability_update
    public :: theta_scale_conditional
    public :: categorical_measurement_error_update
    public :: structural_matrix
    public :: structural_transform
    public :: structural_gaussian_loglik
    public :: structural_parameter_random_walk_update

contains

    pure elemental real(dp) function optimal_acceptance_ratio(dimension) result(value)
        integer, intent(in) :: dimension !! Positive dimension of the adaptive multivariate Metropolis proposal.

        if (dimension < 1) then
            value = 0.44_dp
        else
            value = -real(dimension, dp) / (1.0_dp - 2.75_dp * real(dimension, dp)) - 0.133_dp
        end if
    end function optimal_acceptance_ratio

    pure subroutine adaptive_mh_observe(sample, count, mean_value, proposal_covariance, info)
        real(dp), intent(in) :: sample(:) !! Newly retained liability vector used to update the empirical proposal moments.
        real(dp), intent(inout) :: count !! Number of vectors accumulated before this call; incremented by one.
        real(dp), intent(inout) :: mean_value(:) !! Running proposal mean updated with the upstream online recursion.
        real(dp), intent(inout) :: proposal_covariance(:, :) !! Running proposal covariance accumulator updated in place.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible dimensions or a negative count.
        real(dp), allocatable :: old_mean(:)
        integer :: i
        integer :: j
        integer :: n

        info = 0
        n = size(sample)
        if (n < 1 .or. size(mean_value) /= n .or. size(proposal_covariance, 1) /= n .or. &
            size(proposal_covariance, 2) /= n .or. count < 0.0_dp) then
            info = 1
            return
        end if

        count = count + 1.0_dp
        old_mean = mean_value
        proposal_covariance = proposal_covariance * ((count - 1.0_dp) / count)
        do j = 1, n
            do i = 1, n
                proposal_covariance(i, j) = proposal_covariance(i, j) + old_mean(i) * old_mean(j)
            end do
        end do
        mean_value = ((count - 1.0_dp) * old_mean + sample) / count
        do j = 1, n
            do i = 1, n
                proposal_covariance(i, j) = proposal_covariance(i, j) - &
                    ((count + 1.0_dp) / count) * mean_value(i) * mean_value(j)
                proposal_covariance(i, j) = proposal_covariance(i, j) + sample(i) * sample(j) / count
            end do
        end do
        do i = 1, n
            proposal_covariance(i, i) = proposal_covariance(i, i) + 0.001_dp / count
        end do
    end subroutine adaptive_mh_observe

    pure subroutine adaptive_mh_decay(accepted, accepted_weight, attempted_weight)
        logical, intent(in) :: accepted !! Whether the most recent Metropolis proposal was accepted.
        real(dp), intent(inout) :: accepted_weight !! Exponentially decayed weighted acceptance count.
        real(dp), intent(inout) :: attempted_weight !! Exponentially decayed weighted proposal count.

        accepted_weight = 0.9_dp * accepted_weight
        attempted_weight = 0.9_dp * attempted_weight + 1.0_dp
        if (accepted) accepted_weight = accepted_weight + 1.0_dp
    end subroutine adaptive_mh_decay

    pure subroutine adaptive_mh_finalize(proposal_covariance, mean_value, count, accepted_weight, attempted_weight, &
                                         scale, info)
        real(dp), intent(inout) :: proposal_covariance(:, :) !! Empirical proposal covariance rescaled in place.
        real(dp), intent(in) :: mean_value(:) !! Current running proposal mean.
        real(dp), intent(in) :: count !! Positive number of vectors accumulated in the proposal moments.
        real(dp), intent(inout) :: accepted_weight !! Weighted accepted proposals; reset to zero after scaling.
        real(dp), intent(inout) :: attempted_weight !! Weighted attempted proposals; reset to zero after scaling.
        real(dp), intent(out) :: scale !! MCMCglmm adaptive scale factor based on the current acceptance ratio.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible dimensions or nonpositive count.
        integer :: i
        integer :: j
        integer :: n
        real(dp) :: target

        info = 0
        n = size(mean_value)
        if (n < 1 .or. size(proposal_covariance, 1) /= n .or. size(proposal_covariance, 2) /= n .or. count <= 0.0_dp) then
            scale = 1.0_dp
            info = 1
            return
        end if
        do j = 1, n
            do i = 1, n
                proposal_covariance(i, j) = proposal_covariance(i, j) + mean_value(i) * mean_value(j) / count
            end do
        end do
        target = optimal_acceptance_ratio(n)
        if (attempted_weight > 0.0_dp) then
            scale = 2.0_dp ** (accepted_weight / attempted_weight - target)
        else
            scale = 1.0_dp
        end if
        proposal_covariance = scale * proposal_covariance
        accepted_weight = 0.0_dp
        attempted_weight = 0.0_dp
    end subroutine adaptive_mh_finalize

    pure subroutine binary_slice_liability_update(state, family, y, trials, mean_value, sd_value, current, &
                                                   limit, sample, log_likelihood, lower, upper, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the exponential slice and truncated-normal draws.
        integer, intent(in) :: family !! Native family code; supported slice branches are 3, 14, and 22.
        real(dp), intent(in) :: y !! Binary response; family 14 uses encoded categories 1 and 2.
        real(dp), intent(in) :: trials !! Binomial trial count used only for native family 22.
        real(dp), intent(in) :: mean_value !! Conditional Gaussian liability mean before observation truncation.
        real(dp), intent(in) :: sd_value !! Positive conditional Gaussian liability standard deviation.
        real(dp), intent(in) :: current !! Current scalar liability used to define the slice height.
        real(dp), intent(in) :: limit !! Positive finite truncation guard corresponding to upstream logitt/probitt.
        real(dp), intent(out) :: sample !! Updated liability draw from the exact upstream univariate slice interval.
        real(dp), intent(out) :: log_likelihood !! Observation log likelihood at the current liability before slicing.
        real(dp), intent(out) :: lower !! Lower endpoint of the generated slice interval.
        real(dp), intent(out) :: upper !! Upper endpoint of the generated slice interval.
        integer, intent(out) :: info !! Zero on success; nonzero for unsupported family or invalid scale/trial inputs.
        real(dp) :: exponential_draw
        real(dp) :: probability
        real(dp) :: slice_log_height
        real(dp) :: transformed

        info = 0
        if (sd_value <= 0.0_dp .or. limit <= 0.0_dp) then
            sample = current
            log_likelihood = -huge(1.0_dp)
            lower = -limit
            upper = limit
            info = 1
            return
        end if
        call rng_exponential(state, 1.0_dp, exponential_draw)

        select case (family)
        case (3)
            if (y > 0.5_dp) then
                log_likelihood = current - log1pexp(current)
                slice_log_height = log_likelihood - exponential_draw
                probability = min(max(exp(slice_log_height), tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
                lower = log(probability) - log(1.0_dp - probability)
                upper = limit
            else
                log_likelihood = -log1pexp(current)
                slice_log_height = log_likelihood - exponential_draw
                probability = min(max(exp(slice_log_height), tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
                lower = -limit
                upper = log(1.0_dp - probability) - log(probability)
            end if
        case (14)
            if (y > 1.5_dp) then
                log_likelihood = log(max(normal_cdf(current), tiny(1.0_dp)))
                slice_log_height = log_likelihood - exponential_draw
                probability = min(max(exp(slice_log_height), tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
                lower = max(normal_quantile(probability), -limit)
                upper = limit
            else
                log_likelihood = log(max(normal_cdf(-current), tiny(1.0_dp)))
                slice_log_height = log_likelihood - exponential_draw
                probability = min(max(exp(slice_log_height), tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
                lower = -limit
                upper = min(-normal_quantile(probability), limit)
            end if
        case (22)
            if (trials <= 0.0_dp) then
                sample = current
                log_likelihood = -huge(1.0_dp)
                lower = -limit
                upper = limit
                info = 2
                return
            end if
            if (y > 0.5_dp) then
                transformed = max(1.0_dp - (1.0_dp - logistic(current)) ** trials, tiny(1.0_dp))
                log_likelihood = log(transformed)
                slice_log_height = log_likelihood - exponential_draw
                lower = -limit
                transformed = exp(slice_log_height)
                upper = log(max(1.0_dp - (1.0_dp - logistic(transformed)) ** trials, tiny(1.0_dp)))
                upper = min(upper, limit)
            else
                log_likelihood = trials * log(max(1.0_dp - logistic(current), tiny(1.0_dp)))
                slice_log_height = log_likelihood - exponential_draw
                transformed = exp(slice_log_height)
                lower = trials * log(max(1.0_dp - logistic(transformed), tiny(1.0_dp)))
                lower = max(lower, -limit)
                upper = limit
            end if
        case default
            sample = current
            log_likelihood = -huge(1.0_dp)
            lower = -limit
            upper = limit
            info = 3
            return
        end select

        if (lower >= upper) then
            sample = current
            info = 4
            return
        end if
        call truncated_normal_sample(state, mean_value, sd_value, lower, upper, sample, info)
    end subroutine binary_slice_liability_update

    pure subroutine theta_scale_conditional(state, current_scale, deviation, scale_predictor, residual_variance, &
                                            prior_mean, prior_precision, sample, conditional_mean, &
                                            conditional_variance, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Gaussian theta-scale draw.
        real(dp), intent(in) :: current_scale !! Current theta-scale value already represented in the supplied predictor.
        real(dp), intent(in) :: deviation(:) !! Current liability/response minus the full predictor using current_scale.
        real(dp), intent(in) :: scale_predictor(:) !! Contribution Wscale times location that is multiplied by theta_scale.
        real(dp), intent(in) :: residual_variance !! Positive homogeneous residual variance for scaled observations.
        real(dp), intent(in) :: prior_mean !! Gaussian prior mean for theta_scale.
        real(dp), intent(in) :: prior_precision !! Nonnegative Gaussian prior precision, the inverse prior variance.
        real(dp), intent(out) :: sample !! Gaussian full-conditional draw for theta_scale.
        real(dp), intent(out) :: conditional_mean !! Full-conditional Gaussian mean.
        real(dp), intent(out) :: conditional_variance !! Full-conditional Gaussian variance.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible vectors or invalid variances.
        real(dp) :: information
        real(dp) :: z

        info = 0
        if (size(deviation) /= size(scale_predictor) .or. size(deviation) < 1 .or. residual_variance <= 0.0_dp .or. &
            prior_precision < 0.0_dp) then
            sample = current_scale
            conditional_mean = current_scale
            conditional_variance = 0.0_dp
            info = 1
            return
        end if
        information = dot_product(scale_predictor, scale_predictor)
        if (information + prior_precision <= 0.0_dp) then
            sample = current_scale
            conditional_mean = current_scale
            conditional_variance = 0.0_dp
            info = 2
            return
        end if
        conditional_mean = (dot_product(scale_predictor, deviation) + information * current_scale + &
            prior_mean * prior_precision) / (information + prior_precision)
        conditional_variance = residual_variance / (information + prior_precision)
        call rng_normal(state, z)
        sample = conditional_mean + sqrt(conditional_variance) * z
    end subroutine theta_scale_conditional

    pure subroutine categorical_measurement_error_update(state, liability, base_mean, category_effect, &
                                                         residual_covariance, prior_probability, group, category, &
                                                         posterior_probability, adjusted_mean, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed when latent measurement-error categories are drawn.
        real(dp), intent(in) :: liability(:, :) !! n by t current latent/observed response vectors used in the residual likelihood.
        real(dp), intent(in) :: base_mean(:, :) !! n by t predictor with the current measurement-error contribution removed.
        real(dp), intent(in) :: category_effect(:, :, :) !! n by t by c candidate predictor contributions for each category.
        real(dp), intent(in) :: residual_covariance(:, :) !! t by t residual covariance used by candidate Gaussian likelihoods.
        real(dp), intent(in) :: prior_probability(:, :) !! g by c prior category probabilities for each latent outcome group.
        integer, intent(in) :: group(:) !! Length-n one-based latent outcome group informed by each response row.
        integer, allocatable, intent(out) :: category(:) !! Allocated length-g sampled one-based category assignments.
        real(dp), allocatable, intent(out) :: posterior_probability(:, :) !! Allocated g by c normalized posterior probabilities.
        real(dp), allocatable, intent(out) :: adjusted_mean(:, :) !! Allocated predictor after inserting the sampled
            !! category effects.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shapes, probabilities, groups, or covariance failure.
        real(dp), allocatable :: log_weight(:)
        real(dp) :: log_density
        real(dp) :: max_log_weight
        real(dp) :: total
        real(dp) :: u
        integer :: c
        integer :: g
        integer :: i
        integer :: k
        integer :: n
        integer :: t

        info = 0
        n = size(liability, 1)
        t = size(liability, 2)
        g = size(prior_probability, 1)
        c = size(prior_probability, 2)
        if (n < 1 .or. t < 1 .or. g < 1 .or. c < 1 .or. size(base_mean, 1) /= n .or. &
            size(base_mean, 2) /= t .or. size(category_effect, 1) /= n .or. size(category_effect, 2) /= t .or. &
            size(category_effect, 3) /= c .or. size(residual_covariance, 1) /= t .or. &
            size(residual_covariance, 2) /= t .or. size(group) /= n .or. any(group < 1) .or. any(group > g) .or. &
            any(prior_probability < 0.0_dp) .or. any(sum(prior_probability, dim=2) <= 0.0_dp)) then
            allocate(category(0), posterior_probability(0, 0), adjusted_mean(0, 0))
            info = 1
            return
        end if

        allocate(category(g), posterior_probability(g, c), adjusted_mean(n, t), log_weight(c))
        posterior_probability = 0.0_dp
        do k = 1, g
            do i = 1, c
                if (prior_probability(k, i) > 0.0_dp) then
                    log_weight(i) = log(prior_probability(k, i))
                else
                    log_weight(i) = -huge(1.0_dp)
                end if
            end do
            do i = 1, n
                if (group(i) /= k) cycle
                do c = 1, size(prior_probability, 2)
                    if (log_weight(c) <= -0.5_dp * huge(1.0_dp)) cycle
                    call mvn_log_density(liability(i, :), base_mean(i, :) + category_effect(i, :, c), &
                        residual_covariance, log_density, info)
                    if (info /= 0) return
                    log_weight(c) = log_weight(c) + log_density
                end do
            end do
            max_log_weight = maxval(log_weight)
            posterior_probability(k, :) = exp(log_weight - max_log_weight)
            total = sum(posterior_probability(k, :))
            if (total <= 0.0_dp) then
                info = 2
                return
            end if
            posterior_probability(k, :) = posterior_probability(k, :) / total
            call rng_uniform(state, u)
            total = 0.0_dp
            category(k) = size(prior_probability, 2)
            do c = 1, size(prior_probability, 2)
                total = total + posterior_probability(k, c)
                if (u <= total) then
                    category(k) = c
                    exit
                end if
            end do
        end do
        adjusted_mean = base_mean
        do i = 1, n
            adjusted_mean(i, :) = adjusted_mean(i, :) + category_effect(i, :, category(group(i)))
        end do
    end subroutine categorical_measurement_error_update

    pure subroutine structural_matrix(basis, parameter, lambda_matrix, info)
        real(dp), intent(in) :: basis(:, :, :) !! t by t by p structural design basis with one matrix per structural parameter.
        real(dp), intent(in) :: parameter(:) !! Length-p structural coefficients multiplying the supplied basis matrices.
        real(dp), allocatable, intent(out) :: lambda_matrix(:, :) !! Allocated matrix I minus the weighted structural design.
        integer, intent(out) :: info !! Zero on success; nonzero for nonsquare basis matrices or incompatible parameter count.
        integer :: i
        integer :: n

        info = 0
        n = size(basis, 1)
        if (n < 1 .or. size(basis, 2) /= n .or. size(basis, 3) /= size(parameter)) then
            allocate(lambda_matrix(0, 0))
            info = 1
            return
        end if
        allocate(lambda_matrix(n, n))
        lambda_matrix = 0.0_dp
        do i = 1, n
            lambda_matrix(i, i) = 1.0_dp
        end do
        do i = 1, size(parameter)
            lambda_matrix = lambda_matrix - parameter(i) * basis(:, :, i)
        end do
    end subroutine structural_matrix

    pure subroutine dense_log_abs_det(matrix_value, log_abs_det, sign_value, info)
        real(dp), intent(in) :: matrix_value(:, :) !! Square dense matrix whose determinant sign and log magnitude are requested.
        real(dp), intent(out) :: log_abs_det !! Natural logarithm of the absolute determinant.
        integer, intent(out) :: sign_value !! Determinant sign as minus one or plus one; zero when singular.
        integer, intent(out) :: info !! Zero on success; nonzero for nonsquare or numerically singular input.
        real(dp), allocatable :: work(:, :)
        real(dp), allocatable :: row_tmp(:)
        real(dp) :: pivot_value
        real(dp) :: factor
        integer :: i
        integer :: j
        integer :: n
        integer :: pivot

        n = size(matrix_value, 1)
        if (n < 1 .or. size(matrix_value, 2) /= n) then
            log_abs_det = -huge(1.0_dp)
            sign_value = 0
            info = 1
            return
        end if
        work = matrix_value
        sign_value = 1
        log_abs_det = 0.0_dp
        info = 0
        allocate(row_tmp(n))
        do j = 1, n
            pivot = j - 1 + maxloc(abs(work(j:n, j)), dim=1)
            pivot_value = work(pivot, j)
            if (abs(pivot_value) <= tiny(1.0_dp)) then
                log_abs_det = -huge(1.0_dp)
                sign_value = 0
                info = 2
                return
            end if
            if (pivot /= j) then
                row_tmp = work(j, :)
                work(j, :) = work(pivot, :)
                work(pivot, :) = row_tmp
                sign_value = -sign_value
            end if
            if (work(j, j) < 0.0_dp) sign_value = -sign_value
            log_abs_det = log_abs_det + log(abs(work(j, j)))
            do i = j + 1, n
                factor = work(i, j) / work(j, j)
                work(i, j:n) = work(i, j:n) - factor * work(j, j:n)
            end do
        end do
    end subroutine dense_log_abs_det

    pure subroutine structural_transform(y, basis, parameter, transformed, log_jacobian, determinant_sign, info)
        real(dp), intent(in) :: y(:, :) !! n by t original response/liability matrix before structural transformation.
        real(dp), intent(in) :: basis(:, :, :) !! t by t by p structural design basis matrices.
        real(dp), intent(in) :: parameter(:) !! Length-p structural coefficients defining Lambda = I - sum(parameter*basis).
        real(dp), allocatable, intent(out) :: transformed(:, :) !! Allocated n by t matrix whose rows are Lambda times
            !! y-row columns.
        real(dp), intent(out) :: log_jacobian !! n times log absolute determinant of Lambda.
        integer, intent(out) :: determinant_sign !! Sign of determinant Lambda, used by upstream to reject sign-changing proposals.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shapes or a singular structural matrix.
        real(dp), allocatable :: lambda_matrix(:, :)
        real(dp) :: logdet
        integer :: sign_value

        call structural_matrix(basis, parameter, lambda_matrix, info)
        if (info /= 0 .or. size(y, 2) /= size(lambda_matrix, 1)) then
            allocate(transformed(0, 0))
            log_jacobian = -huge(1.0_dp)
            determinant_sign = 0
            if (info == 0) info = 1
            return
        end if
        call dense_log_abs_det(lambda_matrix, logdet, sign_value, info)
        if (info /= 0) then
            allocate(transformed(0, 0))
            log_jacobian = -huge(1.0_dp)
            determinant_sign = 0
            return
        end if
        transformed = transpose(matmul(lambda_matrix, transpose(y)))
        log_jacobian = real(size(y, 1), dp) * logdet
        determinant_sign = sign_value
    end subroutine structural_transform

    pure subroutine structural_gaussian_loglik(y, mean_value, residual_covariance, basis, parameter, &
                                               log_likelihood, determinant_sign, info)
        real(dp), intent(in) :: y(:, :) !! n by t original Gaussian responses or liabilities.
        real(dp), intent(in) :: mean_value(:, :) !! n by t mean on the transformed Lambda-y scale.
        real(dp), intent(in) :: residual_covariance(:, :) !! t by t Gaussian residual covariance on the transformed scale.
        real(dp), intent(in) :: basis(:, :, :) !! t by t by p structural basis matrices.
        real(dp), intent(in) :: parameter(:) !! Length-p structural coefficients.
        real(dp), intent(out) :: log_likelihood !! Gaussian structural-model log likelihood including the Jacobian.
        integer, intent(out) :: determinant_sign !! Sign of determinant Lambda.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, singular Lambda, or covariance failure.
        real(dp), allocatable :: transformed(:, :)
        real(dp) :: log_jacobian
        real(dp) :: row_log_density
        integer :: i

        if (size(mean_value, 1) /= size(y, 1) .or. size(mean_value, 2) /= size(y, 2) .or. &
            size(residual_covariance, 1) /= size(y, 2) .or. size(residual_covariance, 2) /= size(y, 2)) then
            log_likelihood = -huge(1.0_dp)
            determinant_sign = 0
            info = 1
            return
        end if
        call structural_transform(y, basis, parameter, transformed, log_jacobian, determinant_sign, info)
        if (info /= 0) then
            log_likelihood = -huge(1.0_dp)
            return
        end if
        log_likelihood = log_jacobian
        do i = 1, size(y, 1)
            call mvn_log_density(transformed(i, :), mean_value(i, :), residual_covariance, row_log_density, info)
            if (info /= 0) then
                log_likelihood = -huge(1.0_dp)
                return
            end if
            log_likelihood = log_likelihood + row_log_density
        end do
    end subroutine structural_gaussian_loglik

    pure subroutine structural_parameter_random_walk_update(state, y, mean_value, residual_covariance, basis, &
                                                            prior_mean, prior_precision, proposal_sd, parameter, &
                                                            accepted, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by Gaussian proposals and the MH uniform draw.
        real(dp), intent(in) :: y(:, :) !! n by t original Gaussian responses/liabilities.
        real(dp), intent(in) :: mean_value(:, :) !! n by t mean on the transformed Lambda-y scale.
        real(dp), intent(in) :: residual_covariance(:, :) !! t by t transformed-scale residual covariance.
        real(dp), intent(in) :: basis(:, :, :) !! t by t by p structural basis matrices.
        real(dp), intent(in) :: prior_mean(:) !! Length-p Gaussian prior mean for structural coefficients.
        real(dp), intent(in) :: prior_precision(:, :) !! p by p Gaussian prior precision matrix.
        real(dp), intent(in) :: proposal_sd(:) !! Positive componentwise random-walk proposal standard deviations.
        real(dp), intent(inout) :: parameter(:) !! Current structural coefficients, replaced by an accepted proposal.
        logical, intent(out) :: accepted !! True when the proposed structural vector is accepted.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or likelihood failure.
        real(dp), allocatable :: current_deviation(:)
        real(dp) :: current_log_density
        real(dp) :: current_prior
        integer :: current_sign
        integer :: i
        real(dp), allocatable :: proposal(:)
        real(dp), allocatable :: proposal_deviation(:)
        real(dp) :: proposal_log_density
        real(dp) :: proposal_prior
        integer :: proposal_sign
        real(dp) :: u
        real(dp) :: z

        info = 0
        accepted = .false.
        if (size(parameter) < 1 .or. size(prior_mean) /= size(parameter) .or. &
            size(prior_precision, 1) /= size(parameter) .or. size(prior_precision, 2) /= size(parameter) .or. &
            size(proposal_sd) /= size(parameter) .or. any(proposal_sd <= 0.0_dp)) then
            info = 1
            return
        end if
        call structural_gaussian_loglik(y, mean_value, residual_covariance, basis, parameter, current_log_density, &
            current_sign, info)
        if (info /= 0) return
        current_deviation = parameter - prior_mean
        current_prior = -0.5_dp * dot_product(current_deviation, matmul(prior_precision, current_deviation))
        proposal = parameter
        do i = 1, size(parameter)
            call rng_normal(state, z)
            proposal(i) = proposal(i) + proposal_sd(i) * z
        end do
        call structural_gaussian_loglik(y, mean_value, residual_covariance, basis, proposal, proposal_log_density, &
            proposal_sign, info)
        if (info /= 0) return
        if (proposal_sign /= current_sign) return
        proposal_deviation = proposal - prior_mean
        proposal_prior = -0.5_dp * dot_product(proposal_deviation, matmul(prior_precision, proposal_deviation))
        call rng_uniform(state, u)
        if (log(u) < proposal_log_density + proposal_prior - current_log_density - current_prior) then
            parameter = proposal
            accepted = .true.
        end if
    end subroutine structural_parameter_random_walk_update

end module mcmcglmm_engine_features
