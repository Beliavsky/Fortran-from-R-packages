! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 parameter-expansion updates; see NOTICE.md and upstream/.
module mcmcglmm_px_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_precision
    use mcmcglmm_distributions, only : inverse_wishart_sample
    implicit none
    private

    type, public :: gaussian_px_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: alpha(:, :)
        real(dp), allocatable :: log_likelihood(:)
    end type gaussian_px_mcmc_result

    public :: gaussian_parameter_expanded_mcmc
    public :: px_alpha_conditional

contains

    pure integer function px_joint_index(trait, effect, effects_per_trait) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed or random coefficient index within the trait block.
        integer, intent(in) :: effects_per_trait !! Number of fixed plus random coefficients in each trait block.

        index_value = (trait - 1) * effects_per_trait + effect
    end function px_joint_index

    pure integer function px_beta_index(trait, effect, fixed_effects) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-effect index within the trait block.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect coefficients per trait.

        index_value = (trait - 1) * fixed_effects + effect
    end function px_beta_index

    pure real(dp) function px_design_cross(x, z, alpha, trait_left, effect_left, trait_right, effect_right) result(value)
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared by all traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled random-effect design shared by all traits.
        real(dp), intent(in) :: alpha(:) !! Trait-specific redundant scale parameters multiplying random-effect designs.
        integer, intent(in) :: trait_left !! One-based trait index of the left design column.
        integer, intent(in) :: effect_left !! One-based fixed-then-random column index for the left trait.
        integer, intent(in) :: trait_right !! One-based trait index of the right design column.
        integer, intent(in) :: effect_right !! One-based fixed-then-random column index for the right trait.
        integer :: p

        p = size(x, 2)
        if (effect_left <= p .and. effect_right <= p) then
            value = dot_product(x(:, effect_left), x(:, effect_right))
        else if (effect_left <= p) then
            value = alpha(trait_right) * dot_product(x(:, effect_left), z(:, effect_right - p))
        else if (effect_right <= p) then
            value = alpha(trait_left) * dot_product(z(:, effect_left - p), x(:, effect_right))
        else
            value = alpha(trait_left) * alpha(trait_right) * &
                dot_product(z(:, effect_left - p), z(:, effect_right - p))
        end if
    end function px_design_cross

    pure real(dp) function px_design_response(x, z, alpha, trait, effect, response) result(value)
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared by all traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled random-effect design shared by all traits.
        real(dp), intent(in) :: alpha(:) !! Trait-specific redundant scales multiplying random-effect designs.
        integer, intent(in) :: trait !! One-based trait index owning this design column.
        integer, intent(in) :: effect !! One-based fixed-then-random coefficient index.
        real(dp), intent(in) :: response(:) !! n-vector response column entering the cross-trait precision product.
        integer :: p

        p = size(x, 2)
        if (effect <= p) then
            value = dot_product(x(:, effect), response)
        else
            value = alpha(trait) * dot_product(z(:, effect - p), response)
        end if
    end function px_design_response

    pure subroutine px_unpack_coefficients(packed, fixed_effects, random_levels, traits, beta, random_effects)
        real(dp), intent(in) :: packed(:) !! Trait-major vector containing fixed then unscaled random coefficients.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect coefficients per trait.
        integer, intent(in) :: random_levels !! Number of random-effect coefficients per trait.
        integer, intent(in) :: traits !! Number of response traits.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated p by traits fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated q by traits unscaled random effects.
        integer :: a
        integer :: k
        integer :: m

        m = fixed_effects + random_levels
        allocate(beta(fixed_effects, traits), random_effects(random_levels, traits))
        do a = 1, traits
            do k = 1, fixed_effects
                beta(k, a) = packed(px_joint_index(a, k, m))
            end do
            do k = 1, random_levels
                random_effects(k, a) = packed(px_joint_index(a, fixed_effects + k, m))
            end do
        end do
    end subroutine px_unpack_coefficients

    pure subroutine px_coefficient_conditional(y, x, z, a_inverse, g_working, r_matrix, alpha, beta_prior_mean, &
                                               beta_prior_precision, state, beta, random_effects, info)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled random-effect design shared across traits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q row precision for the unscaled random effects.
        real(dp), intent(in) :: g_working(:, :) !! traits by traits working covariance of unscaled random effects.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: alpha(:) !! Trait-specific redundant scales multiplying the random-effect design.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the joint Gaussian coefficient draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated sampled unscaled random-effect matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or SPD linear-algebra failure.
        real(dp), allocatable :: beta_mean_vector(:)
        real(dp), allocatable :: beta_rhs(:)
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: g_inverse(:, :)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: k
        integer :: m
        integer :: p
        integer :: q
        integer :: traits

        info = 0
        p = size(x, 2)
        q = size(z, 2)
        traits = size(y, 2)
        m = p + q
        if (size(y, 1) /= size(x, 1) .or. size(y, 1) /= size(z, 1) .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(g_working, 1) /= traits .or. size(g_working, 2) /= traits .or. &
            size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits .or. size(alpha) /= traits .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= traits .or. &
            size(beta_prior_precision, 1) /= p * traits .or. size(beta_prior_precision, 2) /= p * traits) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if
        call inverse_matrix(g_working, g_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if

        allocate(precision(m * traits, m * traits), rhs(m * traits), beta_mean_vector(p * traits))
        precision = 0.0_dp
        rhs = 0.0_dp
        do a = 1, traits
            do k = 1, p
                beta_mean_vector(px_beta_index(a, k, p)) = beta_prior_mean(k, a)
            end do
        end do
        beta_rhs = matmul(beta_prior_precision, beta_mean_vector)

        do a = 1, traits
            do b = 1, traits
                do j = 1, m
                    do i = 1, m
                        precision(px_joint_index(a, i, m), px_joint_index(b, j, m)) = &
                            precision(px_joint_index(a, i, m), px_joint_index(b, j, m)) + &
                            r_inverse(a, b) * px_design_cross(x, z, alpha, a, i, b, j)
                    end do
                end do
                do j = 1, q
                    do i = 1, q
                        precision(px_joint_index(a, p + i, m), px_joint_index(b, p + j, m)) = &
                            precision(px_joint_index(a, p + i, m), px_joint_index(b, p + j, m)) + &
                            g_inverse(a, b) * a_inverse(i, j)
                    end do
                end do
            end do
            do i = 1, m
                do b = 1, traits
                    rhs(px_joint_index(a, i, m)) = rhs(px_joint_index(a, i, m)) + &
                        r_inverse(a, b) * px_design_response(x, z, alpha, a, i, y(:, b))
                end do
            end do
            do i = 1, p
                do k = 1, p
                    precision(px_joint_index(a, i, m), px_joint_index(a, k, m)) = &
                        precision(px_joint_index(a, i, m), px_joint_index(a, k, m)) + &
                        beta_prior_precision(px_beta_index(a, i, p), px_beta_index(a, k, p))
                end do
                rhs(px_joint_index(a, i, m)) = rhs(px_joint_index(a, i, m)) + beta_rhs(px_beta_index(a, i, p))
            end do
        end do
        do a = 1, traits
            do b = 1, traits
                if (a == b) cycle
                do j = 1, p
                    do i = 1, p
                        precision(px_joint_index(a, i, m), px_joint_index(b, j, m)) = &
                            precision(px_joint_index(a, i, m), px_joint_index(b, j, m)) + &
                            beta_prior_precision(px_beta_index(a, i, p), px_beta_index(b, j, p))
                    end do
                end do
            end do
        end do

        call sample_mvn_precision(state, rhs, precision, packed, conditional_mean, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call px_unpack_coefficients(packed, p, q, traits, beta, random_effects)
    end subroutine px_coefficient_conditional

    pure subroutine px_alpha_conditional(residual_base, random_prediction, r_matrix, prior_mean, prior_precision, &
                                         state, alpha, conditional_mean, info)
        real(dp), intent(in) :: residual_base(:, :) !! n by traits matrix y-X beta before parameter-expanded random effects.
        real(dp), intent(in) :: random_prediction(:, :) !! n by traits unscaled Z u working random-effect predictions.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: prior_mean(:) !! traits-vector Gaussian prior mean for redundant scale parameters.
        real(dp), intent(in) :: prior_precision(:, :) !! traits by traits Gaussian prior precision for redundant scales.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the redundant-scale Gaussian draw.
        real(dp), allocatable, intent(out) :: alpha(:) !! Allocated sampled redundant scale vector.
        real(dp), allocatable, intent(out) :: conditional_mean(:) !! Allocated full-conditional mean of alpha.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or SPD linear-algebra failure.
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        integer :: a
        integer :: b
        integer :: traits

        info = 0
        traits = size(residual_base, 2)
        if (size(random_prediction, 1) /= size(residual_base, 1) .or. &
            size(random_prediction, 2) /= traits .or. size(r_matrix, 1) /= traits .or. &
            size(r_matrix, 2) /= traits .or. size(prior_mean) /= traits .or. &
            size(prior_precision, 1) /= traits .or. size(prior_precision, 2) /= traits) then
            allocate(alpha(0), conditional_mean(0))
            info = 1
            return
        end if
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(alpha(0), conditional_mean(0))
            return
        end if
        precision = prior_precision
        rhs = matmul(prior_precision, prior_mean)
        do a = 1, traits
            do b = 1, traits
                precision(a, b) = precision(a, b) + &
                    r_inverse(a, b) * dot_product(random_prediction(:, a), random_prediction(:, b))
                rhs(a) = rhs(a) + r_inverse(a, b) * &
                    dot_product(random_prediction(:, a), residual_base(:, b))
            end do
        end do
        call sample_mvn_precision(state, rhs, precision, alpha, conditional_mean, info)
    end subroutine px_alpha_conditional

    pure subroutine expanded_effects(random_effects, g_working, alpha, scaled_effects, expanded_g)
        real(dp), intent(in) :: random_effects(:, :) !! q by traits unscaled working random effects.
        real(dp), intent(in) :: g_working(:, :) !! traits by traits working covariance of unscaled random effects.
        real(dp), intent(in) :: alpha(:) !! Trait-specific redundant scale parameters.
        real(dp), allocatable, intent(out) :: scaled_effects(:, :) !! Allocated random effects on the model/output scale.
        real(dp), allocatable, intent(out) :: expanded_g(:, :) !! Allocated covariance on the model/output scale.
        integer :: a
        integer :: b
        integer :: traits

        traits = size(alpha)
        allocate(scaled_effects(size(random_effects, 1), traits), expanded_g(traits, traits))
        do a = 1, traits
            scaled_effects(:, a) = alpha(a) * random_effects(:, a)
            do b = 1, traits
                expanded_g(a, b) = alpha(a) * g_working(a, b) * alpha(b)
            end do
        end do
    end subroutine expanded_effects

    pure subroutine gaussian_px_log_likelihood(y, x, z, beta, scaled_effects, r_matrix, log_likelihood, info)
        real(dp), intent(in) :: y(:, :) !! n by traits observed Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: beta(:, :) !! p by traits fixed-effect matrix.
        real(dp), intent(in) :: scaled_effects(:, :) !! q by traits parameter-expanded random effects on model scale.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(out) :: log_likelihood !! Sum of row-wise multivariate-normal log densities.
        integer, intent(out) :: info !! Zero on success; nonzero on covariance or shape failure.
        real(dp), allocatable :: mean_row(:)
        real(dp), allocatable :: zero_mean(:)
        real(dp) :: row_log_density
        integer :: i
        integer :: traits

        traits = size(y, 2)
        if (size(x, 1) /= size(y, 1) .or. size(z, 1) /= size(y, 1) .or. size(beta, 1) /= size(x, 2) .or. &
            size(beta, 2) /= traits .or. size(scaled_effects, 1) /= size(z, 2) .or. &
            size(scaled_effects, 2) /= traits .or. size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits) then
            log_likelihood = -huge(1.0_dp)
            info = 1
            return
        end if
        allocate(mean_row(traits), zero_mean(traits))
        zero_mean = 0.0_dp
        log_likelihood = 0.0_dp
        do i = 1, size(y, 1)
            mean_row = matmul(x(i, :), beta) + matmul(z(i, :), scaled_effects)
            call mvn_log_density(y(i, :) - mean_row, zero_mean, r_matrix, row_log_density, info)
            if (info /= 0) then
                log_likelihood = -huge(1.0_dp)
                return
            end if
            log_likelihood = log_likelihood + row_log_density
        end do
    end subroutine gaussian_px_log_likelihood

    pure subroutine gaussian_parameter_expanded_mcmc(y, x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
                                                      g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, &
                                                      alpha_prior_mean, alpha_prior_precision, update_r, &
                                                      iterations, burn, thin, state, result, info)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled random-effect design shared across traits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q row precision for random effects, such as pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :) !! traits by traits IW prior scale for the working G covariance.
        real(dp), intent(in) :: g_prior_df !! Working G inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! traits by traits IW prior scale for residual covariance R.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: alpha_prior_mean(:) !! traits-vector Gaussian prior mean corresponding to alpha.mu.
        real(dp), intent(in) :: alpha_prior_precision(:, :) !! traits-square Gaussian prior precision inverse(alpha.V).
        logical, intent(in) :: update_r !! If true, sample residual covariance; otherwise retain its initialized value.
        integer, intent(in) :: iterations !! Total Gibbs iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by coefficient, covariance, and alpha draws.
        type(gaussian_px_mcmc_result), intent(out) :: result !! Retained parameter-expanded posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: alpha_mean(:)
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: expanded_g(:, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: g_working(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: random_prediction(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: residual_base(:, :)
        real(dp), allocatable :: scaled_effects(:, :)
        real(dp) :: log_likelihood
        integer :: a
        integer :: iteration
        integer :: n
        integer :: nsave
        integer :: p
        integer :: q
        integer :: save_index
        integer :: traits

        info = 0
        n = size(y, 1)
        traits = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        if (n < 1 .or. traits < 1 .or. p < 1 .or. q < 1 .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= traits .or. &
            size(beta_prior_precision, 1) /= p * traits .or. size(beta_prior_precision, 2) /= p * traits .or. &
            size(g_prior_scale, 1) /= traits .or. size(g_prior_scale, 2) /= traits .or. &
            size(r_prior_scale, 1) /= traits .or. size(r_prior_scale, 2) /= traits .or. &
            size(alpha_prior_mean) /= traits .or. size(alpha_prior_precision, 1) /= traits .or. &
            size(alpha_prior_precision, 2) /= traits) then
            info = 1
            return
        end if
        if (iterations <= burn .or. burn < 0 .or. thin < 1 .or. &
            g_prior_df <= real(traits - 1, dp) .or. r_prior_df <= real(traits - 1, dp)) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(result%beta(p, traits, nsave), result%random_effects(q, traits, nsave))
        allocate(result%g(traits, traits, nsave), result%r(traits, traits, nsave))
        allocate(result%alpha(traits, nsave), result%log_likelihood(nsave))
        beta = beta_prior_mean
        allocate(random_effects(q, traits), alpha(traits))
        random_effects = 0.0_dp
        alpha = 1.0_dp
        g_working = g_prior_scale / max(g_prior_df - real(traits + 1, dp), 1.0_dp)
        r_matrix = r_prior_scale / max(r_prior_df - real(traits + 1, dp), 1.0_dp)
        save_index = 0

        do iteration = 1, iterations
            call px_coefficient_conditional(y, x, z, a_inverse, g_working, r_matrix, alpha, beta_prior_mean, &
                beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            g_scale_post = g_prior_scale + matmul(transpose(random_effects), matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_post, g_prior_df + real(q, dp), g_working, info)
            if (info /= 0) return

            residual_base = y - matmul(x, beta)
            random_prediction = matmul(z, random_effects)
            allocate(residual(n, traits))
            do a = 1, traits
                residual(:, a) = residual_base(:, a) - alpha(a) * random_prediction(:, a)
            end do
            if (update_r) then
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if
            deallocate(residual)

            call px_alpha_conditional(residual_base, random_prediction, r_matrix, alpha_prior_mean, &
                alpha_prior_precision, state, alpha, alpha_mean, info)
            if (info /= 0) return

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                call expanded_effects(random_effects, g_working, alpha, scaled_effects, expanded_g)
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = scaled_effects
                result%g(:, :, save_index) = expanded_g
                result%r(:, :, save_index) = r_matrix
                result%alpha(:, save_index) = alpha
                call gaussian_px_log_likelihood(y, x, z, beta, scaled_effects, r_matrix, log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine gaussian_parameter_expanded_mcmc

end module mcmcglmm_px_sampler
