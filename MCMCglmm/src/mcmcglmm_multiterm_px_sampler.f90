! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 parameter-expansion and multi-G updates; see NOTICE.md and upstream/.
module mcmcglmm_multiterm_px_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state, rng_normal, rng_uniform
    use mcmcglmm_matrix, only : sample_mvn_precision, mvn_log_density
    use mcmcglmm_distributions, only : inverse_wishart_sample
    use mcmcglmm_multiterm_sampler, only : multi_term_gaussian_loglik
    use mcmcglmm_families, only : scalar_family_loglik
    implicit none
    private

    type, public :: multi_term_gaussian_px_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: alpha(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
    end type multi_term_gaussian_px_mcmc_result

    type, public :: multi_term_family_px_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: alpha(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:, :)
        real(dp) :: acceptance_rate = 0.0_dp
    end type multi_term_family_px_mcmc_result

    public :: multi_term_px_alpha_conditional
    public :: multi_term_gaussian_parameter_expanded_mcmc
    public :: heterogeneous_multi_term_parameter_expanded_mcmc

contains

    pure integer function pxmt_joint_index(trait, effect, effects_per_trait) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-then-random coefficient index within a trait block.
        integer, intent(in) :: effects_per_trait !! Number of fixed plus random coefficients per trait.

        index_value = (trait - 1) * effects_per_trait + effect
    end function pxmt_joint_index

    pure integer function pxmt_beta_index(trait, effect, fixed_effects) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-effect index within a trait block.
        integer, intent(in) :: fixed_effects !! Number of fixed effects per trait.

        index_value = (trait - 1) * fixed_effects + effect
    end function pxmt_beta_index

    pure integer function pxmt_alpha_index(trait, term, traits) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: term !! One-based random-effect covariance term.
        integer, intent(in) :: traits !! Number of response traits.

        index_value = (term - 1) * traits + trait
    end function pxmt_alpha_index

    pure real(dp) function scaled_design_cross(x, z, random_term, alpha, trait_left, effect_left, &
                                               trait_right, effect_right) result(value)
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term labels for random columns.
        real(dp), intent(in) :: alpha(:, :) !! traits by n_term redundant scale matrix.
        integer, intent(in) :: trait_left !! One-based trait owning the left design column.
        integer, intent(in) :: effect_left !! One-based fixed-then-random left effect index.
        integer, intent(in) :: trait_right !! One-based trait owning the right design column.
        integer, intent(in) :: effect_right !! One-based fixed-then-random right effect index.
        real(dp) :: left_scale
        integer :: p
        real(dp) :: right_scale

        p = size(x, 2)
        left_scale = 1.0_dp
        right_scale = 1.0_dp
        if (effect_left > p) left_scale = alpha(trait_left, random_term(effect_left - p))
        if (effect_right > p) right_scale = alpha(trait_right, random_term(effect_right - p))
        if (effect_left <= p .and. effect_right <= p) then
            value = dot_product(x(:, effect_left), x(:, effect_right))
        else if (effect_left <= p) then
            value = right_scale * dot_product(x(:, effect_left), z(:, effect_right - p))
        else if (effect_right <= p) then
            value = left_scale * dot_product(z(:, effect_left - p), x(:, effect_right))
        else
            value = left_scale * right_scale * dot_product(z(:, effect_left - p), z(:, effect_right - p))
        end if
    end function scaled_design_cross

    pure real(dp) function scaled_design_response(x, z, random_term, alpha, trait, effect, response) result(value)
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term labels for random columns.
        real(dp), intent(in) :: alpha(:, :) !! traits by n_term redundant scale matrix.
        integer, intent(in) :: trait !! One-based trait owning the design column.
        integer, intent(in) :: effect !! One-based fixed-then-random effect index.
        real(dp), intent(in) :: response(:) !! n-vector response column entering the residual precision product.
        integer :: p

        p = size(x, 2)
        if (effect <= p) then
            value = dot_product(x(:, effect), response)
        else
            value = alpha(trait, random_term(effect - p)) * dot_product(z(:, effect - p), response)
        end if
    end function scaled_design_response

    pure subroutine multi_term_px_coefficient_conditional(y, x, z, random_term, a_inverse, g_working, r_matrix, &
                                                          alpha, beta_prior_mean, beta_prior_precision, state, &
                                                          beta, random_effects, info)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term label for each random coefficient.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random effects.
        real(dp), intent(in) :: g_working(:, :, :) !! traits by traits by n_term working G covariances.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: alpha(:, :) !! traits by n_term redundant random-effect scales.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Gaussian coefficient draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated unscaled q by traits random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or SPD failure.
        real(dp), allocatable :: beta_mean_vector(:)
        real(dp), allocatable :: beta_rhs(:)
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: g_inverse(:, :, :)
        real(dp), allocatable :: inverse_block(:, :)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        integer :: a
        integer :: b
        integer :: effect_left
        integer :: effect_right
        integer :: m
        integer :: nterm
        integer :: p
        integer :: q
        integer :: term
        integer :: traits

        info = 0
        p = size(x, 2)
        q = size(z, 2)
        traits = size(y, 2)
        nterm = size(g_working, 3)
        m = p + q
        if (size(y, 1) /= size(x, 1) .or. size(y, 1) /= size(z, 1) .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(g_working, 1) /= traits .or. &
            size(g_working, 2) /= traits .or. size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits .or. &
            size(alpha, 1) /= traits .or. size(alpha, 2) /= nterm .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= traits .or. size(beta_prior_precision, 1) /= p * traits .or. &
            size(beta_prior_precision, 2) /= p * traits .or. any(random_term < 1) .or. &
            any(random_term > nterm)) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if

        allocate(g_inverse(traits, traits, nterm))
        do term = 1, nterm
            call inverse_matrix(g_working(:, :, term), inverse_block, info)
            if (info /= 0) then
                allocate(beta(0, 0), random_effects(0, 0))
                return
            end if
            g_inverse(:, :, term) = inverse_block
        end do
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if

        allocate(beta_mean_vector(p * traits), precision(m * traits, m * traits), rhs(m * traits))
        do a = 1, traits
            do effect_left = 1, p
                beta_mean_vector(pxmt_beta_index(a, effect_left, p)) = beta_prior_mean(effect_left, a)
            end do
        end do
        beta_rhs = matmul(beta_prior_precision, beta_mean_vector)
        precision = 0.0_dp
        rhs = 0.0_dp
        do a = 1, traits
            do b = 1, traits
                do effect_right = 1, m
                    do effect_left = 1, m
                        precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(b, effect_right, m)) = &
                            precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(b, effect_right, m)) + &
                            r_inverse(a, b) * scaled_design_cross(x, z, random_term, alpha, a, effect_left, &
                                b, effect_right)
                    end do
                end do
                do effect_right = 1, q
                    do effect_left = 1, q
                        if (random_term(effect_left) /= random_term(effect_right)) cycle
                        term = random_term(effect_left)
                        precision(pxmt_joint_index(a, p + effect_left, m), &
                            pxmt_joint_index(b, p + effect_right, m)) = &
                            precision(pxmt_joint_index(a, p + effect_left, m), &
                            pxmt_joint_index(b, p + effect_right, m)) + &
                            g_inverse(a, b, term) * a_inverse(effect_left, effect_right)
                    end do
                end do
            end do
            do effect_left = 1, m
                do b = 1, traits
                    rhs(pxmt_joint_index(a, effect_left, m)) = rhs(pxmt_joint_index(a, effect_left, m)) + &
                        r_inverse(a, b) * scaled_design_response(x, z, random_term, alpha, a, effect_left, y(:, b))
                end do
            end do
            do effect_left = 1, p
                do effect_right = 1, p
                    precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(a, effect_right, m)) = &
                        precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(a, effect_right, m)) + &
                        beta_prior_precision(pxmt_beta_index(a, effect_left, p), pxmt_beta_index(a, effect_right, p))
                end do
                rhs(pxmt_joint_index(a, effect_left, m)) = rhs(pxmt_joint_index(a, effect_left, m)) + &
                    beta_rhs(pxmt_beta_index(a, effect_left, p))
            end do
        end do
        do a = 1, traits
            do b = 1, traits
                if (a == b) cycle
                do effect_right = 1, p
                    do effect_left = 1, p
                        precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(b, effect_right, m)) = &
                            precision(pxmt_joint_index(a, effect_left, m), pxmt_joint_index(b, effect_right, m)) + &
                            beta_prior_precision(pxmt_beta_index(a, effect_left, p), pxmt_beta_index(b, effect_right, p))
                    end do
                end do
            end do
        end do

        call sample_mvn_precision(state, rhs, precision, packed, conditional_mean, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        allocate(beta(p, traits), random_effects(q, traits))
        do a = 1, traits
            do effect_left = 1, p
                beta(effect_left, a) = packed(pxmt_joint_index(a, effect_left, m))
            end do
            do effect_left = 1, q
                random_effects(effect_left, a) = packed(pxmt_joint_index(a, p + effect_left, m))
            end do
        end do
    end subroutine multi_term_px_coefficient_conditional

    pure subroutine multi_term_px_predictions(z, random_term, random_effects, predictions)
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design matrix.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term labels.
        real(dp), intent(in) :: random_effects(:, :) !! q by traits unscaled random effects.
        real(dp), allocatable, intent(out) :: predictions(:, :, :) !! Allocated n by traits by n_term term predictions.
        integer :: j
        integer :: term

        allocate(predictions(size(z, 1), size(random_effects, 2), maxval(random_term)))
        predictions = 0.0_dp
        do j = 1, size(z, 2)
            term = random_term(j)
            predictions(:, :, term) = predictions(:, :, term) + &
                spread(z(:, j), 2, size(random_effects, 2)) * spread(random_effects(j, :), 1, size(z, 1))
        end do
    end subroutine multi_term_px_predictions

    pure subroutine multi_term_px_alpha_conditional(residual_base, term_predictions, r_matrix, prior_mean, &
                                                    prior_precision, state, alpha, conditional_mean, info)
        real(dp), intent(in) :: residual_base(:, :) !! n by traits response after subtracting fixed effects.
        real(dp), intent(in) :: term_predictions(:, :, :) !! n by traits by n_term unscaled random predictions.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: prior_mean(:, :) !! traits by n_term Gaussian prior mean for redundant scales.
        real(dp), intent(in) :: prior_precision(:, :) !! traits*n_term square Gaussian prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the redundant-scale Gaussian draw.
        real(dp), allocatable, intent(out) :: alpha(:, :) !! Allocated sampled traits by n_term redundant scales.
        real(dp), allocatable, intent(out) :: conditional_mean(:, :) !! Allocated full-conditional mean matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or SPD failure.
        real(dp), allocatable :: mean_vector(:)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: packed_mean(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        integer :: a
        integer :: b
        integer :: nterm
        integer :: s
        integer :: t
        integer :: traits

        info = 0
        traits = size(residual_base, 2)
        nterm = size(term_predictions, 3)
        if (size(term_predictions, 1) /= size(residual_base, 1) .or. size(term_predictions, 2) /= traits .or. &
            size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits .or. &
            size(prior_mean, 1) /= traits .or. size(prior_mean, 2) /= nterm .or. &
            size(prior_precision, 1) /= traits * nterm .or. size(prior_precision, 2) /= traits * nterm) then
            allocate(alpha(0, 0), conditional_mean(0, 0))
            info = 1
            return
        end if
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(alpha(0, 0), conditional_mean(0, 0))
            return
        end if
        allocate(mean_vector(traits * nterm), precision(traits * nterm, traits * nterm), rhs(traits * nterm))
        do t = 1, nterm
            do a = 1, traits
                mean_vector(pxmt_alpha_index(a, t, traits)) = prior_mean(a, t)
            end do
        end do
        precision = prior_precision
        rhs = matmul(prior_precision, mean_vector)
        do t = 1, nterm
            do a = 1, traits
                do b = 1, traits
                    rhs(pxmt_alpha_index(a, t, traits)) = rhs(pxmt_alpha_index(a, t, traits)) + &
                        r_inverse(a, b) * dot_product(term_predictions(:, a, t), residual_base(:, b))
                    do s = 1, nterm
                        precision(pxmt_alpha_index(a, t, traits), pxmt_alpha_index(b, s, traits)) = &
                            precision(pxmt_alpha_index(a, t, traits), pxmt_alpha_index(b, s, traits)) + &
                            r_inverse(a, b) * &
                            dot_product(term_predictions(:, a, t), term_predictions(:, b, s))
                    end do
                end do
            end do
        end do
        call sample_mvn_precision(state, rhs, precision, packed, packed_mean, info)
        if (info /= 0) then
            allocate(alpha(0, 0), conditional_mean(0, 0))
            return
        end if
        allocate(alpha(traits, nterm), conditional_mean(traits, nterm))
        do t = 1, nterm
            do a = 1, traits
                alpha(a, t) = packed(pxmt_alpha_index(a, t, traits))
                conditional_mean(a, t) = packed_mean(pxmt_alpha_index(a, t, traits))
            end do
        end do
    end subroutine multi_term_px_alpha_conditional

    pure subroutine scale_multi_term_effects(random_term, random_effects, g_working, alpha, scaled_effects, expanded_g)
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term label for each random coefficient.
        real(dp), intent(in) :: random_effects(:, :) !! q by traits unscaled working random effects.
        real(dp), intent(in) :: g_working(:, :, :) !! traits by traits by n_term working covariance blocks.
        real(dp), intent(in) :: alpha(:, :) !! traits by n_term redundant scale matrix.
        real(dp), allocatable, intent(out) :: scaled_effects(:, :) !! Allocated q by traits output-scale random effects.
        real(dp), allocatable, intent(out) :: expanded_g(:, :, :) !! Allocated output-scale covariance blocks.
        integer :: a
        integer :: b
        integer :: j
        integer :: term
        integer :: traits

        traits = size(random_effects, 2)
        allocate(scaled_effects(size(random_effects, 1), traits), expanded_g(traits, traits, size(g_working, 3)))
        do j = 1, size(random_effects, 1)
            term = random_term(j)
            do a = 1, traits
                scaled_effects(j, a) = alpha(a, term) * random_effects(j, a)
            end do
        end do
        do term = 1, size(g_working, 3)
            do b = 1, traits
                do a = 1, traits
                    expanded_g(a, b, term) = alpha(a, term) * g_working(a, b, term) * alpha(b, term)
                end do
            end do
        end do
    end subroutine scale_multi_term_effects

    pure subroutine multi_term_gaussian_parameter_expanded_mcmc(y, x, z, random_term, a_inverse, beta_prior_mean, &
                                                                beta_prior_precision, g_prior_scale, g_prior_df, &
                                                                r_prior_scale, r_prior_df, alpha_prior_mean, &
                                                                alpha_prior_precision, update_r, iterations, burn, &
                                                                thin, state, result, info)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term label for each random coefficient.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random effects.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! traits by traits by n_term working-G prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term working-G prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! traits by traits residual covariance prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: alpha_prior_mean(:, :) !! traits by n_term redundant-scale Gaussian prior means.
        real(dp), intent(in) :: alpha_prior_precision(:, :) !! traits*n_term square redundant-scale prior precision.
        logical, intent(in) :: update_r !! If true, update residual covariance; otherwise retain its initialized value.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by coefficient, covariance, and alpha draws.
        type(multi_term_gaussian_px_mcmc_result), intent(out) :: result !! Retained multi-G parameter-expanded samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or numerical failure.
        real(dp), allocatable :: alpha(:, :)
        real(dp), allocatable :: alpha_mean(:, :)
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: expanded_g(:, :, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: g_working(:, :, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: residual_base(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: scaled_effects(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        real(dp), allocatable :: term_predictions(:, :, :)
        integer :: i
        integer :: iteration
        real(dp) :: log_likelihood
        integer :: n
        integer :: nsave
        integer :: nterm
        integer :: p
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: save_index
        integer :: term
        integer :: traits

        info = 0
        n = size(y, 1)
        traits = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (n < 1 .or. traits < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. size(x, 1) /= n .or. &
            size(z, 1) /= n .or. size(random_term) /= q .or. size(a_inverse, 1) /= q .or. &
            size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= traits .or. size(beta_prior_precision, 1) /= p * traits .or. &
            size(beta_prior_precision, 2) /= p * traits .or. size(g_prior_scale, 1) /= traits .or. &
            size(g_prior_scale, 2) /= traits .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= traits .or. size(r_prior_scale, 2) /= traits .or. &
            size(alpha_prior_mean, 1) /= traits .or. size(alpha_prior_mean, 2) /= nterm .or. &
            size(alpha_prior_precision, 1) /= traits * nterm .or. &
            size(alpha_prior_precision, 2) /= traits * nterm .or. any(random_term < 1) .or. &
            any(random_term > nterm)) then
            info = 1
            return
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(traits - 1, dp)) then
                info = 2
                return
            end if
        end do
        if (r_prior_df <= real(traits - 1, dp) .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(result%beta(p, traits, nsave), result%random_effects(q, traits, nsave))
        allocate(result%g(traits, traits, nterm, nsave), result%r(traits, traits, nsave))
        allocate(result%alpha(traits, nterm, nsave), result%log_likelihood(nsave))
        allocate(g_working(traits, traits, nterm), random_effects(q, traits))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        alpha = alpha_prior_mean
        do term = 1, nterm
            g_working(:, :, term) = g_prior_scale(:, :, term) / &
                max(g_prior_df(term) - real(traits + 1, dp), 1.0_dp)
        end do
        r_matrix = r_prior_scale / max(r_prior_df - real(traits + 1, dp), 1.0_dp)
        save_index = 0

        do iteration = 1, iterations
            call multi_term_px_coefficient_conditional(y, x, z, random_term, a_inverse, g_working, r_matrix, alpha, &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            do term = 1, nterm
                qterm = count(random_term == term)
                allocate(term_indices(qterm))
                i = 0
                do random_index = 1, q
                    if (random_term(random_index) /= term) cycle
                    i = i + 1
                    term_indices(i) = random_index
                end do
                term_effects = random_effects(term_indices, :)
                term_a_inverse = a_inverse(term_indices, term_indices)
                g_scale_post = g_prior_scale(:, :, term) + &
                    matmul(transpose(term_effects), matmul(term_a_inverse, term_effects))
                call inverse_wishart_sample(state, g_scale_post, g_prior_df(term) + real(qterm, dp), sampled_g, info)
                deallocate(term_indices)
                if (info /= 0) return
                g_working(:, :, term) = sampled_g
            end do

            residual_base = y - matmul(x, beta)
            call multi_term_px_predictions(z, random_term, random_effects, term_predictions)
            call multi_term_px_alpha_conditional(residual_base, term_predictions, r_matrix, alpha_prior_mean, &
                alpha_prior_precision, state, alpha, alpha_mean, info)
            if (info /= 0) return
            call scale_multi_term_effects(random_term, random_effects, g_working, alpha, scaled_effects, expanded_g)

            if (update_r) then
                residual = y - matmul(x, beta) - matmul(z, scaled_effects)
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = scaled_effects
                result%g(:, :, :, save_index) = expanded_g
                result%r(:, :, save_index) = r_matrix
                result%alpha(:, :, save_index) = alpha
                call multi_term_gaussian_loglik(y, x, z, beta, scaled_effects, r_matrix, log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine multi_term_gaussian_parameter_expanded_mcmc


    pure logical function px_scalar_family_supported(family) result(supported)
        integer, intent(in) :: family !! Native scalar family code checked for heterogeneous parameter expansion.

        select case (family)
        case (1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 22, 23, 24)
            supported = .true.
        case default
            supported = .false.
        end select
    end function px_scalar_family_supported

    pure subroutine heterogeneous_multi_term_parameter_expanded_mcmc(family, y, additional, additional2, x, z, &
                                                                       random_term, a_inverse, beta_prior_mean, &
                                                                       beta_prior_precision, g_prior_scale, &
                                                                       g_prior_df, r_prior_scale, r_prior_df, &
                                                                       alpha_prior_mean, alpha_prior_precision, &
                                                                       update_r, proposal_scale, iterations, burn, &
                                                                       thin, state, result, info, observed)
        integer, intent(in) :: family(:) !! Length-t native scalar family codes for the heterogeneous response traits.
        real(dp), intent(in) :: y(:, :) !! n by t response matrix after native-family preprocessing.
        real(dp), intent(in) :: additional(:, :) !! n by t first family-specific auxiliary values such as trials or bounds.
        real(dp), intent(in) :: additional2(:, :) !! n by t second family-specific auxiliary values such as degrees of freedom.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-term label for each random coefficient.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random effects.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by t Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*t square Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! t by t by n_term working-G inverse-Wishart prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term working-G inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! t by t latent residual inverse-Wishart prior scale.
        real(dp), intent(in) :: r_prior_df !! Latent residual inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: alpha_prior_mean(:, :) !! t by n_term Gaussian redundant-scale prior means.
        real(dp), intent(in) :: alpha_prior_precision(:, :) !! t*n_term square redundant-scale prior precision.
        logical, intent(in) :: update_r !! If true, sample the full latent residual covariance; otherwise keep it fixed.
        real(dp), intent(in) :: proposal_scale !! Positive one-coordinate latent random-walk multiplier.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by latent, coefficient, covariance, and alpha draws.
        type(multi_term_family_px_mcmc_result), intent(out) :: result !! Retained output-scale heterogeneous PX samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        logical, optional, intent(in) :: observed(:, :) !! Optional n by t response mask; false entries omit observation likelihood.
        real(dp), allocatable :: alpha(:, :)
        real(dp), allocatable :: alpha_mean(:, :)
        real(dp), allocatable :: beta(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: expanded_g(:, :, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: g_working(:, :, :)
        real(dp), allocatable :: liability(:, :)
        real(dp) :: log_likelihood
        logical, allocatable :: observed_mask(:, :)
        real(dp) :: proposal
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp), allocatable :: proposal_row(:)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: residual_base(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: scaled_effects(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        real(dp), allocatable :: term_predictions(:, :, :)
        real(dp) :: u
        real(dp) :: zdraw
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: j
        integer :: n
        integer :: nsave
        integer :: nterm
        integer :: p
        integer :: proposals
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: save_index
        integer :: term
        integer :: traits

        info = 0
        n = size(y, 1)
        traits = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (n < 1 .or. traits < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. size(family) /= traits .or. &
            size(additional, 1) /= n .or. size(additional, 2) /= traits .or. &
            size(additional2, 1) /= n .or. size(additional2, 2) /= traits .or. size(x, 1) /= n .or. &
            size(z, 1) /= n .or. size(random_term) /= q .or. size(a_inverse, 1) /= q .or. &
            size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= traits .or. size(beta_prior_precision, 1) /= p * traits .or. &
            size(beta_prior_precision, 2) /= p * traits .or. size(g_prior_scale, 1) /= traits .or. &
            size(g_prior_scale, 2) /= traits .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= traits .or. size(r_prior_scale, 2) /= traits .or. &
            size(alpha_prior_mean, 1) /= traits .or. size(alpha_prior_mean, 2) /= nterm .or. &
            size(alpha_prior_precision, 1) /= traits * nterm .or. &
            size(alpha_prior_precision, 2) /= traits * nterm .or. any(random_term < 1) .or. &
            any(random_term > nterm)) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed, 1) /= n .or. size(observed, 2) /= traits) then
                info = 1
                return
            end if
        end if
        do j = 1, traits
            if (.not. px_scalar_family_supported(family(j))) then
                info = 2
                return
            end if
        end do
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(traits - 1, dp)) then
                info = 3
                return
            end if
        end do
        if (r_prior_df <= real(traits - 1, dp) .or. proposal_scale <= 0.0_dp .or. iterations <= burn .or. &
            burn < 0 .or. thin < 1) then
            info = 3
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, traits, nsave), result%random_effects(q, traits, nsave))
        allocate(result%g(traits, traits, nterm, nsave), result%r(traits, traits, nsave))
        allocate(result%alpha(traits, nterm, nsave), result%log_likelihood(nsave))
        allocate(result%last_liability(n, traits), g_working(traits, traits, nterm))
        allocate(random_effects(q, traits), liability(n, traits), proposal_row(traits), observed_mask(n, traits))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        alpha = alpha_prior_mean
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        do term = 1, nterm
            g_working(:, :, term) = g_prior_scale(:, :, term) / &
                max(g_prior_df(term) - real(traits + 1, dp), 1.0_dp)
        end do
        r_matrix = r_prior_scale / max(r_prior_df - real(traits + 1, dp), 1.0_dp)
        call scale_multi_term_effects(random_term, random_effects, g_working, alpha, scaled_effects, expanded_g)
        eta = matmul(x, beta) + matmul(z, scaled_effects)
        liability = eta
        do j = 1, traits
            if (family(j) == 1) then
                do i = 1, n
                    if (observed_mask(i, j)) liability(i, j) = y(i, j)
                end do
            else if (family(j) == 6) then
                do i = 1, n
                    if (observed_mask(i, j)) liability(i, j) = min(max(liability(i, j), y(i, j)), additional(i, j))
                end do
            end if
        end do
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            call scale_multi_term_effects(random_term, random_effects, g_working, alpha, scaled_effects, expanded_g)
            eta = matmul(x, beta) + matmul(z, scaled_effects)
            do i = 1, n
                do j = 1, traits
                    if (family(j) == 1 .and. observed_mask(i, j)) cycle
                    proposal_row = liability(i, :)
                    call rng_normal(state, zdraw)
                    proposal = liability(i, j) + proposal_scale * sqrt(r_matrix(j, j)) * zdraw
                    proposal_row(j) = proposal
                    call mvn_log_density(liability(i, :), eta(i, :), r_matrix, current_log_prior, info)
                    if (info /= 0) return
                    call mvn_log_density(proposal_row, eta(i, :), r_matrix, proposal_log_prior, info)
                    if (info /= 0) return
                    current_log_target = current_log_prior
                    proposal_log_target = proposal_log_prior
                    if (observed_mask(i, j)) then
                        current_log_target = current_log_target + scalar_family_loglik(family(j), y(i, j), &
                            liability(i, j), additional(i, j), additional2(i, j))
                        proposal_log_target = proposal_log_target + scalar_family_loglik(family(j), y(i, j), &
                            proposal, additional(i, j), additional2(i, j))
                    end if
                    call rng_uniform(state, u)
                    proposals = proposals + 1
                    if (log(u) < proposal_log_target - current_log_target) then
                        liability(i, j) = proposal
                        accepted = accepted + 1
                    end if
                end do
            end do

            call multi_term_px_coefficient_conditional(liability, x, z, random_term, a_inverse, g_working, r_matrix, &
                alpha, beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            do term = 1, nterm
                qterm = count(random_term == term)
                allocate(term_indices(qterm))
                i = 0
                do random_index = 1, q
                    if (random_term(random_index) /= term) cycle
                    i = i + 1
                    term_indices(i) = random_index
                end do
                term_effects = random_effects(term_indices, :)
                term_a_inverse = a_inverse(term_indices, term_indices)
                g_scale_post = g_prior_scale(:, :, term) + &
                    matmul(transpose(term_effects), matmul(term_a_inverse, term_effects))
                call inverse_wishart_sample(state, g_scale_post, g_prior_df(term) + real(qterm, dp), sampled_g, info)
                deallocate(term_indices)
                if (info /= 0) return
                g_working(:, :, term) = sampled_g
            end do

            residual_base = liability - matmul(x, beta)
            call multi_term_px_predictions(z, random_term, random_effects, term_predictions)
            call multi_term_px_alpha_conditional(residual_base, term_predictions, r_matrix, alpha_prior_mean, &
                alpha_prior_precision, state, alpha, alpha_mean, info)
            if (info /= 0) return
            call scale_multi_term_effects(random_term, random_effects, g_working, alpha, scaled_effects, expanded_g)
            if (update_r) then
                residual = liability - matmul(x, beta) - matmul(z, scaled_effects)
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = scaled_effects
                result%g(:, :, :, save_index) = expanded_g
                result%r(:, :, save_index) = r_matrix
                result%alpha(:, :, save_index) = alpha
                eta = matmul(x, beta) + matmul(z, scaled_effects)
                log_likelihood = 0.0_dp
                do i = 1, n
                    call mvn_log_density(liability(i, :), eta(i, :), r_matrix, current_log_prior, info)
                    if (info /= 0) return
                    log_likelihood = log_likelihood + current_log_prior
                    do j = 1, traits
                        if (family(j) /= 1 .and. observed_mask(i, j)) then
                            log_likelihood = log_likelihood + scalar_family_loglik(family(j), y(i, j), &
                                liability(i, j), additional(i, j), additional2(i, j))
                        end if
                    end do
                end do
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine heterogeneous_multi_term_parameter_expanded_mcmc

end module mcmcglmm_multiterm_px_sampler
