! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 theta_scale and multi-G Gaussian update blocks; see NOTICE.md and upstream/.
module mcmcglmm_theta_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state, rng_normal
    use mcmcglmm_distributions, only : inverse_wishart_sample
    use mcmcglmm_multiterm_sampler, only : multi_term_coefficient_conditional, multi_term_gaussian_loglik
    implicit none
    private

    type, public :: theta_scale_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: theta_scale(:)
        real(dp), allocatable :: log_likelihood(:)
    end type theta_scale_mcmc_result

    public :: theta_scale_multivariate_conditional
    public :: theta_scale_gaussian_mixed_mcmc

contains

    pure subroutine theta_scale_multivariate_conditional(state, current_scale, residual, scale_predictor, &
                                                         residual_covariance, prior_mean, prior_precision, sample, &
                                                         conditional_mean, conditional_variance, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the scalar Gaussian theta draw.
        real(dp), intent(in) :: current_scale !! Current theta value represented in residual and scale_predictor.
        real(dp), intent(in) :: residual(:, :) !! n by t response minus the full predictor at current_scale.
        real(dp), intent(in) :: scale_predictor(:, :) !! n by t contribution multiplied by theta_scale.
        real(dp), intent(in) :: residual_covariance(:, :) !! t by t positive-definite residual covariance.
        real(dp), intent(in) :: prior_mean !! Gaussian prior mean for theta_scale.
        real(dp), intent(in) :: prior_precision !! Nonnegative Gaussian prior precision for theta_scale.
        real(dp), intent(out) :: sample !! Full-conditional Gaussian theta_scale draw.
        real(dp), intent(out) :: conditional_mean !! Full-conditional Gaussian mean.
        real(dp), intent(out) :: conditional_variance !! Full-conditional Gaussian variance.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or covariance inversion failure.
        real(dp) :: information
        real(dp), allocatable :: inverse_r(:, :)
        real(dp) :: numerator
        real(dp) :: z
        integer :: i
        integer :: t

        info = 0
        t = size(residual, 2)
        if (size(residual, 1) < 1 .or. t < 1 .or. size(scale_predictor, 1) /= size(residual, 1) .or. &
            size(scale_predictor, 2) /= t .or. size(residual_covariance, 1) /= t .or. &
            size(residual_covariance, 2) /= t .or. prior_precision < 0.0_dp) then
            sample = current_scale
            conditional_mean = current_scale
            conditional_variance = 0.0_dp
            info = 1
            return
        end if
        call inverse_matrix(residual_covariance, inverse_r, info)
        if (info /= 0) then
            sample = current_scale
            conditional_mean = current_scale
            conditional_variance = 0.0_dp
            return
        end if
        information = prior_precision
        numerator = prior_mean * prior_precision
        do i = 1, size(residual, 1)
            information = information + dot_product(scale_predictor(i, :), matmul(inverse_r, scale_predictor(i, :)))
            numerator = numerator + dot_product(scale_predictor(i, :), &
                matmul(inverse_r, residual(i, :) + current_scale * scale_predictor(i, :)))
        end do
        if (information <= 0.0_dp) then
            sample = current_scale
            conditional_mean = current_scale
            conditional_variance = 0.0_dp
            info = 2
            return
        end if
        conditional_mean = numerator / information
        conditional_variance = 1.0_dp / information
        call rng_normal(state, z)
        sample = conditional_mean + sqrt(conditional_variance) * z
    end subroutine theta_scale_multivariate_conditional

    pure subroutine theta_scale_gaussian_mixed_mcmc(y, x, z, x_scale, z_scale, random_term, a_inverse, &
                                                     beta_prior_mean, beta_prior_precision, g_prior_scale, &
                                                     g_prior_df, r_prior_scale, r_prior_df, theta_prior_mean, &
                                                     theta_prior_precision, iterations, burn, thin, state, result, &
                                                     info, update_g, update_r, initial_g, initial_r, initial_theta)
        real(dp), intent(in) :: y(:, :) !! n by t Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p unscaled fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q unscaled concatenated random-effect design matrix.
        real(dp), intent(in) :: x_scale(:, :) !! n by p subset of fixed design entries multiplied by theta_scale.
        real(dp), intent(in) :: z_scale(:, :) !! n by q subset of random design entries multiplied by theta_scale.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect precision matrix, block structured by random_term.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by t Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*t square Gaussian prior precision for fixed effects.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! t by t by n_term inverse-Wishart G prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart G prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! t by t inverse-Wishart residual prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: theta_prior_mean !! Gaussian prior mean for theta_scale.
        real(dp), intent(in) :: theta_prior_precision !! Nonnegative Gaussian prior precision for theta_scale.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by all Gaussian and covariance updates.
        type(theta_scale_mcmc_result), intent(out) :: result !! Retained theta-scale Gaussian posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term flag indicating G matrices to update.
        logical, optional, intent(in) :: update_r !! Optional residual-covariance update flag; defaults to true.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional t by t by n_term initial or fixed G covariances.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional t by t initial or fixed residual covariance.
        real(dp), optional, intent(in) :: initial_theta !! Optional initial theta_scale value; defaults to one.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: effective_x(:, :)
        real(dp), allocatable :: effective_z(:, :)
        real(dp), allocatable :: g_matrix(:, :, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: scale_predictor(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        logical, allocatable :: update_g_mask(:)
        logical :: update_r_value
        real(dp) :: conditional_mean
        real(dp) :: conditional_variance
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
        integer :: t
        integer :: term
        real(dp) :: theta_value
        real(dp) :: theta_sample

        info = 0
        n = size(y, 1)
        t = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (n < 1 .or. t < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. size(x, 1) /= n .or. &
            size(z, 1) /= n .or. size(x_scale, 1) /= n .or. size(x_scale, 2) /= p .or. &
            size(z_scale, 1) /= n .or. size(z_scale, 2) /= q .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= t .or. size(beta_prior_precision, 1) /= p * t .or. &
            size(beta_prior_precision, 2) /= p * t .or. size(g_prior_scale, 1) /= t .or. &
            size(g_prior_scale, 2) /= t .or. size(g_prior_df) /= nterm .or. size(r_prior_scale, 1) /= t .or. &
            size(r_prior_scale, 2) /= t .or. any(random_term < 1) .or. any(random_term > nterm) .or. &
            theta_prior_precision < 0.0_dp) then
            info = 1
            return
        end if
        if (present(update_g)) then
            if (size(update_g) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_g)) then
            if (size(initial_g, 1) /= t .or. size(initial_g, 2) /= t .or. size(initial_g, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_r)) then
            if (size(initial_r, 1) /= t .or. size(initial_r, 2) /= t) then
                info = 1
                return
            end if
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(t - 1, dp)) then
                info = 2
                return
            end if
        end do
        if (r_prior_df <= real(t - 1, dp) .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(update_g_mask(nterm))
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        update_r_value = .true.
        if (present(update_r)) update_r_value = update_r
        theta_value = 1.0_dp
        if (present(initial_theta)) theta_value = initial_theta

        allocate(result%beta(p, t, nsave), result%random_effects(q, t, nsave))
        allocate(result%g(t, t, nterm, nsave), result%r(t, t, nsave), result%theta_scale(nsave))
        allocate(result%log_likelihood(nsave), g_matrix(t, t, nterm), random_effects(q, t))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        do term = 1, nterm
            if (present(initial_g)) then
                g_matrix(:, :, term) = initial_g(:, :, term)
            else
                g_matrix(:, :, term) = g_prior_scale(:, :, term) / max(g_prior_df(term) - real(t + 1, dp), 1.0_dp)
            end if
        end do
        if (present(initial_r)) then
            r_matrix = initial_r
        else
            r_matrix = r_prior_scale / max(r_prior_df - real(t + 1, dp), 1.0_dp)
        end if
        save_index = 0

        do iteration = 1, iterations
            effective_x = x + (theta_value - 1.0_dp) * x_scale
            effective_z = z + (theta_value - 1.0_dp) * z_scale
            call multi_term_coefficient_conditional(y, effective_x, effective_z, random_term, a_inverse, g_matrix, &
                r_matrix, beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            do term = 1, nterm
                if (.not. update_g_mask(term)) cycle
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
                g_matrix(:, :, term) = sampled_g
            end do

            effective_x = x + (theta_value - 1.0_dp) * x_scale
            effective_z = z + (theta_value - 1.0_dp) * z_scale
            residual = y - matmul(effective_x, beta) - matmul(effective_z, random_effects)
            if (update_r_value) then
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            scale_predictor = matmul(x_scale, beta) + matmul(z_scale, random_effects)
            call theta_scale_multivariate_conditional(state, theta_value, residual, scale_predictor, r_matrix, &
                theta_prior_mean, theta_prior_precision, theta_sample, conditional_mean, conditional_variance, info)
            if (info /= 0) return
            theta_value = theta_sample

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                result%theta_scale(save_index) = theta_value
                effective_x = x + (theta_value - 1.0_dp) * x_scale
                effective_z = z + (theta_value - 1.0_dp) * z_scale
                call multi_term_gaussian_loglik(y, effective_x, effective_z, beta, random_effects, r_matrix, &
                    log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine theta_scale_gaussian_mixed_mcmc

end module mcmcglmm_theta_sampler
