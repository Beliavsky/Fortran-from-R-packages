! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 path/link transformation and mixed-model update blocks; see NOTICE.md and upstream/.
module mcmcglmm_structural_sampler
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_distributions, only : inverse_wishart_sample
    use mcmcglmm_multiterm_sampler, only : multi_term_coefficient_conditional
    use mcmcglmm_engine_features, only : structural_transform, structural_gaussian_loglik, &
        structural_parameter_random_walk_update
    implicit none
    private

    type, public :: structural_gaussian_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: structural_parameter(:, :)
        real(dp), allocatable :: log_likelihood(:)
        real(dp) :: structural_acceptance_rate = 0.0_dp
    end type structural_gaussian_mcmc_result

    public :: structural_gaussian_multi_term_mcmc

contains

    pure subroutine structural_gaussian_multi_term_mcmc(y, x, z, random_term, a_inverse, beta_prior_mean, &
                                                         beta_prior_precision, g_prior_scale, g_prior_df, &
                                                         r_prior_scale, r_prior_df, basis, structural_prior_mean, &
                                                         structural_prior_precision, structural_proposal_sd, &
                                                         iterations, burn, thin, state, result, info, update_g, &
                                                         update_r, initial_g, initial_r, initial_structural)
        real(dp), intent(in) :: y(:, :) !! n by t original Gaussian responses before Lambda structural transformation.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design on the transformed response scale.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design on the transformed response scale.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision spanning all random terms.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by t Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! pt by pt Gaussian prior precision in trait-major order.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! t by t by n_term inverse-Wishart random-effect prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart degrees of freedom for G blocks.
        real(dp), intent(in) :: r_prior_scale(:, :) !! t by t inverse-Wishart residual prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: basis(:, :, :) !! t by t by n_lambda path/SIR basis defining Lambda = I - sum(lambda B).
        real(dp), intent(in) :: structural_prior_mean(:) !! Length-n_lambda Gaussian prior mean for structural coefficients.
        real(dp), intent(in) :: structural_prior_precision(:, :) !! n_lambda square Gaussian prior precision matrix.
        real(dp), intent(in) :: structural_proposal_sd(:) !! Positive componentwise structural random-walk proposal SDs.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by all coefficient, covariance, and MH updates.
        type(structural_gaussian_mcmc_result), intent(out) :: result !! Retained path-model posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or a numerical failure.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term flag; false keeps that G covariance fixed.
        logical, optional, intent(in) :: update_r !! Optional residual-covariance update flag; defaults to true.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional t by t by n_term initial or fixed G covariances.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional t by t initial or fixed residual covariance.
        real(dp), optional, intent(in) :: initial_structural(:) !! Optional starting structural coefficients.
        logical :: accepted
        real(dp), allocatable :: beta(:, :)
        integer :: determinant_sign
        real(dp), allocatable :: g_matrix(:, :, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp) :: log_jacobian
        real(dp) :: log_likelihood
        real(dp), allocatable :: mean_value(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: structural_parameter(:)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        real(dp), allocatable :: transformed(:, :)
        logical, allocatable :: update_g_mask(:)
        logical :: update_r_value
        integer :: accepted_count
        integer :: i
        integer :: iteration
        integer :: n
        integer :: nlambda
        integer :: nsave
        integer :: nterm
        integer :: p
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: save_index
        integer :: t
        integer :: term

        info = 0
        n = size(y, 1)
        t = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        nlambda = size(basis, 3)
        if (n < 1 .or. t < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. nlambda < 1 .or. &
            size(x, 1) /= n .or. size(z, 1) /= n .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= t .or. &
            size(beta_prior_precision, 1) /= p * t .or. size(beta_prior_precision, 2) /= p * t .or. &
            size(g_prior_scale, 1) /= t .or. size(g_prior_scale, 2) /= t .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= t .or. size(r_prior_scale, 2) /= t .or. &
            size(basis, 1) /= t .or. size(basis, 2) /= t .or. size(structural_prior_mean) /= nlambda .or. &
            size(structural_prior_precision, 1) /= nlambda .or. size(structural_prior_precision, 2) /= nlambda .or. &
            size(structural_proposal_sd) /= nlambda .or. any(structural_proposal_sd <= 0.0_dp) .or. &
            any(random_term < 1) .or. any(random_term > nterm)) then
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
        if (present(initial_structural)) then
            if (size(initial_structural) /= nlambda) then
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

        allocate(update_g_mask(nterm), structural_parameter(nlambda), g_matrix(t, t, nterm))
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        update_r_value = .true.
        if (present(update_r)) update_r_value = update_r
        structural_parameter = structural_prior_mean
        if (present(initial_structural)) structural_parameter = initial_structural
        do term = 1, nterm
            if (present(initial_g)) then
                g_matrix(:, :, term) = initial_g(:, :, term)
            else
                g_matrix(:, :, term) = g_prior_scale(:, :, term) / &
                    max(g_prior_df(term) - real(t + 1, dp), 1.0_dp)
            end if
        end do
        if (present(initial_r)) then
            r_matrix = initial_r
        else
            r_matrix = r_prior_scale / max(r_prior_df - real(t + 1, dp), 1.0_dp)
        end if
        beta = beta_prior_mean
        allocate(random_effects(q, t))
        random_effects = 0.0_dp
        allocate(result%beta(p, t, nsave), result%random_effects(q, t, nsave))
        allocate(result%g(t, t, nterm, nsave), result%r(t, t, nsave))
        allocate(result%structural_parameter(nlambda, nsave), result%log_likelihood(nsave))
        accepted_count = 0
        save_index = 0

        do iteration = 1, iterations
            call structural_transform(y, basis, structural_parameter, transformed, log_jacobian, determinant_sign, info)
            if (info /= 0 .or. determinant_sign == 0) then
                if (info == 0) info = 4
                return
            end if
            call multi_term_coefficient_conditional(transformed, x, z, random_term, a_inverse, g_matrix, r_matrix, &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
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

            mean_value = matmul(x, beta) + matmul(z, random_effects)
            if (update_r_value) then
                residual = transformed - mean_value
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if
            call structural_parameter_random_walk_update(state, y, mean_value, r_matrix, basis, &
                structural_prior_mean, structural_prior_precision, structural_proposal_sd, structural_parameter, &
                accepted, info)
            if (info /= 0) return
            if (accepted) accepted_count = accepted_count + 1

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                result%structural_parameter(:, save_index) = structural_parameter
                call structural_gaussian_loglik(y, mean_value, r_matrix, basis, structural_parameter, &
                    log_likelihood, determinant_sign, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%structural_acceptance_rate = real(accepted_count, dp) / real(iterations, dp)
    end subroutine structural_gaussian_multi_term_mcmc

end module mcmcglmm_structural_sampler
