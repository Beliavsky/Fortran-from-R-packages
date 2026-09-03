! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 grouped/two-process and multi-G update blocks; see NOTICE.md and upstream/.
module mcmcglmm_grouped_multiterm_sampler
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state, rng_uniform
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_covariance
    use mcmcglmm_covariance, only : covariance_update_dispatch
    use mcmcglmm_families, only : two_part_family_loglik, multinomial_log_kernel, ztmb_log_kernel, &
        ztmultinomial_log_kernel
    use mcmcglmm_multiterm_sampler, only : multi_term_coefficient_conditional
    implicit none
    private

    type, public :: grouped_multi_term_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:, :)
        real(dp) :: acceptance_rate = 0.0_dp
    end type grouped_multi_term_mcmc_result

    public :: two_part_multi_term_mixed_mcmc
    public :: multinomial_multi_term_mixed_mcmc

contains

    pure real(dp) function grouped_loglik_multi(family, response, eta) result(value)
        integer, intent(in) :: family !! Grouped family code: 3 for multinomial, 26 for ztmb, or 27 for ztmultinomial.
        integer, intent(in) :: response(:) !! Counts or binary category indicators for one observational group.
        real(dp), intent(in) :: eta(:) !! Latent logits for non-reference categories or Bernoulli dimensions.

        select case (family)
        case (3)
            value = multinomial_log_kernel(response, eta)
        case (26)
            value = ztmb_log_kernel(response, eta)
        case (27)
            value = ztmultinomial_log_kernel(response, eta)
        case default
            value = -huge(1.0_dp)
        end select
    end function grouped_loglik_multi

    pure subroutine two_part_multi_term_mixed_mcmc(family, y, trials, x, z, random_term, a_inverse, &
                                                    beta_prior_mean, beta_prior_precision, g_prior_scale, &
                                                    g_prior_df, r_prior_scale, r_prior_df, update_r, &
                                                    proposal_scale, iterations, burn, thin, state, result, info, &
                                                    observed, update_g, initial_g, initial_r, g_update_mode, &
                                                    r_update_mode, g_split, r_split, g_fixed_covariance, &
                                                    r_fixed_covariance, g_ante_order, r_ante_order, &
                                                    g_ante_common_beta, r_ante_common_beta, &
                                                    g_ante_common_variance, r_ante_common_variance)
        integer, intent(in) :: family !! Native two-process family code: 11, 15, 18, 19, or 25.
        real(dp), intent(in) :: y(:) !! Original response before main/zero-process latent expansion.
        integer, intent(in) :: trials(:) !! Binomial trials for families 19/25; positive placeholders otherwise.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across the two latent processes.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across latent processes.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision spanning every random term.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by 2 Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! 2p by 2p Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! 2 by 2 by n_term inverse-Wishart G prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart G prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! 2 by 2 inverse-Wishart residual prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        logical, intent(in) :: update_r !! If true, update the two-process latent residual covariance.
        real(dp), intent(in) :: proposal_scale !! Positive multiplier for row-wise latent Gaussian random walks.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by latent, coefficient, and covariance updates.
        type(grouped_multi_term_mcmc_result), intent(out) :: result !! Retained multi-G two-process posterior draws.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or a numerical failure.
        logical, optional, intent(in) :: observed(:) !! Optional length-n mask; false rows are latent-model imputations.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term mask; false keeps that G block fixed.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional 2 by 2 by n_term initial or fixed G covariances.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional 2 by 2 initial or fixed latent residual covariance.
        integer, optional, intent(in) :: g_update_mode(:) !! Optional per-G covariance mode using update codes zero through six.
        integer, optional, intent(in) :: r_update_mode !! Optional residual covariance mode using update codes zero through six.
        integer, optional, intent(in) :: g_split(:) !! Optional leading dimensions used by G covariance modes 2, 4, and 6.
        integer, optional, intent(in) :: r_split !! Optional leading dimension used by residual covariance modes 2, 4, and 6.
        real(dp), optional, intent(in) :: g_fixed_covariance(:, :, :) !! Full matrices supplying fixed lower-right G blocks.
        real(dp), optional, intent(in) :: r_fixed_covariance(:, :) !! Full matrix supplying a fixed lower-right residual block.
        integer, optional, intent(in) :: g_ante_order(:) !! Optional per-G antedependence lag orders; defaults to one.
        integer, optional, intent(in) :: r_ante_order !! Optional residual antedependence lag order; defaults to one.
        logical, optional, intent(in) :: g_ante_common_beta(:) !! Optional per-G flags sharing one coefficient per ante lag.
        logical, optional, intent(in) :: r_ante_common_beta !! Optional residual flag sharing one coefficient per ante lag.
        logical, optional, intent(in) :: g_ante_common_variance(:) !! Optional per-G flags sharing one ante innovation variance.
        logical, optional, intent(in) :: r_ante_common_variance !! Optional residual flag sharing one ante innovation variance.
        real(dp), allocatable :: beta(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: g_matrix(:, :, :)
        real(dp), allocatable :: fixed_block_in(:, :)
        real(dp), allocatable :: fixed_block_out(:, :)
        integer, allocatable :: g_ante_order_value(:)
        logical, allocatable :: g_ante_common_beta_value(:)
        logical, allocatable :: g_ante_common_variance_value(:)
        integer, allocatable :: g_mode(:)
        integer, allocatable :: g_split_value(:)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: liability(:, :)
        real(dp) :: log_likelihood
        logical, allocatable :: observed_mask(:)
        real(dp), allocatable :: proposal(:)
        real(dp), allocatable :: proposal_covariance(:, :)
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        logical, allocatable :: update_g_mask(:)
        real(dp) :: u
        integer :: accepted
        integer :: i
        integer :: iteration
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
        integer :: r_ante_order_value
        logical :: r_ante_common_beta_value
        logical :: r_ante_common_variance_value
        integer :: r_mode
        integer :: r_split_value
        logical :: covariance_accepted

        info = 0
        n = size(y)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (.not. any(family == [11, 15, 18, 19, 25]) .or. n < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. &
            size(trials) /= n .or. size(x, 1) /= n .or. size(z, 1) /= n .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= 2 .or. size(beta_prior_precision, 1) /= 2 * p .or. &
            size(beta_prior_precision, 2) /= 2 * p .or. size(g_prior_scale, 1) /= 2 .or. &
            size(g_prior_scale, 2) /= 2 .or. size(g_prior_df) /= nterm .or. size(r_prior_scale, 1) /= 2 .or. &
            size(r_prior_scale, 2) /= 2 .or. any(random_term < 1) .or. any(random_term > nterm)) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
        end if
        if (present(update_g)) then
            if (size(update_g) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_g)) then
            if (size(initial_g, 1) /= 2 .or. size(initial_g, 2) /= 2 .or. size(initial_g, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_r)) then
            if (size(initial_r, 1) /= 2 .or. size(initial_r, 2) /= 2) then
                info = 1
                return
            end if
        end if
        if (present(g_update_mode)) then
            if (size(g_update_mode) /= nterm .or. any(.not. valid_covariance_mode(g_update_mode))) then
                info = 1
                return
            end if
        end if
        if (present(r_update_mode)) then
            if (.not. valid_covariance_mode(r_update_mode)) then
                info = 1
                return
            end if
        end if
        if (present(g_split)) then
            if (size(g_split) /= nterm .or. any(g_split < 0) .or. any(g_split > 2)) then
                info = 1
                return
            end if
        end if
        if (present(r_split)) then
            if (r_split < 0 .or. r_split > 2) then
                info = 1
                return
            end if
        end if
        if (present(g_fixed_covariance)) then
            if (size(g_fixed_covariance, 1) /= 2 .or. size(g_fixed_covariance, 2) /= 2 .or. &
                size(g_fixed_covariance, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(r_fixed_covariance)) then
            if (size(r_fixed_covariance, 1) /= 2 .or. size(r_fixed_covariance, 2) /= 2) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_order)) then
            if (size(g_ante_order) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_beta)) then
            if (size(g_ante_common_beta) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_variance)) then
            if (size(g_ante_common_variance) /= nterm) then
                info = 1
                return
            end if
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= 1.0_dp) then
                info = 2
                return
            end if
        end do
        if (r_prior_df <= 1.0_dp .or. proposal_scale <= 0.0_dp .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(observed_mask(n), update_g_mask(nterm), g_matrix(2, 2, nterm))
        allocate(g_mode(nterm), g_split_value(nterm), g_ante_order_value(nterm))
        allocate(g_ante_common_beta_value(nterm), g_ante_common_variance_value(nterm))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        g_mode = 1
        if (present(g_update_mode)) g_mode = g_update_mode
        where (.not. update_g_mask) g_mode = 0
        g_split_value = 2
        if (present(g_split)) g_split_value = g_split
        r_mode = 1
        if (present(r_update_mode)) r_mode = r_update_mode
        if (.not. update_r) r_mode = 0
        r_split_value = 2
        if (present(r_split)) r_split_value = r_split
        g_ante_order_value = 1
        if (present(g_ante_order)) g_ante_order_value = g_ante_order
        g_ante_common_beta_value = .false.
        if (present(g_ante_common_beta)) g_ante_common_beta_value = g_ante_common_beta
        g_ante_common_variance_value = .false.
        if (present(g_ante_common_variance)) g_ante_common_variance_value = g_ante_common_variance
        r_ante_order_value = 1
        if (present(r_ante_order)) r_ante_order_value = r_ante_order
        r_ante_common_beta_value = .false.
        if (present(r_ante_common_beta)) r_ante_common_beta_value = r_ante_common_beta
        r_ante_common_variance_value = .false.
        if (present(r_ante_common_variance)) r_ante_common_variance_value = r_ante_common_variance
        do term = 1, nterm
            if (g_mode(term) == 5) then
                if (g_ante_order_value(term) /= 1) then
                    info = 5
                    return
                end if
            end if
        end do
        if (r_mode == 5 .and. r_ante_order_value /= 1) then
            info = 5
            return
        end if
        do term = 1, nterm
            if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. .not. present(g_fixed_covariance)) then
                info = 1
                return
            end if
        end do
        if ((r_mode == 2 .or. r_mode == 4) .and. .not. present(r_fixed_covariance)) then
            info = 1
            return
        end if
        allocate(result%beta(p, 2, nsave), result%random_effects(q, 2, nsave))
        allocate(result%g(2, 2, nterm, nsave), result%r(2, 2, nsave))
        allocate(result%log_likelihood(nsave), result%last_liability(n, 2))
        beta = beta_prior_mean
        allocate(random_effects(q, 2), liability(n, 2))
        random_effects = 0.0_dp
        do term = 1, nterm
            if (present(initial_g)) then
                g_matrix(:, :, term) = initial_g(:, :, term)
            else
                g_matrix(:, :, term) = g_prior_scale(:, :, term) / max(g_prior_df(term) - 3.0_dp, 1.0_dp)
            end if
        end do
        if (present(initial_r)) then
            r_matrix = initial_r
        else
            r_matrix = r_prior_scale / max(r_prior_df - 3.0_dp, 1.0_dp)
        end if
        eta = matmul(x, beta) + matmul(z, random_effects)
        liability = eta
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            eta = matmul(x, beta) + matmul(z, random_effects)
            proposal_covariance = proposal_scale * proposal_scale * r_matrix
            do i = 1, n
                call sample_mvn_covariance(state, liability(i, :), proposal_covariance, proposal, info)
                if (info /= 0) return
                call mvn_log_density(liability(i, :), eta(i, :), r_matrix, current_log_prior, info)
                if (info /= 0) return
                call mvn_log_density(proposal, eta(i, :), r_matrix, proposal_log_prior, info)
                if (info /= 0) return
                current_log_target = current_log_prior
                proposal_log_target = proposal_log_prior
                if (observed_mask(i)) then
                    current_log_target = current_log_target + &
                        two_part_family_loglik(family, y(i), liability(i, 1), liability(i, 2), trials(i))
                    proposal_log_target = proposal_log_target + &
                        two_part_family_loglik(family, y(i), proposal(1), proposal(2), trials(i))
                end if
                call rng_uniform(state, u)
                proposals = proposals + 1
                if (log(u) < proposal_log_target - current_log_target) then
                    liability(i, :) = proposal
                    accepted = accepted + 1
                end if
            end do

            call multi_term_coefficient_conditional(liability, x, z, random_term, a_inverse, g_matrix, r_matrix, &
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
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. present(g_fixed_covariance)) then
                    fixed_block_in = g_fixed_covariance(g_split_value(term) + 1:2, g_split_value(term) + 1:2, term)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (g_mode(term) == 5) then
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, covariance_accepted, info, &
                        ante_location=term_effects, ante_lag_order=g_ante_order_value(term), &
                        ante_common_beta=g_ante_common_beta_value(term), &
                        ante_common_variance=g_ante_common_variance_value(term), ante_a_inverse=term_a_inverse)
                else
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, covariance_accepted, info)
                end if
                deallocate(term_indices)
                if (info /= 0) return
                g_matrix(:, :, term) = sampled_g
            end do
            if (r_mode /= 0) then
                eta = matmul(x, beta) + matmul(z, random_effects)
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((r_mode == 2 .or. r_mode == 4) .and. present(r_fixed_covariance)) then
                    fixed_block_in = r_fixed_covariance(r_split_value + 1:2, r_split_value + 1:2)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (r_mode == 5) then
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, &
                        covariance_accepted, info, ante_location=residual, ante_lag_order=r_ante_order_value, &
                        ante_common_beta=r_ante_common_beta_value, ante_common_variance=r_ante_common_variance_value)
                else
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, &
                        covariance_accepted, info)
                end if
                if (info /= 0) return
                r_matrix = sampled_g
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                log_likelihood = 0.0_dp
                do i = 1, n
                    if (observed_mask(i)) then
                        log_likelihood = log_likelihood + &
                            two_part_family_loglik(family, y(i), liability(i, 1), liability(i, 2), trials(i))
                    end if
                end do
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine two_part_multi_term_mixed_mcmc

    pure subroutine multinomial_multi_term_mixed_mcmc(family, response, x, z, random_term, a_inverse, &
                                                       beta_prior_mean, beta_prior_precision, g_prior_scale, &
                                                       g_prior_df, r_prior_scale, r_prior_df, update_r, &
                                                       proposal_scale, iterations, burn, thin, state, result, info, &
                                                       observed, update_g, initial_g, initial_r, g_update_mode, &
                                                       r_update_mode, g_split, r_split, g_fixed_covariance, &
                                                       r_fixed_covariance, g_ante_order, r_ante_order, &
                                                       g_ante_common_beta, r_ante_common_beta, &
                                                       g_ante_common_variance, r_ante_common_variance)
        integer, intent(in) :: family !! Grouped native family code: 3, 26, or 27.
        integer, intent(in) :: response(:, :) !! n by category count/indicator matrix in the native grouped layout.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across latent logits.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across latent logits.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision spanning every random term.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by n_latent Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*n_latent square Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! n_latent square by n_term inverse-Wishart G prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart G prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! n_latent square inverse-Wishart latent residual prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        logical, intent(in) :: update_r !! If true, update the grouped latent residual covariance.
        real(dp), intent(in) :: proposal_scale !! Positive multiplier for row-wise multivariate latent random walks.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by latent, coefficient, and covariance updates.
        type(grouped_multi_term_mcmc_result), intent(out) :: result !! Retained multi-G grouped-family posterior draws.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or a numerical failure.
        logical, optional, intent(in) :: observed(:) !! Optional length-n mask; false groups omit likelihood and are imputed.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term mask; false keeps that G block fixed.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional latent by latent by n_term initial/fixed G blocks.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional latent by latent initial/fixed residual covariance.
        integer, optional, intent(in) :: g_update_mode(:) !! Optional per-G covariance mode using update codes zero through six.
        integer, optional, intent(in) :: r_update_mode !! Optional residual covariance mode using update codes zero through six.
        integer, optional, intent(in) :: g_split(:) !! Optional leading dimensions used by G covariance modes 2, 4, and 6.
        integer, optional, intent(in) :: r_split !! Optional leading dimension used by residual covariance modes 2, 4, and 6.
        real(dp), optional, intent(in) :: g_fixed_covariance(:, :, :) !! Full matrices supplying fixed lower-right G blocks.
        real(dp), optional, intent(in) :: r_fixed_covariance(:, :) !! Full matrix supplying a fixed lower-right residual block.
        integer, optional, intent(in) :: g_ante_order(:) !! Optional per-G antedependence lag orders; defaults to one.
        integer, optional, intent(in) :: r_ante_order !! Optional residual antedependence lag order; defaults to one.
        logical, optional, intent(in) :: g_ante_common_beta(:) !! Optional per-G flags sharing one coefficient per ante lag.
        logical, optional, intent(in) :: r_ante_common_beta !! Optional residual flag sharing one coefficient per ante lag.
        logical, optional, intent(in) :: g_ante_common_variance(:) !! Optional per-G flags sharing one ante innovation variance.
        logical, optional, intent(in) :: r_ante_common_variance !! Optional residual flag sharing one ante innovation variance.
        real(dp), allocatable :: beta(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: g_matrix(:, :, :)
        real(dp), allocatable :: fixed_block_in(:, :)
        real(dp), allocatable :: fixed_block_out(:, :)
        integer, allocatable :: g_ante_order_value(:)
        logical, allocatable :: g_ante_common_beta_value(:)
        logical, allocatable :: g_ante_common_variance_value(:)
        integer, allocatable :: g_mode(:)
        integer, allocatable :: g_split_value(:)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: liability(:, :)
        real(dp) :: log_likelihood
        logical, allocatable :: observed_mask(:)
        real(dp), allocatable :: proposal(:)
        real(dp), allocatable :: proposal_covariance(:, :)
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        logical, allocatable :: update_g_mask(:)
        real(dp) :: u
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: n
        integer :: n_latent
        integer :: nsave
        integer :: nterm
        integer :: p
        integer :: proposals
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: save_index
        integer :: term
        integer :: r_ante_order_value
        logical :: r_ante_common_beta_value
        logical :: r_ante_common_variance_value
        integer :: r_mode
        integer :: r_split_value
        logical :: covariance_accepted

        info = 0
        n = size(response, 1)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (family == 26) then
            n_latent = size(response, 2)
        else
            n_latent = size(response, 2) - 1
        end if
        if (.not. any(family == [3, 26, 27]) .or. n < 1 .or. n_latent < 1 .or. p < 1 .or. q < 1 .or. &
            nterm < 1 .or. size(x, 1) /= n .or. size(z, 1) /= n .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= n_latent .or. size(beta_prior_precision, 1) /= p * n_latent .or. &
            size(beta_prior_precision, 2) /= p * n_latent .or. size(g_prior_scale, 1) /= n_latent .or. &
            size(g_prior_scale, 2) /= n_latent .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= n_latent .or. size(r_prior_scale, 2) /= n_latent .or. &
            any(random_term < 1) .or. any(random_term > nterm)) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
        end if
        if (present(update_g)) then
            if (size(update_g) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_g)) then
            if (size(initial_g, 1) /= n_latent .or. size(initial_g, 2) /= n_latent .or. &
                size(initial_g, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_r)) then
            if (size(initial_r, 1) /= n_latent .or. size(initial_r, 2) /= n_latent) then
                info = 1
                return
            end if
        end if
        if (present(g_update_mode)) then
            if (size(g_update_mode) /= nterm .or. any(.not. valid_covariance_mode(g_update_mode))) then
                info = 1
                return
            end if
        end if
        if (present(r_update_mode)) then
            if (.not. valid_covariance_mode(r_update_mode)) then
                info = 1
                return
            end if
        end if
        if (present(g_split)) then
            if (size(g_split) /= nterm .or. any(g_split < 0) .or. any(g_split > n_latent)) then
                info = 1
                return
            end if
        end if
        if (present(r_split)) then
            if (r_split < 0 .or. r_split > n_latent) then
                info = 1
                return
            end if
        end if
        if (present(g_fixed_covariance)) then
            if (size(g_fixed_covariance, 1) /= n_latent .or. size(g_fixed_covariance, 2) /= n_latent .or. &
                size(g_fixed_covariance, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(r_fixed_covariance)) then
            if (size(r_fixed_covariance, 1) /= n_latent .or. size(r_fixed_covariance, 2) /= n_latent) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_order)) then
            if (size(g_ante_order) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_beta)) then
            if (size(g_ante_common_beta) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_variance)) then
            if (size(g_ante_common_variance) /= nterm) then
                info = 1
                return
            end if
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(n_latent - 1, dp)) then
                info = 2
                return
            end if
        end do
        if (r_prior_df <= real(n_latent - 1, dp) .or. proposal_scale <= 0.0_dp .or. &
            iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        allocate(observed_mask(n))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        do i = 1, n
            if (.not. observed_mask(i)) cycle
            if (any(response(i, :) < 0) .or. sum(response(i, :)) < 1) then
                info = 3
                return
            end if
            if (family == 26 .and. any(response(i, :) > 1)) then
                info = 3
                return
            end if
        end do
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(update_g_mask(nterm), g_matrix(n_latent, n_latent, nterm))
        allocate(g_mode(nterm), g_split_value(nterm), g_ante_order_value(nterm))
        allocate(g_ante_common_beta_value(nterm), g_ante_common_variance_value(nterm))
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        g_mode = 1
        if (present(g_update_mode)) g_mode = g_update_mode
        where (.not. update_g_mask) g_mode = 0
        g_split_value = n_latent
        if (present(g_split)) g_split_value = g_split
        r_mode = 1
        if (present(r_update_mode)) r_mode = r_update_mode
        if (.not. update_r) r_mode = 0
        r_split_value = n_latent
        if (present(r_split)) r_split_value = r_split
        g_ante_order_value = 1
        if (present(g_ante_order)) g_ante_order_value = g_ante_order
        g_ante_common_beta_value = .false.
        if (present(g_ante_common_beta)) g_ante_common_beta_value = g_ante_common_beta
        g_ante_common_variance_value = .false.
        if (present(g_ante_common_variance)) g_ante_common_variance_value = g_ante_common_variance
        r_ante_order_value = 1
        if (present(r_ante_order)) r_ante_order_value = r_ante_order
        r_ante_common_beta_value = .false.
        if (present(r_ante_common_beta)) r_ante_common_beta_value = r_ante_common_beta
        r_ante_common_variance_value = .false.
        if (present(r_ante_common_variance)) r_ante_common_variance_value = r_ante_common_variance
        do term = 1, nterm
            if (g_mode(term) == 5) then
                if (n_latent < 2 .or. g_ante_order_value(term) < 1 .or. g_ante_order_value(term) >= n_latent) then
                    info = 5
                    return
                end if
            end if
        end do
        if (r_mode == 5) then
            if (n_latent < 2 .or. r_ante_order_value < 1 .or. r_ante_order_value >= n_latent) then
                info = 5
                return
            end if
        end if
        do term = 1, nterm
            if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. .not. present(g_fixed_covariance)) then
                info = 1
                return
            end if
        end do
        if ((r_mode == 2 .or. r_mode == 4) .and. .not. present(r_fixed_covariance)) then
            info = 1
            return
        end if
        allocate(result%beta(p, n_latent, nsave), result%random_effects(q, n_latent, nsave))
        allocate(result%g(n_latent, n_latent, nterm, nsave), result%r(n_latent, n_latent, nsave))
        allocate(result%log_likelihood(nsave), result%last_liability(n, n_latent))
        beta = beta_prior_mean
        allocate(random_effects(q, n_latent), liability(n, n_latent))
        random_effects = 0.0_dp
        do term = 1, nterm
            if (present(initial_g)) then
                g_matrix(:, :, term) = initial_g(:, :, term)
            else
                g_matrix(:, :, term) = g_prior_scale(:, :, term) / &
                    max(g_prior_df(term) - real(n_latent + 1, dp), 1.0_dp)
            end if
        end do
        if (present(initial_r)) then
            r_matrix = initial_r
        else
            r_matrix = r_prior_scale / max(r_prior_df - real(n_latent + 1, dp), 1.0_dp)
        end if
        eta = matmul(x, beta) + matmul(z, random_effects)
        liability = eta
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            eta = matmul(x, beta) + matmul(z, random_effects)
            proposal_covariance = proposal_scale * proposal_scale * r_matrix
            do i = 1, n
                call sample_mvn_covariance(state, liability(i, :), proposal_covariance, proposal, info)
                if (info /= 0) return
                call mvn_log_density(liability(i, :), eta(i, :), r_matrix, current_log_prior, info)
                if (info /= 0) return
                call mvn_log_density(proposal, eta(i, :), r_matrix, proposal_log_prior, info)
                if (info /= 0) return
                current_log_target = current_log_prior
                proposal_log_target = proposal_log_prior
                if (observed_mask(i)) then
                    current_log_target = current_log_target + grouped_loglik_multi(family, response(i, :), liability(i, :))
                    proposal_log_target = proposal_log_target + grouped_loglik_multi(family, response(i, :), proposal)
                end if
                call rng_uniform(state, u)
                proposals = proposals + 1
                if (log(u) < proposal_log_target - current_log_target) then
                    liability(i, :) = proposal
                    accepted = accepted + 1
                end if
            end do

            call multi_term_coefficient_conditional(liability, x, z, random_term, a_inverse, g_matrix, r_matrix, &
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
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. present(g_fixed_covariance)) then
                    fixed_block_in = g_fixed_covariance(g_split_value(term) + 1:n_latent, &
                        g_split_value(term) + 1:n_latent, term)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (g_mode(term) == 5) then
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, covariance_accepted, info, &
                        ante_location=term_effects, ante_lag_order=g_ante_order_value(term), &
                        ante_common_beta=g_ante_common_beta_value(term), &
                        ante_common_variance=g_ante_common_variance_value(term), ante_a_inverse=term_a_inverse)
                else
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, covariance_accepted, info)
                end if
                deallocate(term_indices)
                if (info /= 0) return
                g_matrix(:, :, term) = sampled_g
            end do
            if (r_mode /= 0) then
                eta = matmul(x, beta) + matmul(z, random_effects)
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((r_mode == 2 .or. r_mode == 4) .and. present(r_fixed_covariance)) then
                    fixed_block_in = r_fixed_covariance(r_split_value + 1:n_latent, r_split_value + 1:n_latent)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (r_mode == 5) then
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, &
                        covariance_accepted, info, ante_location=residual, ante_lag_order=r_ante_order_value, &
                        ante_common_beta=r_ante_common_beta_value, ante_common_variance=r_ante_common_variance_value)
                else
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, &
                        covariance_accepted, info)
                end if
                if (info /= 0) return
                r_matrix = sampled_g
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                log_likelihood = 0.0_dp
                do i = 1, n
                    if (observed_mask(i)) then
                        log_likelihood = log_likelihood + grouped_loglik_multi(family, response(i, :), liability(i, :))
                    end if
                end do
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine multinomial_multi_term_mixed_mcmc

    pure elemental logical function valid_covariance_mode(mode) result(valid)
        integer, intent(in) :: mode !! Candidate MCMCglmm covariance update code.

        valid = any(mode == [0, 1, 2, 3, 4, 5, 6])
    end function valid_covariance_mode

end module mcmcglmm_grouped_multiterm_sampler
