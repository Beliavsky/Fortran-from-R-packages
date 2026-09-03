! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 mixed-family and multi-G conditionals; see NOTICE.md and upstream/.
module mcmcglmm_unified_sampler
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state, rng_normal, rng_uniform
    use mcmcglmm_matrix, only : mvn_log_density
    use mcmcglmm_covariance, only : covariance_update_dispatch
    use mcmcglmm_families, only : scalar_family_loglik
    use mcmcglmm_multiterm_sampler, only : multi_term_coefficient_conditional, &
        multi_term_coefficient_conditional_sparse
    use mcmcglmm_engine_features, only : categorical_measurement_error_update
    use mcmcglmm_sparse, only : mcmcglmm_sparse_matrix, sparse_matmul_matrix, sparse_validate
    use mcmcglmm_sparse_factorization, only : sparse_precision_cache
    implicit none
    private

    type, public :: unified_family_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:, :)
        integer, allocatable :: measurement_category(:)
        real(dp), allocatable :: measurement_probability(:, :)
        real(dp) :: acceptance_rate = 0.0_dp
        integer :: symbolic_analyses = 0
    end type unified_family_mcmc_result

    public :: heterogeneous_multi_term_mixed_mcmc

contains

    pure logical function unified_scalar_family_supported(family) result(supported)
        integer, intent(in) :: family !! Native scalar family code being checked for the unified dense engine.

        select case (family)
        case (1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 22, 23, 24)
            supported = .true.
        case default
            supported = .false.
        end select
    end function unified_scalar_family_supported

    pure subroutine heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, &
                                                          a_inverse, beta_prior_mean, beta_prior_precision, &
                                                          g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, &
                                                          update_r, proposal_scale, iterations, burn, thin, state, &
                                                          result, info, observed, update_g, initial_g, initial_r, &
                                                          g_update_mode, r_update_mode, g_split, r_split, &
                                                          g_fixed_covariance, r_fixed_covariance, &
                                                          g_ante_order, r_ante_order, g_ante_common_beta, &
                                                          r_ante_common_beta, g_ante_common_variance, &
                                                          r_ante_common_variance, measurement_category_effect, &
                                                          measurement_prior_probability, measurement_group, &
                                                          sparse_x, sparse_z)
        integer, intent(in) :: family(:) !! Length-t native family code for each heterogeneous scalar-response trait.
        real(dp), intent(in) :: y(:, :) !! n by t responses after the corresponding native-family preprocessing.
        real(dp), intent(in) :: additional(:, :) !! n by t first family-specific auxiliary values such as trials or bounds.
        real(dp), intent(in) :: additional2(:, :) !! n by t second family-specific auxiliary values such as degrees of freedom.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across traits.
        integer, intent(in) :: random_term(:) !! Length-q one-based covariance-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q row precision spanning all random-effect levels.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by t Gaussian fixed-effect prior mean matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! pt by pt fixed-effect prior precision in trait-major order.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! t by t by n_term inverse-Wishart scales for random effects.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart degrees of freedom for random-effect covariances.
        real(dp), intent(in) :: r_prior_scale(:, :) !! t by t inverse-Wishart scale for latent residual covariance.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample the full latent residual covariance; otherwise keep it fixed.
        real(dp), intent(in) :: proposal_scale !! Positive multiplier for one-coordinate latent random-walk proposals.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by latent, coefficient, and covariance draws.
        type(unified_family_mcmc_result), intent(out) :: result !! Retained mixed-family, multi-G posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        logical, optional, intent(in) :: observed(:, :) !! Optional n by t mask; false responses are latent-model imputations.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term mask; false keeps that G block fixed.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional t by t by n_term initial or fixed G covariances.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional t by t initial or fixed latent residual covariance.
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
        real(dp), optional, intent(in) :: measurement_category_effect(:, :, :) !! n by t by c candidate fixed predictor offsets.
        real(dp), optional, intent(in) :: measurement_prior_probability(:, :) !! g by c prior category probabilities.
        integer, optional, intent(in) :: measurement_group(:) !! Length-n one-based measurement-error group labels.
        type(mcmcglmm_sparse_matrix), optional, intent(in) :: sparse_x !! Optional CSR replacement for dense x.
        type(mcmcglmm_sparse_matrix), optional, intent(in) :: sparse_z !! Optional CSR replacement for dense z.
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
        real(dp), allocatable :: measurement_adjusted_mean(:, :)
        integer, allocatable :: measurement_category(:)
        real(dp), allocatable :: measurement_offset(:, :)
        real(dp), allocatable :: measurement_probability(:, :)
        real(dp) :: log_likelihood
        logical, allocatable :: observed_mask(:, :)
        logical, allocatable :: update_g_mask(:)
        real(dp) :: proposal
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp), allocatable :: proposal_row(:)
        real(dp), allocatable :: random_effects(:, :)
        integer :: random_index
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
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
        integer :: save_index
        integer :: t
        integer :: term
        integer :: r_ante_order_value
        logical :: r_ante_common_beta_value
        logical :: r_ante_common_variance_value
        integer :: r_mode
        integer :: r_split_value
        logical :: covariance_accepted
        logical :: measurement_enabled
        logical :: use_sparse_design
        type(sparse_precision_cache) :: coefficient_factor_cache

        info = 0
        n = size(y, 1)
        t = size(y, 2)
        use_sparse_design = present(sparse_x) .and. present(sparse_z)
        if (present(sparse_x) .neqv. present(sparse_z)) then
            info = 1
            return
        end if
        if (use_sparse_design) then
            call sparse_validate(sparse_x, info)
            if (info /= 0) return
            call sparse_validate(sparse_z, info)
            if (info /= 0 .or. sparse_x%nrow /= n .or. sparse_z%nrow /= n) then
                info = 1
                return
            end if
            p = sparse_x%ncol
            q = sparse_z%ncol
        else
            p = size(x, 2)
            q = size(z, 2)
        end if
        nterm = size(g_prior_scale, 3)
        if (n < 1 .or. t < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. size(family) /= t .or. &
            size(additional, 1) /= n .or. size(additional, 2) /= t .or. &
            size(additional2, 1) /= n .or. size(additional2, 2) /= t .or. &
            size(random_term) /= q .or. size(a_inverse, 1) /= q .or. &
            size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= t .or. &
            size(beta_prior_precision, 1) /= p * t .or. size(beta_prior_precision, 2) /= p * t .or. &
            size(g_prior_scale, 1) /= t .or. size(g_prior_scale, 2) /= t .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= t .or. size(r_prior_scale, 2) /= t .or. any(random_term < 1) .or. &
            any(random_term > nterm)) then
            info = 1
            return
        end if
        if (.not. use_sparse_design) then
            if (size(x, 1) /= n .or. size(z, 1) /= n) then
                info = 1
                return
            end if
        end if
        if (present(observed)) then
            if (size(observed, 1) /= n .or. size(observed, 2) /= t) then
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
            if (size(g_split) /= nterm .or. any(g_split < 0) .or. any(g_split > t)) then
                info = 1
                return
            end if
        end if
        if (present(r_split)) then
            if (r_split < 0 .or. r_split > t) then
                info = 1
                return
            end if
        end if
        if (present(g_fixed_covariance)) then
            if (size(g_fixed_covariance, 1) /= t .or. size(g_fixed_covariance, 2) /= t .or. &
                size(g_fixed_covariance, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(r_fixed_covariance)) then
            if (size(r_fixed_covariance, 1) /= t .or. size(r_fixed_covariance, 2) /= t) then
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
        measurement_enabled = present(measurement_category_effect) .and. &
            present(measurement_prior_probability) .and. present(measurement_group)
        if ((present(measurement_category_effect) .or. present(measurement_prior_probability) .or. &
            present(measurement_group)) .and. .not. measurement_enabled) then
            info = 1
            return
        end if
        if (measurement_enabled) then
            if (size(measurement_category_effect, 1) /= n .or. size(measurement_category_effect, 2) /= t .or. &
                size(measurement_category_effect, 3) /= size(measurement_prior_probability, 2) .or. &
                size(measurement_group) /= n) then
                info = 1
                return
            end if
        end if
        do j = 1, t
            if (.not. unified_scalar_family_supported(family(j))) then
                info = 2
                return
            end if
        end do
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(t - 1, dp)) then
                info = 3
                return
            end if
        end do
        if (r_prior_df <= real(t - 1, dp) .or. proposal_scale <= 0.0_dp .or. iterations <= burn .or. &
            burn < 0 .or. thin < 1) then
            info = 3
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, t, nsave), result%random_effects(q, t, nsave))
        allocate(result%g(t, t, nterm, nsave), result%r(t, t, nsave), result%log_likelihood(nsave))
        allocate(result%last_liability(n, t), beta(p, t), random_effects(q, t), liability(n, t), proposal_row(t))
        allocate(measurement_offset(n, t), source=0.0_dp)
        allocate(observed_mask(n, t), update_g_mask(nterm), g_matrix(t, t, nterm))
        allocate(g_mode(nterm), g_split_value(nterm), g_ante_order_value(nterm))
        allocate(g_ante_common_beta_value(nterm), g_ante_common_variance_value(nterm))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        g_mode = 1
        if (present(g_update_mode)) g_mode = g_update_mode
        where (.not. update_g_mask) g_mode = 0
        g_split_value = t
        if (present(g_split)) g_split_value = g_split
        r_mode = 1
        if (present(r_update_mode)) r_mode = r_update_mode
        if (.not. update_r) r_mode = 0
        r_split_value = t
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
                if (t < 2 .or. g_ante_order_value(term) < 1 .or. g_ante_order_value(term) >= t) then
                    info = 5
                    return
                end if
            end if
        end do
        if (r_mode == 5) then
            if (t < 2 .or. r_ante_order_value < 1 .or. r_ante_order_value >= t) then
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
        beta = beta_prior_mean
        random_effects = 0.0_dp
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
        call unified_linear_predictor(x, z, beta, random_effects, eta, info, sparse_x, sparse_z)
        if (info /= 0) return
        liability = eta
        do j = 1, t
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
            call unified_linear_predictor(x, z, beta, random_effects, eta, info, sparse_x, sparse_z)
            if (info /= 0) return
            if (measurement_enabled) then
                call categorical_measurement_error_update(state, liability, eta, measurement_category_effect, &
                    r_matrix, measurement_prior_probability, measurement_group, measurement_category, &
                    measurement_probability, measurement_adjusted_mean, info)
                if (info /= 0) return
                measurement_offset = measurement_adjusted_mean - eta
                eta = measurement_adjusted_mean
            end if
            do i = 1, n
                do j = 1, t
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

            if (use_sparse_design) then
                call multi_term_coefficient_conditional_sparse(liability - measurement_offset, sparse_x, sparse_z, &
                    random_term, a_inverse, g_matrix, r_matrix, beta_prior_mean, beta_prior_precision, state, beta, &
                    random_effects, info, factor_cache=coefficient_factor_cache)
            else
                call multi_term_coefficient_conditional(liability - measurement_offset, x, z, random_term, &
                    a_inverse, g_matrix, r_matrix, beta_prior_mean, beta_prior_precision, state, beta, random_effects, &
                    info)
            end if
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
                    fixed_block_in = g_fixed_covariance(g_split_value(term) + 1:t, g_split_value(term) + 1:t, term)
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
                call unified_linear_predictor(x, z, beta, random_effects, eta, info, sparse_x, sparse_z)
                if (info /= 0) return
                eta = eta + measurement_offset
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((r_mode == 2 .or. r_mode == 4) .and. present(r_fixed_covariance)) then
                    fixed_block_in = r_fixed_covariance(r_split_value + 1:t, r_split_value + 1:t)
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
                call unified_linear_predictor(x, z, beta, random_effects, eta, info, sparse_x, sparse_z)
                if (info /= 0) return
                eta = eta + measurement_offset
                log_likelihood = 0.0_dp
                do i = 1, n
                    call mvn_log_density(liability(i, :), eta(i, :), r_matrix, current_log_prior, info)
                    if (info /= 0) return
                    log_likelihood = log_likelihood + current_log_prior
                    do j = 1, t
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
        if (measurement_enabled) then
            result%measurement_category = measurement_category
            result%measurement_probability = measurement_probability
        end if
        if (use_sparse_design) result%symbolic_analyses = coefficient_factor_cache%analysis_count
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine heterogeneous_multi_term_mixed_mcmc

    pure subroutine unified_linear_predictor(x, z, beta, random_effects, eta, info, sparse_x, sparse_z)
        !! Form X*beta + Z*u from either dense designs or a matched pair of CSR designs.
        real(dp), intent(in) :: x(:, :) !! Dense fixed-effect design, used when sparse_x is absent.
        real(dp), intent(in) :: z(:, :) !! Dense random-effect design, used when sparse_z is absent.
        real(dp), intent(in) :: beta(:, :) !! Fixed-effect coefficient matrix.
        real(dp), intent(in) :: random_effects(:, :) !! Random-effect coefficient matrix.
        real(dp), allocatable, intent(out) :: eta(:, :) !! Allocated linear predictor.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid design storage or shapes.
        type(mcmcglmm_sparse_matrix), optional, intent(in) :: sparse_x !! Optional CSR fixed-effect design.
        type(mcmcglmm_sparse_matrix), optional, intent(in) :: sparse_z !! Optional CSR random-effect design.
        real(dp), allocatable :: fixed_predictor(:, :)
        real(dp), allocatable :: random_predictor(:, :)

        info = 0
        if (present(sparse_x) .and. present(sparse_z)) then
            call sparse_matmul_matrix(sparse_x, beta, fixed_predictor, info)
            if (info /= 0) return
            call sparse_matmul_matrix(sparse_z, random_effects, random_predictor, info)
            if (info /= 0 .or. any(shape(fixed_predictor) /= shape(random_predictor))) then
                if (info == 0) info = 1
                return
            end if
            eta = fixed_predictor + random_predictor
        else if (.not. present(sparse_x) .and. .not. present(sparse_z)) then
            if (size(x, 2) /= size(beta, 1) .or. size(z, 2) /= size(random_effects, 1) .or. &
                size(x, 1) /= size(z, 1)) then
                info = 1
                return
            end if
            eta = matmul(x, beta) + matmul(z, random_effects)
        else
            info = 1
        end if
    end subroutine unified_linear_predictor

    pure elemental logical function valid_covariance_mode(mode) result(valid)
        integer, intent(in) :: mode !! Candidate MCMCglmm covariance update code.

        valid = any(mode == [0, 1, 2, 3, 4, 5, 6])
    end function valid_covariance_mode

end module mcmcglmm_unified_sampler
