! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_family_sampler
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state, rng_normal, rng_uniform
    use mcmcglmm_math, only : normal_cdf, normal_logpdf
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_covariance
    use mcmcglmm_distributions, only : inverse_wishart_sample, truncated_normal_sample
    use mcmcglmm_families, only : scalar_family_loglik, two_part_family_loglik, &
        multinomial_log_kernel, ordinal_probit_loglik, ztmb_log_kernel, ztmultinomial_log_kernel
    use mcmcglmm_sampler, only : coefficient_conditional
    use mcmcglmm_engine_features, only : adaptive_mh_decay, adaptive_mh_finalize, adaptive_mh_observe, &
        binary_slice_liability_update
    implicit none
    private

    type, public :: family_mcmc_result
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:)
        real(dp) :: acceptance_rate = 0.0_dp
        real(dp) :: final_proposal_sd = 0.0_dp
    end type family_mcmc_result

    type, public :: threshold_cutpoint_mcmc_result
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: cutpoints(:, :)
        real(dp), allocatable :: last_liability(:)
        real(dp) :: final_proposal_sd = 0.0_dp
        real(dp) :: cutpoint_acceptance_rate = 0.0_dp
    end type threshold_cutpoint_mcmc_result

    type, public :: multivariate_family_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:, :)
        real(dp) :: acceptance_rate = 0.0_dp
    end type multivariate_family_mcmc_result

    type, public :: ordinal_native_mcmc_result
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: cutpoints(:, :)
        real(dp), allocatable :: last_liability(:)
        real(dp) :: liability_acceptance_rate = 0.0_dp
        real(dp) :: cutpoint_acceptance_rate = 0.0_dp
        real(dp) :: final_cutpoint_proposal_sd = 0.0_dp
    end type ordinal_native_mcmc_result

    public :: latent_family_mixed_mcmc
    public :: heterogeneous_family_mixed_mcmc
    public :: two_part_mixed_mcmc
    public :: multinomial_family_mixed_mcmc
    public :: threshold_cutpoint_mixed_mcmc
    public :: ordinal_native_mixed_mcmc

contains

    pure logical function scalar_family_supported(family) result(supported)
        integer, intent(in) :: family !! Native scalar family code being checked for dense latent-MH support.

        select case (family)
        case (2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 22, 23, 24)
            supported = .true.
        case default
            supported = .false.
        end select
    end function scalar_family_supported

    pure subroutine update_scalar_coefficients(liability, x, z, a_inverse, g_value, r_value, beta_prior_mean, &
                                               beta_prior_precision, state, beta, random_effects, info)
        real(dp), intent(in) :: liability(:) !! Current n-vector of latent Gaussian link values.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision matrix.
        real(dp), intent(in) :: g_value !! Positive scalar random-effect variance.
        real(dp), intent(in) :: r_value !! Positive scalar residual/liability variance.
        real(dp), intent(in) :: beta_prior_mean(:) !! p-vector fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p by p fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Gaussian coefficient update.
        real(dp), allocatable, intent(out) :: beta(:) !! Allocated p-vector fixed-effect draw.
        real(dp), allocatable, intent(out) :: random_effects(:) !! Allocated q-vector random-effect draw.
        integer, intent(out) :: info !! Zero on success; nonzero if the shared Gaussian conditional fails.
        real(dp), allocatable :: beta_matrix(:, :)
        real(dp) :: beta_mean_matrix(size(beta_prior_mean), 1)
        real(dp) :: g_matrix(1, 1)
        real(dp) :: r_matrix(1, 1)
        real(dp), allocatable :: random_matrix(:, :)

        beta_mean_matrix(:, 1) = beta_prior_mean
        g_matrix(1, 1) = g_value
        r_matrix(1, 1) = r_value
        call coefficient_conditional(reshape(liability, [size(liability), 1]), x, z, a_inverse, g_matrix, r_matrix, &
            beta_mean_matrix, beta_prior_precision, state, beta_matrix, random_matrix, info)
        if (info /= 0) then
            allocate(beta(0), random_effects(0))
            return
        end if
        beta = beta_matrix(:, 1)
        random_effects = random_matrix(:, 1)
    end subroutine update_scalar_coefficients

    pure subroutine latent_family_mixed_mcmc(family, y, additional, additional2, x, z, a_inverse, beta_prior_mean, &
                                             beta_prior_precision, g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, &
                                             update_r, proposal_scale, iterations, burn, thin, state, result, info, &
                                             slice_sampling, adaptive_mh, slice_limit)
        integer, intent(in) :: family !! Active scalar native family code supported by the dense latent-Gaussian engine.
        real(dp), intent(in) :: y(:) !! Response vector after native-family preprocessing.
        real(dp), intent(in) :: additional(:) !! Family-specific first auxiliary vector; ignored when the family does not use it.
        real(dp), intent(in) :: additional2(:) !! Family-specific second auxiliary vector; ignored when unused.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision, typically pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:) !! p-vector Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p by p Gaussian prior precision for fixed effects.
        real(dp), intent(in) :: g_prior_scale !! Positive scalar inverse-Wishart scale for the random-effect variance.
        real(dp), intent(in) :: g_prior_df !! Positive random-effect inverse-Wishart degrees of freedom.
        real(dp), intent(in) :: r_prior_scale !! Positive scalar inverse-Wishart scale for the latent residual variance.
        real(dp), intent(in) :: r_prior_df !! Positive residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample latent residual variance; otherwise keep its initialized value fixed.
        real(dp), intent(in) :: proposal_scale !! Positive random-walk scale relative to the current latent residual SD.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by latent, coefficient, and variance updates.
        type(family_mcmc_result), intent(out) :: result !! Retained samples, final liabilities, and aggregate MH acceptance rate.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        logical, optional, intent(in) :: slice_sampling !! Use upstream univariate slice updates when family/trials permit.
        logical, optional, intent(in) :: adaptive_mh !! Adapt scalar MH variance during burn-in with upstream recursions.
        real(dp), optional, intent(in) :: slice_limit !! Positive finite liability bound used by slice updates; default 1e30.
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: liability(:)
        real(dp), allocatable :: random_effects(:)
        real(dp) :: current_log_target
        real(dp) :: eta
        real(dp), allocatable :: g_matrix(:, :)
        real(dp) :: g_scale_matrix(1, 1)
        real(dp) :: log_likelihood
        real(dp) :: proposal
        real(dp) :: proposal_log_target
        real(dp), allocatable :: r_matrix(:, :)
        real(dp) :: r_scale_matrix(1, 1)
        real(dp) :: residual
        real(dp) :: sd_value
        real(dp) :: u
        real(dp) :: zdraw
        real(dp) :: accepted_weight
        real(dp) :: adaptive_count
        real(dp) :: adaptive_scale
        real(dp) :: attempted_weight
        real(dp) :: proposal_covariance(1, 1)
        real(dp) :: proposal_mean(1)
        real(dp) :: sample_vector(1)
        real(dp) :: slice_lower
        real(dp) :: slice_log_likelihood
        real(dp) :: slice_upper
        real(dp) :: slice_limit_value
        logical :: accepted_this
        logical :: adaptive_value
        logical :: slice_value
        logical :: use_slice
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: n
        integer :: nsave
        integer :: p
        integer :: proposals
        integer :: q
        integer :: save_index

        info = 0
        n = size(y)
        p = size(x, 2)
        q = size(z, 2)
        if (.not. scalar_family_supported(family) .or. n < 1 .or. p < 1 .or. q < 1 .or. &
            size(additional) /= n .or. size(additional2) /= n .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean) /= p .or. &
            size(beta_prior_precision, 1) /= p .or. size(beta_prior_precision, 2) /= p) then
            info = 1
            return
        end if
        if (g_prior_scale <= 0.0_dp .or. g_prior_df <= 0.0_dp .or. r_prior_scale <= 0.0_dp .or. &
            r_prior_df <= 0.0_dp .or. proposal_scale <= 0.0_dp .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if
        slice_value = .false.
        if (present(slice_sampling)) slice_value = slice_sampling
        adaptive_value = .false.
        if (present(adaptive_mh)) adaptive_value = adaptive_mh
        slice_limit_value = 1.0e30_dp
        if (present(slice_limit)) slice_limit_value = slice_limit
        if (slice_limit_value <= 0.0_dp) then
            info = 2
            return
        end if
        use_slice = slice_value .and. (family == 22 .or. (family == 3 .and. all(nint(additional) == 1)))

        allocate(result%beta(p, nsave), result%random_effects(q, nsave), result%g(nsave), result%r(nsave))
        allocate(result%log_likelihood(nsave), result%last_liability(n), liability(n), beta(p), random_effects(q))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        allocate(g_matrix(1, 1), r_matrix(1, 1))
        g_matrix(1, 1) = g_prior_scale / max(g_prior_df - 2.0_dp, 1.0_dp)
        r_matrix(1, 1) = r_prior_scale / max(r_prior_df - 2.0_dp, 1.0_dp)
        liability = matmul(x, beta)
        do i = 1, n
            if (family == 6) liability(i) = min(max(liability(i), y(i)), additional(i))
        end do
        accepted_weight = 0.0_dp
        attempted_weight = 0.0_dp
        adaptive_count = 0.0_dp
        proposal_mean = 0.0_dp
        proposal_covariance(1, 1) = proposal_scale * proposal_scale * r_matrix(1, 1)
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            sd_value = sqrt(r_matrix(1, 1))
            do i = 1, n
                eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
                if (family == 6) then
                    call truncated_normal_sample(state, eta, sd_value, y(i), additional(i), liability(i), info)
                    if (info /= 0) return
                else if (use_slice) then
                    call binary_slice_liability_update(state, family, y(i), additional(i), eta, sd_value, liability(i), &
                        slice_limit_value, liability(i), slice_log_likelihood, slice_lower, slice_upper, info)
                    if (info /= 0) return
                else
                    call rng_normal(state, zdraw)
                    if (adaptive_value .and. iteration <= burn) then
                        proposal = liability(i) + sqrt(max(proposal_covariance(1, 1), tiny(1.0_dp))) * zdraw
                    else
                        proposal = liability(i) + proposal_scale * sd_value * zdraw
                    end if
                    current_log_target = scalar_family_loglik(family, y(i), liability(i), additional(i), additional2(i)) + &
                        normal_logpdf(liability(i), eta, sd_value)
                    proposal_log_target = scalar_family_loglik(family, y(i), proposal, additional(i), additional2(i)) + &
                        normal_logpdf(proposal, eta, sd_value)
                    call rng_uniform(state, u)
                    proposals = proposals + 1
                    accepted_this = log(u) < proposal_log_target - current_log_target
                    if (accepted_this) then
                        liability(i) = proposal
                        accepted = accepted + 1
                    end if
                    if (adaptive_value .and. iteration <= burn) then
                        call adaptive_mh_decay(accepted_this, accepted_weight, attempted_weight)
                        sample_vector(1) = liability(i)
                        call adaptive_mh_observe(sample_vector, adaptive_count, proposal_mean, proposal_covariance, info)
                        if (info /= 0) return
                    end if
                end if
            end do
            if (adaptive_value .and. .not. use_slice .and. iteration <= burn .and. adaptive_count > 0.0_dp) then
                call adaptive_mh_finalize(proposal_covariance, proposal_mean, adaptive_count, accepted_weight, &
                    attempted_weight, adaptive_scale, info)
                if (info /= 0) return
            end if

            call update_scalar_coefficients(liability, x, z, a_inverse, g_matrix(1, 1), r_matrix(1, 1), &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            g_scale_matrix(1, 1) = g_prior_scale + dot_product(random_effects, matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_matrix, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return

            if (update_r) then
                r_scale_matrix(1, 1) = r_prior_scale
                do i = 1, n
                    eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
                    residual = liability(i) - eta
                    r_scale_matrix(1, 1) = r_scale_matrix(1, 1) + residual * residual
                end do
                call inverse_wishart_sample(state, r_scale_matrix, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, save_index) = beta
                result%random_effects(:, save_index) = random_effects
                result%g(save_index) = g_matrix(1, 1)
                result%r(save_index) = r_matrix(1, 1)
                log_likelihood = 0.0_dp
                do i = 1, n
                    log_likelihood = log_likelihood + &
                        scalar_family_loglik(family, y(i), liability(i), additional(i), additional2(i))
                end do
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
        if (adaptive_value .and. .not. use_slice) then
            result%final_proposal_sd = sqrt(max(proposal_covariance(1, 1), 0.0_dp))
        else
            result%final_proposal_sd = proposal_scale * sqrt(r_matrix(1, 1))
        end if
    end subroutine latent_family_mixed_mcmc

    pure logical function heterogeneous_family_supported(family) result(supported)
        integer, intent(in) :: family !! Native family code being checked for the heterogeneous dense engine.

        select case (family)
        case (1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 22, 23, 24)
            supported = .true.
        case default
            supported = .false.
        end select
    end function heterogeneous_family_supported

    pure subroutine heterogeneous_family_mixed_mcmc(family, y, additional, additional2, x, z, a_inverse, &
                                                      beta_prior_mean, beta_prior_precision, g_prior_scale, &
                                                      g_prior_df, r_prior_scale, r_prior_df, update_r, proposal_scale, &
                                                      iterations, burn, thin, state, result, info, observed)
        integer, intent(in) :: family(:) !! Length-t native family codes for heterogeneous scalar-response traits.
        real(dp), intent(in) :: y(:, :) !! n by t responses after native-family preprocessing.
        real(dp), intent(in) :: additional(:, :) !! n by t first auxiliary values such as trials, bounds, or scales.
        real(dp), intent(in) :: additional2(:, :) !! n by t second auxiliary values such as Student-t degrees of freedom.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design shared across traits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by t fixed-effect Gaussian prior means.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! pt by pt fixed-effect prior precision in trait-major order.
        real(dp), intent(in) :: g_prior_scale(:, :) !! t by t inverse-Wishart scale for random-effect covariance.
        real(dp), intent(in) :: g_prior_df !! Random-effect inverse-Wishart degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! t by t inverse-Wishart scale for residual/liability covariance.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample the full residual/liability covariance matrix.
        real(dp), intent(in) :: proposal_scale !! Positive scalar multiplier for one-coordinate liability random walks.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by liability, coefficient, and covariance draws.
        type(multivariate_family_mcmc_result), intent(out) :: result !! Retained heterogeneous posterior draws and liabilities.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or numerical failure.
        logical, optional, intent(in) :: observed(:, :) !! Optional n by t mask; false entries are treated as missing responses.
        real(dp), allocatable :: beta(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: liability(:, :)
        real(dp) :: log_likelihood
        real(dp) :: proposal
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp), allocatable :: proposal_row(:)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp) :: u
        real(dp) :: zdraw
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: j
        integer :: n
        integer :: nsave
        integer :: p
        integer :: proposals
        integer :: q
        integer :: save_index
        integer :: t
        logical, allocatable :: observed_mask(:, :)

        info = 0
        n = size(y, 1)
        t = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        if (n < 1 .or. t < 1 .or. p < 1 .or. q < 1 .or. size(family) /= t .or. &
            size(additional, 1) /= n .or. size(additional, 2) /= t .or. &
            size(additional2, 1) /= n .or. size(additional2, 2) /= t .or. &
            size(x, 1) /= n .or. size(z, 1) /= n .or. size(a_inverse, 1) /= q .or. &
            size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= t .or. size(beta_prior_precision, 1) /= p * t .or. &
            size(beta_prior_precision, 2) /= p * t .or. size(g_prior_scale, 1) /= t .or. &
            size(g_prior_scale, 2) /= t .or. size(r_prior_scale, 1) /= t .or. size(r_prior_scale, 2) /= t) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed, 1) /= n .or. size(observed, 2) /= t) then
                info = 1
                return
            end if
        end if
        do j = 1, t
            if (.not. heterogeneous_family_supported(family(j))) then
                info = 2
                return
            end if
        end do
        if (g_prior_df <= real(t - 1, dp) .or. r_prior_df <= real(t - 1, dp) .or. proposal_scale <= 0.0_dp .or. &
            iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 3
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, t, nsave), result%random_effects(q, t, nsave))
        allocate(result%g(t, t, nsave), result%r(t, t, nsave), result%log_likelihood(nsave))
        allocate(result%last_liability(n, t), beta(p, t), random_effects(q, t), liability(n, t), proposal_row(t))
        allocate(observed_mask(n, t))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        beta = beta_prior_mean
        random_effects = 0.0_dp
        g_matrix = g_prior_scale / max(g_prior_df - real(t + 1, dp), 1.0_dp)
        r_matrix = r_prior_scale / max(r_prior_df - real(t + 1, dp), 1.0_dp)
        eta = matmul(x, beta) + matmul(z, random_effects)
        liability = eta
        do j = 1, t
            if (family(j) == 1) then
                do i = 1, n
                    if (observed_mask(i, j)) liability(i, j) = y(i, j)
                end do
            else if (family(j) == 6) then
                do i = 1, n
                    if (observed_mask(i, j)) then
                        liability(i, j) = min(max(liability(i, j), y(i, j)), additional(i, j))
                    end if
                end do
            end if
        end do
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            eta = matmul(x, beta) + matmul(z, random_effects)
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

            call coefficient_conditional(liability, x, z, a_inverse, g_matrix, r_matrix, beta_prior_mean, &
                beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            g_scale_post = g_prior_scale + matmul(transpose(random_effects), matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_post, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return
            if (update_r) then
                eta = matmul(x, beta) + matmul(z, random_effects)
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                eta = matmul(x, beta) + matmul(z, random_effects)
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
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine heterogeneous_family_mixed_mcmc

    pure subroutine two_part_mixed_mcmc(family, y, trials, x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
                                        g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, update_r, proposal_scale, &
                                        iterations, burn, thin, state, result, info, observed)
        integer, intent(in) :: family !! Native two-process family code: 11, 15, 18, 19, or 25.
        real(dp), intent(in) :: y(:) !! Original response before MCMCglmm's two-trait zero-process expansion.
        integer, intent(in) :: trials(:) !! Binomial trial counts for families 19/25; positive placeholders for Poisson families.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared by main and zero processes.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design shared by main and zero processes.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by 2 fixed-effect prior means for main and zero processes.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! 2p by 2p fixed-effect prior precision in trait-major order.
        real(dp), intent(in) :: g_prior_scale(:, :) !! 2 by 2 inverse-Wishart scale for random-effect covariance.
        real(dp), intent(in) :: g_prior_df !! Random-effect inverse-Wishart degrees of freedom, required to exceed one.
        real(dp), intent(in) :: r_prior_scale(:, :) !! 2 by 2 inverse-Wishart scale for latent residual covariance.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample the two-process latent residual covariance.
        real(dp), intent(in) :: proposal_scale !! Positive scale multiplying residual covariance for row-wise MH proposals.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by MH, coefficient, and covariance updates.
        type(multivariate_family_mcmc_result), intent(out) :: result !! Retained two-process posterior samples and liabilities.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/settings or numerical failure.
        logical, optional, intent(in) :: observed(:) !! Optional length-n mask; false rows omit the likelihood and are imputed.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: liability(:, :)
        real(dp), allocatable :: proposal(:)
        real(dp), allocatable :: proposal_covariance(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp) :: log_likelihood
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp) :: u
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: n
        integer :: nsave
        integer :: p
        integer :: proposals
        integer :: q
        integer :: save_index
        logical, allocatable :: observed_mask(:)

        info = 0
        n = size(y)
        p = size(x, 2)
        q = size(z, 2)
        if (.not. any(family == [11, 15, 18, 19, 25]) .or. n < 1 .or. p < 1 .or. q < 1 .or. &
            size(trials) /= n .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= 2 .or. &
            size(beta_prior_precision, 1) /= 2 * p .or. size(beta_prior_precision, 2) /= 2 * p .or. &
            size(g_prior_scale, 1) /= 2 .or. size(g_prior_scale, 2) /= 2 .or. &
            size(r_prior_scale, 1) /= 2 .or. size(r_prior_scale, 2) /= 2) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
        end if
        if (g_prior_df <= 1.0_dp .or. r_prior_df <= 1.0_dp .or. proposal_scale <= 0.0_dp .or. &
            iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(result%beta(p, 2, nsave), result%random_effects(q, 2, nsave))
        allocate(result%g(2, 2, nsave), result%r(2, 2, nsave), result%log_likelihood(nsave))
        allocate(result%last_liability(n, 2), liability(n, 2), random_effects(q, 2), beta(p, 2))
        allocate(observed_mask(n))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        beta = beta_prior_mean
        random_effects = 0.0_dp
        g_matrix = g_prior_scale / max(g_prior_df - 3.0_dp, 1.0_dp)
        r_matrix = r_prior_scale / max(r_prior_df - 3.0_dp, 1.0_dp)
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

            call coefficient_conditional(liability, x, z, a_inverse, g_matrix, r_matrix, beta_prior_mean, &
                beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            g_scale_post = g_prior_scale + matmul(transpose(random_effects), matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_post, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return
            if (update_r) then
                eta = matmul(x, beta) + matmul(z, random_effects)
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, save_index) = g_matrix
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
    end subroutine two_part_mixed_mcmc

    pure real(dp) function grouped_loglik(family, response, eta) result(value)
        integer, intent(in) :: family !! Grouped family code: 3 for multinomial, 26 for ztmb, or 27 for ztmultinomial.
        integer, intent(in) :: response(:) !! Counts or binary category indicators for one observational group.
        real(dp), intent(in) :: eta(:) !! Current latent logits for the group's non-reference or Bernoulli dimensions.

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
    end function grouped_loglik

    pure subroutine multinomial_family_mixed_mcmc(family, response, x, z, a_inverse, beta_prior_mean, &
                                                  beta_prior_precision, g_prior_scale, g_prior_df, r_prior_scale, &
                                                  r_prior_df, update_r, proposal_scale, iterations, burn, thin, &
                                                  state, result, info, observed)
        integer, intent(in) :: family !! Grouped native family code: 3, 26, or 27.
        integer, intent(in) :: response(:, :) !! n by category count/indicator matrix; see API coverage for each family layout.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across latent logits.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design shared across latent logits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by n_latent fixed-effect prior mean matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*n_latent square fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :) !! n_latent square inverse-Wishart scale for random effects.
        real(dp), intent(in) :: g_prior_df !! Random-effect inverse-Wishart degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! n_latent square inverse-Wishart scale for latent residuals.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample latent residual covariance; otherwise keep it fixed.
        real(dp), intent(in) :: proposal_scale !! Positive scale for row-wise multivariate random-walk proposals.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by MH, coefficient, and covariance updates.
        type(multivariate_family_mcmc_result), intent(out) :: result !! Retained grouped-family posterior samples and liabilities.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/settings or numerical failure.
        logical, optional, intent(in) :: observed(:) !! Optional length-n mask; false groups are imputed from the latent model.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: eta(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: liability(:, :)
        real(dp), allocatable :: proposal(:)
        real(dp), allocatable :: proposal_covariance(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp) :: current_log_prior
        real(dp) :: current_log_target
        real(dp) :: log_likelihood
        real(dp) :: proposal_log_prior
        real(dp) :: proposal_log_target
        real(dp) :: u
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: n
        integer :: n_latent
        integer :: nsave
        integer :: p
        integer :: proposals
        integer :: q
        integer :: save_index
        logical, allocatable :: observed_mask(:)

        info = 0
        n = size(response, 1)
        p = size(x, 2)
        q = size(z, 2)
        if (family == 26) then
            n_latent = size(response, 2)
        else
            n_latent = size(response, 2) - 1
        end if
        if (.not. any(family == [3, 26, 27]) .or. n < 1 .or. n_latent < 1 .or. p < 1 .or. q < 1 .or. &
            size(x, 1) /= n .or. size(z, 1) /= n .or. size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= n_latent .or. &
            size(beta_prior_precision, 1) /= p * n_latent .or. &
            size(beta_prior_precision, 2) /= p * n_latent .or. &
            size(g_prior_scale, 1) /= n_latent .or. size(g_prior_scale, 2) /= n_latent .or. &
            size(r_prior_scale, 1) /= n_latent .or. size(r_prior_scale, 2) /= n_latent) then
            info = 1
            return
        end if
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
        end if
        if (g_prior_df <= real(n_latent - 1, dp) .or. &
            r_prior_df <= real(n_latent - 1, dp) .or. proposal_scale <= 0.0_dp .or. &
            iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        allocate(observed_mask(n))
        observed_mask = .true.
        if (present(observed)) observed_mask = observed
        do i = 1, n
            if (.not. observed_mask(i)) cycle
            if (any(response(i, :) < 0)) then
                info = 3
                return
            end if
            if (family == 26 .and. any(response(i, :) > 1)) then
                info = 3
                return
            end if
            if (sum(response(i, :)) < 1) then
                info = 3
                return
            end if
        end do
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, n_latent, nsave), result%random_effects(q, n_latent, nsave))
        allocate(result%g(n_latent, n_latent, nsave), result%r(n_latent, n_latent, nsave))
        allocate(result%log_likelihood(nsave), result%last_liability(n, n_latent))
        allocate(beta(p, n_latent), random_effects(q, n_latent), liability(n, n_latent))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        g_matrix = g_prior_scale / max(g_prior_df - real(n_latent + 1, dp), 1.0_dp)
        r_matrix = r_prior_scale / max(r_prior_df - real(n_latent + 1, dp), 1.0_dp)
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
                    current_log_target = current_log_target + grouped_loglik(family, response(i, :), liability(i, :))
                    proposal_log_target = proposal_log_target + grouped_loglik(family, response(i, :), proposal)
                end if
                call rng_uniform(state, u)
                proposals = proposals + 1
                if (log(u) < proposal_log_target - current_log_target) then
                    liability(i, :) = proposal
                    accepted = accepted + 1
                end if
            end do

            call coefficient_conditional(liability, x, z, a_inverse, g_matrix, r_matrix, beta_prior_mean, &
                beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            g_scale_post = g_prior_scale + matmul(transpose(random_effects), matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_post, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return
            if (update_r) then
                eta = matmul(x, beta) + matmul(z, random_effects)
                residual = liability - eta
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                log_likelihood = 0.0_dp
                do i = 1, n
                    if (observed_mask(i)) then
                        log_likelihood = log_likelihood + grouped_loglik(family, response(i, :), liability(i, :))
                    end if
                end do
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
        if (proposals > 0) result%acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine multinomial_family_mixed_mcmc


    pure real(dp) function ordinal_liability_loglik(y_category, cutpoints, liability, observed) result(value)
        integer, intent(in) :: y_category(:) !! One-based ordinal categories; ignored where observed is false.
        real(dp), intent(in) :: cutpoints(:) !! Increasing ordered-probit category boundaries.
        real(dp), intent(in) :: liability(:) !! Current latent link values before the unit-variance probit observation layer.
        logical, optional, intent(in) :: observed(:) !! Optional mask selecting responses that contribute an ordinal likelihood.
        integer :: i

        value = 0.0_dp
        do i = 1, size(y_category)
            if (present(observed)) then
                if (.not. observed(i)) cycle
            end if
            value = value + ordinal_probit_loglik(y_category(i), liability(i), cutpoints)
        end do
    end function ordinal_liability_loglik

    pure subroutine ordinal_native_mixed_mcmc(y_category, initial_cutpoints, cutpoint_proposal_sd, adapt_cutpoints, &
                                               x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
                                               g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, update_r, &
                                               liability_proposal_scale, iterations, burn, thin, state, result, info, observed, &
                                               slice_sampling, slice_limit)
        integer, intent(in) :: y_category(:) !! One-based ordinal categories for each observation.
        real(dp), intent(in) :: initial_cutpoints(:) !! Increasing boundaries; first two and final values remain fixed.
        real(dp), intent(in) :: cutpoint_proposal_sd !! Positive initial SD for ordered cutpoint random walks.
        logical, intent(in) :: adapt_cutpoints !! If true, adapt the cutpoint proposal during burn-in toward 0.44 acceptance.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision, typically pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:) !! p-vector Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p by p Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale !! Positive inverse-Wishart scale for the scalar random-effect variance.
        real(dp), intent(in) :: g_prior_df !! Positive inverse-Wishart degrees of freedom for random-effect variance.
        real(dp), intent(in) :: r_prior_scale !! Positive inverse-Wishart scale for the latent-link residual variance.
        real(dp), intent(in) :: r_prior_df !! Positive residual inverse-Wishart degrees of freedom when update_r is true.
        logical, intent(in) :: update_r !! If true, sample the latent-link residual variance; otherwise keep it fixed.
        real(dp), intent(in) :: liability_proposal_scale !! Positive random-walk scale relative to latent residual SD.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded and eligible for cutpoint adaptation.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by all Metropolis and Gibbs updates.
        type(ordinal_native_mcmc_result), intent(out) :: result !! Retained family-14 posterior draws and diagnostics.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or numerical failure.
        logical, optional, intent(in) :: observed(:) !! Optional mask; false responses are imputed without likelihood contribution.
        logical, optional, intent(in) :: slice_sampling !! Use native binary family-14 slice updates when there are two categories.
        real(dp), optional, intent(in) :: slice_limit !! Positive finite liability guard for binary slice updates.
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: cutpoints(:)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp) :: g_scale_matrix(1, 1)
        real(dp), allocatable :: liability(:)
        logical, allocatable :: observed_mask(:)
        real(dp), allocatable :: proposed_cutpoints(:)
        real(dp), allocatable :: random_effects(:)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp) :: r_scale_matrix(1, 1)
        real(dp) :: acceptance_weight
        real(dp) :: accepted_weight
        real(dp) :: current_cutpoint_loglik
        real(dp) :: current_log_target
        real(dp) :: eta
        real(dp) :: log_hastings
        real(dp) :: proposal
        real(dp) :: proposal_cutpoint_loglik
        real(dp) :: proposal_log_target
        real(dp) :: proposal_sd
        real(dp) :: residual
        real(dp) :: sd_value
        real(dp) :: slice_log_likelihood
        real(dp) :: slice_lower
        real(dp) :: slice_upper
        real(dp) :: slice_limit_value
        logical :: use_slice
        real(dp) :: u
        real(dp) :: zdraw
        integer :: accepted_cutpoints
        integer :: accepted_liabilities
        integer :: i
        integer :: iteration
        integer :: n
        integer :: ncut
        integer :: nsave
        integer :: p
        integer :: proposed_cutpoints_count
        integer :: proposed_liabilities
        integer :: q
        integer :: save_index

        info = 0
        n = size(y_category)
        ncut = size(initial_cutpoints)
        p = size(x, 2)
        q = size(z, 2)
        if (n < 1 .or. ncut < 3 .or. p < 1 .or. q < 1 .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean) /= p .or. &
            size(beta_prior_precision, 1) /= p .or. size(beta_prior_precision, 2) /= p) then
            info = 1
            return
        end if
        allocate(observed_mask(n))
        observed_mask = .true.
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
            observed_mask = observed
        end if
        do i = 1, n
            if (observed_mask(i)) then
                if (y_category(i) < 1 .or. y_category(i) >= ncut) then
                    info = 2
                    return
                end if
            end if
        end do
        if (cutpoint_proposal_sd <= 0.0_dp .or. &
            liability_proposal_scale <= 0.0_dp .or. g_prior_scale <= 0.0_dp .or. g_prior_df <= 0.0_dp .or. &
            r_prior_scale <= 0.0_dp .or. r_prior_df <= 0.0_dp .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        do i = 2, ncut
            if (initial_cutpoints(i) <= initial_cutpoints(i - 1)) then
                info = 3
                return
            end if
        end do
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, nsave), result%random_effects(q, nsave), result%g(nsave), result%r(nsave))
        allocate(result%log_likelihood(nsave), result%cutpoints(ncut, nsave), result%last_liability(n))
        allocate(beta(p), random_effects(q), liability(n), cutpoints(ncut), g_matrix(1, 1), r_matrix(1, 1))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        cutpoints = initial_cutpoints
        g_matrix(1, 1) = g_prior_scale / max(g_prior_df - 2.0_dp, 1.0_dp)
        r_matrix(1, 1) = r_prior_scale / max(r_prior_df - 2.0_dp, 1.0_dp)
        liability = matmul(x, beta)
        proposal_sd = cutpoint_proposal_sd
        acceptance_weight = 1.0_dp
        accepted_weight = 1.0_dp
        accepted_cutpoints = 0
        accepted_liabilities = 0
        proposed_cutpoints_count = 0
        proposed_liabilities = 0
        save_index = 0
        use_slice = .false.
        if (present(slice_sampling)) use_slice = slice_sampling
        use_slice = use_slice .and. ncut == 3
        slice_limit_value = 20.0_dp
        if (present(slice_limit)) slice_limit_value = slice_limit
        if (use_slice .and. slice_limit_value <= 0.0_dp) then
            info = 2
            return
        end if

        do iteration = 1, iterations
            if (ncut > 3) then
                call propose_cutpoints(state, cutpoints, proposal_sd, proposed_cutpoints, log_hastings, info)
                if (info /= 0) return
                current_cutpoint_loglik = ordinal_liability_loglik(y_category, cutpoints, liability, observed_mask)
                proposal_cutpoint_loglik = ordinal_liability_loglik(y_category, proposed_cutpoints, liability, observed_mask)
                call rng_uniform(state, u)
                proposed_cutpoints_count = proposed_cutpoints_count + 1
                acceptance_weight = 0.9_dp * acceptance_weight + 1.0_dp
                accepted_weight = 0.9_dp * accepted_weight
                if (log(u) < proposal_cutpoint_loglik - current_cutpoint_loglik + log_hastings) then
                    cutpoints = proposed_cutpoints
                    accepted_cutpoints = accepted_cutpoints + 1
                    accepted_weight = accepted_weight + 1.0_dp
                end if
                if (adapt_cutpoints .and. iteration <= burn) then
                    proposal_sd = proposal_sd * 2.0_dp ** (accepted_weight / acceptance_weight - 0.44_dp)
                end if
            end if

            sd_value = sqrt(r_matrix(1, 1))
            do i = 1, n
                eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
                if (.not. observed_mask(i)) then
                    call rng_normal(state, zdraw)
                    liability(i) = eta + sd_value * zdraw
                    cycle
                end if
                if (use_slice) then
                    call binary_slice_liability_update(state, 14, real(y_category(i), dp), 1.0_dp, eta, sd_value, &
                        liability(i), slice_limit_value, liability(i), slice_log_likelihood, slice_lower, slice_upper, info)
                    if (info /= 0) return
                    proposed_liabilities = proposed_liabilities + 1
                    accepted_liabilities = accepted_liabilities + 1
                else
                    call rng_normal(state, zdraw)
                    proposal = liability(i) + liability_proposal_scale * sd_value * zdraw
                    current_log_target = ordinal_probit_loglik(y_category(i), liability(i), cutpoints) + &
                        normal_logpdf(liability(i), eta, sd_value)
                    proposal_log_target = ordinal_probit_loglik(y_category(i), proposal, cutpoints) + &
                        normal_logpdf(proposal, eta, sd_value)
                    call rng_uniform(state, u)
                    proposed_liabilities = proposed_liabilities + 1
                    if (log(u) < proposal_log_target - current_log_target) then
                        liability(i) = proposal
                        accepted_liabilities = accepted_liabilities + 1
                    end if
                end if
            end do

            call update_scalar_coefficients(liability, x, z, a_inverse, g_matrix(1, 1), r_matrix(1, 1), &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            g_scale_matrix(1, 1) = g_prior_scale + dot_product(random_effects, matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_matrix, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return

            if (update_r) then
                r_scale_matrix(1, 1) = r_prior_scale
                do i = 1, n
                    eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
                    residual = liability(i) - eta
                    r_scale_matrix(1, 1) = r_scale_matrix(1, 1) + residual * residual
                end do
                call inverse_wishart_sample(state, r_scale_matrix, r_prior_df + real(n, dp), r_matrix, info)
                if (info /= 0) return
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, save_index) = beta
                result%random_effects(:, save_index) = random_effects
                result%g(save_index) = g_matrix(1, 1)
                result%r(save_index) = r_matrix(1, 1)
                result%cutpoints(:, save_index) = cutpoints
                result%log_likelihood(save_index) = ordinal_liability_loglik(y_category, cutpoints, liability, observed_mask)
            end if
        end do
        result%last_liability = liability
        result%final_cutpoint_proposal_sd = proposal_sd
        if (proposed_liabilities > 0) then
            result%liability_acceptance_rate = real(accepted_liabilities, dp) / real(proposed_liabilities, dp)
        end if
        if (proposed_cutpoints_count > 0) then
            result%cutpoint_acceptance_rate = real(accepted_cutpoints, dp) / real(proposed_cutpoints_count, dp)
        end if
    end subroutine ordinal_native_mixed_mcmc

    pure real(dp) function threshold_observed_loglik(y_category, cutpoints, x, z, beta, random_effects, observed) result(value)
        integer, intent(in) :: y_category(:) !! One-based threshold categories; ignored where observed is false.
        real(dp), intent(in) :: cutpoints(:) !! Increasing latent-normal category boundaries.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: beta(:) !! Current p-vector fixed effects.
        real(dp), intent(in) :: random_effects(:) !! Current q-vector random effects.
        logical, optional, intent(in) :: observed(:) !! Optional mask selecting responses that contribute category likelihood.
        real(dp) :: eta
        real(dp) :: probability
        integer :: i

        value = 0.0_dp
        do i = 1, size(y_category)
            if (present(observed)) then
                if (.not. observed(i)) cycle
            end if
            eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
            probability = normal_cdf(cutpoints(y_category(i) + 1) - eta) - &
                normal_cdf(cutpoints(y_category(i)) - eta)
            value = value + log(max(probability, tiny(1.0_dp)))
        end do
    end function threshold_observed_loglik

    pure real(dp) function truncated_normal_normalizer(center, sd_value, lower, upper) result(value)
        real(dp), intent(in) :: center !! Mean of the normal proposal before truncation.
        real(dp), intent(in) :: sd_value !! Positive proposal standard deviation.
        real(dp), intent(in) :: lower !! Lower truncation boundary.
        real(dp), intent(in) :: upper !! Upper truncation boundary.

        value = normal_cdf((upper - center) / sd_value) - normal_cdf((lower - center) / sd_value)
        value = max(value, tiny(1.0_dp))
    end function truncated_normal_normalizer

    pure subroutine propose_cutpoints(state, old_cutpoints, proposal_sd, new_cutpoints, log_hastings, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by sequential truncated-normal proposals.
        real(dp), intent(in) :: old_cutpoints(:) !! Current cutpoints with the first two and final boundaries fixed.
        real(dp), intent(in) :: proposal_sd !! Positive random-walk proposal SD for free cutpoints.
        real(dp), allocatable, intent(out) :: new_cutpoints(:) !! Proposed ordered cutpoint vector.
        real(dp), intent(out) :: log_hastings !! Log reverse-to-forward proposal-density ratio.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid cutpoints or a failed truncated-normal draw.
        real(dp) :: forward_normalizer
        real(dp) :: reverse_normalizer
        integer :: j
        integer :: ncut

        info = 0
        ncut = size(old_cutpoints)
        allocate(new_cutpoints(ncut))
        new_cutpoints = old_cutpoints
        log_hastings = 0.0_dp
        if (proposal_sd <= 0.0_dp .or. ncut < 3) then
            info = 1
            return
        end if
        do j = 2, ncut
            if (old_cutpoints(j) <= old_cutpoints(j - 1)) then
                info = 2
                return
            end if
        end do
        if (ncut <= 3) return

        do j = 3, ncut - 1
            call truncated_normal_sample(state, old_cutpoints(j), proposal_sd, new_cutpoints(j - 1), &
                old_cutpoints(j + 1), new_cutpoints(j), info)
            if (info /= 0) return
        end do
        do j = 3, ncut - 1
            forward_normalizer = truncated_normal_normalizer(old_cutpoints(j), proposal_sd, &
                new_cutpoints(j - 1), old_cutpoints(j + 1))
            reverse_normalizer = truncated_normal_normalizer(new_cutpoints(j), proposal_sd, &
                old_cutpoints(j - 1), new_cutpoints(j + 1))
            log_hastings = log_hastings + log(forward_normalizer) - log(reverse_normalizer)
        end do
    end subroutine propose_cutpoints

    pure subroutine threshold_cutpoint_mixed_mcmc(y_category, initial_cutpoints, cutpoint_proposal_sd, adapt_cutpoints, &
                                                  x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
                                                  g_prior_scale, g_prior_df, iterations, burn, thin, state, result, info, observed)
        integer, intent(in) :: y_category(:) !! One-based threshold categories for each observation.
        real(dp), intent(in) :: initial_cutpoints(:) !! Increasing boundaries; first two and final values remain fixed.
        real(dp), intent(in) :: cutpoint_proposal_sd !! Positive initial SD for MCMCglmm-style cutpoint random walks.
        logical, intent(in) :: adapt_cutpoints !! If true, adapt proposal SD during burn-in toward 0.44 acceptance.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision, typically pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:) !! p-vector Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p by p Gaussian prior precision for fixed effects.
        real(dp), intent(in) :: g_prior_scale !! Positive inverse-Wishart scale for the scalar random-effect variance.
        real(dp), intent(in) :: g_prior_df !! Positive inverse-Wishart degrees of freedom for random-effect variance.
        integer, intent(in) :: iterations !! Total Gibbs/Metropolis iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded and eligible for adaptation.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by cutpoint, liability, coefficient, and variance draws.
        type(threshold_cutpoint_mcmc_result), intent(out) :: result !! Retained threshold posterior draws and cutpoints.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or a numerical sampling failure.
        logical, optional, intent(in) :: observed(:) !! Optional mask; false responses have unconstrained latent liabilities.
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: cutpoints(:)
        real(dp), allocatable :: liability(:)
        logical, allocatable :: observed_mask(:)
        real(dp), allocatable :: proposed_cutpoints(:)
        real(dp), allocatable :: random_effects(:)
        real(dp) :: acceptance_weight
        real(dp) :: accepted_weight
        real(dp) :: current_loglik
        real(dp) :: eta
        real(dp), allocatable :: g_matrix(:, :)
        real(dp) :: g_scale_matrix(1, 1)
        real(dp) :: log_hastings
        real(dp) :: proposal_loglik
        real(dp) :: proposal_sd
        real(dp) :: r_matrix(1, 1)
        real(dp) :: u
        integer :: accepted
        integer :: i
        integer :: iteration
        integer :: n
        integer :: ncut
        integer :: nsave
        integer :: p
        integer :: proposals
        integer :: q
        integer :: save_index

        info = 0
        n = size(y_category)
        ncut = size(initial_cutpoints)
        p = size(x, 2)
        q = size(z, 2)
        if (n < 1 .or. ncut < 3 .or. p < 1 .or. q < 1 .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean) /= p .or. &
            size(beta_prior_precision, 1) /= p .or. size(beta_prior_precision, 2) /= p) then
            info = 1
            return
        end if
        allocate(observed_mask(n))
        observed_mask = .true.
        if (present(observed)) then
            if (size(observed) /= n) then
                info = 1
                return
            end if
            observed_mask = observed
        end if
        do i = 1, n
            if (observed_mask(i)) then
                if (y_category(i) < 1 .or. y_category(i) >= ncut) then
                    info = 2
                    return
                end if
            end if
        end do
        if (cutpoint_proposal_sd <= 0.0_dp .or. &
            g_prior_scale <= 0.0_dp .or. g_prior_df <= 0.0_dp .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        do i = 2, ncut
            if (initial_cutpoints(i) <= initial_cutpoints(i - 1)) then
                info = 3
                return
            end if
        end do
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if

        allocate(result%beta(p, nsave), result%random_effects(q, nsave), result%g(nsave))
        allocate(result%log_likelihood(nsave), result%cutpoints(ncut, nsave), result%last_liability(n))
        allocate(beta(p), random_effects(q), liability(n), cutpoints(ncut), g_matrix(1, 1))
        beta = beta_prior_mean
        random_effects = 0.0_dp
        cutpoints = initial_cutpoints
        g_matrix(1, 1) = g_prior_scale / max(g_prior_df - 2.0_dp, 1.0_dp)
        r_matrix(1, 1) = 1.0_dp
        proposal_sd = cutpoint_proposal_sd
        acceptance_weight = 1.0_dp
        accepted_weight = 1.0_dp
        accepted = 0
        proposals = 0
        save_index = 0

        do iteration = 1, iterations
            if (ncut > 3) then
                call propose_cutpoints(state, cutpoints, proposal_sd, proposed_cutpoints, log_hastings, info)
                if (info /= 0) return
                current_loglik = threshold_observed_loglik(y_category, cutpoints, x, z, beta, random_effects, observed_mask)
                proposal_loglik = threshold_observed_loglik(y_category, proposed_cutpoints, x, z, beta, &
                    random_effects, observed_mask)
                call rng_uniform(state, u)
                proposals = proposals + 1
                acceptance_weight = 0.9_dp * acceptance_weight + 1.0_dp
                accepted_weight = 0.9_dp * accepted_weight
                if (log(u) < proposal_loglik - current_loglik + log_hastings) then
                    cutpoints = proposed_cutpoints
                    accepted = accepted + 1
                    accepted_weight = accepted_weight + 1.0_dp
                end if
                if (adapt_cutpoints .and. iteration <= burn) then
                    proposal_sd = proposal_sd * 2.0_dp ** (accepted_weight / acceptance_weight - 0.44_dp)
                end if
            end if

            do i = 1, n
                eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
                if (observed_mask(i)) then
                    call truncated_normal_sample(state, eta, 1.0_dp, cutpoints(y_category(i)), &
                        cutpoints(y_category(i) + 1), liability(i), info)
                    if (info /= 0) return
                else
                    call rng_normal(state, u)
                    liability(i) = eta + u
                end if
            end do
            call update_scalar_coefficients(liability, x, z, a_inverse, g_matrix(1, 1), r_matrix(1, 1), &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            g_scale_matrix(1, 1) = g_prior_scale + dot_product(random_effects, matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_matrix, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, save_index) = beta
                result%random_effects(:, save_index) = random_effects
                result%g(save_index) = g_matrix(1, 1)
                result%cutpoints(:, save_index) = cutpoints
                result%log_likelihood(save_index) = &
                    threshold_observed_loglik(y_category, cutpoints, x, z, beta, random_effects, observed_mask)
            end if
        end do
        result%last_liability = liability
        result%final_proposal_sd = proposal_sd
        if (proposals > 0) result%cutpoint_acceptance_rate = real(accepted, dp) / real(proposals, dp)
    end subroutine threshold_cutpoint_mixed_mcmc

end module mcmcglmm_family_sampler
