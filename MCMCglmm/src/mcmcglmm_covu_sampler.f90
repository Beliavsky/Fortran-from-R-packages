! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 covu joint G-R updates; see NOTICE.md and upstream/.
module mcmcglmm_covu_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_precision
    use mcmcglmm_joint_gr, only : joint_gr_covariance_update, joint_gr_decompose
    implicit none
    private

    type, public :: covu_gaussian_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: joint_covariance(:, :, :)
        real(dp), allocatable :: g(:, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: regression(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
    end type covu_gaussian_mcmc_result

    public :: covu_coefficient_conditional
    public :: covu_gaussian_loglik
    public :: covu_gaussian_mixed_mcmc

contains

    pure integer function covu_beta_index(trait, effect, fixed_effects) result(index_value)
        integer, intent(in) :: trait !! One-based response-trait index.
        integer, intent(in) :: effect !! One-based fixed-effect index within the response trait.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect coefficients for each response trait.

        index_value = (trait - 1) * fixed_effects + effect
    end function covu_beta_index

    pure integer function covu_random_index(level, coordinate, fixed_size, random_coordinates) result(index_value)
        integer, intent(in) :: level !! One-based iid covu level, corresponding to one residual level.
        integer, intent(in) :: coordinate !! One-based covu random-effect coordinate within the level.
        integer, intent(in) :: fixed_size !! Total number of packed fixed-effect coefficients.
        integer, intent(in) :: random_coordinates !! Number of covu coordinates at each level.

        index_value = fixed_size + (level - 1) * random_coordinates + coordinate
    end function covu_random_index

    pure subroutine covu_coefficient_conditional(y, x, random_loading, g_matrix, r_matrix, regression, &
                                                 beta_prior_mean, beta_prior_precision, state, beta, &
                                                 random_effects, info)
        real(dp), intent(in) :: y(:, :) !! n by r Gaussian responses, with one row for each covu/residual level.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix shared across response traits.
        real(dp), intent(in) :: random_loading(:, :) !! r by g direct loading of each covu random coordinate on each response.
        real(dp), intent(in) :: g_matrix(:, :) !! g by g marginal covariance of the covu random-effect coordinates.
        real(dp), intent(in) :: r_matrix(:, :) !! r by r residual covariance conditional on the covu random effects.
        real(dp), intent(in) :: regression(:, :) !! r by g residual-on-random regression from the joint covariance.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by r Gaussian fixed-effect prior mean matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! pr by pr fixed-effect prior precision in trait-major order.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the joint coefficient Gaussian draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated p by r sampled fixed-effect coefficient matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated n by g sampled iid covu random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes or a linear-algebra failure.
        real(dp), allocatable :: beta_mean_vector(:)
        real(dp), allocatable :: beta_rhs(:)
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: effective_loading(:, :)
        real(dp), allocatable :: g_inverse(:, :)
        real(dp), allocatable :: h(:, :)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: weighted_h(:, :)
        real(dp), allocatable :: weighted_y(:)
        integer :: a
        integer :: b
        integer :: fixed_size
        integer :: g
        integer :: i
        integer :: j
        integer :: k
        integer :: m
        integer :: n
        integer :: p
        integer :: r

        info = 0
        n = size(y, 1)
        r = size(y, 2)
        p = size(x, 2)
        g = size(random_loading, 2)
        fixed_size = p * r
        m = fixed_size + n * g
        if (n < 1 .or. r < 1 .or. p < 1 .or. g < 1 .or. size(x, 1) /= n .or. &
            size(random_loading, 1) /= r .or. size(g_matrix, 1) /= g .or. size(g_matrix, 2) /= g .or. &
            size(r_matrix, 1) /= r .or. size(r_matrix, 2) /= r .or. size(regression, 1) /= r .or. &
            size(regression, 2) /= g .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= r .or. size(beta_prior_precision, 1) /= fixed_size .or. &
            size(beta_prior_precision, 2) /= fixed_size) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if

        call inverse_matrix(g_matrix, g_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        effective_loading = random_loading + regression
        allocate(precision(m, m), rhs(m), beta_mean_vector(fixed_size))
        precision = 0.0_dp
        rhs = 0.0_dp
        do a = 1, r
            do j = 1, p
                beta_mean_vector(covu_beta_index(a, j, p)) = beta_prior_mean(j, a)
            end do
        end do
        beta_rhs = matmul(beta_prior_precision, beta_mean_vector)
        precision(1:fixed_size, 1:fixed_size) = beta_prior_precision
        rhs(1:fixed_size) = beta_rhs

        do i = 1, n
            allocate(h(r, m))
            h = 0.0_dp
            do a = 1, r
                do j = 1, p
                    h(a, covu_beta_index(a, j, p)) = x(i, j)
                end do
                do k = 1, g
                    h(a, covu_random_index(i, k, fixed_size, g)) = effective_loading(a, k)
                end do
            end do
            weighted_h = matmul(r_inverse, h)
            weighted_y = matmul(r_inverse, y(i, :))
            precision = precision + matmul(transpose(h), weighted_h)
            rhs = rhs + matmul(transpose(h), weighted_y)
            deallocate(h)
            do a = 1, g
                do b = 1, g
                    precision(covu_random_index(i, a, fixed_size, g), &
                        covu_random_index(i, b, fixed_size, g)) = &
                        precision(covu_random_index(i, a, fixed_size, g), &
                        covu_random_index(i, b, fixed_size, g)) + g_inverse(a, b)
                end do
            end do
        end do

        call sample_mvn_precision(state, rhs, precision, packed, conditional_mean, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        allocate(beta(p, r), random_effects(n, g))
        do a = 1, r
            do j = 1, p
                beta(j, a) = packed(covu_beta_index(a, j, p))
            end do
        end do
        do i = 1, n
            do k = 1, g
                random_effects(i, k) = packed(covu_random_index(i, k, fixed_size, g))
            end do
        end do
    end subroutine covu_coefficient_conditional

    pure subroutine covu_gaussian_loglik(y, x, random_loading, beta, random_effects, r_matrix, regression, &
                                         log_likelihood, info)
        real(dp), intent(in) :: y(:, :) !! n by r Gaussian responses, one row for each covu/residual level.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: random_loading(:, :) !! r by g direct random-effect loading matrix.
        real(dp), intent(in) :: beta(:, :) !! p by r fixed-effect coefficient matrix.
        real(dp), intent(in) :: random_effects(:, :) !! n by g sampled covu random effects.
        real(dp), intent(in) :: r_matrix(:, :) !! r by r residual covariance conditional on random effects.
        real(dp), intent(in) :: regression(:, :) !! r by g residual-on-random regression beta_rr.
        real(dp), intent(out) :: log_likelihood !! Conditional Gaussian log likelihood summed over residual levels.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes or covariance failure.
        real(dp), allocatable :: effective_loading(:, :)
        real(dp), allocatable :: mean_row(:)
        real(dp) :: row_log_density
        integer :: i
        integer :: n
        integer :: r

        info = 0
        n = size(y, 1)
        r = size(y, 2)
        if (size(x, 1) /= n .or. size(beta, 1) /= size(x, 2) .or. size(beta, 2) /= r .or. &
            size(random_effects, 1) /= n .or. size(random_effects, 2) /= size(random_loading, 2) .or. &
            size(random_loading, 1) /= r .or. size(r_matrix, 1) /= r .or. size(r_matrix, 2) /= r .or. &
            size(regression, 1) /= r .or. size(regression, 2) /= size(random_loading, 2)) then
            log_likelihood = -huge(1.0_dp)
            info = 1
            return
        end if
        effective_loading = random_loading + regression
        log_likelihood = 0.0_dp
        do i = 1, n
            mean_row = matmul(x(i, :), beta) + matmul(effective_loading, random_effects(i, :))
            call mvn_log_density(y(i, :), mean_row, r_matrix, row_log_density, info)
            if (info /= 0) then
                log_likelihood = -huge(1.0_dp)
                return
            end if
            log_likelihood = log_likelihood + row_log_density
        end do
    end subroutine covu_gaussian_loglik

    pure subroutine covu_gaussian_mixed_mcmc(y, x, random_loading, beta_prior_mean, beta_prior_precision, &
                                             joint_prior_scale, joint_prior_df, iterations, burn, thin, state, &
                                             result, info, initial_joint_covariance, update_mode, fixed_block)
        real(dp), intent(in) :: y(:, :) !! n by r Gaussian responses with one row per joint random/residual level.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across response traits.
        real(dp), intent(in) :: random_loading(:, :) !! r by g direct loading of covu random effects onto responses.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by r Gaussian fixed-effect prior mean matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! pr by pr fixed-effect prior precision in trait-major order.
        real(dp), intent(in) :: joint_prior_scale(:, :) !! (g+r) square prior scale for the joint G-R covariance.
        real(dp), intent(in) :: joint_prior_df !! Prior degrees of freedom for the joint covariance update.
        integer, intent(in) :: iterations !! Total MCMC iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by coefficients and joint covariance draws.
        type(covu_gaussian_mcmc_result), intent(out) :: result !! Retained Gaussian covu posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or numerical failure.
        real(dp), optional, intent(in) :: initial_joint_covariance(:, :) !! Optional starting joint covariance matrix.
        integer, optional, intent(in) :: update_mode !! Joint covariance update code 0, 1, 2, 3, 4, or 6; defaults to one.
        real(dp), optional, intent(in) :: fixed_block(:, :) !! Fixed lower-right block required by update modes 2 or 4.
        logical :: accepted
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: conditional_covariance(:, :)
        real(dp), allocatable :: fixed_block_value(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: joint_covariance(:, :)
        real(dp), allocatable :: joint_data(:, :)
        real(dp), allocatable :: posterior_sum(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: raw_residual(:, :)
        real(dp), allocatable :: regression(:, :)
        real(dp), allocatable :: updated_conditional(:, :)
        real(dp), allocatable :: updated_joint(:, :)
        real(dp), allocatable :: updated_regression(:, :)
        real(dp) :: log_likelihood
        integer :: dim_joint
        integer :: g
        integer :: iteration
        integer :: mode
        integer :: n
        integer :: nsave
        integer :: p
        integer :: r
        integer :: save_index

        info = 0
        n = size(y, 1)
        r = size(y, 2)
        p = size(x, 2)
        g = size(random_loading, 2)
        dim_joint = g + r
        if (n < 1 .or. r < 1 .or. p < 1 .or. g < 1 .or. size(x, 1) /= n .or. &
            size(random_loading, 1) /= r .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= r .or. size(beta_prior_precision, 1) /= p * r .or. &
            size(beta_prior_precision, 2) /= p * r .or. size(joint_prior_scale, 1) /= dim_joint .or. &
            size(joint_prior_scale, 2) /= dim_joint .or. joint_prior_df <= real(dim_joint - 1, dp) .or. &
            iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 1
            return
        end if
        mode = 1
        if (present(update_mode)) mode = update_mode
        if (.not. any(mode == [0, 1, 2, 3, 4, 6])) then
            info = 1
            return
        end if
        if (present(initial_joint_covariance)) then
            if (size(initial_joint_covariance, 1) /= dim_joint .or. &
                size(initial_joint_covariance, 2) /= dim_joint) then
                info = 1
                return
            end if
            joint_covariance = initial_joint_covariance
        else
            joint_covariance = joint_prior_scale / max(joint_prior_df - real(dim_joint + 1, dp), 1.0_dp)
        end if
        if (present(fixed_block)) then
            fixed_block_value = fixed_block
        else
            allocate(fixed_block_value(0, 0))
        end if
        if ((mode == 2 .or. mode == 4) .and. &
            (size(fixed_block_value, 1) /= r .or. size(fixed_block_value, 2) /= r)) then
            info = 1
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 1
            return
        end if

        call joint_gr_decompose(joint_covariance, g, g_matrix, regression, conditional_covariance, info)
        if (info /= 0) return
        allocate(result%beta(p, r, nsave), result%random_effects(n, g, nsave))
        allocate(result%joint_covariance(dim_joint, dim_joint, nsave), result%g(g, g, nsave))
        allocate(result%r(r, r, nsave), result%regression(r, g, nsave), result%log_likelihood(nsave))
        save_index = 0

        do iteration = 1, iterations
            call covu_coefficient_conditional(y, x, random_loading, g_matrix, conditional_covariance, regression, &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            raw_residual = y - matmul(x, beta) - matmul(random_effects, transpose(random_loading))
            allocate(joint_data(n, dim_joint))
            joint_data(:, 1:g) = random_effects
            joint_data(:, g + 1:dim_joint) = raw_residual
            posterior_sum = joint_prior_scale + matmul(transpose(joint_data), joint_data)
            deallocate(joint_data)
            call joint_gr_covariance_update(state, mode, posterior_sum, real(n, dp), joint_prior_df, &
                joint_prior_scale, joint_covariance, g, fixed_block_value, updated_joint, updated_conditional, &
                updated_regression, accepted, info)
            if (info /= 0) return
            joint_covariance = updated_joint
            conditional_covariance = updated_conditional
            regression = updated_regression
            g_matrix = joint_covariance(1:g, 1:g)

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%joint_covariance(:, :, save_index) = joint_covariance
                result%g(:, :, save_index) = g_matrix
                result%r(:, :, save_index) = conditional_covariance
                result%regression(:, :, save_index) = regression
                call covu_gaussian_loglik(y, x, random_loading, beta, random_effects, conditional_covariance, &
                    regression, log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine covu_gaussian_mixed_mcmc

end module mcmcglmm_covu_sampler
