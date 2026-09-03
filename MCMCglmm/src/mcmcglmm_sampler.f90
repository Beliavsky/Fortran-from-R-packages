! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_precision
    use mcmcglmm_distributions, only : inverse_wishart_sample, truncated_normal_sample
    use mcmcglmm_math, only : normal_cdf
    implicit none
    private

    type, public :: ordinal_mcmc_result
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: log_likelihood(:)
        real(dp), allocatable :: last_liability(:)
    end type ordinal_mcmc_result

    type, public :: gaussian_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
    end type gaussian_mcmc_result

    public :: coefficient_conditional
    public :: gaussian_mixed_mcmc
    public :: ordinal_probit_mixed_mcmc

contains

    pure integer function joint_index(trait, effect, effects_per_trait) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based effect index within a trait block.
        integer, intent(in) :: effects_per_trait !! Number of fixed plus random coefficients per trait.

        index_value = (trait - 1) * effects_per_trait + effect
    end function joint_index

    pure integer function beta_index(trait, effect, fixed_effects) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-effect index within a trait.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect columns.

        index_value = (trait - 1) * fixed_effects + effect
    end function beta_index

    pure subroutine unpack_coefficients(packed, fixed_effects, random_levels, traits, beta, random_effects)
        real(dp), intent(in) :: packed(:) !! Trait-major vector containing fixed then random effects for each trait.
        integer, intent(in) :: fixed_effects !! Number p of fixed-effect coefficients per trait.
        integer, intent(in) :: random_levels !! Number q of random-effect levels per trait.
        integer, intent(in) :: traits !! Number of response traits.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated p by traits fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated q by traits random-effect matrix.
        integer :: a
        integer :: k
        integer :: m

        m = fixed_effects + random_levels
        allocate(beta(fixed_effects, traits), random_effects(random_levels, traits))
        do a = 1, traits
            do k = 1, fixed_effects
                beta(k, a) = packed(joint_index(a, k, m))
            end do
            do k = 1, random_levels
                random_effects(k, a) = packed(joint_index(a, fixed_effects + k, m))
            end do
        end do
    end subroutine unpack_coefficients

    pure subroutine gaussian_log_likelihood(y, x, z, beta, random_effects, r_matrix, log_likelihood, info)
        real(dp), intent(in) :: y(:, :) !! n by n_trait response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: beta(:, :) !! p by n_trait fixed-effect matrix.
        real(dp), intent(in) :: random_effects(:, :) !! q by n_trait random-effect matrix.
        real(dp), intent(in) :: r_matrix(:, :) !! n_trait by n_trait residual covariance matrix.
        real(dp), intent(out) :: log_likelihood !! Sum of row-wise multivariate-normal log densities.
        integer, intent(out) :: info !! Zero on success; nonzero on shape mismatch or covariance failure.
        real(dp), allocatable :: mean_row(:)
        real(dp), allocatable :: zero_mean(:)
        real(dp) :: row_log_density
        integer :: i
        integer :: n
        integer :: traits

        info = 0
        n = size(y, 1)
        traits = size(y, 2)
        if (size(x, 1) /= n .or. size(z, 1) /= n .or. size(beta, 1) /= size(x, 2) .or. &
            size(beta, 2) /= traits .or. size(random_effects, 1) /= size(z, 2) .or. &
            size(random_effects, 2) /= traits .or. size(r_matrix, 1) /= traits .or. &
            size(r_matrix, 2) /= traits) then
            log_likelihood = -huge(1.0_dp)
            info = 1
            return
        end if
        allocate(mean_row(traits), zero_mean(traits))
        zero_mean = 0.0_dp
        log_likelihood = 0.0_dp
        do i = 1, n
            mean_row = matmul(x(i, :), beta) + matmul(z(i, :), random_effects)
            call mvn_log_density(y(i, :) - mean_row, zero_mean, r_matrix, row_log_density, info)
            if (info /= 0) then
                log_likelihood = -huge(1.0_dp)
                return
            end if
            log_likelihood = log_likelihood + row_log_density
        end do
    end subroutine gaussian_log_likelihood

    pure subroutine coefficient_conditional(y, x, z, a_inverse, g_matrix, r_matrix, beta_prior_mean, &
                                            beta_prior_precision, state, beta, random_effects, info)
        real(dp), intent(in) :: y(:, :) !! n by n_trait Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix shared across traits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q row precision for random effects, such as pedigree A^{-1}.
        real(dp), intent(in) :: g_matrix(:, :) !! n_trait by n_trait random-effect covariance G.
        real(dp), intent(in) :: r_matrix(:, :) !! n_trait by n_trait residual covariance R.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by n_trait Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*n_trait square prior precision in trait-major ordering.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the joint Gaussian coefficient draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated sampled random-effect matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for shape or SPD linear-algebra failure.
        real(dp), allocatable :: beta_mean_vector(:)
        real(dp), allocatable :: beta_rhs(:)
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: g_inverse(:, :)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: w(:, :)
        real(dp), allocatable :: wtw(:, :)
        real(dp), allocatable :: wty(:, :)
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
            size(g_matrix, 1) /= traits .or. size(g_matrix, 2) /= traits .or. &
            size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits .or. &
            size(beta_prior_mean, 1) /= p .or. size(beta_prior_mean, 2) /= traits .or. &
            size(beta_prior_precision, 1) /= p * traits .or. size(beta_prior_precision, 2) /= p * traits) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if

        allocate(w(size(y, 1), m))
        w(:, 1:p) = x
        w(:, p + 1:m) = z
        wtw = matmul(transpose(w), w)
        wty = matmul(transpose(w), y)
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

        allocate(precision(m * traits, m * traits), rhs(m * traits), beta_mean_vector(p * traits))
        precision = 0.0_dp
        rhs = 0.0_dp
        do a = 1, traits
            do k = 1, p
                beta_mean_vector(beta_index(a, k, p)) = beta_prior_mean(k, a)
            end do
        end do
        beta_rhs = matmul(beta_prior_precision, beta_mean_vector)

        do a = 1, traits
            do b = 1, traits
                do j = 1, m
                    do i = 1, m
                        precision(joint_index(a, i, m), joint_index(b, j, m)) = &
                            precision(joint_index(a, i, m), joint_index(b, j, m)) + r_inverse(a, b) * wtw(i, j)
                    end do
                end do
                do j = 1, q
                    do i = 1, q
                        precision(joint_index(a, p + i, m), joint_index(b, p + j, m)) = &
                            precision(joint_index(a, p + i, m), joint_index(b, p + j, m)) + &
                            g_inverse(a, b) * a_inverse(i, j)
                    end do
                end do
            end do
            do i = 1, m
                do b = 1, traits
                    rhs(joint_index(a, i, m)) = rhs(joint_index(a, i, m)) + r_inverse(a, b) * wty(i, b)
                end do
            end do
        end do

        do a = 1, traits
            do b = 1, traits
                do j = 1, p
                    do i = 1, p
                        precision(joint_index(a, i, m), joint_index(b, j, m)) = &
                            precision(joint_index(a, i, m), joint_index(b, j, m)) + &
                            beta_prior_precision(beta_index(a, i, p), beta_index(b, j, p))
                    end do
                end do
            end do
            do i = 1, p
                rhs(joint_index(a, i, m)) = rhs(joint_index(a, i, m)) + beta_rhs(beta_index(a, i, p))
            end do
        end do

        call sample_mvn_precision(state, rhs, precision, packed, conditional_mean, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call unpack_coefficients(packed, p, q, traits, beta, random_effects)
    end subroutine coefficient_conditional

    pure subroutine gaussian_mixed_mcmc(y, x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
                                        g_prior_scale, g_prior_df, r_prior_scale, r_prior_df, iterations, burn, thin, &
                                        state, result, info)
        real(dp), intent(in) :: y(:, :) !! n by n_trait Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix shared across traits.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision, typically pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by n_trait fixed-effect Gaussian prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*n_trait square fixed-effect prior precision in trait-major ordering.
        real(dp), intent(in) :: g_prior_scale(:, :) !! n_trait by n_trait standard inverse-Wishart scale matrix for G.
        real(dp), intent(in) :: g_prior_df !! G inverse-Wishart prior degrees of freedom, required to exceed n_trait minus one.
        real(dp), intent(in) :: r_prior_scale(:, :) !! n_trait by n_trait standard inverse-Wishart scale matrix for R.
        real(dp), intent(in) :: r_prior_df !! R inverse-Wishart prior degrees of freedom, required to exceed n_trait minus one.
        integer, intent(in) :: iterations !! Total Gibbs iterations, including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded before thinning.
        integer, intent(in) :: thin !! Positive thinning interval applied after burn-in.
        type(rng_state), intent(inout) :: state !! Deterministic generator state advanced through the Gibbs chain.
        type(gaussian_mcmc_result), intent(out) :: result !! Posterior samples retained after burn-in and thinning.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/settings or a failed matrix operation.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp) :: log_likelihood
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
            size(r_prior_scale, 1) /= traits .or. size(r_prior_scale, 2) /= traits) then
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
        allocate(result%g(traits, traits, nsave), result%r(traits, traits, nsave), result%log_likelihood(nsave))
        beta = beta_prior_mean
        allocate(random_effects(q, traits))
        random_effects = 0.0_dp
        g_matrix = g_prior_scale / max(g_prior_df - real(traits + 1, dp), 1.0_dp)
        r_matrix = r_prior_scale / max(r_prior_df - real(traits + 1, dp), 1.0_dp)
        save_index = 0

        do iteration = 1, iterations
            call coefficient_conditional(y, x, z, a_inverse, g_matrix, r_matrix, beta_prior_mean, &
                                         beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            g_scale_post = g_prior_scale + matmul(transpose(random_effects), matmul(a_inverse, random_effects))
            call inverse_wishart_sample(state, g_scale_post, g_prior_df + real(q, dp), g_matrix, info)
            if (info /= 0) return

            residual = y - matmul(x, beta) - matmul(z, random_effects)
            r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
            call inverse_wishart_sample(state, r_scale_post, r_prior_df + real(n, dp), r_matrix, info)
            if (info /= 0) return

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                call gaussian_log_likelihood(y, x, z, beta, random_effects, r_matrix, log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine gaussian_mixed_mcmc


    pure subroutine ordinal_log_likelihood(y_category, cutpoints, x, z, beta, random_effects, log_likelihood, info)
        integer, intent(in) :: y_category(:) !! One-based ordered response categories in the range 1 to n_cutpoints minus one.
        real(dp), intent(in) :: cutpoints(:) !! Strictly increasing category boundaries of length n_category plus one.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: beta(:) !! p-vector of fixed effects.
        real(dp), intent(in) :: random_effects(:) !! q-vector of random effects.
        real(dp), intent(out) :: log_likelihood !! Sum of ordered-probit log category probabilities.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid categories, cutpoints, or design shapes.
        real(dp) :: eta
        real(dp) :: probability
        integer :: i
        integer :: n

        info = 0
        n = size(y_category)
        if (size(x, 1) /= n .or. size(z, 1) /= n .or. size(beta) /= size(x, 2) .or. &
            size(random_effects) /= size(z, 2) .or. size(cutpoints) < 2) then
            log_likelihood = -huge(1.0_dp)
            info = 1
            return
        end if
        if (any(y_category < 1) .or. any(y_category >= size(cutpoints))) then
            log_likelihood = -huge(1.0_dp)
            info = 2
            return
        end if
        do i = 2, size(cutpoints)
            if (cutpoints(i) <= cutpoints(i - 1)) then
                log_likelihood = -huge(1.0_dp)
                info = 3
                return
            end if
        end do
        log_likelihood = 0.0_dp
        do i = 1, n
            eta = dot_product(x(i, :), beta) + dot_product(z(i, :), random_effects)
            probability = normal_cdf(cutpoints(y_category(i) + 1) - eta) - &
                normal_cdf(cutpoints(y_category(i)) - eta)
            log_likelihood = log_likelihood + log(max(probability, tiny(1.0_dp)))
        end do
    end subroutine ordinal_log_likelihood

    pure subroutine ordinal_probit_mixed_mcmc(y_category, cutpoints, x, z, a_inverse, beta_prior_mean, &
                                              beta_prior_precision, g_prior_scale, g_prior_df, iterations, burn, thin, &
                                              state, result, info)
        integer, intent(in) :: y_category(:) !! One-based ordered response categories for each observation.
        real(dp), intent(in) :: cutpoints(:) !! Increasing latent-normal category bounds; large endpoints may represent infinity.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q random-effect design matrix.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-effect row precision, typically pedigree A^{-1}.
        real(dp), intent(in) :: beta_prior_mean(:) !! p-vector Gaussian prior mean for fixed effects.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p by p Gaussian prior precision for fixed effects.
        real(dp), intent(in) :: g_prior_scale !! Positive scalar inverse-Wishart/inverse-gamma scale for random-effect variance.
        real(dp), intent(in) :: g_prior_df !! Positive prior degrees of freedom for random-effect variance.
        integer, intent(in) :: iterations !! Total latent-variable Gibbs iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by liabilities, coefficients, and variance draws.
        type(ordinal_mcmc_result), intent(out) :: result !! Retained ordered-probit posterior draws and final liability state.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/settings or a numerical failure.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: beta_mean_matrix(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: g_scale_matrix(:, :)
        real(dp), allocatable :: liability(:)
        real(dp), allocatable :: random_effects(:, :)
        real(dp) :: log_likelihood
        real(dp) :: mean_value
        real(dp) :: r_matrix(1, 1)
        integer :: i
        integer :: iteration
        integer :: n
        integer :: nsave
        integer :: p
        integer :: q
        integer :: save_index

        info = 0
        n = size(y_category)
        p = size(x, 2)
        q = size(z, 2)
        if (n < 1 .or. p < 1 .or. q < 1 .or. size(x, 1) /= n .or. size(z, 1) /= n .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(beta_prior_mean) /= p .or. &
            size(beta_prior_precision, 1) /= p .or. size(beta_prior_precision, 2) /= p .or. &
            size(cutpoints) < 2 .or. g_prior_scale <= 0.0_dp .or. g_prior_df <= 0.0_dp) then
            info = 1
            return
        end if
        if (any(y_category < 1) .or. any(y_category >= size(cutpoints)) .or. iterations <= burn .or. &
            burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        do i = 2, size(cutpoints)
            if (cutpoints(i) <= cutpoints(i - 1)) then
                info = 3
                return
            end if
        end do
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 4
            return
        end if
        allocate(result%beta(p, nsave), result%random_effects(q, nsave), result%g(nsave), result%log_likelihood(nsave))
        allocate(result%last_liability(n), liability(n), beta_mean_matrix(p, 1), random_effects(q, 1), g_matrix(1, 1))
        beta_mean_matrix(:, 1) = beta_prior_mean
        allocate(beta(p, 1))
        beta(:, 1) = beta_prior_mean
        random_effects = 0.0_dp
        g_matrix(1, 1) = g_prior_scale / max(g_prior_df - 2.0_dp, 1.0_dp)
        r_matrix(1, 1) = 1.0_dp
        liability = 0.0_dp
        save_index = 0

        do iteration = 1, iterations
            do i = 1, n
                mean_value = dot_product(x(i, :), beta(:, 1)) + dot_product(z(i, :), random_effects(:, 1))
                call truncated_normal_sample(state, mean_value, 1.0_dp, cutpoints(y_category(i)), &
                    cutpoints(y_category(i) + 1), liability(i), info)
                if (info /= 0) return
            end do
            call coefficient_conditional(reshape(liability, [n, 1]), x, z, a_inverse, g_matrix, r_matrix, &
                beta_mean_matrix, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return
            allocate(g_scale_matrix(1, 1))
            g_scale_matrix(1, 1) = g_prior_scale + dot_product(random_effects(:, 1), &
                matmul(a_inverse, random_effects(:, 1)))
            call inverse_wishart_sample(state, g_scale_matrix, g_prior_df + real(q, dp), g_matrix, info)
            if (allocated(g_scale_matrix)) deallocate(g_scale_matrix)
            if (info /= 0) return

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, save_index) = beta(:, 1)
                result%random_effects(:, save_index) = random_effects(:, 1)
                result%g(save_index) = g_matrix(1, 1)
                call ordinal_log_likelihood(y_category, cutpoints, x, z, beta(:, 1), random_effects(:, 1), &
                    log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
        result%last_liability = liability
    end subroutine ordinal_probit_mixed_mcmc

end module mcmcglmm_sampler
