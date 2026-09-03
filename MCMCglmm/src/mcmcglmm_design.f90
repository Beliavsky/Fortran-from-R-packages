! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_design
    use r_kinds, only : dp
    use r_linalg, only : solve_spd
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_math, only : logistic
    use mcmcglmm_matrix, only : kronecker_product, mvn_log_density, sample_mvn_covariance
    implicit none
    private

    public :: path_matrix
    public :: sir_matrix
    public :: multiple_membership_design
    public :: gelman_prior_design
    public :: random_effect_covariance
    public :: d_divergence_mc

contains

    pure subroutine path_matrix(cause, effect, dimension, matrix_psi, info)
        integer, intent(in) :: cause(:) !! One-based cause indices for directed paths.
        integer, intent(in) :: effect(:) !! One-based effect indices paired elementwise with cause.
        integer, intent(in) :: dimension !! Positive square path-matrix order k.
        real(dp), allocatable, intent(out) :: matrix_psi(:, :) !! k by k matrix with one at each (effect,cause) path.
        integer, intent(out) :: info !! Zero on success; nonzero for size mismatch, invalid dimension, or out-of-range path indices.
        integer :: i

        info = 0
        if (size(cause) /= size(effect) .or. dimension < 1) then
            allocate(matrix_psi(0, 0))
            info = 1
            return
        end if
        if (any(cause < 1) .or. any(cause > dimension) .or. any(effect < 1) .or. any(effect > dimension)) then
            allocate(matrix_psi(0, 0))
            info = 2
            return
        end if
        allocate(matrix_psi(dimension, dimension))
        matrix_psi = 0.0_dp
        do i = 1, size(cause)
            matrix_psi(effect(i), cause(i)) = 1.0_dp
        end do
    end subroutine path_matrix

    pure subroutine sir_matrix(x_cause, x_effect, zero_diagonal, interaction, info)
        real(dp), intent(in) :: x_cause(:, :) !! n by p numerical cause-side design used by sir's interaction map.
        real(dp), intent(in) :: x_effect(:, :) !! n by p effect-side design paired columnwise with x_cause.
        logical, intent(in) :: zero_diagonal !! If true, set self-interaction diagonal entries to zero after multiplication.
        real(dp), allocatable, intent(out) :: interaction(:, :) !! Allocated n by n matrix X_cause times transpose(X_effect).
        integer, intent(out) :: info !! Zero on success; nonzero when row or column dimensions are incompatible.
        integer :: i
        integer :: n

        info = 0
        n = size(x_cause, 1)
        if (n < 1 .or. size(x_effect, 1) /= n .or. size(x_cause, 2) /= size(x_effect, 2)) then
            allocate(interaction(0, 0))
            info = 1
            return
        end if
        interaction = matmul(x_cause, transpose(x_effect))
        if (zero_diagonal) then
            do i = 1, n
                interaction(i, i) = 0.0_dp
            end do
        end if
    end subroutine sir_matrix

    pure subroutine multiple_membership_design(parts, combined, active_columns)
        real(dp), intent(in) :: parts(:, :, :) !! Common n by p design matrices, one membership component per plane.
        real(dp), allocatable, intent(out) :: combined(:, :) !! Sum of component designs with original column order.
        integer, allocatable, intent(out) :: active_columns(:) !! One-based combined-design columns containing nonzero entries.
        logical, allocatable :: active(:)
        integer :: nactive
        integer :: j

        combined = sum(parts, dim=3)
        allocate(active(size(combined, 2)))
        do j = 1, size(combined, 2)
            active(j) = any(abs(combined(:, j)) > 0.0_dp)
        end do
        nactive = count(active)
        allocate(active_columns(nactive))
        nactive = 0
        do j = 1, size(active)
            if (active(j)) then
                nactive = nactive + 1
                active_columns(nactive) = j
            end if
        end do
    end subroutine multiple_membership_design

    pure subroutine gelman_prior_design(x_original, x_scaled, scale_value, intercept_scale, covariance, info)
        real(dp), intent(in) :: x_original(:, :) !! Original X1 matrix mapping standardized coefficients to the fitted design.
        real(dp), intent(in) :: x_scaled(:, :) !! Standardized model matrix X2 constructed from scaled continuous/binary predictors.
        real(dp), intent(in) :: scale_value !! Prior standard deviation assigned to non-intercept standardized coefficients.
        real(dp), intent(in) :: intercept_scale !! Prior standard deviation assigned to the first standardized coefficient.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Fixed-effect prior covariance P I P^T from gelman.prior.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible design dimensions or singular X1^T X1.
        real(dp), allocatable :: cross(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: mapping(:, :)
        real(dp), allocatable :: prior_scaled(:, :)
        integer :: i
        integer :: p

        info = 0
        if (size(x_original, 1) /= size(x_scaled, 1) .or. size(x_original, 2) /= size(x_scaled, 2) .or. &
            scale_value < 0.0_dp .or. intercept_scale < 0.0_dp) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        p = size(x_original, 2)
        gram = matmul(transpose(x_original), x_original)
        cross = matmul(transpose(x_original), x_scaled)
        allocate(mapping(p, p))
        call solve_spd(gram, cross, mapping, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        allocate(prior_scaled(p, p))
        prior_scaled = 0.0_dp
        do i = 1, p
            prior_scaled(i, i) = scale_value ** 2
        end do
        if (p > 0) prior_scaled(1, 1) = intercept_scale ** 2
        covariance = matmul(mapping, matmul(prior_scaled, transpose(mapping)))
    end subroutine gelman_prior_design

    pure subroutine random_effect_covariance(design, trait_covariance, level_covariance, covariance, info)
        real(dp), intent(in) :: design(:, :) !! Random-effect design Z with columns ordered by trait blocks over levels.
        real(dp), intent(in) :: trait_covariance(:, :) !! Trait covariance matrix Vtmp.
        real(dp), intent(in) :: level_covariance(:, :) !! Level covariance; identity for independent levels or inverse(ginverse).
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated marginal covariance Z (trait_cov tensor level_cov) Z^T.
        integer, intent(out) :: info !! Zero on success; nonzero for nonsquare factors or incompatible design column count.
        real(dp), allocatable :: random_covariance(:, :)
        integer :: levels
        integer :: traits

        info = 0
        traits = size(trait_covariance, 1)
        levels = size(level_covariance, 1)
        if (size(trait_covariance, 2) /= traits .or. size(level_covariance, 2) /= levels .or. &
            size(design, 2) /= traits * levels) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        call kronecker_product(trait_covariance, level_covariance, random_covariance)
        covariance = matmul(design, matmul(random_covariance, transpose(design)))
    end subroutine random_effect_covariance

    pure subroutine d_divergence_mc(state, covariance_a, covariance_b, n_sample, divergence, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by Monte Carlo draws from the first covariance matrix.
        real(dp), intent(in) :: covariance_a(:, :) !! First positive-definite zero-mean Gaussian covariance matrix CA.
        real(dp), intent(in) :: covariance_b(:, :) !! Second matching positive-definite zero-mean Gaussian covariance CB.
        integer, intent(in) :: n_sample !! Positive Monte Carlo sample count corresponding to Ddivergence's n argument.
        real(dp), intent(out) :: divergence !! Monte Carlo divergence sqrt(1 - 2 mean(g/(f+g))).
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions/sample count or a failed Gaussian operation.
        real(dp), allocatable :: draw(:)
        real(dp), allocatable :: zero_mean(:)
        real(dp) :: log_f
        real(dp) :: log_g
        real(dp) :: ratio_sum
        integer :: i
        integer :: n

        info = 0
        n = size(covariance_a, 1)
        if (size(covariance_a, 2) /= n .or. size(covariance_b, 1) /= n .or. size(covariance_b, 2) /= n .or. &
            n < 1 .or. n_sample < 1) then
            divergence = 0.0_dp
            info = 1
            return
        end if
        allocate(zero_mean(n))
        zero_mean = 0.0_dp
        ratio_sum = 0.0_dp
        do i = 1, n_sample
            call sample_mvn_covariance(state, zero_mean, covariance_a, draw, info)
            if (info /= 0) then
                divergence = 0.0_dp
                return
            end if
            call mvn_log_density(draw, zero_mean, covariance_a, log_f, info)
            if (info /= 0) then
                divergence = 0.0_dp
                return
            end if
            call mvn_log_density(draw, zero_mean, covariance_b, log_g, info)
            if (info /= 0) then
                divergence = 0.0_dp
                return
            end if
            ratio_sum = ratio_sum + logistic(log_g - log_f)
        end do
        divergence = sqrt(max(0.0_dp, 1.0_dp - 2.0_dp * ratio_sum / real(n_sample, dp)))
    end subroutine d_divergence_mc

end module mcmcglmm_design
