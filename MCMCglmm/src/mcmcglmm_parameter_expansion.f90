! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_parameter_expansion
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : sample_mvn_precision
    implicit none
    private

    public :: parameter_expansion_conditional
    public :: apply_parameter_expansion
    public :: expanded_covariance

contains

    pure subroutine parameter_expansion_conditional(state, residual, expansion_design, observation_precision, &
                                                     prior_mean, prior_precision, alpha, mean_alpha, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Gaussian redundant-parameter draw.
        real(dp), intent(in) :: residual(:) !! Working response after subtracting all non-expanded model contributions.
        real(dp), intent(in) :: expansion_design(:, :) !! n by a design whose columns are unscaled random-effect contributions.
        real(dp), intent(in) :: observation_precision(:, :) !! n by n Gaussian working-response precision matrix.
        real(dp), intent(in) :: prior_mean(:) !! a-vector prior mean corresponding to upstream `alpha.mu`.
        real(dp), intent(in) :: prior_precision(:, :) !! a by a prior precision, inverse of upstream `alpha.V`.
        real(dp), allocatable, intent(out) :: alpha(:) !! Allocated redundant working-parameter draw.
        real(dp), allocatable, intent(out) :: mean_alpha(:) !! Allocated conditional posterior mean of the redundant parameters.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes or a failed SPD solve/sample.
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: weighted_design(:, :)
        real(dp), allocatable :: weighted_residual(:)
        integer :: a
        integer :: n

        info = 0
        n = size(residual)
        a = size(prior_mean)
        if (n < 1 .or. a < 1 .or. size(expansion_design, 1) /= n .or. size(expansion_design, 2) /= a .or. &
            size(observation_precision, 1) /= n .or. size(observation_precision, 2) /= n .or. &
            size(prior_precision, 1) /= a .or. size(prior_precision, 2) /= a) then
            allocate(alpha(0), mean_alpha(0))
            info = 1
            return
        end if

        weighted_design = matmul(observation_precision, expansion_design)
        weighted_residual = matmul(observation_precision, residual)
        precision = prior_precision + matmul(transpose(expansion_design), weighted_design)
        rhs = matmul(prior_precision, prior_mean) + matmul(transpose(expansion_design), weighted_residual)
        call sample_mvn_precision(state, rhs, precision, alpha, mean_alpha, info)
    end subroutine parameter_expansion_conditional

    pure subroutine apply_parameter_expansion(unscaled_effects, alpha, scaled_effects, info)
        real(dp), intent(in) :: unscaled_effects(:, :) !! q by a working random-effect matrix before redundant scaling.
        real(dp), intent(in) :: alpha(:) !! a-vector redundant scale parameters.
        real(dp), allocatable, intent(out) :: scaled_effects(:, :) !! Allocated q by a posterior-scale random-effect matrix.
        integer, intent(out) :: info !! Zero on success; one when the trait/scale dimensions do not agree.
        integer :: j

        info = 0
        if (size(unscaled_effects, 2) /= size(alpha)) then
            allocate(scaled_effects(0, 0))
            info = 1
            return
        end if
        allocate(scaled_effects(size(unscaled_effects, 1), size(unscaled_effects, 2)))
        do j = 1, size(alpha)
            scaled_effects(:, j) = unscaled_effects(:, j) * alpha(j)
        end do
    end subroutine apply_parameter_expansion

    pure subroutine expanded_covariance(working_covariance, alpha, covariance, info)
        real(dp), intent(in) :: working_covariance(:, :) !! a by a covariance matrix on the redundant working scale.
        real(dp), intent(in) :: alpha(:) !! a-vector redundant scale parameters.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated posterior-scale covariance diag(alpha) G diag(alpha).
        integer, intent(out) :: info !! Zero on success; one when dimensions do not agree.
        integer :: i
        integer :: j
        integer :: a

        info = 0
        a = size(alpha)
        if (size(working_covariance, 1) /= a .or. size(working_covariance, 2) /= a) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        allocate(covariance(a, a))
        do j = 1, a
            do i = 1, a
                covariance(i, j) = working_covariance(i, j) * alpha(i) * alpha(j)
            end do
        end do
    end subroutine expanded_covariance

end module mcmcglmm_parameter_expansion
