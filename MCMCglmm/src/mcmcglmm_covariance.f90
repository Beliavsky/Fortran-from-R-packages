! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 covariance-update kernels; see NOTICE.md and upstream/.
module mcmcglmm_covariance
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix, spd_inverse_logdet
    use mcmcglmm_rng, only : rng_state, rng_uniform
    use mcmcglmm_matrix, only : covariance_to_correlation, kronecker_product, sample_mvn_precision
    use mcmcglmm_distributions, only : inverse_wishart_sample
    use mcmcglmm_ante, only : ante_covariance_samples
    implicit none
    private

    public :: conditioned_covariance_update
    public :: correlation_structure_update
    public :: correlation_submatrix_update
    public :: identity_direct_sum_update
    public :: covariance_update_dispatch

contains

    pure subroutine conditioned_inverse_wishart_precision(state, precision, degrees_freedom, split, fixed_block, &
                                                          covariance, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by inverse-Wishart and Gaussian draws.
        real(dp), intent(in) :: precision(:, :) !! Full positive-definite precision matrix used by the conditional IW law.
        real(dp), intent(in) :: degrees_freedom !! Inverse-Wishart degrees of freedom for the full covariance matrix.
        integer, intent(in) :: split !! Number of leading coordinates sampled before the fixed lower-right block.
        real(dp), intent(in) :: fixed_block(:, :) !! Positive-definite lower-right covariance block held fixed.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated conditional inverse-Wishart covariance draw.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shapes or a failed matrix operation.
        real(dp), allocatable :: a11(:, :)
        real(dp), allocatable :: a11_inverse(:, :)
        real(dp), allocatable :: a12(:, :)
        real(dp), allocatable :: a22(:, :)
        real(dp), allocatable :: half(:, :)
        real(dp), allocatable :: half_fixed(:, :)
        real(dp), allocatable :: iw11(:, :)
        real(dp), allocatable :: noise(:)
        real(dp), allocatable :: noise_mean(:)
        real(dp), allocatable :: noise_precision(:, :)
        real(dp), allocatable :: sampled_mean(:)
        real(dp), allocatable :: schur(:, :)
        real(dp), allocatable :: schur_inverse(:, :)
        real(dp), allocatable :: t1(:, :)
        real(dp), allocatable :: t1_inverse(:, :)
        integer :: n
        integer :: n_fixed

        info = 0
        n = size(precision, 1)
        n_fixed = n - split
        if (size(precision, 2) /= n .or. split < 1 .or. split >= n .or. &
            size(fixed_block, 1) /= n_fixed .or. size(fixed_block, 2) /= n_fixed .or. &
            degrees_freedom <= real(n - 1, dp)) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if

        a11 = precision(1:split, 1:split)
        a12 = precision(1:split, split + 1:n)
        a22 = precision(split + 1:n, split + 1:n)
        call inverse_matrix(a11, a11_inverse, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        schur = a22 - matmul(transpose(a12), matmul(a11_inverse, a12))
        call inverse_matrix(schur, schur_inverse, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if

        call inverse_wishart_sample(state, a11_inverse, degrees_freedom, t1_inverse, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call inverse_matrix(t1_inverse, t1, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if

        half = matmul(a11_inverse, a12)
        call kronecker_product(t1, schur_inverse, noise_precision)
        allocate(noise_mean(split * n_fixed))
        noise_mean = 0.0_dp
        call sample_mvn_precision(state, noise_mean, noise_precision, noise, sampled_mean, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        half = -(half + reshape(noise, [split, n_fixed]))
        half_fixed = matmul(fixed_block, transpose(half))
        iw11 = t1_inverse + matmul(half, half_fixed)

        allocate(covariance(n, n))
        covariance = 0.0_dp
        covariance(1:split, 1:split) = iw11
        covariance(1:split, split + 1:n) = transpose(half_fixed)
        covariance(split + 1:n, 1:split) = half_fixed
        covariance(split + 1:n, split + 1:n) = fixed_block
    end subroutine conditioned_inverse_wishart_precision

    pure subroutine conditioned_covariance_update(state, posterior_sum, degrees_freedom, fix_index, fixed_block, &
                                                  covariance, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the conditional inverse-Wishart draw.
        real(dp), intent(in) :: posterior_sum(:, :) !! Cross-product plus prior scale matrix before inversion.
        real(dp), intent(in) :: degrees_freedom !! Full posterior inverse-Wishart degrees of freedom.
        integer, intent(in) :: fix_index !! One-based first coordinate of the covariance block held fixed.
        real(dp), intent(in) :: fixed_block(:, :) !! Lower-right covariance block beginning at fix_index.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated covariance draw satisfying the constraint.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid input or a failed matrix operation.
        real(dp), allocatable :: precision(:, :)
        integer :: n

        n = size(posterior_sum, 1)
        if (size(posterior_sum, 2) /= n .or. fix_index < 1 .or. fix_index > n) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        if (fix_index == 1) then
            if (size(fixed_block, 1) /= n .or. size(fixed_block, 2) /= n) then
                allocate(covariance(0, 0))
                info = 1
                return
            end if
            allocate(covariance(n, n))
            covariance = fixed_block
            info = 0
            return
        end if
        call inverse_matrix(posterior_sum, precision, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call conditioned_inverse_wishart_precision(state, precision, degrees_freedom, fix_index - 1, fixed_block, &
            covariance, info)
    end subroutine conditioned_covariance_update

    pure subroutine identity_direct_sum_update(state, posterior_sum, degrees_freedom, split, covariance, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the leading-block inverse-Wishart draw.
        real(dp), intent(in) :: posterior_sum(:, :) !! Cross-product plus prior scale matrix before inversion.
        real(dp), intent(in) :: degrees_freedom !! Posterior degrees of freedom used for the sampled leading block.
        integer, intent(in) :: split !! Number of leading coordinates sampled; remaining coordinates form an identity block.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated direct-sum covariance draw.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid input or a matrix/sampling failure.
        real(dp), allocatable :: leading_precision(:, :)
        real(dp), allocatable :: leading_scale(:, :)
        real(dp), allocatable :: leading_sample(:, :)
        real(dp), allocatable :: precision(:, :)
        integer :: i
        integer :: n

        info = 0
        n = size(posterior_sum, 1)
        if (size(posterior_sum, 2) /= n .or. split < 1 .or. split > n .or. &
            degrees_freedom <= real(split - 1, dp)) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        call inverse_matrix(posterior_sum, precision, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        leading_precision = precision(1:split, 1:split)
        call inverse_matrix(leading_precision, leading_scale, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call inverse_wishart_sample(state, leading_scale, degrees_freedom, leading_sample, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        allocate(covariance(n, n))
        covariance = 0.0_dp
        covariance(1:split, 1:split) = leading_sample
        do i = split + 1, n
            covariance(i, i) = 1.0_dp
        end do
    end subroutine identity_direct_sum_update

    pure subroutine correlation_structure_update(state, posterior_sum, sample_degrees_freedom, prior_degrees_freedom, &
                                                 prior_scale, old_covariance, covariance, accepted, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the proposal draw and MH uniform variate.
        real(dp), intent(in) :: posterior_sum(:, :) !! Cross-product plus prior scale matrix used by the correlation update.
        real(dp), intent(in) :: sample_degrees_freedom !! Effective sample degrees of freedom used in the IW proposal.
        real(dp), intent(in) :: prior_degrees_freedom !! Prior degrees of freedom entering the Liu-Daniels MH ratio.
        real(dp), intent(in) :: prior_scale(:, :) !! Prior scale; its diagonal fixes marginal variances of the output.
        real(dp), intent(in) :: old_covariance(:, :) !! Current covariance/correlation matrix before the update.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated accepted covariance with fixed marginal variances.
        logical, intent(out) :: accepted !! True when the proposed correlation matrix is accepted.
        integer, intent(out) :: info !! Zero on success; nonzero for shape, SPD, or sampling failure.
        real(dp), allocatable :: candidate(:, :)
        real(dp), allocatable :: candidate_correlation(:, :)
        real(dp), allocatable :: candidate_inverse(:, :)
        real(dp), allocatable :: crossproduct(:, :)
        real(dp), allocatable :: crossproduct_correlation(:, :)
        real(dp), allocatable :: old_correlation(:, :)
        real(dp), allocatable :: old_inverse(:, :)
        real(dp) :: log_alpha
        real(dp) :: logdet_candidate
        real(dp) :: logdet_old
        real(dp) :: u
        integer :: i
        integer :: j
        integer :: n

        info = 0
        accepted = .false.
        n = size(posterior_sum, 1)
        if (size(posterior_sum, 2) /= n .or. size(prior_scale, 1) /= n .or. size(prior_scale, 2) /= n .or. &
            size(old_covariance, 1) /= n .or. size(old_covariance, 2) /= n .or. &
            sample_degrees_freedom <= real(n - 1, dp) .or. prior_degrees_freedom < 0.0_dp) then
            allocate(covariance(0, 0))
            info = 1
            return
        end if
        do i = 1, n
            if (prior_scale(i, i) <= 0.0_dp) then
                allocate(covariance(0, 0))
                info = 1
                return
            end if
        end do

        crossproduct = posterior_sum - prior_scale
        call covariance_to_correlation(crossproduct, crossproduct_correlation, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call inverse_wishart_sample(state, crossproduct_correlation, sample_degrees_freedom, candidate, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call covariance_to_correlation(candidate, candidate_correlation, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call covariance_to_correlation(old_covariance, old_correlation, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call spd_inverse_logdet(old_correlation, old_inverse, logdet_old, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call spd_inverse_logdet(candidate_correlation, candidate_inverse, logdet_candidate, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if

        log_alpha = logdet_old - logdet_candidate
        do i = 1, n
            log_alpha = log_alpha + log(old_inverse(i, i)) - log(candidate_inverse(i, i))
        end do
        log_alpha = 0.5_dp * prior_degrees_freedom * log_alpha
        call rng_uniform(state, u)
        if (logdet_candidate >= log(1.0e-7_dp) .and. log(u) < log_alpha) then
            old_correlation = candidate_correlation
            accepted = .true.
        end if

        allocate(covariance(n, n))
        do j = 1, n
            do i = 1, n
                covariance(i, j) = old_correlation(i, j) * sqrt(prior_scale(i, i) * prior_scale(j, j))
            end do
        end do
    end subroutine correlation_structure_update

    pure subroutine correlation_submatrix_update(state, posterior_sum, sample_degrees_freedom, prior_degrees_freedom, &
                                                 prior_scale, split, old_fixed_block, covariance, fixed_block, &
                                                 accepted, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by correlation and conditional-IW draws.
        real(dp), intent(in) :: posterior_sum(:, :) !! Cross-product plus prior scale matrix for the full covariance.
        real(dp), intent(in) :: sample_degrees_freedom !! Effective sample degrees of freedom for the full covariance update.
        real(dp), intent(in) :: prior_degrees_freedom !! Prior degrees of freedom used by the correlation-block MH ratio.
        real(dp), intent(in) :: prior_scale(:, :) !! Full prior scale matrix; lower-right diagonal fixes block variances.
        integer, intent(in) :: split !! Number of leading covariance coordinates outside the correlation submatrix.
        real(dp), intent(in) :: old_fixed_block(:, :) !! Current lower-right correlation-constrained covariance block.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated full covariance after the constrained update.
        real(dp), allocatable, intent(out) :: fixed_block(:, :) !! Allocated updated lower-right constrained block.
        logical, intent(out) :: accepted !! True when the lower-right correlation proposal is accepted.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid input or a numerical failure.
        real(dp), allocatable :: posterior_precision(:, :)
        integer :: n
        integer :: n_fixed

        info = 0
        accepted = .false.
        n = size(posterior_sum, 1)
        n_fixed = n - split
        if (size(posterior_sum, 2) /= n .or. size(prior_scale, 1) /= n .or. size(prior_scale, 2) /= n .or. &
            split < 1 .or. split >= n .or. size(old_fixed_block, 1) /= n_fixed .or. &
            size(old_fixed_block, 2) /= n_fixed .or. sample_degrees_freedom <= real(n - 1, dp)) then
            allocate(covariance(0, 0), fixed_block(0, 0))
            info = 1
            return
        end if

        call correlation_structure_update(state, posterior_sum(split + 1:n, split + 1:n), &
            sample_degrees_freedom - real(split, dp), prior_degrees_freedom, &
            prior_scale(split + 1:n, split + 1:n), old_fixed_block, fixed_block, accepted, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call inverse_matrix(posterior_sum, posterior_precision, info)
        if (info /= 0) then
            allocate(covariance(0, 0))
            return
        end if
        call conditioned_inverse_wishart_precision(state, posterior_precision, sample_degrees_freedom, split, &
            fixed_block, covariance, info)
    end subroutine correlation_submatrix_update


    pure subroutine covariance_update_dispatch(state, update_code, posterior_sum, sample_degrees_freedom, &
                                               prior_degrees_freedom, prior_scale, old_covariance, split, &
                                               fixed_block_in, covariance, fixed_block_out, accepted, info, &
                                               ante_location, ante_lag_order, ante_common_beta, &
                                               ante_common_variance, ante_a_inverse)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the requested stochastic covariance update.
        integer, intent(in) :: update_code !! MCMCglmm covariance update code; code 5 uses the antedependence sampler.
        real(dp), intent(in) :: posterior_sum(:, :) !! Data cross-product plus prior scale matrix for this covariance block.
        real(dp), intent(in) :: sample_degrees_freedom !! Effective data degrees of freedom, normally the number of levels.
        real(dp), intent(in) :: prior_degrees_freedom !! Prior degrees of freedom supplied by the covariance prior.
        real(dp), intent(in) :: prior_scale(:, :) !! Prior scale matrix, including fixed marginal variances for code 3 or 4.
        real(dp), intent(in) :: old_covariance(:, :) !! Current covariance, retained by code 0 and used by correlation MH updates.
        integer, intent(in) :: split !! Leading dimension for codes 2, 4, and 6; ignored by codes 0, 1, and 3.
        real(dp), intent(in) :: fixed_block_in(:, :) !! Fixed lower-right block for codes 2/4; may be zero-size otherwise.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated covariance after applying the selected update.
        real(dp), allocatable, intent(out) :: fixed_block_out(:, :) !! Updated constrained block for code 4, or a copy for code 2.
        logical, intent(out) :: accepted !! Correlation-MH acceptance flag for codes 3/4; true for direct Gibbs/fixed updates.
        integer, intent(out) :: info !! Zero on success; code 5 needs antedependence locations and a valid lag order.
        real(dp), optional, intent(in) :: ante_location(:, :) !! Level-by-trait locations used by antedependence code 5.
        integer, optional, intent(in) :: ante_lag_order !! Positive antedependence lag order required for update code 5.
        logical, optional, intent(in) :: ante_common_beta !! If true, share one antedependence coefficient per lag.
        logical, optional, intent(in) :: ante_common_variance !! If true, share one innovation variance across traits.
        real(dp), optional, intent(in) :: ante_a_inverse(:, :) !! Optional level precision matrix for antedependence code 5.
        real(dp), allocatable :: ante_samples(:, :, :)
        logical :: common_beta_value
        logical :: common_variance_value
        integer :: n
        real(dp) :: posterior_degrees_freedom

        info = 0
        accepted = .true.
        n = size(posterior_sum, 1)
        if (n < 1 .or. size(posterior_sum, 2) /= n .or. size(prior_scale, 1) /= n .or. &
            size(prior_scale, 2) /= n .or. size(old_covariance, 1) /= n .or. size(old_covariance, 2) /= n .or. &
            sample_degrees_freedom < 0.0_dp .or. prior_degrees_freedom < 0.0_dp) then
            allocate(covariance(0, 0), fixed_block_out(0, 0))
            info = 1
            return
        end if
        posterior_degrees_freedom = sample_degrees_freedom + prior_degrees_freedom

        select case (update_code)
        case (0)
            covariance = old_covariance
            allocate(fixed_block_out(0, 0))
        case (1)
            call inverse_wishart_sample(state, posterior_sum, posterior_degrees_freedom, covariance, info)
            allocate(fixed_block_out(0, 0))
        case (2)
            if (split < 0 .or. split >= n .or. size(fixed_block_in, 1) /= n - split .or. &
                size(fixed_block_in, 2) /= n - split) then
                allocate(covariance(0, 0), fixed_block_out(0, 0))
                info = 2
                return
            end if
            if (split == 0) then
                covariance = fixed_block_in
            else
                call conditioned_covariance_update(state, posterior_sum, posterior_degrees_freedom, split + 1, &
                    fixed_block_in, covariance, info)
            end if
            fixed_block_out = fixed_block_in
        case (3)
            call correlation_structure_update(state, posterior_sum, sample_degrees_freedom, prior_degrees_freedom, &
                prior_scale, old_covariance, covariance, accepted, info)
            allocate(fixed_block_out(0, 0))
        case (4)
            if (split < 1 .or. split >= n .or. size(fixed_block_in, 1) /= n - split .or. &
                size(fixed_block_in, 2) /= n - split) then
                allocate(covariance(0, 0), fixed_block_out(0, 0))
                info = 2
                return
            end if
            call correlation_submatrix_update(state, posterior_sum, sample_degrees_freedom, prior_degrees_freedom, &
                prior_scale, split, fixed_block_in, covariance, fixed_block_out, accepted, info)
        case (5)
            if (.not. present(ante_location) .or. .not. present(ante_lag_order)) then
                allocate(covariance(0, 0), fixed_block_out(0, 0))
                accepted = .false.
                info = 5
                return
            end if
            common_beta_value = .false.
            if (present(ante_common_beta)) common_beta_value = ante_common_beta
            common_variance_value = .false.
            if (present(ante_common_variance)) common_variance_value = ante_common_variance
            if (present(ante_a_inverse)) then
                call ante_covariance_samples(state, ante_location, ante_lag_order, 1, common_beta_value, &
                    common_variance_value, ante_samples, info, ante_a_inverse)
            else
                call ante_covariance_samples(state, ante_location, ante_lag_order, 1, common_beta_value, &
                    common_variance_value, ante_samples, info)
            end if
            if (info /= 0) then
                allocate(covariance(0, 0), fixed_block_out(0, 0))
                accepted = .false.
                return
            end if
            covariance = ante_samples(:, :, 1)
            allocate(fixed_block_out(0, 0))
        case (6)
            if (split < 1 .or. split > n) then
                allocate(covariance(0, 0), fixed_block_out(0, 0))
                info = 2
                return
            end if
            call identity_direct_sum_update(state, posterior_sum, posterior_degrees_freedom, split, covariance, info)
            allocate(fixed_block_out(0, 0))
        case default
            allocate(covariance(0, 0), fixed_block_out(0, 0))
            accepted = .false.
            info = 4
        end select
    end subroutine covariance_update_dispatch

end module mcmcglmm_covariance
