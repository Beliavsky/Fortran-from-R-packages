! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_ante
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state, rng_chisq
    use mcmcglmm_matrix, only : sample_mvn_precision
    implicit none
    private

    public :: ante_covariance_samples

contains

    pure subroutine ante_design(location, lag_order, common_beta, design, info)
        real(dp), intent(in) :: location(:, :) !! Level-by-trait matrix used as the antedependence response and lagged predictors.
        integer, intent(in) :: lag_order !! Positive antedependence order strictly less than the number of traits.
        logical, intent(in) :: common_beta !! True for one coefficient per lag; false for trait-specific lag coefficients.
        real(dp), allocatable, intent(out) :: design(:, :) !! Trait-major antedependence design with one row per level/trait.
        integer, intent(out) :: info !! Zero on success; nonzero for an invalid lag order or empty location matrix.
        integer :: column
        integer :: dim_g
        integer :: j
        integer :: k
        integer :: levels
        integer :: nbeta
        integer :: row_start

        info = 0
        levels = size(location, 1)
        dim_g = size(location, 2)
        if (levels < 1 .or. dim_g < 2 .or. lag_order < 1 .or. lag_order >= dim_g) then
            allocate(design(0, 0))
            info = 1
            return
        end if
        if (common_beta) then
            nbeta = lag_order
        else
            nbeta = dim_g * lag_order - lag_order * (lag_order + 1) / 2
        end if
        allocate(design(levels * dim_g, nbeta))
        design = 0.0_dp
        column = 0
        if (common_beta) then
            do k = 1, lag_order
                column = k
                do j = k + 1, dim_g
                    row_start = (j - 1) * levels
                    design(row_start + 1:row_start + levels, column) = location(:, j - k)
                end do
            end do
        else
            do k = 1, lag_order
                do j = k + 1, dim_g
                    column = column + 1
                    row_start = (j - 1) * levels
                    design(row_start + 1:row_start + levels, column) = location(:, j - k)
                end do
            end do
        end if
    end subroutine ante_design

    pure subroutine ante_residual_crossproduct(residual, a_inverse, crossproduct, info)
        real(dp), intent(in) :: residual(:, :) !! Level-by-trait antedependence residual matrix.
        real(dp), intent(in), optional :: a_inverse(:, :) !! Optional level precision such as A^{-1}; identity when absent.
        real(dp), allocatable, intent(out) :: crossproduct(:, :) !! Allocated trait cross-product residual^T A^{-1} residual.
        integer, intent(out) :: info !! Zero on success; nonzero when an optional row precision has incompatible dimensions.

        info = 0
        if (present(a_inverse)) then
            if (size(a_inverse, 1) /= size(residual, 1) .or. size(a_inverse, 2) /= size(residual, 1)) then
                allocate(crossproduct(0, 0))
                info = 1
                return
            end if
            crossproduct = matmul(transpose(residual), matmul(a_inverse, residual))
        else
            crossproduct = matmul(transpose(residual), residual)
        end if
    end subroutine ante_residual_crossproduct

    pure subroutine ante_precision(ivar, levels, a_inverse, precision, info)
        real(dp), intent(in) :: ivar(:) !! Positive innovation variances, one per trait.
        integer, intent(in) :: levels !! Number of levels/rows per trait block.
        real(dp), intent(in), optional :: a_inverse(:, :) !! Optional level precision matrix; identity is used when absent.
        real(dp), allocatable, intent(out) :: precision(:, :) !! Block residual precision diag(1/ivar) tensor A^{-1}.
        integer, intent(out) :: info !! Zero on success; nonzero for nonpositive variance or incompatible A^{-1} dimensions.
        integer :: i
        integer :: row_start

        info = 0
        if (any(ivar <= 0.0_dp)) then
            allocate(precision(0, 0))
            info = 1
            return
        end if
        if (present(a_inverse)) then
            if (size(a_inverse, 1) /= levels .or. size(a_inverse, 2) /= levels) then
                allocate(precision(0, 0))
                info = 2
                return
            end if
        end if
        allocate(precision(levels * size(ivar), levels * size(ivar)))
        precision = 0.0_dp
        do i = 1, size(ivar)
            row_start = (i - 1) * levels
            if (present(a_inverse)) then
                precision(row_start + 1:row_start + levels, row_start + 1:row_start + levels) = a_inverse / ivar(i)
            else
                precision(row_start + 1:row_start + levels, row_start + 1:row_start + levels) = 0.0_dp
                precision(row_start + 1:row_start + levels, row_start + 1:row_start + levels) = &
                    identity_matrix(levels) / ivar(i)
            end if
        end do
    end subroutine ante_precision

    pure function identity_matrix(n) result(matrix_value)
        integer, intent(in) :: n !! Nonnegative order of the identity matrix.
        real(dp) :: matrix_value(n, n)
        integer :: i

        matrix_value = 0.0_dp
        do i = 1, n
            matrix_value(i, i) = 1.0_dp
        end do
    end function identity_matrix

    pure subroutine build_ante_cholesky_inverse(beta, ivar, lag_order, common_beta, factor, info)
        real(dp), intent(in) :: beta(:) !! Sampled antedependence regression coefficients in lag-major ordering.
        real(dp), intent(in) :: ivar(:) !! Positive innovation variances by trait.
        integer, intent(in) :: lag_order !! Antedependence lag order.
        logical, intent(in) :: common_beta !! True for one beta per lag; false for trait-specific lag coefficients.
        real(dp), allocatable, intent(out) :: factor(:, :) !! Lower triangular factor L satisfying G^{-1} = L^T L.
        integer, intent(out) :: info !! Zero on success; nonzero for inconsistent beta length or nonpositive innovation variances.
        integer :: count
        integer :: dim_g
        integer :: i
        integer :: j

        info = 0
        dim_g = size(ivar)
        if (any(ivar <= 0.0_dp)) then
            allocate(factor(0, 0))
            info = 1
            return
        end if
        if (common_beta) then
            if (size(beta) /= lag_order) then
                allocate(factor(0, 0))
                info = 2
                return
            end if
        else
            if (size(beta) /= dim_g * lag_order - lag_order * (lag_order + 1) / 2) then
                allocate(factor(0, 0))
                info = 3
                return
            end if
        end if
        allocate(factor(dim_g, dim_g))
        factor = 0.0_dp
        do i = 1, dim_g
            factor(i, i) = 1.0_dp / sqrt(ivar(i))
        end do
        if (common_beta) then
            do j = 1, lag_order
                do i = 1, dim_g - j
                    factor(i + j, i) = -beta(j) / sqrt(ivar(i + j))
                end do
            end do
        else
            count = 0
            do j = 1, lag_order
                do i = 1, dim_g - j
                    count = count + 1
                    factor(i + j, i) = -beta(count) / sqrt(ivar(i + j))
                end do
            end do
        end if
    end subroutine build_ante_cholesky_inverse

    pure subroutine ante_covariance_samples(state, location, lag_order, n_sample, common_beta, common_variance, &
                                            samples, info, a_inverse)
        type(rng_state), intent(inout) :: state !! Generator state consumed by antedependence coefficient and variance sampling.
        real(dp), intent(in) :: location(:, :) !! Level-by-trait matrix corresponding to rante's y argument.
        integer, intent(in) :: lag_order !! Positive antedependence order less than the number of traits.
        integer, intent(in) :: n_sample !! Positive number of covariance matrices to sample.
        logical, intent(in) :: common_beta !! Match rante(cbeta=TRUE) when true; otherwise use trait-specific lag coefficients.
        logical, intent(in) :: common_variance !! True matches rante(cvar=TRUE); false samples variance by trait.
        real(dp), allocatable, intent(out) :: samples(:, :, :) !! Allocated trait covariance samples, dim by dim by n_sample.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or sampling/matrix failure.
        real(dp), intent(in), optional :: a_inverse(:, :) !! Optional level precision matrix corresponding to rante's Ainv argument.
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: beta_mean(:)
        real(dp), allocatable :: beta_prior_precision(:, :)
        real(dp), allocatable :: crossproduct(:, :)
        real(dp), allocatable :: design(:, :)
        real(dp), allocatable :: factor(:, :)
        real(dp), allocatable :: g_inverse(:, :)
        real(dp), allocatable :: g_matrix(:, :)
        real(dp), allocatable :: ivar(:)
        real(dp), allocatable :: location_vector(:)
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: residual_precision(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp) :: chi
        integer :: dim_g
        integer :: i
        integer :: levels
        integer :: nbeta
        integer :: sample_index

        info = 0
        levels = size(location, 1)
        dim_g = size(location, 2)
        if (levels < 1 .or. dim_g < 2 .or. n_sample < 1 .or. lag_order < 1 .or. lag_order >= dim_g) then
            allocate(samples(0, 0, 0))
            info = 1
            return
        end if
        if (present(a_inverse)) then
            if (size(a_inverse, 1) /= levels .or. size(a_inverse, 2) /= levels) then
                allocate(samples(0, 0, 0))
                info = 2
                return
            end if
        end if
        call ante_design(location, lag_order, common_beta, design, info)
        if (info /= 0) then
            allocate(samples(0, 0, 0))
            return
        end if
        nbeta = size(design, 2)
        allocate(samples(dim_g, dim_g, n_sample), beta_prior_precision(nbeta, nbeta), ivar(dim_g))
        allocate(beta_mean(nbeta))
        beta_prior_precision = 0.0_dp
        beta_mean = 0.0_dp
        do i = 1, nbeta
            beta_prior_precision(i, i) = 1.0e-9_dp
        end do
        ivar = 1.0_dp
        location_vector = reshape(location, [levels * dim_g])

        do sample_index = 1, n_sample
            if (present(a_inverse)) then
                call ante_precision(ivar, levels, a_inverse, residual_precision, info)
            else
                call ante_precision(ivar, levels, precision=residual_precision, info=info)
            end if
            if (info /= 0) return
            precision = matmul(transpose(design), matmul(residual_precision, design)) + beta_prior_precision
            rhs = matmul(transpose(design), matmul(residual_precision, location_vector))
            call sample_mvn_precision(state, rhs, precision, beta, beta_mean, info)
            if (info /= 0) return
            residual = location - reshape(matmul(design, beta), [levels, dim_g])
            if (present(a_inverse)) then
                call ante_residual_crossproduct(residual, a_inverse, crossproduct, info)
            else
                call ante_residual_crossproduct(residual, crossproduct=crossproduct, info=info)
            end if
            if (info /= 0) return
            if (common_variance) then
                call rng_chisq(state, real(levels * dim_g, dp), chi)
                ivar = sum([(crossproduct(i, i), i=1, dim_g)]) / chi
            else
                do i = 1, dim_g
                    call rng_chisq(state, real(levels, dp), chi)
                    ivar(i) = crossproduct(i, i) / chi
                end do
            end if
            call build_ante_cholesky_inverse(beta, ivar, lag_order, common_beta, factor, info)
            if (info /= 0) return
            g_inverse = matmul(transpose(factor), factor)
            call inverse_matrix(g_inverse, g_matrix, info)
            if (info /= 0) return
            samples(:, :, sample_index) = g_matrix
        end do
    end subroutine ante_covariance_samples

end module mcmcglmm_ante
