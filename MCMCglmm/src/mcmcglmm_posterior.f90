! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_posterior
    use r_kinds, only : dp
    use r_linalg, only : cholesky_factor, inverse_matrix, symmetric_eigen
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mcmcglmm_matrix, only : covariance_to_correlation
    implicit none
    private

    public :: posterior_correlations
    public :: posterior_inverses
    public :: posterior_eigenvalues
    public :: posterior_modes
    public :: ante_parameters

contains


    pure subroutine posterior_modes(samples, adjust, modes, info)
        real(dp), intent(in) :: samples(:, :) !! Posterior draws stored as iterations by variables, matching as.matrix(mcmc).
        real(dp), optional, intent(in) :: adjust !! Optional positive bw.nrd0 multiplier; defaults to MCMCglmm posterior.mode's 0.1.
        real(dp), allocatable, intent(out) :: modes(:) !! One KDE grid mode for each posterior variable column.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid sample size, adjustment, or bandwidth.
        real(dp), allocatable :: sorted(:)
        real(dp) :: adjust_value
        real(dp) :: bandwidth
        real(dp) :: density_value
        real(dp) :: from_value
        real(dp) :: grid_value
        real(dp) :: highest_density
        real(dp) :: hi
        real(dp) :: iqr
        real(dp) :: lo
        real(dp) :: mean_value
        real(dp) :: step
        real(dp) :: to_value
        integer, parameter :: grid_points = 512
        integer :: i
        integer :: j
        integer :: n
        integer :: variable

        info = 0
        n = size(samples, 1)
        adjust_value = 0.1_dp
        if (present(adjust)) adjust_value = adjust
        if (n < 2 .or. size(samples, 2) < 1 .or. adjust_value <= 0.0_dp) then
            allocate(modes(0))
            info = 1
            return
        end if
        allocate(modes(size(samples, 2)), sorted(n))
        do variable = 1, size(samples, 2)
            sorted = samples(:, variable)
            if (any(.not. ieee_is_finite(sorted))) then
                modes = 0.0_dp
                info = 2
                return
            end if
            call sort_increasing(sorted)
            mean_value = sum(sorted) / real(n, dp)
            hi = sqrt(sum((sorted - mean_value) ** 2) / real(n - 1, dp))
            iqr = quantile_type7(sorted, 0.75_dp) - quantile_type7(sorted, 0.25_dp)
            lo = min(hi, iqr / 1.34_dp)
            if (lo <= 0.0_dp) lo = hi
            if (lo <= 0.0_dp) lo = abs(sorted(1))
            if (lo <= 0.0_dp) lo = 1.0_dp
            bandwidth = adjust_value * 0.9_dp * lo * real(n, dp) ** (-0.2_dp)
            if (.not. ieee_is_finite(bandwidth) .or. bandwidth <= 0.0_dp) then
                modes = 0.0_dp
                info = 2
                return
            end if
            from_value = sorted(1) - 3.0_dp * bandwidth
            to_value = sorted(n) + 3.0_dp * bandwidth
            step = (to_value - from_value) / real(grid_points - 1, dp)
            highest_density = -1.0_dp
            modes(variable) = from_value
            do i = 1, grid_points
                grid_value = from_value + real(i - 1, dp) * step
                density_value = 0.0_dp
                do j = 1, n
                    density_value = density_value + exp(-0.5_dp * ((grid_value - sorted(j)) / bandwidth) ** 2)
                end do
                if (density_value > highest_density) then
                    highest_density = density_value
                    modes(variable) = grid_value
                end if
            end do
        end do
    end subroutine posterior_modes

    pure subroutine sort_increasing(values)
        real(dp), intent(inout) :: values(:) !! Finite real values sorted increasingly in place by insertion sort.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_increasing

    pure real(dp) function quantile_type7(sorted, probability) result(value)
        real(dp), intent(in) :: sorted(:) !! Increasing finite values used for R type-7 interpolation.
        real(dp), intent(in) :: probability !! Quantile probability in the closed interval [0,1].
        real(dp) :: fraction
        real(dp) :: index_value
        integer :: lower_index

        if (probability <= 0.0_dp) then
            value = sorted(1)
            return
        end if
        if (probability >= 1.0_dp) then
            value = sorted(size(sorted))
            return
        end if
        index_value = 1.0_dp + real(size(sorted) - 1, dp) * probability
        lower_index = int(floor(index_value))
        fraction = index_value - real(lower_index, dp)
        if (lower_index >= size(sorted)) then
            value = sorted(size(sorted))
        else
            value = (1.0_dp - fraction) * sorted(lower_index) + fraction * sorted(lower_index + 1)
        end if
    end function quantile_type7

    pure subroutine posterior_correlations(samples, correlations, info)
        real(dp), intent(in) :: samples(:, :, :) !! Sequence of square covariance matrices stored as n by n by n_sample.
        real(dp), allocatable, intent(out) :: correlations(:, :, :) !! Covariance-to-correlation transform for every sample.
        integer, intent(out) :: info !! Zero on success; otherwise the first covariance conversion failure code.
        real(dp), allocatable :: current(:, :)
        integer :: k
        integer :: n

        info = 0
        n = size(samples, 1)
        if (size(samples, 2) /= n) then
            allocate(correlations(0, 0, 0))
            info = 1
            return
        end if
        allocate(correlations(n, n, size(samples, 3)))
        do k = 1, size(samples, 3)
            call covariance_to_correlation(samples(:, :, k), current, info)
            if (info /= 0) then
                correlations = 0.0_dp
                return
            end if
            correlations(:, :, k) = current
        end do
    end subroutine posterior_correlations

    pure subroutine posterior_inverses(samples, inverses, info)
        real(dp), intent(in) :: samples(:, :, :) !! Sequence of nonsingular square matrices stored as n by n by n_sample.
        real(dp), allocatable, intent(out) :: inverses(:, :, :) !! Allocated inverse matrix for every sample plane.
        integer, intent(out) :: info !! Zero on success; otherwise the first inversion failure code.
        real(dp), allocatable :: current(:, :)
        integer :: k
        integer :: n

        info = 0
        n = size(samples, 1)
        if (size(samples, 2) /= n) then
            allocate(inverses(0, 0, 0))
            info = 1
            return
        end if
        allocate(inverses(n, n, size(samples, 3)))
        do k = 1, size(samples, 3)
            call inverse_matrix(samples(:, :, k), current, info)
            if (info /= 0) then
                inverses = 0.0_dp
                return
            end if
            inverses(:, :, k) = current
        end do
    end subroutine posterior_inverses

    subroutine posterior_eigenvalues(samples, eigenvalues, info)
        real(dp), intent(in) :: samples(:, :, :) !! Sequence of symmetric square matrices stored as n by n by n_sample.
        real(dp), allocatable, intent(out) :: eigenvalues(:, :) !! n by sample eigenvalues ordered largest to smallest.
        integer, intent(out) :: info !! Zero on success; otherwise the first symmetric eigensolver failure code.
        real(dp), allocatable :: current_values(:)
        real(dp), allocatable :: vectors(:, :)
        integer :: k
        integer :: n

        info = 0
        n = size(samples, 1)
        if (size(samples, 2) /= n) then
            allocate(eigenvalues(0, 0))
            info = 1
            return
        end if
        allocate(eigenvalues(n, size(samples, 3)))
        do k = 1, size(samples, 3)
            call symmetric_eigen(samples(:, :, k), current_values, vectors, info, descending=.true.)
            if (info /= 0) then
                eigenvalues = 0.0_dp
                return
            end if
            eigenvalues(:, k) = current_values
        end do
    end subroutine posterior_eigenvalues

    pure subroutine ante_parameters(covariance, lag_order, innovation_variances, coefficients, info)
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite covariance matrix V.
        integer, intent(in) :: lag_order !! Positive antedependence order k, strictly less than matrix dimension.
        real(dp), allocatable, intent(out) :: innovation_variances(:) !! Squared Cholesky diagonals used by posterior.ante.
        real(dp), allocatable, intent(out) :: coefficients(:, :) !! Coefficients by lag/start coordinate; unused cells are zero.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape/order or Cholesky/inversion failure.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: diagonal_matrix(:, :)
        real(dp), allocatable :: inverse_factor(:, :)
        real(dp), allocatable :: upper_factor(:, :)
        integer :: i
        integer :: lag
        integer :: n

        info = 0
        n = size(covariance, 1)
        if (size(covariance, 2) /= n .or. lag_order < 1 .or. lag_order >= n) then
            allocate(innovation_variances(0), coefficients(0, 0))
            info = 1
            return
        end if
        call cholesky_factor(covariance, upper_factor, info, upper=.true.)
        if (info /= 0) then
            allocate(innovation_variances(0), coefficients(0, 0))
            return
        end if
        call inverse_matrix(upper_factor, inverse_factor, info)
        if (info /= 0) then
            allocate(innovation_variances(0), coefficients(0, 0))
            return
        end if
        allocate(diagonal_matrix(n, n), innovation_variances(n), coefficients(lag_order, n - 1))
        diagonal_matrix = 0.0_dp
        coefficients = 0.0_dp
        do i = 1, n
            diagonal_matrix(i, i) = upper_factor(i, i)
            innovation_variances(i) = upper_factor(i, i) ** 2
        end do
        beta = transpose(matmul(inverse_factor, diagonal_matrix))
        do lag = 1, lag_order
            do i = 1, n - lag
                coefficients(lag, i) = -beta(i + lag, i)
            end do
        end do
    end subroutine ante_parameters

end module mcmcglmm_posterior
