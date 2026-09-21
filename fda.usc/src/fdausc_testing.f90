module fdausc_testing
    use r_kinds, only : dp
    use fdausc_kernels, only : ker_norm
    use fdausc_integration, only : integrate_curve
    implicit none
    private

    public :: flm_f_statistic, dfv_statistics, mmd_rbf_statistic
    public :: distribution_equality_statistic, mean_difference_statistic, covariance_hs_statistic
    public :: fanova_mean_distance_statistic

contains

    pure real(dp) function flm_f_statistic(t, curves, y) result(statistic)
        real(dp), intent(in) :: t(:) !! Increasing functional grid shared by all predictor curves.
        real(dp), intent(in) :: curves(:, :) !! Functional predictors, one curve per response observation.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with curve rows.
        real(dp), allocatable :: mean_curve(:)
        real(dp), allocatable :: weighted_mean(:)
        real(dp) :: ymean
        integer :: i

        allocate(mean_curve(size(curves, 2)), weighted_mean(size(curves, 2)))
        mean_curve = sum(curves, dim=1)/real(size(curves, 1), dp)
        ymean = sum(y)/real(size(y), dp)
        weighted_mean = 0.0_dp
        do i = 1, size(curves, 1)
            weighted_mean = weighted_mean + (curves(i, :) - mean_curve)*(y(i) - ymean)
        end do
        weighted_mean = weighted_mean/real(size(curves, 1), dp)
        statistic = sqrt(max(0.0_dp, integrate_curve(t, weighted_mean*weighted_mean, 2)))
    end function flm_f_statistic

    pure subroutine dfv_statistics(distances, y, bandwidths, statistics, weights)
        real(dp), intent(in) :: distances(:, :) !! Square functional distance matrix aligned with response y.
        real(dp), intent(in) :: y(:) !! Scalar response vector entering the Delsol-Ferraty-Vieu statistic.
        real(dp), intent(in) :: bandwidths(:) !! Positive bandwidth values at which the statistic is evaluated.
        real(dp), intent(out) :: statistics(:) !! DFV statistic for every bandwidth.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights, default one.
        real(dp) :: fitted
        real(dp) :: term
        integer :: h
        integer :: i
        integer :: j

        do h = 1, size(bandwidths)
            statistics(h) = 0.0_dp
            do i = 1, size(y)
                fitted = 0.0_dp
                do j = 1, size(y)
                    fitted = fitted + 2.0_dp*ker_norm(abs(distances(i, j)/bandwidths(h)))*y(j)
                end do
                term = fitted*fitted
                if (present(weights)) term = term*weights(i)
                statistics(h) = statistics(h) + term
            end do
            statistics(h) = statistics(h)/real(size(y), dp)
        end do
    end subroutine dfv_statistics

    pure real(dp) function mmd_rbf_statistic(dx, dy, dxy, bandwidth, scale_by_n) result(statistic)
        real(dp), intent(in) :: dx(:, :) !! Pairwise distance matrix within the first sample.
        real(dp), intent(in) :: dy(:, :) !! Pairwise distance matrix within the second sample.
        real(dp), intent(in) :: dxy(:, :) !! Distances from first-sample rows to second-sample columns.
        real(dp), intent(in) :: bandwidth !! Positive RBF-kernel bandwidth.
        logical, intent(in), optional :: scale_by_n !! If true multiply MMD squared by first-sample size as in MMDA.test.
        real(dp) :: kx
        real(dp) :: ky
        real(dp) :: kxy
        logical :: scaled

        kx = sum(exp(-0.5_dp*(dx/bandwidth)**2))/real(size(dx), dp)
        ky = sum(exp(-0.5_dp*(dy/bandwidth)**2))/real(size(dy), dp)
        kxy = sum(exp(-0.5_dp*(dxy/bandwidth)**2))/real(size(dxy), dp)
        statistic = kx + ky - 2.0_dp*kxy
        scaled = .false.
        if (present(scale_by_n)) scaled = scale_by_n
        if (scaled) statistic = real(size(dx, 1), dp)*statistic
    end function mmd_rbf_statistic

    pure subroutine distribution_equality_statistic(dx, dy, dxy, statistic, studentized)
        real(dp), intent(in) :: dx(:, :) !! Pairwise distance matrix within the first sample.
        real(dp), intent(in) :: dy(:, :) !! Pairwise distance matrix within the second sample.
        real(dp), intent(in) :: dxy(:, :) !! Cross-sample distance matrix.
        real(dp), intent(out) :: statistic !! Exchangeable energy-style equality-of-distribution statistic Tn.
        real(dp), intent(out) :: studentized !! Studentized statistic Tn divided by pooled within-sample scale.
        integer :: n
        integer :: m
        real(dp) :: vx
        real(dp) :: vy
        real(dp) :: denom

        n = size(dx, 1)
        m = size(dy, 1)
        statistic = real(n, dp)*(2.0_dp*sum(dxy)/real(n*m, dp) &
            - sum(dx)/real(n*n, dp) - sum(dy)/real(m*m, dp))
        vx = sum(dx)/(2.0_dp*real(n*n, dp))
        vy = sum(dy)/(2.0_dp*real(m*m, dp))
        denom = vx + real(n, dp)/real(m, dp)*vy
        if (abs(denom) <= tiny(1.0_dp)) then
            studentized = 0.0_dp
        else
            studentized = statistic/denom
        end if
    end subroutine distribution_equality_statistic

    pure real(dp) function mean_difference_statistic(t, x, y) result(statistic)
        real(dp), intent(in) :: t(:) !! Increasing functional grid shared by the two samples.
        real(dp), intent(in) :: x(:, :) !! First functional sample, curves in rows.
        real(dp), intent(in) :: y(:, :) !! Second functional sample, curves in rows.
        real(dp), allocatable :: mx(:)
        real(dp), allocatable :: my(:)
        integer :: n
        integer :: m

        n = size(x, 1)
        m = size(y, 1)
        mx = sum(x, dim=1)/real(n, dp)
        my = sum(y, dim=1)/real(m, dp)
        statistic = real(n*m, dp)/real(n + m, dp)*integrate_curve(t, (mx - my)**2, 2)
    end function mean_difference_statistic

    pure real(dp) function covariance_hs_statistic(x, y) result(statistic)
        real(dp), intent(in) :: x(:, :) !! First functional sample on a common discretization grid, curves in rows.
        real(dp), intent(in) :: y(:, :) !! Second functional sample on the same grid, curves in rows.
        real(dp), allocatable :: mx(:)
        real(dp), allocatable :: my(:)
        real(dp), allocatable :: sx(:, :)
        real(dp), allocatable :: sy(:, :)
        real(dp), allocatable :: d(:)
        integer :: i
        integer :: j
        integer :: k

        mx = sum(x, dim=1)/real(size(x, 1), dp)
        my = sum(y, dim=1)/real(size(y, 1), dp)
        allocate(sx(size(x, 2), size(x, 2)), sy(size(y, 2), size(y, 2)), d(size(x, 2)))
        sx = 0.0_dp
        sy = 0.0_dp
        do i = 1, size(x, 1)
            d = x(i, :) - mx
            do j = 1, size(d)
                do k = 1, size(d)
                    sx(j, k) = sx(j, k) + d(j)*d(k)
                end do
            end do
        end do
        do i = 1, size(y, 1)
            d = y(i, :) - my
            do j = 1, size(d)
                do k = 1, size(d)
                    sy(j, k) = sy(j, k) + d(j)*d(k)
                end do
            end do
        end do
        sx = sx/real(max(1, size(x, 1) - 1), dp)
        sy = sy/real(max(1, size(y, 1) - 1), dp)
        statistic = sum((sx - sy)**2)
    end function covariance_hs_statistic

    pure real(dp) function fanova_mean_distance_statistic(t, curves, groups) result(statistic)
        real(dp), intent(in) :: t(:) !! Increasing functional grid shared by all curves.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        integer, intent(in) :: groups(:) !! Positive group labels aligned with curve rows.
        real(dp), allocatable :: means(:, :)
        real(dp), allocatable :: diff(:)
        integer, allocatable :: counts(:)
        integer :: gmax
        integer :: g
        integer :: h
        integer :: i

        gmax = maxval(groups)
        allocate(means(gmax, size(curves, 2)), counts(gmax), diff(size(curves, 2)))
        means = 0.0_dp
        counts = 0
        do i = 1, size(curves, 1)
            if (groups(i) >= 1 .and. groups(i) <= gmax) then
                means(groups(i), :) = means(groups(i), :) + curves(i, :)
                counts(groups(i)) = counts(groups(i)) + 1
            end if
        end do
        do g = 1, gmax
            if (counts(g) > 0) means(g, :) = means(g, :)/real(counts(g), dp)
        end do
        statistic = 0.0_dp
        do g = 1, gmax - 1
            if (counts(g) == 0) cycle
            do h = g + 1, gmax
                if (counts(h) == 0) cycle
                diff = means(g, :) - means(h, :)
                statistic = statistic + integrate_curve(t, diff*diff, 2)
            end do
        end do
    end function fanova_mean_distance_statistic

end module fdausc_testing
