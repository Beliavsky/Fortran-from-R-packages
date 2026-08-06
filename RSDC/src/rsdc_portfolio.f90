! SPDX-License-Identifier: GPL-3.0-only
module rsdc_portfolio
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_portfolio_result
    use rsdc_parameters, only: correlations_to_matrix
    use rsdc_linalg, only: inverse_spd, is_positive_definite, standard_deviation
    implicit none
    private
    public :: rsdc_minvar, rsdc_maxdiv

contains

    subroutine rsdc_minvar(sigma, predicted_corr, returns_matrix, result, long_only, lagged, ok)
        real(dp), intent(in) :: sigma(:, :), predicted_corr(:, :), returns_matrix(:, :)
        type(rsdc_portfolio_result), intent(out) :: result
        logical, intent(in), optional :: long_only, lagged
        logical, intent(out), optional :: ok
        logical :: lo, lag, good
        real(dp), allocatable :: covariance(:, :)
        integer :: t, nobs, k
        lo = .true.; if (present(long_only)) lo = long_only
        lag = .false.; if (present(lagged)) lag = lagged
        nobs = size(sigma, 1); k = size(sigma, 2)
        good = size(predicted_corr, 1) == nobs .and. size(predicted_corr, 2) == k * (k - 1) / 2
        good = good .and. all(shape(returns_matrix) == [nobs, k])
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(result%weights(nobs, k), result%returns(nobs), covariance(k, k))
        result%n_fallback = 0
        do t = 1, nobs
            call build_covariance(sigma(t, :), predicted_corr(t, :), covariance, good)
            if (good) call minimum_variance_weights(covariance, result%weights(t, :), lo, good)
            if (.not. good) then
                result%weights(t, :) = 1.0_dp / real(k, dp)
                result%n_fallback = result%n_fallback + 1
                good = .true.
            end if
        end do
        call portfolio_returns(returns_matrix, result%weights, lag, result%returns)
        call summarize_returns(result)
        if (present(ok)) ok = .true.
    end subroutine rsdc_minvar

    subroutine rsdc_maxdiv(sigma, predicted_corr, returns_matrix, result, long_only, lagged, ok)
        real(dp), intent(in) :: sigma(:, :), predicted_corr(:, :), returns_matrix(:, :)
        type(rsdc_portfolio_result), intent(out) :: result
        logical, intent(in), optional :: long_only, lagged
        logical, intent(out), optional :: ok
        logical :: lo, lag, good
        real(dp), allocatable :: covariance(:, :)
        integer :: t, nobs, k
        lo = .true.; if (present(long_only)) lo = long_only
        lag = .false.; if (present(lagged)) lag = lagged
        nobs = size(sigma, 1); k = size(sigma, 2)
        good = size(predicted_corr, 1) == nobs .and. size(predicted_corr, 2) == k * (k - 1) / 2
        good = good .and. all(shape(returns_matrix) == [nobs, k])
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(result%weights(nobs, k), result%returns(nobs), result%diversification_ratios(nobs), covariance(k, k))
        result%n_fallback = 0
        do t = 1, nobs
            call build_covariance(sigma(t, :), predicted_corr(t, :), covariance, good)
            if (good) call maximum_diversification_weights(covariance, sigma(t, :), result%weights(t, :), lo, good)
            if (.not. good) then
                result%weights(t, :) = 1.0_dp / real(k, dp)
                result%n_fallback = result%n_fallback + 1
                good = .true.
            end if
            result%diversification_ratios(t) = dot_product(result%weights(t, :), sigma(t, :)) / &
                sqrt(max(dot_product(result%weights(t, :), matmul(covariance, result%weights(t, :))), tiny(1.0_dp)))
        end do
        call portfolio_returns(returns_matrix, result%weights, lag, result%returns)
        call summarize_returns(result)
        result%mean_diversification = sum(result%diversification_ratios) / real(nobs, dp)
        if (present(ok)) ok = .true.
    end subroutine rsdc_maxdiv

    subroutine build_covariance(sd, rho, covariance, ok)
        real(dp), intent(in) :: sd(:), rho(:)
        real(dp), intent(out) :: covariance(:, :)
        logical, intent(out) :: ok
        real(dp) :: r(size(sd), size(sd))
        integer :: i, j, k
        k = size(sd)
        ok = all(sd > 0.0_dp) .and. all(abs(rho) < 1.0_dp)
        if (.not. ok) return
        call correlations_to_matrix(rho, k, r)
        ok = is_positive_definite(r, 1.0e-10_dp)
        if (.not. ok) return
        do i = 1, k
            do j = 1, k
                covariance(i, j) = sd(i) * r(i, j) * sd(j)
            end do
        end do
    end subroutine build_covariance

    subroutine minimum_variance_weights(covariance, w, long_only, ok)
        real(dp), intent(in) :: covariance(:, :)
        real(dp), intent(out) :: w(:)
        logical, intent(in) :: long_only
        logical, intent(out) :: ok
        real(dp) :: target(size(w))
        target = 1.0_dp
        call constrained_inverse_weights(covariance, target, w, long_only, ok)
    end subroutine minimum_variance_weights

    subroutine maximum_diversification_weights(covariance, sd, w, long_only, ok)
        real(dp), intent(in) :: covariance(:, :), sd(:)
        real(dp), intent(out) :: w(:)
        logical, intent(in) :: long_only
        logical, intent(out) :: ok
        if (long_only) then
            call constrained_inverse_weights(covariance, sd, w, .true., ok)
        else
            call bounded_maxdiv_weights(covariance, sd, w, ok)
        end if
    end subroutine maximum_diversification_weights

    subroutine bounded_maxdiv_weights(covariance, sd, w, ok)
        real(dp), intent(in) :: covariance(:, :), sd(:)
        real(dp), intent(out) :: w(:)
        logical, intent(out) :: ok
        real(dp) :: grad(size(w)), trial(size(w)), sw(size(w))
        real(dp) :: q, numerator, denominator, f, ftrial, step
        integer :: iter
        w = 1.0_dp / real(size(w), dp)
        f = maxdiv_objective(w)
        do iter = 1, 1000
            q = dot_product(w, matmul(covariance, w))
            numerator = dot_product(w, sd)
            denominator = sqrt(max(q, tiny(1.0_dp)))
            grad = -sd / denominator + numerator * matmul(covariance, w) / max(q * denominator, tiny(1.0_dp))
            step = 0.2_dp
            do
                trial = w - step * grad
                call project_box_sum(trial, sw)
                ftrial = maxdiv_objective(sw)
                if (ftrial <= f .or. step < 1.0e-10_dp) exit
                step = 0.5_dp * step
            end do
            if (maxval(abs(sw - w)) < 1.0e-10_dp) exit
            w = sw
            f = ftrial
        end do
        ok = all(w >= -1.0_dp - 1.0e-10_dp) .and. all(w <= 1.0_dp + 1.0e-10_dp) .and. &
             abs(sum(w) - 1.0_dp) < 1.0e-8_dp
    contains
        real(dp) function maxdiv_objective(v) result(value)
            real(dp), intent(in) :: v(:)
            value = -dot_product(v, sd) / sqrt(max(dot_product(v, matmul(covariance, v)), tiny(1.0_dp)))
        end function maxdiv_objective
    end subroutine bounded_maxdiv_weights

    subroutine project_box_sum(v, w)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        real(dp) :: lo, hi, mid, s
        integer :: iter
        lo = minval(v) - 2.0_dp
        hi = maxval(v) + 2.0_dp
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            s = sum(max(min(v - mid, 1.0_dp), -1.0_dp))
            if (s > 1.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        w = max(min(v - 0.5_dp * (lo + hi), 1.0_dp), -1.0_dp)
    end subroutine project_box_sum

    subroutine constrained_inverse_weights(covariance, target, w, long_only, ok)
        real(dp), intent(in) :: covariance(:, :), target(:)
        real(dp), intent(out) :: w(:)
        logical, intent(in) :: long_only
        logical, intent(out) :: ok
        integer, allocatable :: active(:)
        real(dp), allocatable :: sub(:, :), inv(:, :), v(:)
        integer :: n, m, i, j, worst, newm
        real(dp) :: denom
        n = size(w)
        allocate(active(n)); active = [(i, i=1,n)]
        w = 0.0_dp
        do
            m = size(active)
            allocate(sub(m, m), inv(m, m), v(m))
            do i = 1, m
                do j = 1, m
                    sub(i, j) = covariance(active(i), active(j))
                end do
            end do
            call inverse_spd(sub, inv, ok)
            if (.not. ok) return
            v = matmul(inv, target(active))
            denom = sum(v)
            if (abs(denom) <= tiny(1.0_dp)) then
                ok = .false.; return
            end if
            v = v / denom
            if (.not. long_only .or. all(v >= -1.0e-12_dp)) then
                w = 0.0_dp
                do i = 1, m
                    w(active(i)) = max(v(i), 0.0_dp)
                end do
                if (long_only) w = w / sum(w)
                ok = .true.
                return
            end if
            worst = minloc(v, dim=1)
            if (m <= 1) then
                ok = .false.; return
            end if
            newm = m - 1
            active = [active(1:worst - 1), active(worst + 1:m)]
            deallocate(sub, inv, v)
        end do
    end subroutine constrained_inverse_weights

    subroutine portfolio_returns(y, weights, lagged, r)
        real(dp), intent(in) :: y(:, :), weights(:, :)
        logical, intent(in) :: lagged
        real(dp), intent(out) :: r(:)
        integer :: t
        if (lagged) then
            r(1) = 0.0_dp
            do t = 2, size(y, 1)
                r(t) = dot_product(y(t, :), weights(t - 1, :))
            end do
        else
            do t = 1, size(y, 1)
                r(t) = dot_product(y(t, :), weights(t, :))
            end do
        end if
    end subroutine portfolio_returns

    subroutine summarize_returns(result)
        type(rsdc_portfolio_result), intent(inout) :: result
        result%mean_return = sum(result%returns) / real(size(result%returns), dp)
        result%realized_volatility = standard_deviation(result%returns)
    end subroutine summarize_returns
end module rsdc_portfolio
