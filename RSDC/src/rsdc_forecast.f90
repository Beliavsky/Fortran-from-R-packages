! SPDX-License-Identifier: GPL-3.0-only
module rsdc_forecast
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_model, rsdc_filter_result, rsdc_forecast_result
    use rsdc_types, only: rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc_filter, only: rsdc_hamilton
    use rsdc_parameters, only: transition_at, correlations_to_matrix, lower_tri_size
    implicit none
    private
    public :: rsdc_forecast_path, rsdc_forecast_ahead

contains

    subroutine rsdc_forecast_path(model, residuals, sigma, result, x, use_smoothed, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :), sigma(:, :)
        type(rsdc_forecast_result), intent(out) :: result
        real(dp), intent(in), optional :: x(:, :)
        logical, intent(in), optional :: use_smoothed
        logical, intent(out), optional :: ok
        type(rsdc_filter_result) :: filt
        real(dp), allocatable :: weights(:, :), r(:, :), d(:, :)
        logical :: smooth, good
        integer :: t, s, i, k, c, nobs, npar
        k = model%n_series; c = lower_tri_size(k); nobs = size(residuals, 1)
        smooth = .true.; if (present(use_smoothed)) smooth = use_smoothed
        good = size(residuals, 2) == k .and. size(sigma, 1) == nobs .and. size(sigma, 2) == k
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(result%predicted_correlations(nobs, c), result%covariance(k, k, nobs))
        if (model%method == rsdc_const) then
            result%predicted_correlations = spread(model%correlations(1, :), 1, nobs)
            allocate(result%regime_probabilities(nobs, 1)); result%regime_probabilities = 1.0_dp
        else
            if (model%method == rsdc_tvtp) then
                if (.not. present(x)) then
                    if (present(ok)) ok = .false.
                    return
                end if
                call rsdc_hamilton(residuals, model%correlations, filt, x=x, beta=model%beta)
            else
                call rsdc_hamilton(residuals, model%correlations, filt, pmat=model%transition_matrix)
            end if
            if (.not. filt%ok) then
                if (present(ok)) ok = .false.
                return
            end if
            allocate(weights(model%n_regimes, nobs))
            if (smooth) then
                weights = filt%smoothed
            else
                weights = filt%filtered
            end if
            allocate(result%regime_probabilities(nobs, model%n_regimes))
            result%regime_probabilities = transpose(weights)
            result%predicted_correlations = 0.0_dp
            do t = 1, nobs
                do s = 1, model%n_regimes
                    result%predicted_correlations(t, :) = result%predicted_correlations(t, :) + &
                        weights(s, t) * model%correlations(s, :)
                end do
            end do
        end if
        allocate(r(k, k), d(k, k))
        do t = 1, nobs
            call correlations_to_matrix(result%predicted_correlations(t, :), k, r)
            d = 0.0_dp
            do i = 1, k
                d(i, i) = sigma(t, i)
            end do
            result%covariance(:, :, t) = matmul(d, matmul(r, d))
        end do
        select case (model%method)
        case (rsdc_const)
            npar = c
        case (rsdc_nox)
            npar = model%n_regimes * (model%n_regimes - 1) + model%n_regimes * c
        case default
            if (model%n_regimes == 2) then
                npar = model%n_regimes * model%n_covariates + model%n_regimes * c
            else
                npar = model%n_regimes * (model%n_regimes - 1) * model%n_covariates + &
                    model%n_regimes * c
            end if
        end select
        result%bic = log(real(nobs, dp)) * real(npar, dp) - 2.0_dp * model%log_likelihood
        if (present(ok)) ok = .true.
    end subroutine rsdc_forecast_path

    subroutine rsdc_forecast_ahead(model, residuals, horizon, regime_probs, predicted_corr, x, x_future, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        integer, intent(in) :: horizon
        real(dp), allocatable, intent(out) :: regime_probs(:, :), predicted_corr(:, :)
        real(dp), intent(in), optional :: x(:, :), x_future(:, :)
        logical, intent(out), optional :: ok
        type(rsdc_filter_result) :: filt
        real(dp), allocatable :: pi_k(:), pmat(:, :), xlast(:)
        integer :: h, n, c
        logical :: good
        n = model%n_regimes; c = size(model%correlations, 2)
        good = horizon >= 1
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(regime_probs(horizon, n), predicted_corr(horizon, c), pi_k(n), pmat(n, n))
        if (model%method == rsdc_const) then
            pi_k = 1.0_dp
        else if (model%method == rsdc_nox) then
            call rsdc_hamilton(residuals, model%correlations, filt, pmat=model%transition_matrix)
            good = filt%ok
            if (good) pi_k = filt%filtered(:, size(residuals, 1))
        else
            if (.not. present(x)) then
                good = .false.
            else
                call rsdc_hamilton(residuals, model%correlations, filt, x=x, beta=model%beta)
                good = filt%ok
                if (good) pi_k = filt%filtered(:, size(residuals, 1))
            end if
        end if
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        do h = 1, horizon
            select case (model%method)
            case (rsdc_const)
                pmat(1, 1) = 1.0_dp
            case (rsdc_nox)
                pmat = model%transition_matrix
            case (rsdc_tvtp)
                if (present(x_future)) then
                    if (size(x_future, 1) < horizon) then
                        if (present(ok)) ok = .false.
                        return
                    end if
                    call transition_at(model%beta, x_future(h, :), pmat)
                else
                    allocate(xlast(size(x, 2)))
                    xlast = x(size(x, 1), :)
                    call transition_at(model%beta, xlast, pmat)
                    deallocate(xlast)
                end if
            end select
            pi_k = matmul(pi_k, pmat)
            regime_probs(h, :) = pi_k
            predicted_corr(h, :) = matmul(pi_k, model%correlations)
        end do
        if (present(ok)) ok = .true.
    end subroutine rsdc_forecast_ahead
end module rsdc_forecast
