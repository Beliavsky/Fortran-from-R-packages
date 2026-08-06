! SPDX-License-Identifier: GPL-3.0-only
module rsdc_bootstrap
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_model, rsdc_control, rsdc_bootstrap_result, rsdc_simulation_result
    use rsdc_types, only: rsdc_tvtp
    use rsdc_simulation, only: rsdc_simulate, rsdc_simulate_fixed
    use rsdc_estimation, only: rsdc_estimate
    use rsdc_linalg, only: standard_deviation, quantile_type7, mean_columns
    implicit none
    private
    public :: rsdc_parametric_bootstrap

contains

    subroutine rsdc_parametric_bootstrap(model, residuals, b, result, x, seed, control, level, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        integer, intent(in) :: b
        type(rsdc_bootstrap_result), intent(out) :: result
        real(dp), intent(in), optional :: x(:, :)
        integer, intent(in), optional :: seed
        type(rsdc_control), intent(in), optional :: control
        real(dp), intent(in), optional :: level
        logical, intent(out), optional :: ok
        type(rsdc_control) :: cfg
        type(rsdc_simulation_result) :: sim
        type(rsdc_model) :: fit
        real(dp), allocatable :: mu(:, :), draws(:, :), ysim(:, :), work(:), center(:), scale(:)
        real(dp) :: lev, alpha
        integer :: i, j, base_seed, np, nobs, n, k, success
        logical :: good
        cfg = rsdc_control()
        if (present(control)) cfg = control
        if (allocated(cfg%start)) deallocate(cfg%start)
        allocate(cfg%start(size(model%parameters)))
        cfg%start = model%parameters
        cfg%compute_vcov = .false.
        lev = 0.95_dp
        if (present(level)) lev = level
        base_seed = 12345
        if (present(seed)) base_seed = seed
        nobs = size(residuals, 1)
        n = model%n_regimes
        k = model%n_series
        np = size(model%parameters)
        good = b >= 2 .and. lev > 0.0_dp .and. lev < 1.0_dp
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(mu(n, k), draws(b, np), ysim(nobs, k), center(k), scale(k))
        mu = 0.0_dp
        draws = 0.0_dp
        success = 0
        do i = 1, b
            if (model%method == rsdc_tvtp) then
                if (.not. present(x)) cycle
                call rsdc_simulate(nobs, x, model%beta, mu, model%covariance, sim, base_seed + i, good)
            else
                call rsdc_simulate_fixed(nobs, model%transition_matrix, mu, model%covariance, sim, &
                    base_seed + i, good)
            end if
            if (.not. good) cycle
            ysim = sim%observations
            center = mean_columns(ysim)
            do j = 1, k
                ysim(:, j) = ysim(:, j) - center(j)
                scale(j) = standard_deviation(ysim(:, j))
                if (scale(j) > 0.0_dp) ysim(:, j) = ysim(:, j) / scale(j)
            end do
            if (model%method == rsdc_tvtp) then
                call rsdc_estimate(model%method, ysim, n, fit, x=x, control=cfg, ok=good)
            else
                call rsdc_estimate(model%method, ysim, n, fit, control=cfg, ok=good)
            end if
            if (good .and. size(fit%parameters) == np) then
                success = success + 1
                draws(success, :) = fit%parameters
            end if
        end do
        result%n_success = success
        result%level = lev
        allocate(result%standard_errors(np), result%covariance(np, np), result%confidence_interval(np, 2))
        if (success >= 2) then
            allocate(result%parameter_draws(success, np), work(success))
            result%parameter_draws = draws(1:success, :)
            result%covariance = 0.0_dp
            do i = 1, np
                result%standard_errors(i) = standard_deviation(result%parameter_draws(:, i))
            end do
            do i = 1, success
                result%covariance = result%covariance + outer_product( &
                    result%parameter_draws(i, :) - sum(result%parameter_draws, dim=1) / real(success, dp))
            end do
            result%covariance = result%covariance / real(success - 1, dp)
            alpha = (1.0_dp - lev) / 2.0_dp
            do j = 1, np
                work = result%parameter_draws(:, j)
                result%confidence_interval(j, 1) = quantile_type7(work, alpha)
                result%confidence_interval(j, 2) = quantile_type7(work, 1.0_dp - alpha)
            end do
            good = .true.
        else
            allocate(result%parameter_draws(0, np))
            result%standard_errors = 0.0_dp
            result%covariance = 0.0_dp
            result%confidence_interval = 0.0_dp
            good = .false.
        end if
        if (present(ok)) ok = good
    end subroutine rsdc_parametric_bootstrap

    pure function outer_product(v) result(a)
        real(dp), intent(in) :: v(:)
        real(dp) :: a(size(v), size(v))
        a = spread(v, 2, size(v)) * spread(v, 1, size(v))
    end function outer_product
end module rsdc_bootstrap
