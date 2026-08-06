! SPDX-License-Identifier: GPL-3.0-only
module rsdc_bands
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_model, rsdc_bands_result, rsdc_filter_result
    use rsdc_types, only: rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc_parameters, only: unpack_natural_parameters
    use rsdc_filter, only: rsdc_hamilton
    use rsdc_random, only: seed_rsdc, random_multivariate_normal
    use rsdc_linalg, only: quantile_type7
    implicit none
    private
    public :: rsdc_corr_bands

contains

    subroutine rsdc_corr_bands(model, residuals, b, bands, x, level, seed, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        integer, intent(in) :: b
        type(rsdc_bands_result), intent(out) :: bands
        real(dp), intent(in), optional :: x(:, :)
        real(dp), intent(in), optional :: level
        integer, intent(in), optional :: seed
        logical, intent(out), optional :: ok
        real(dp), allocatable :: point(:, :), paths(:, :, :), draw(:), work(:)
        real(dp) :: lev, a
        integer :: i, t, c, used, base_seed
        logical :: good
        lev = 0.95_dp; if (present(level)) lev = level
        base_seed = 123; if (present(seed)) base_seed = seed
        good = allocated(model%vcov) .and. b >= 2 .and. lev > 0.0_dp .and. lev < 1.0_dp
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        call path_from_parameters(model%parameters, point, good)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        t = size(point, 1); c = size(point, 2)
        allocate(paths(t, c, b), draw(size(model%parameters)), work(b))
        call seed_rsdc(base_seed)
        used = 0
        do i = 1, b
            call random_multivariate_normal(model%parameters, model%vcov, draw, good)
            if (.not. good) exit
            call add_path(draw, good)
        end do
        bands%n_used = used; bands%level = lev
        allocate(bands%fit(t, c), bands%lower(t, c), bands%upper(t, c))
        bands%fit = point
        if (used >= 2) then
            a = (1.0_dp - lev) / 2.0_dp
            do i = 1, c
                do t = 1, size(point, 1)
                    work(1:used) = paths(t, i, 1:used)
                    bands%lower(t, i) = quantile_type7(work(1:used), a)
                    bands%upper(t, i) = quantile_type7(work(1:used), 1.0_dp - a)
                end do
            end do
            good = .true.
        else
            bands%lower = point; bands%upper = point; good = .false.
        end if
        if (present(ok)) ok = good

    contains
        subroutine add_path(parameters, local_ok)
            real(dp), intent(in) :: parameters(:)
            logical, intent(out) :: local_ok
            real(dp), allocatable :: one(:, :)
            call path_from_parameters(parameters, one, local_ok)
            if (local_ok) then
                used = used + 1
                paths(:, :, used) = one
            end if
        end subroutine add_path

        subroutine path_from_parameters(parameters, path, local_ok)
            real(dp), intent(in) :: parameters(:)
            real(dp), allocatable, intent(out) :: path(:, :)
            logical, intent(out) :: local_ok
            real(dp), allocatable :: beta(:, :), rho(:, :), pmat(:, :)
            type(rsdc_filter_result) :: filt
            integer :: s, tt
            call unpack_natural_parameters(parameters, model%method, model%n_regimes, model%n_series, &
                model%n_covariates, beta, rho, pmat, local_ok)
            if (.not. local_ok .or. any(abs(rho) >= 1.0_dp)) return
            if (model%method == rsdc_const) then
                allocate(path(size(residuals, 1), size(rho, 2)))
                path = spread(rho(1, :), 1, size(residuals, 1))
                return
            else if (model%method == rsdc_nox) then
                call rsdc_hamilton(residuals, rho, filt, pmat=pmat)
            else
                if (.not. present(x)) then
                    local_ok = .false.; return
                end if
                call rsdc_hamilton(residuals, rho, filt, x=x, beta=beta)
            end if
            local_ok = filt%ok
            if (.not. local_ok) return
            allocate(path(size(residuals, 1), size(rho, 2))); path = 0.0_dp
            do tt = 1, size(residuals, 1)
                do s = 1, model%n_regimes
                    path(tt, :) = path(tt, :) + filt%smoothed(s, tt) * rho(s, :)
                end do
            end do
        end subroutine path_from_parameters
    end subroutine rsdc_corr_bands
end module rsdc_bands
