! SPDX-License-Identifier: GPL-3.0-only
module rsdc_starts
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_starts_result, rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc_parameters, only: lower_tri_size, expected_parameter_count, matrix_to_correlations
    use rsdc_linalg, only: sample_correlation, quantile_type7
    use rsdc_likelihood, only: rsdc_negative_log_likelihood
    implicit none
    private
    public :: rsdc_make_starts

contains

    subroutine rsdc_make_starts(residuals, n_regimes, method, result, x, window, n_starts, stay, ok)
        real(dp), intent(in) :: residuals(:, :)
        integer, intent(in) :: n_regimes, method
        type(rsdc_starts_result), intent(out) :: result
        real(dp), intent(in), optional :: x(:, :)
        integer, intent(in), optional :: window, n_starts
        real(dp), intent(in), optional :: stay
        logical, intent(out), optional :: ok
        integer :: nobs, k, n, p, c, w, ns, t, s, b, nh, np, i, n_group
        real(dp) :: stay0
        real(dp), allocatable :: indicator(:), cuts(:), group_data(:, :), cor(:, :), rho(:)
        real(dp), allocatable :: base_rho(:, :), params(:), shr(:)
        logical :: good
        nobs = size(residuals, 1); k = size(residuals, 2)
        n = merge(1, n_regimes, method == rsdc_const)
        p = 0; if (present(x)) p = size(x, 2)
        w = min(126, nobs); if (present(window)) w = min(max(window, 3), nobs)
        ns = 5; if (present(n_starts)) ns = max(1, n_starts)
        stay0 = 0.95_dp; if (present(stay)) stay0 = stay
        c = lower_tri_size(k); np = expected_parameter_count(method, n, k, p)
        good = nobs >= 3 .and. k >= 2 .and. stay0 > 0.0_dp .and. stay0 < 1.0_dp
        if (method == rsdc_tvtp) good = good .and. present(x)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(indicator(nobs), result%regime_split(nobs), cuts(max(n - 1, 0)))
        do t = 1, nobs
            call mean_pairwise_window(residuals(max(1, t - w + 1):t, :), indicator(t))
        end do
        if (n == 1) then
            result%regime_split = 1
        else
            do s = 1, n - 1
                cuts(s) = quantile_type7(indicator, real(s, dp) / real(n, dp))
            end do
            do t = 1, nobs
                result%regime_split(t) = 1
                do s = 1, n - 1
                    if (indicator(t) > cuts(s)) result%regime_split(t) = s + 1
                end do
            end do
        end if
        allocate(base_rho(n, c), cor(k, k), rho(c))
        do s = 1, n
            n_group = count(result%regime_split == s)
            if (n_group >= max(3, k + 1)) then
                allocate(group_data(n_group, k)); b = 0
                do t = 1, nobs
                    if (result%regime_split(t) == s) then
                        b = b + 1; group_data(b, :) = residuals(t, :)
                    end if
                end do
                call sample_correlation(group_data, cor)
                deallocate(group_data)
            else
                call sample_correlation(residuals, cor)
            end if
            call matrix_to_correlations(cor, rho)
            base_rho(s, :) = rho
        end do
        allocate(result%starts(ns, np), result%initial_log_likelihood(ns), result%shrinkage(ns), shr(ns))
        if (ns == 1) then
            shr(1) = 1.0_dp
        else
            do i = 1, ns
                shr(i) = 1.0_dp - 0.55_dp * real(i - 1, dp) / real(ns - 1, dp)
            end do
        end if
        result%shrinkage = shr
        nh = np - n * c
        allocate(params(np))
        do i = 1, ns
            params = 0.0_dp
            if (method == rsdc_nox) then
                if (n == 2) then
                    params(1:2) = stay0
                else
                    do s = 1, n
                        params((s - 1) * (n - 1) + 1:s * (n - 1)) = (1.0_dp - stay0) / real(n - 1, dp)
                        if (s < n) params((s - 1) * (n - 1) + s) = stay0
                    end do
                end if
            else if (method == rsdc_tvtp) then
                if (n == 2) then
                    do s = 1, n
                        params((s - 1) * p + 1) = log(stay0 / (1.0_dp - stay0))
                    end do
                else
                    do s = 1, n
                        if (s < n) then
                            params((s - 1) * (n - 1) * p + (s - 1) * p + 1) = &
                                log(stay0 * real(n - 1, dp) / (1.0_dp - stay0))
                        else
                            do b = 1, n - 1
                                params((s - 1) * (n - 1) * p + (b - 1) * p + 1) = &
                                    log((1.0_dp - stay0) / (stay0 * real(n - 1, dp)))
                            end do
                        end if
                    end do
                end if
            end if
            params(nh + 1:) = reshape(transpose(shr(i) * base_rho), [n * c])
            result%starts(i, :) = params
            if (method == rsdc_tvtp) then
                result%initial_log_likelihood(i) = -rsdc_negative_log_likelihood(params, residuals, method, n, x)
            else
                result%initial_log_likelihood(i) = -rsdc_negative_log_likelihood(params, residuals, method, n)
            end if
        end do
        if (present(ok)) ok = .true.
    end subroutine rsdc_make_starts

    subroutine mean_pairwise_window(y, value)
        real(dp), intent(in) :: y(:, :)
        real(dp), intent(out) :: value
        real(dp), allocatable :: cor(:, :)
        integer :: i, j, k, q
        k = size(y, 2); allocate(cor(k, k)); call sample_correlation(y, cor)
        value = 0.0_dp; q = 0
        do j = 1, k - 1
            do i = j + 1, k
                value = value + cor(i, j); q = q + 1
            end do
        end do
        value = value / real(max(q, 1), dp)
    end subroutine mean_pairwise_window
end module rsdc_starts
