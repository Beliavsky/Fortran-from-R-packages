! SPDX-License-Identifier: GPL-3.0-only
module rsdc_filter
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use rsdc_kinds, only: dp, rsdc_pi
    use rsdc_types, only: rsdc_filter_result
    use rsdc_linalg, only: inverse_spd, logdet_spd, is_positive_definite
    use rsdc_parameters, only: correlations_to_matrix, transition_at
    implicit none
    private
    public :: rsdc_hamilton, rsdc_log_likelihood

contains

    subroutine rsdc_hamilton(y, rho, result, pmat, x, beta, xi_init)
        real(dp), intent(in) :: y(:, :)
        real(dp), intent(in) :: rho(:, :)
        type(rsdc_filter_result), intent(out) :: result
        real(dp), intent(in), optional :: pmat(:, :)
        real(dp), intent(in), optional :: x(:, :), beta(:, :), xi_init(:)

        integer :: t, s, nobs, k, n, i
        real(dp), allocatable :: logdens(:, :), sigma(:, :, :), inv(:, :)
        real(dp), allocatable :: xi(:), pred(:), w(:), ld(:), ratio(:), temp(:)
        real(dp), allocatable :: pt(:, :), rr(:, :)
        real(dp) :: c, sw, logdet, quad, denom, ridge
        logical :: ok, tvtp

        nobs = size(y, 1)
        k = size(y, 2)
        n = size(rho, 1)
        result%ok = .false.
        result%log_likelihood = -huge(1.0_dp)
        if (nobs < 1 .or. k < 1 .or. size(rho, 2) /= k * (k - 1) / 2) return
        tvtp = present(x) .and. present(beta)
        if (present(x) .neqv. present(beta)) return
        if (tvtp) then
            if (size(x, 1) /= nobs .or. size(beta, 1) /= n) return
        else
            if (.not. present(pmat)) return
            if (size(pmat, 1) /= n .or. size(pmat, 2) /= n) return
        end if

        allocate(logdens(n, nobs), sigma(k, k, n), inv(k, k), rr(k, k))
        ridge = 1.0e-8_dp
        do s = 1, n
            call correlations_to_matrix(rho(s, :), k, rr)
            if (.not. is_positive_definite(rr, 1.0e-8_dp)) return
            sigma(:, :, s) = rr
            do i = 1, k
                rr(i, i) = rr(i, i) + ridge
            end do
            call inverse_spd(rr, inv, ok)
            if (.not. ok) return
            logdet = logdet_spd(rr, ok)
            if (.not. ok) return
            do t = 1, nobs
                quad = dot_product(y(t, :), matmul(inv, y(t, :)))
                logdens(s, t) = -0.5_dp * (logdet + quad + real(k, dp) * log(2.0_dp * rsdc_pi))
            end do
        end do

        allocate(result%filtered(n, nobs), result%predicted(n, nobs))
        allocate(result%smoothed(n, nobs), result%loglik_t(nobs))
        allocate(xi(n), pred(n), w(n), ld(n), ratio(n), temp(n), pt(n, n))
        if (present(xi_init)) then
            if (size(xi_init) /= n .or. any(xi_init < 0.0_dp) .or. sum(xi_init) <= 0.0_dp) return
            xi = xi_init / sum(xi_init)
        else
            xi = 1.0_dp / real(n, dp)
        end if
        result%log_likelihood = 0.0_dp
        do t = 1, nobs
            if (tvtp) then
                call transition_at(beta, x(t, :), pt)
            else
                pt = pmat
            end if
            pred = matmul(transpose(pt), xi)
            result%predicted(:, t) = pred
            ld = logdens(:, t)
            c = maxval(ld)
            w = pred * exp(ld - c)
            sw = sum(w)
            if (.not. ieee_is_finite(sw) .or. sw <= 0.0_dp) return
            result%filtered(:, t) = w / sw
            result%loglik_t(t) = log(sw) + c
            result%log_likelihood = result%log_likelihood + result%loglik_t(t)
            xi = result%filtered(:, t)
        end do

        result%smoothed(:, nobs) = result%filtered(:, nobs)
        do t = nobs - 1, 1, -1
            if (tvtp) then
                call transition_at(beta, x(t + 1, :), pt)
            else
                pt = pmat
            end if
            do i = 1, n
                denom = max(result%predicted(i, t + 1), epsilon(1.0_dp))
                ratio(i) = result%smoothed(i, t + 1) / denom
            end do
            temp = result%filtered(:, t) * matmul(pt, ratio)
            sw = sum(temp)
            if (.not. ieee_is_finite(sw) .or. sw <= epsilon(1.0_dp)) return
            result%smoothed(:, t) = temp / sw
        end do
        result%ok = .true.
    end subroutine rsdc_hamilton

    real(dp) function rsdc_log_likelihood(y, rho, pmat, x, beta, xi_init) result(ll)
        real(dp), intent(in) :: y(:, :), rho(:, :)
        real(dp), intent(in), optional :: pmat(:, :), x(:, :), beta(:, :), xi_init(:)
        type(rsdc_filter_result) :: out
        if (present(x) .and. present(beta)) then
            if (present(xi_init)) then
                call rsdc_hamilton(y, rho, out, x=x, beta=beta, xi_init=xi_init)
            else
                call rsdc_hamilton(y, rho, out, x=x, beta=beta)
            end if
        else if (present(pmat)) then
            if (present(xi_init)) then
                call rsdc_hamilton(y, rho, out, pmat=pmat, xi_init=xi_init)
            else
                call rsdc_hamilton(y, rho, out, pmat=pmat)
            end if
        else
            ll = -huge(1.0_dp)
            return
        end if
        ll = out%log_likelihood
    end function rsdc_log_likelihood
end module rsdc_filter
