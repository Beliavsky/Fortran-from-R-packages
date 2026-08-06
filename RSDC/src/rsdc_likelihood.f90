! SPDX-License-Identifier: GPL-3.0-only
module rsdc_likelihood
    use rsdc_kinds, only: dp, rsdc_pi
    use rsdc_types, only: rsdc_const, rsdc_nox, rsdc_tvtp, rsdc_filter_result
    use rsdc_parameters, only: expected_parameter_count, unpack_natural_parameters
    use rsdc_parameters, only: correlations_to_matrix
    use rsdc_linalg, only: inverse_spd, logdet_spd, is_positive_definite
    use rsdc_filter, only: rsdc_hamilton
    implicit none
    private
    public :: rsdc_negative_log_likelihood, rsdc_loglik_contributions

contains

    real(dp) function rsdc_negative_log_likelihood(parameters, y, method, n, x) result(nll)
        real(dp), intent(in) :: parameters(:), y(:, :)
        integer, intent(in) :: method, n
        real(dp), intent(in), optional :: x(:, :)
        real(dp), allocatable :: beta(:, :), rho(:, :), pmat(:, :)
        integer :: p
        logical :: ok
        type(rsdc_filter_result) :: out
        nll = 1.0e10_dp
        if (size(y, 2) < 2 .or. any(abs(parameters) >= huge(1.0_dp))) return
        p = 0
        if (present(x)) p = size(x, 2)
        if (size(parameters) /= expected_parameter_count(method, n, size(y, 2), p)) return
        call unpack_natural_parameters(parameters, method, n, size(y, 2), p, beta, rho, pmat, ok)
        if (.not. ok .or. any(abs(rho) >= 1.0_dp)) return
        select case (method)
        case (rsdc_const)
            allocate(out%loglik_t(size(y, 1)))
            call constant_contributions(y, rho(1, :), out%loglik_t, ok)
            if (.not. ok) return
            nll = -sum(out%loglik_t)
        case (rsdc_nox)
            call rsdc_hamilton(y, rho, out, pmat=pmat)
            if (out%ok) nll = -out%log_likelihood
        case (rsdc_tvtp)
            if (.not. present(x)) return
            call rsdc_hamilton(y, rho, out, x=x, beta=beta)
            if (out%ok) nll = -out%log_likelihood
        end select
    end function rsdc_negative_log_likelihood

    subroutine rsdc_loglik_contributions(parameters, y, method, n, values, ok, x)
        real(dp), intent(in) :: parameters(:), y(:, :)
        integer, intent(in) :: method, n
        real(dp), intent(out) :: values(:)
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: x(:, :)
        real(dp), allocatable :: beta(:, :), rho(:, :), pmat(:, :)
        integer :: p
        type(rsdc_filter_result) :: out
        p = 0
        if (present(x)) p = size(x, 2)
        call unpack_natural_parameters(parameters, method, n, size(y, 2), p, beta, rho, pmat, ok)
        if (.not. ok) return
        select case (method)
        case (rsdc_const)
            call constant_contributions(y, rho(1, :), values, ok)
        case (rsdc_nox)
            call rsdc_hamilton(y, rho, out, pmat=pmat)
            ok = out%ok
            if (ok) values = out%loglik_t
        case (rsdc_tvtp)
            if (.not. present(x)) then
                ok = .false.; return
            end if
            call rsdc_hamilton(y, rho, out, x=x, beta=beta)
            ok = out%ok
            if (ok) values = out%loglik_t
        end select
    end subroutine rsdc_loglik_contributions

    subroutine constant_contributions(y, rho, values, ok)
        real(dp), intent(in) :: y(:, :), rho(:)
        real(dp), intent(out) :: values(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: r(:, :), inv(:, :)
        real(dp) :: ld, quad
        integer :: i, k, t
        k = size(y, 2)
        allocate(r(k, k), inv(k, k))
        call correlations_to_matrix(rho, k, r)
        if (.not. is_positive_definite(r, 1.0e-8_dp)) then
            ok = .false.; return
        end if
        do i = 1, k
            r(i, i) = r(i, i) + 1.0e-8_dp
        end do
        call inverse_spd(r, inv, ok)
        if (.not. ok) return
        ld = logdet_spd(r, ok)
        if (.not. ok) return
        do t = 1, size(y, 1)
            quad = dot_product(y(t, :), matmul(inv, y(t, :)))
            values(t) = -0.5_dp * (ld + quad + real(k, dp) * log(2.0_dp * rsdc_pi))
        end do
    end subroutine constant_contributions
end module rsdc_likelihood
