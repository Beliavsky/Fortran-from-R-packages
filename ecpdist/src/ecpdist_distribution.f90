! SPDX-License-Identifier: GPL-3.0-only
module ecpdist_distribution
    use ecpdist_kinds, only: dp
    use ecpdist_math, only: expm1_stable, log_abs_expm1, logsumexp2, nan_dp, inf_dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: decp, pecp, qecp, secp, hecp, recp
    public :: ecp_cumhaz, ecp_valid_parameters

contains

    pure logical function ecp_valid_parameters(lambda, gamma, phi) result(ok)
        real(dp), intent(in) :: lambda, gamma, phi
        ok = ieee_is_finite(lambda) .and. ieee_is_finite(gamma) .and. &
            ieee_is_finite(phi) .and. lambda > 0.0_dp .and. gamma > 0.0_dp .and. &
            abs(phi) > tiny(1.0_dp)
    end function ecp_valid_parameters

    pure real(dp) function log_survival_ecp(x, lambda, gamma, phi) result(ls)
        real(dp), intent(in) :: x, lambda, gamma, phi
        real(dp) :: xg, log_a

        if (.not. ecp_valid_parameters(lambda, gamma, phi) .or. x < 0.0_dp) then
            ls = nan_dp()
            return
        end if
        if (.not. ieee_is_finite(x)) then
            ls = -inf_dp()
            return
        end if
        if (x <= 0.0_dp) then
            ls = 0.0_dp
            return
        end if

        xg = x**gamma
        if (.not. ieee_is_finite(xg) .or. xg > 700.0_dp) then
            ls = -inf_dp()
            return
        end if
        log_a = lambda*(1.0_dp - exp(xg))
        if (log_a < log(tiny(1.0_dp))) then
            ls = -inf_dp()
            return
        end if
        ls = log_abs_expm1(-phi*exp(log_a)) - log_abs_expm1(-phi)
        if (ls > 0.0_dp .and. ls < 100.0_dp*epsilon(1.0_dp)) ls = 0.0_dp
    end function log_survival_ecp

    elemental real(dp) function decp(x, lambda, gamma, phi, log_value) result(ans)
        real(dp), intent(in) :: x, lambda, gamma, phi
        logical, intent(in), optional :: log_value
        logical :: want_log
        real(dp) :: xg, log_a, log_pdf

        want_log = .false.
        if (present(log_value)) want_log = log_value

        if (.not. ecp_valid_parameters(lambda, gamma, phi) .or. x < 0.0_dp) then
            ans = nan_dp()
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (want_log) then
                ans = -inf_dp()
            else
                ans = 0.0_dp
            end if
            return
        end if

        if (x <= 0.0_dp) then
            if (gamma < 1.0_dp) then
                ans = inf_dp()
                if (want_log) ans = inf_dp()
                return
            else if (gamma > 1.0_dp) then
                if (want_log) then
                    ans = -inf_dp()
                else
                    ans = 0.0_dp
                end if
                return
            end if
            xg = 0.0_dp
            log_a = 0.0_dp
            log_pdf = log(lambda*gamma) + log(abs(phi)) - log_abs_expm1(-phi) - phi
        else
            xg = x**gamma
            if (.not. ieee_is_finite(xg) .or. xg > 700.0_dp) then
                if (want_log) then
                    ans = -inf_dp()
                else
                    ans = 0.0_dp
                end if
                return
            end if
            log_a = lambda*(1.0_dp - exp(xg))
            log_pdf = log(lambda*gamma) + log(abs(phi)) - log_abs_expm1(-phi) + &
                (gamma - 1.0_dp)*log(x) - phi*exp(log_a) + log_a + xg
        end if

        if (want_log) then
            ans = log_pdf
        else
            ans = exp(log_pdf)
        end if
    end function decp

    elemental real(dp) function pecp(q, lambda, gamma, phi, lower_tail, log_p) result(ans)
        real(dp), intent(in) :: q, lambda, gamma, phi
        logical, intent(in), optional :: lower_tail, log_p
        logical :: lower, want_log
        real(dp) :: ls, lf

        lower = .true.
        want_log = .false.
        if (present(lower_tail)) lower = lower_tail
        if (present(log_p)) want_log = log_p

        ls = log_survival_ecp(q, lambda, gamma, phi)
        if (.not. ieee_is_finite(ls) .and. ls > 0.0_dp) then
            ans = ls
            return
        end if

        if (lower) then
            if (.not. ieee_is_finite(ls) .and. ls < 0.0_dp) then
                if (want_log) then
                    ans = 0.0_dp
                else
                    ans = 1.0_dp
                end if
            else if (ls >= 0.0_dp) then
                if (want_log) then
                    ans = -inf_dp()
                else
                    ans = 0.0_dp
                end if
            else
                lf = log(-expm1_stable(ls))
                if (want_log) then
                    ans = lf
                else
                    ans = exp(lf)
                end if
            end if
        else
            if (want_log) then
                ans = ls
            else
                ans = exp(ls)
            end if
        end if
    end function pecp

    elemental real(dp) function secp(q, lambda, gamma, phi, lower_tail, cum_haz) result(ans)
        real(dp), intent(in) :: q, lambda, gamma, phi
        logical, intent(in), optional :: lower_tail, cum_haz
        logical :: lower, ch

        lower = .false.
        ch = .false.
        if (present(lower_tail)) lower = lower_tail
        if (present(cum_haz)) ch = cum_haz

        if (ch) then
            if (lower) then
                ans = -pecp(q, lambda, gamma, phi, lower_tail=.true., log_p=.true.)
            else
                ans = -pecp(q, lambda, gamma, phi, lower_tail=.false., log_p=.true.)
            end if
        else
            ans = pecp(q, lambda, gamma, phi, lower_tail=lower, log_p=.false.)
        end if
    end function secp

    elemental real(dp) function ecp_cumhaz(x, lambda, gamma, phi) result(ans)
        real(dp), intent(in) :: x, lambda, gamma, phi
        ans = -pecp(x, lambda, gamma, phi, lower_tail=.false., log_p=.true.)
    end function ecp_cumhaz

    elemental real(dp) function hecp(x, lambda, gamma, phi, log_value) result(ans)
        real(dp), intent(in) :: x, lambda, gamma, phi
        logical, intent(in), optional :: log_value
        logical :: want_log
        real(dp) :: log_h

        want_log = .false.
        if (present(log_value)) want_log = log_value
        log_h = decp(x, lambda, gamma, phi, log_value=.true.) - &
            pecp(x, lambda, gamma, phi, lower_tail=.false., log_p=.true.)
        if (want_log) then
            ans = log_h
        else
            ans = exp(log_h)
        end if
    end function hecp

    elemental real(dp) function qecp(p, lambda, gamma, phi, lower_tail, log_p) result(ans)
        real(dp), intent(in) :: p, lambda, gamma, phi
        logical, intent(in), optional :: lower_tail, log_p
        logical :: lower, p_is_log
        real(dp) :: u, log_u, log_1mu, log_b, a, z

        lower = .true.
        p_is_log = .false.
        if (present(lower_tail)) lower = lower_tail
        if (present(log_p)) p_is_log = log_p

        if (.not. ecp_valid_parameters(lambda, gamma, phi)) then
            ans = nan_dp()
            return
        end if

        if (p_is_log) then
            if (p > 0.0_dp) then
                ans = nan_dp()
                return
            end if
            u = exp(p)
        else
            if (p < 0.0_dp .or. p > 1.0_dp) then
                ans = nan_dp()
                return
            end if
            u = p
        end if
        if (.not. lower) u = 1.0_dp - u

        if (u <= 0.0_dp) then
            ans = 0.0_dp
            return
        end if
        if (u >= 1.0_dp) then
            ans = inf_dp()
            return
        end if

        log_u = log(u)
        log_1mu = log(1.0_dp - u)
        log_b = logsumexp2(log_u, log_1mu - phi)
        a = -log_b/phi

        if (a <= 0.0_dp) then
            ans = inf_dp()
            return
        end if
        if (a > 1.0_dp .and. a < 1.0_dp + 100.0_dp*epsilon(1.0_dp)) a = 1.0_dp
        if (a > 1.0_dp) then
            ans = nan_dp()
            return
        end if

        z = log(1.0_dp - log(a)/lambda)
        if (z < 0.0_dp .and. z > -100.0_dp*epsilon(1.0_dp)) z = 0.0_dp
        if (z < 0.0_dp) then
            ans = nan_dp()
        else
            ans = z**(1.0_dp/gamma)
        end if
    end function qecp

    subroutine recp(n, lambda, gamma, phi, x, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda, gamma, phi
        real(dp), intent(out) :: x(:)
        integer, intent(out), optional :: status
        real(dp), allocatable :: u(:)

        if (present(status)) status = 0
        if (n < 1 .or. size(x) < n .or. .not. ecp_valid_parameters(lambda, gamma, phi)) then
            if (present(status)) status = 1
            if (size(x) > 0) x = nan_dp()
            return
        end if

        allocate(u(n))
        call random_number(u)
        x(1:n) = qecp(u, lambda, gamma, phi)
    end subroutine recp

end module ecpdist_distribution
