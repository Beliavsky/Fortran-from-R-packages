! SPDX-License-Identifier: GPL-2.0-or-later
module bgfd_core
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use adequacy_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: bgfd_e = 1
    integer, parameter, public :: bgfd_ee = 2
    integer, parameter, public :: bgfd_w = 3
    integer, parameter, public :: bgfd_ew = 4
    integer, parameter, public :: bgfd_f = 5
    integer, parameter, public :: bgfd_l = 6
    integer, parameter, public :: bgfd_b = 7
    integer, parameter, public :: bgfd_bx = 8

    public :: bgfd_npar, bgfd_pdf, bgfd_cdf, bgfd_quantile
    public :: bgfd_survival, bgfd_hazard, bgfd_random
    public :: base_pdf, base_cdf, base_quantile

contains

    pure integer function bgfd_npar(family) result(n)
        integer, intent(in) :: family
        select case (family)
        case (bgfd_e, bgfd_bx)
            n = 2
        case (bgfd_ee, bgfd_w, bgfd_f, bgfd_l)
            n = 3
        case (bgfd_ew, bgfd_b)
            n = 4
        case default
            n = 0
        end select
    end function bgfd_npar

    pure function expm1_stable(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (abs(x) < 1.0e-5_dp) then
            y = x * (1.0_dp + x * (0.5_dp + x * (1.0_dp/6.0_dp + x/24.0_dp)))
        else
            y = exp(x) - 1.0_dp
        end if
    end function expm1_stable

    pure function log1p_stable(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (abs(x) < 1.0e-5_dp) then
            y = x * (1.0_dp + x * (-0.5_dp + x * (1.0_dp/3.0_dp - x/4.0_dp)))
        else
            y = log(1.0_dp + x)
        end if
    end function log1p_stable

    pure function log_expm1_pos(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (x <= 0.0_dp) then
            y = -ieee_value(1.0_dp, ieee_positive_inf)
        else if (x > 40.0_dp) then
            y = x + log1p_stable(-exp(-x))
        else
            y = log(expm1_stable(x))
        end if
    end function log_expm1_pos

    pure function log1pexp(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (x > 40.0_dp) then
            y = x + log1p_stable(exp(-x))
        else if (x < -40.0_dp) then
            y = exp(x)
        else
            y = log1p_stable(exp(x))
        end if
    end function log1pexp

    pure function boundary_power_pdf(power, coefficient) result(value)
        real(dp), intent(in) :: power, coefficient
        real(dp) :: value
        real(dp), parameter :: tol = 32.0_dp * epsilon(1.0_dp)
        if (power > tol) then
            value = 0.0_dp
        else if (power < -tol) then
            value = ieee_value(1.0_dp, ieee_positive_inf)
        else
            value = coefficient
        end if
    end function boundary_power_pdf

    pure function base_cdf(family, x, par) result(value)
        integer, intent(in) :: family
        real(dp), intent(in) :: x, par(:)
        real(dp) :: value, z

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (.not. base_params_valid(family, par)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if

        select case (family)
        case (bgfd_e)
            value = -expm1_stable(-par(1)*x)
        case (bgfd_ee)
            value = (-expm1_stable(-par(1)*x))**par(2)
        case (bgfd_w)
            value = -expm1_stable(-par(1)*x**par(2))
        case (bgfd_ew)
            value = (-expm1_stable(-par(1)*x**par(2)))**par(3)
        case (bgfd_f)
            z = (x/par(1))**par(2)
            value = z/(1.0_dp + z)
        case (bgfd_l)
            value = 1.0_dp - (1.0_dp + x/par(1))**(-par(2))
        case (bgfd_b)
            z = (x/par(1))**par(2)
            value = 1.0_dp - (1.0_dp + z)**(-par(3))
        case (bgfd_bx)
            value = (-expm1_stable(-x*x))**par(1)
        case default
            value = ieee_value(1.0_dp, ieee_quiet_nan)
        end select
        value = min(1.0_dp, max(0.0_dp, value))
    end function base_cdf

    pure function base_pdf(family, x, par) result(value)
        integer, intent(in) :: family
        real(dp), intent(in) :: x, par(:)
        real(dp) :: value, z, epart

        if (x < 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (.not. base_params_valid(family, par)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        if (x <= tiny(1.0_dp)) then
            select case (family)
            case (bgfd_e)
                value = par(1)
            case (bgfd_ee)
                value = boundary_power_pdf(par(2)-1.0_dp, par(2)*par(1)**par(2))
            case (bgfd_w)
                value = boundary_power_pdf(par(2)-1.0_dp, par(1)*par(2))
            case (bgfd_ew)
                value = boundary_power_pdf(par(2)*par(3)-1.0_dp, &
                    par(2)*par(3)*par(1)**par(3))
            case (bgfd_f)
                value = boundary_power_pdf(par(2)-1.0_dp, par(2)*par(1)**(-par(2)))
            case (bgfd_l)
                value = par(2)/par(1)
            case (bgfd_b)
                value = boundary_power_pdf(par(2)-1.0_dp, &
                    par(3)*par(2)*par(1)**(-par(2)))
            case (bgfd_bx)
                value = boundary_power_pdf(2.0_dp*par(1)-1.0_dp, 2.0_dp*par(1))
            case default
                value = ieee_value(1.0_dp, ieee_quiet_nan)
            end select
            return
        end if

        select case (family)
        case (bgfd_e)
            value = par(1)*exp(-par(1)*x)
        case (bgfd_ee)
            epart = exp(-par(1)*x)
            value = par(1)*epart*par(2)*(1.0_dp-epart)**(par(2)-1.0_dp)
        case (bgfd_w)
            value = par(1)*par(2)*x**(par(2)-1.0_dp)*exp(-par(1)*x**par(2))
        case (bgfd_ew)
            epart = exp(-par(1)*x**par(2))
            value = par(1)*par(2)*x**(par(2)-1.0_dp)*epart*par(3) * &
                (1.0_dp-epart)**(par(3)-1.0_dp)
        case (bgfd_f)
            z = (x/par(1))**par(2)
            value = par(2)*x**(par(2)-1.0_dp)*par(1)**(-par(2))/(1.0_dp+z)**2
        case (bgfd_l)
            value = (par(2)/par(1))*(1.0_dp+x/par(1))**(-par(2)-1.0_dp)
        case (bgfd_b)
            z = (x/par(1))**par(2)
            value = par(3)*par(2)*x**(par(2)-1.0_dp)*par(1)**(-par(2)) * &
                (1.0_dp+z)**(-par(3)-1.0_dp)
        case (bgfd_bx)
            epart = exp(-x*x)
            value = 2.0_dp*par(1)*x*epart*(1.0_dp-epart)**(par(1)-1.0_dp)
        case default
            value = ieee_value(1.0_dp, ieee_quiet_nan)
        end select
    end function base_pdf

    pure function base_quantile(family, p, par) result(value)
        integer, intent(in) :: family
        real(dp), intent(in) :: p, par(:)
        real(dp) :: value, u

        if (.not. base_params_valid(family, par) .or. p < 0.0_dp .or. p > 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        if (p <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (p >= 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_positive_inf)
            return
        end if

        select case (family)
        case (bgfd_e)
            value = -log1p_stable(-p)/par(1)
        case (bgfd_ee)
            u = p**(1.0_dp/par(2))
            value = -log1p_stable(-u)/par(1)
        case (bgfd_w)
            value = (-log1p_stable(-p)/par(1))**(1.0_dp/par(2))
        case (bgfd_ew)
            u = p**(1.0_dp/par(3))
            value = (-log1p_stable(-u)/par(1))**(1.0_dp/par(2))
        case (bgfd_f)
            value = par(1)*(p/(1.0_dp-p))**(1.0_dp/par(2))
        case (bgfd_l)
            value = par(1)*((1.0_dp-p)**(-1.0_dp/par(2))-1.0_dp)
        case (bgfd_b)
            value = par(1)*((1.0_dp-p)**(-1.0_dp/par(3))-1.0_dp)**(1.0_dp/par(2))
        case (bgfd_bx)
            u = p**(1.0_dp/par(1))
            value = sqrt(-log1p_stable(-u))
        case default
            value = ieee_value(1.0_dp, ieee_quiet_nan)
        end select
    end function base_quantile

    pure function bgfd_cdf(family, complementary, x, par) result(value)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: x, par(:)
        real(dp) :: value, u, lambda, e_lambda, den, z, log_den, log_num
        integer :: np

        np = bgfd_npar(family)
        if (np == 0 .or. size(par) /= np .or. any(par <= 0.0_dp)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        lambda = par(np)
        u = base_cdf(family, x, par(:np-1))
        if (u >= 1.0_dp) then
            value = 1.0_dp
            return
        end if

        e_lambda = exp(lambda)
        if (.not. complementary) then
            den = -expm1_stable(1.0_dp-e_lambda)
            z = -e_lambda*(1.0_dp-exp(-lambda*u))
            value = (-expm1_stable(z))/den
        else
            log_den = log_expm1_pos(e_lambda-1.0_dp)
            log_num = log_expm1_pos(exp(lambda*u)-1.0_dp)
            value = exp(log_num-log_den)
        end if
        value = min(1.0_dp, max(0.0_dp, value))
    end function bgfd_cdf

    pure function bgfd_pdf(family, complementary, x, par, log_value) result(value)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: x, par(:)
        logical, intent(in), optional :: log_value
        real(dp) :: value, u, g, lambda, e_lambda, den, log_den, logf
        integer :: np
        logical :: want_log

        want_log = .false.
        if (present(log_value)) want_log = log_value
        np = bgfd_npar(family)
        if (np == 0 .or. size(par) /= np .or. any(par <= 0.0_dp)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        if (x < 0.0_dp) then
            if (want_log) then
                value = -ieee_value(1.0_dp, ieee_positive_inf)
            else
                value = 0.0_dp
            end if
            return
        end if
        lambda = par(np)
        u = base_cdf(family, x, par(:np-1))
        g = base_pdf(family, x, par(:np-1))
        if (g <= 0.0_dp) then
            if (want_log) then
                value = -ieee_value(1.0_dp, ieee_positive_inf)
            else
                value = 0.0_dp
            end if
            return
        end if

        e_lambda = exp(lambda)
        if (.not. complementary) then
            den = -expm1_stable(1.0_dp-e_lambda)
            logf = log(lambda) + log(g) + lambda*(1.0_dp-u) - &
                e_lambda*(1.0_dp-exp(-lambda*u)) - log(den)
        else
            log_den = log_expm1_pos(e_lambda-1.0_dp)
            logf = log(lambda) + log(g) + lambda*u + exp(lambda*u)-1.0_dp - log_den
        end if
        if (want_log) then
            value = logf
        else
            value = exp(logf)
        end if
    end function bgfd_pdf

    pure function bgfd_quantile(family, complementary, p, par, log_p, lower_tail) result(value)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: p, par(:)
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, prob, u, lambda, e_lambda, den, z, log_den, log_pd, y
        integer :: np
        logical :: lp, lower

        lp = .false.
        lower = .true.
        if (present(log_p)) lp = log_p
        if (present(lower_tail)) lower = lower_tail
        prob = p
        if (lp) prob = exp(prob)
        if (.not. lower) prob = 1.0_dp - prob

        np = bgfd_npar(family)
        if (np == 0 .or. size(par) /= np .or. any(par <= 0.0_dp) .or. &
            prob < 0.0_dp .or. prob > 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        if (prob <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (prob >= 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_positive_inf)
            return
        end if

        lambda = par(np)
        e_lambda = exp(lambda)
        if (.not. complementary) then
            den = -expm1_stable(1.0_dp-e_lambda)
            z = log1p_stable(-prob*den)
            u = -log(1.0_dp + exp(-lambda)*z)/lambda
        else
            log_den = log_expm1_pos(e_lambda-1.0_dp)
            log_pd = log(prob) + log_den
            y = log1pexp(log_pd)
            u = log1p_stable(y)/lambda
        end if
        u = min(1.0_dp, max(0.0_dp, u))
        value = base_quantile(family, u, par(:np-1))
    end function bgfd_quantile

    pure function bgfd_survival(family, complementary, x, par, log_value) result(value)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: x, par(:)
        logical, intent(in), optional :: log_value
        real(dp) :: value, cdf
        logical :: want_log

        want_log = .false.
        if (present(log_value)) want_log = log_value
        cdf = bgfd_cdf(family, complementary, x, par)
        if (want_log) then
            value = log1p_stable(-cdf)
        else
            value = 1.0_dp - cdf
        end if
    end function bgfd_survival

    pure function bgfd_hazard(family, complementary, x, par, log_value) result(value)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: x, par(:)
        logical, intent(in), optional :: log_value
        real(dp) :: value, logh
        logical :: want_log

        want_log = .false.
        if (present(log_value)) want_log = log_value
        logh = bgfd_pdf(family, complementary, x, par, .true.) - &
            bgfd_survival(family, complementary, x, par, .true.)
        if (want_log) then
            value = logh
        else
            value = exp(logh)
        end if
    end function bgfd_hazard

    subroutine bgfd_random(family, complementary, par, x)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        real(dp), intent(in) :: par(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: u(size(x))
        integer :: i

        call random_number(u)
        do i = 1, size(x)
            x(i) = bgfd_quantile(family, complementary, u(i), par)
        end do
    end subroutine bgfd_random

    pure logical function base_params_valid(family, par) result(ok)
        integer, intent(in) :: family
        real(dp), intent(in) :: par(:)
        integer :: need
        select case (family)
        case (bgfd_e, bgfd_bx)
            need = 1
        case (bgfd_ee, bgfd_w, bgfd_f, bgfd_l)
            need = 2
        case (bgfd_ew, bgfd_b)
            need = 3
        case default
            need = -1
        end select
        ok = size(par) == need
        if (ok) ok = all(par > 0.0_dp)
    end function base_params_valid

end module bgfd_core
