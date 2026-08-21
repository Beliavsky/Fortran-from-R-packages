! SPDX-License-Identifier: GPL-2.0-or-later
module bgfd_distributions
    use adequacy_kinds, only: dp
    use bgfd_core, only: bgfd_pdf, bgfd_cdf, bgfd_quantile, bgfd_survival, bgfd_hazard, bgfd_random
    use bgfd_core, only: id_e => bgfd_e, id_ee => bgfd_ee, id_w => bgfd_w, id_ew => bgfd_ew
    use bgfd_core, only: id_f => bgfd_f, id_l => bgfd_l, id_burr => bgfd_b, id_bx => bgfd_bx
    implicit none
    private

    public :: d_bell_e
    public :: p_bell_e
    public :: q_bell_e
    public :: s_bell_e
    public :: h_bell_e
    public :: r_bell_e
    public :: d_bell_ee
    public :: p_bell_ee
    public :: q_bell_ee
    public :: s_bell_ee
    public :: h_bell_ee
    public :: r_bell_ee
    public :: d_bell_w
    public :: p_bell_w
    public :: q_bell_w
    public :: s_bell_w
    public :: h_bell_w
    public :: r_bell_w
    public :: d_bell_ew
    public :: p_bell_ew
    public :: q_bell_ew
    public :: s_bell_ew
    public :: h_bell_ew
    public :: r_bell_ew
    public :: d_bell_f
    public :: p_bell_f
    public :: q_bell_f
    public :: s_bell_f
    public :: h_bell_f
    public :: r_bell_f
    public :: d_bell_l
    public :: p_bell_l
    public :: q_bell_l
    public :: s_bell_l
    public :: h_bell_l
    public :: r_bell_l
    public :: d_bell_b
    public :: p_bell_b
    public :: q_bell_b
    public :: s_bell_b
    public :: h_bell_b
    public :: r_bell_b
    public :: d_bell_bx
    public :: p_bell_bx
    public :: q_bell_bx
    public :: s_bell_bx
    public :: h_bell_bx
    public :: r_bell_bx
    public :: d_cbell_e
    public :: p_cbell_e
    public :: q_cbell_e
    public :: s_cbell_e
    public :: h_cbell_e
    public :: r_cbell_e
    public :: d_cbell_ee
    public :: p_cbell_ee
    public :: q_cbell_ee
    public :: s_cbell_ee
    public :: h_cbell_ee
    public :: r_cbell_ee
    public :: d_cbell_w
    public :: p_cbell_w
    public :: q_cbell_w
    public :: s_cbell_w
    public :: h_cbell_w
    public :: r_cbell_w
    public :: d_cbell_ew
    public :: p_cbell_ew
    public :: q_cbell_ew
    public :: s_cbell_ew
    public :: h_cbell_ew
    public :: r_cbell_ew
    public :: d_cbell_f
    public :: p_cbell_f
    public :: q_cbell_f
    public :: s_cbell_f
    public :: h_cbell_f
    public :: r_cbell_f
    public :: d_cbell_l
    public :: p_cbell_l
    public :: q_cbell_l
    public :: s_cbell_l
    public :: h_cbell_l
    public :: r_cbell_l
    public :: d_cbell_b
    public :: p_cbell_b
    public :: q_cbell_b
    public :: s_cbell_b
    public :: h_cbell_b
    public :: r_cbell_b
    public :: d_cbell_bx
    public :: p_cbell_bx
    public :: q_cbell_bx
    public :: s_cbell_bx
    public :: h_cbell_bx
    public :: r_cbell_bx

contains

    pure function probability_output(cdf, log_p, lower_tail) result(value)
        real(dp), intent(in) :: cdf
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, p
        logical :: lp, lower
        lp = .false.
        lower = .true.
        if (present(log_p)) lp = log_p
        if (present(lower_tail)) lower = lower_tail
        if (lower) then
            p = cdf
        else
            p = 1.0_dp - cdf
        end if
        if (lp) then
            value = log(p)
        else
            value = p
        end if
    end function probability_output

    pure function d_bell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_e, .false., x, [alpha, lambda], log_value)
    end function d_bell_e

    pure function p_bell_e(x, alpha, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_e, .false., x, [alpha, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_e

    pure function q_bell_e(p, alpha, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_e, .false., p, [alpha, lambda], log_p, lower_tail)
    end function q_bell_e

    pure function s_bell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_e, .false., x, [alpha, lambda], log_value)
    end function s_bell_e

    pure function h_bell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_e, .false., x, [alpha, lambda], log_value)
    end function h_bell_e

    subroutine r_bell_e(n, alpha, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_e, .false., [alpha, lambda], x)
    end subroutine r_bell_e

    pure function d_bell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_ee, .false., x, [alpha, beta, lambda], log_value)
    end function d_bell_ee

    pure function p_bell_ee(x, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_ee, .false., x, [alpha, beta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_ee

    pure function q_bell_ee(p, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_ee, .false., p, [alpha, beta, lambda], log_p, lower_tail)
    end function q_bell_ee

    pure function s_bell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_ee, .false., x, [alpha, beta, lambda], log_value)
    end function s_bell_ee

    pure function h_bell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_ee, .false., x, [alpha, beta, lambda], log_value)
    end function h_bell_ee

    subroutine r_bell_ee(n, alpha, beta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_ee, .false., [alpha, beta, lambda], x)
    end subroutine r_bell_ee

    pure function d_bell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_w, .false., x, [alpha, beta, lambda], log_value)
    end function d_bell_w

    pure function p_bell_w(x, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_w, .false., x, [alpha, beta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_w

    pure function q_bell_w(p, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_w, .false., p, [alpha, beta, lambda], log_p, lower_tail)
    end function q_bell_w

    pure function s_bell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_w, .false., x, [alpha, beta, lambda], log_value)
    end function s_bell_w

    pure function h_bell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_w, .false., x, [alpha, beta, lambda], log_value)
    end function h_bell_w

    subroutine r_bell_w(n, alpha, beta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_w, .false., [alpha, beta, lambda], x)
    end subroutine r_bell_w

    pure function d_bell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_ew, .false., x, [alpha, beta, theta, lambda], log_value)
    end function d_bell_ew

    pure function p_bell_ew(x, alpha, beta, theta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_ew, .false., x, [alpha, beta, theta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_ew

    pure function q_bell_ew(p, alpha, beta, theta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_ew, .false., p, [alpha, beta, theta, lambda], log_p, lower_tail)
    end function q_bell_ew

    pure function s_bell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_ew, .false., x, [alpha, beta, theta, lambda], log_value)
    end function s_bell_ew

    pure function h_bell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_ew, .false., x, [alpha, beta, theta, lambda], log_value)
    end function h_bell_ew

    subroutine r_bell_ew(n, alpha, beta, theta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, theta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_ew, .false., [alpha, beta, theta, lambda], x)
    end subroutine r_bell_ew

    pure function d_bell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_f, .false., x, [a, b, lambda], log_value)
    end function d_bell_f

    pure function p_bell_f(x, a, b, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_f, .false., x, [a, b, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_f

    pure function q_bell_f(p, a, b, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_f, .false., p, [a, b, lambda], log_p, lower_tail)
    end function q_bell_f

    pure function s_bell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_f, .false., x, [a, b, lambda], log_value)
    end function s_bell_f

    pure function h_bell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_f, .false., x, [a, b, lambda], log_value)
    end function h_bell_f

    subroutine r_bell_f(n, a, b, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_f, .false., [a, b, lambda], x)
    end subroutine r_bell_f

    pure function d_bell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_l, .false., x, [b, q, lambda], log_value)
    end function d_bell_l

    pure function p_bell_l(x, b, q, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_l, .false., x, [b, q, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_l

    pure function q_bell_l(p, b, q, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_l, .false., p, [b, q, lambda], log_p, lower_tail)
    end function q_bell_l

    pure function s_bell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_l, .false., x, [b, q, lambda], log_value)
    end function s_bell_l

    pure function h_bell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_l, .false., x, [b, q, lambda], log_value)
    end function h_bell_l

    subroutine r_bell_l(n, b, q, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: b, q, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_l, .false., [b, q, lambda], x)
    end subroutine r_bell_l

    pure function d_bell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_burr, .false., x, [a, b, k, lambda], log_value)
    end function d_bell_b

    pure function p_bell_b(x, a, b, k, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_burr, .false., x, [a, b, k, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_b

    pure function q_bell_b(p, a, b, k, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_burr, .false., p, [a, b, k, lambda], log_p, lower_tail)
    end function q_bell_b

    pure function s_bell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_burr, .false., x, [a, b, k, lambda], log_value)
    end function s_bell_b

    pure function h_bell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_burr, .false., x, [a, b, k, lambda], log_value)
    end function h_bell_b

    subroutine r_bell_b(n, a, b, k, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, k, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_burr, .false., [a, b, k, lambda], x)
    end subroutine r_bell_b

    pure function d_bell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_bx, .false., x, [a, lambda], log_value)
    end function d_bell_bx

    pure function p_bell_bx(x, a, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_bx, .false., x, [a, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_bell_bx

    pure function q_bell_bx(p, a, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_bx, .false., p, [a, lambda], log_p, lower_tail)
    end function q_bell_bx

    pure function s_bell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_bx, .false., x, [a, lambda], log_value)
    end function s_bell_bx

    pure function h_bell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_bx, .false., x, [a, lambda], log_value)
    end function h_bell_bx

    subroutine r_bell_bx(n, a, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_bx, .false., [a, lambda], x)
    end subroutine r_bell_bx

    pure function d_cbell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_e, .true., x, [alpha, lambda], log_value)
    end function d_cbell_e

    pure function p_cbell_e(x, alpha, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_e, .true., x, [alpha, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_e

    pure function q_cbell_e(p, alpha, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_e, .true., p, [alpha, lambda], log_p, lower_tail)
    end function q_cbell_e

    pure function s_cbell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_e, .true., x, [alpha, lambda], log_value)
    end function s_cbell_e

    pure function h_cbell_e(x, alpha, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_e, .true., x, [alpha, lambda], log_value)
    end function h_cbell_e

    subroutine r_cbell_e(n, alpha, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_e, .true., [alpha, lambda], x)
    end subroutine r_cbell_e

    pure function d_cbell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_ee, .true., x, [alpha, beta, lambda], log_value)
    end function d_cbell_ee

    pure function p_cbell_ee(x, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_ee, .true., x, [alpha, beta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_ee

    pure function q_cbell_ee(p, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_ee, .true., p, [alpha, beta, lambda], log_p, lower_tail)
    end function q_cbell_ee

    pure function s_cbell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_ee, .true., x, [alpha, beta, lambda], log_value)
    end function s_cbell_ee

    pure function h_cbell_ee(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_ee, .true., x, [alpha, beta, lambda], log_value)
    end function h_cbell_ee

    subroutine r_cbell_ee(n, alpha, beta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_ee, .true., [alpha, beta, lambda], x)
    end subroutine r_cbell_ee

    pure function d_cbell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_w, .true., x, [alpha, beta, lambda], log_value)
    end function d_cbell_w

    pure function p_cbell_w(x, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_w, .true., x, [alpha, beta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_w

    pure function q_cbell_w(p, alpha, beta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_w, .true., p, [alpha, beta, lambda], log_p, lower_tail)
    end function q_cbell_w

    pure function s_cbell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_w, .true., x, [alpha, beta, lambda], log_value)
    end function s_cbell_w

    pure function h_cbell_w(x, alpha, beta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_w, .true., x, [alpha, beta, lambda], log_value)
    end function h_cbell_w

    subroutine r_cbell_w(n, alpha, beta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_w, .true., [alpha, beta, lambda], x)
    end subroutine r_cbell_w

    pure function d_cbell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_ew, .true., x, [alpha, beta, theta, lambda], log_value)
    end function d_cbell_ew

    pure function p_cbell_ew(x, alpha, beta, theta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_ew, .true., x, [alpha, beta, theta, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_ew

    pure function q_cbell_ew(p, alpha, beta, theta, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_ew, .true., p, [alpha, beta, theta, lambda], log_p, lower_tail)
    end function q_cbell_ew

    pure function s_cbell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_ew, .true., x, [alpha, beta, theta, lambda], log_value)
    end function s_cbell_ew

    pure function h_cbell_ew(x, alpha, beta, theta, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: alpha, beta, theta, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_ew, .true., x, [alpha, beta, theta, lambda], log_value)
    end function h_cbell_ew

    subroutine r_cbell_ew(n, alpha, beta, theta, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, theta, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_ew, .true., [alpha, beta, theta, lambda], x)
    end subroutine r_cbell_ew

    pure function d_cbell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_f, .true., x, [a, b, lambda], log_value)
    end function d_cbell_f

    pure function p_cbell_f(x, a, b, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_f, .true., x, [a, b, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_f

    pure function q_cbell_f(p, a, b, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_f, .true., p, [a, b, lambda], log_p, lower_tail)
    end function q_cbell_f

    pure function s_cbell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_f, .true., x, [a, b, lambda], log_value)
    end function s_cbell_f

    pure function h_cbell_f(x, a, b, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_f, .true., x, [a, b, lambda], log_value)
    end function h_cbell_f

    subroutine r_cbell_f(n, a, b, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_f, .true., [a, b, lambda], x)
    end subroutine r_cbell_f

    pure function d_cbell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_l, .true., x, [b, q, lambda], log_value)
    end function d_cbell_l

    pure function p_cbell_l(x, b, q, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_l, .true., x, [b, q, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_l

    pure function q_cbell_l(p, b, q, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_l, .true., p, [b, q, lambda], log_p, lower_tail)
    end function q_cbell_l

    pure function s_cbell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_l, .true., x, [b, q, lambda], log_value)
    end function s_cbell_l

    pure function h_cbell_l(x, b, q, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: b, q, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_l, .true., x, [b, q, lambda], log_value)
    end function h_cbell_l

    subroutine r_cbell_l(n, b, q, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: b, q, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_l, .true., [b, q, lambda], x)
    end subroutine r_cbell_l

    pure function d_cbell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_burr, .true., x, [a, b, k, lambda], log_value)
    end function d_cbell_b

    pure function p_cbell_b(x, a, b, k, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_burr, .true., x, [a, b, k, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_b

    pure function q_cbell_b(p, a, b, k, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_burr, .true., p, [a, b, k, lambda], log_p, lower_tail)
    end function q_cbell_b

    pure function s_cbell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_burr, .true., x, [a, b, k, lambda], log_value)
    end function s_cbell_b

    pure function h_cbell_b(x, a, b, k, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, b, k, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_burr, .true., x, [a, b, k, lambda], log_value)
    end function h_cbell_b

    subroutine r_cbell_b(n, a, b, k, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, k, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_burr, .true., [a, b, k, lambda], x)
    end subroutine r_cbell_b

    pure function d_cbell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_pdf(id_bx, .true., x, [a, lambda], log_value)
    end function d_cbell_bx

    pure function p_cbell_bx(x, a, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value, cdf
        cdf = bgfd_cdf(id_bx, .true., x, [a, lambda])
        value = probability_output(cdf, log_p, lower_tail)
    end function p_cbell_bx

    pure function q_cbell_bx(p, a, lambda, log_p, lower_tail) result(value)
        real(dp), intent(in) :: p
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_p, lower_tail
        real(dp) :: value
        value = bgfd_quantile(id_bx, .true., p, [a, lambda], log_p, lower_tail)
    end function q_cbell_bx

    pure function s_cbell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_survival(id_bx, .true., x, [a, lambda], log_value)
    end function s_cbell_bx

    pure function h_cbell_bx(x, a, lambda, log_value) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a, lambda
        logical, intent(in), optional :: log_value
        real(dp) :: value
        value = bgfd_hazard(id_bx, .true., x, [a, lambda], log_value)
    end function h_cbell_bx

    subroutine r_cbell_bx(n, a, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, lambda
        real(dp), intent(out) :: x(n)
        call bgfd_random(id_bx, .true., [a, lambda], x)
    end subroutine r_cbell_bx

end module bgfd_distributions
