! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_tests
    use hyper2_kinds, only : dp
    use hyper2_types, only : hyper2_model
    use hyper2_likelihood, only : loglik_h2, equalp
    use hyper2_optimize, only : fit_result, maxp
    implicit none
    private

    type, public :: hyper2_test_result
        real(dp) :: statistic = 0.0_dp
        real(dp) :: p_value = 1.0_dp
        integer :: df = 0
        real(dp) :: null_support = -huge(1.0_dp)
        real(dp) :: alternative_support = -huge(1.0_dp)
        real(dp), allocatable :: null_estimate(:)
        real(dp), allocatable :: alternative_estimate(:)
        logical :: converged = .false.
    end type hyper2_test_result

    public :: equalp_test, knownp_test, chisq_survival

contains

    function equalp_test(h, startp) result(out)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in), optional :: startp(:)
        type(hyper2_test_result) :: out
        type(fit_result) :: fit

        if (present(startp)) then
            fit = maxp(h, startp=startp)
        else
            fit = maxp(h)
        end if
        out%alternative_estimate = fit%p
        out%alternative_support = fit%log_likelihood
        out%null_estimate = equalp(h%size())
        out%null_support = loglik_h2(out%null_estimate, h)
        out%statistic = out%alternative_support - out%null_support
        out%df = max(0, h%size()-1)
        out%p_value = chisq_survival(max(0.0_dp, 2.0_dp*out%statistic), out%df)
        out%converged = fit%converged
    end function equalp_test

    function knownp_test(h, p, startp) result(out)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in) :: p(:)
        real(dp), intent(in), optional :: startp(:)
        type(hyper2_test_result) :: out
        type(fit_result) :: fit

        if (size(p) /= h%size()) error stop "knownp_test: p size mismatch"
        if (any(p <= 0.0_dp) .or. abs(sum(p)-1.0_dp) > 1.0e-10_dp) then
            error stop "knownp_test: p must be a positive probability vector"
        end if
        if (present(startp)) then
            fit = maxp(h, startp=startp)
        else
            fit = maxp(h)
        end if
        out%alternative_estimate = fit%p
        out%alternative_support = fit%log_likelihood
        out%null_estimate = p
        out%null_support = loglik_h2(p, h)
        out%statistic = out%alternative_support - out%null_support
        out%df = max(0, h%size()-1)
        out%p_value = chisq_survival(max(0.0_dp, 2.0_dp*out%statistic), out%df)
        out%converged = fit%converged
    end function knownp_test

    real(dp) function chisq_survival(x, df) result(q)
        real(dp), intent(in) :: x
        integer, intent(in) :: df
        if (df <= 0) then
            if (x <= 0.0_dp) then
                q = 1.0_dp
            else
                q = 0.0_dp
            end if
            return
        end if
        q = gammq(0.5_dp*real(df,dp), 0.5_dp*max(0.0_dp,x))
    end function chisq_survival

    real(dp) function gammq(a, x) result(q)
        real(dp), intent(in) :: a, x
        if (a <= 0.0_dp .or. x < 0.0_dp) error stop "gammq: invalid argument"
        if (x <= 0.0_dp) then
            q = 1.0_dp
        else if (x < a + 1.0_dp) then
            q = 1.0_dp - gser(a,x)
        else
            q = gcf(a,x)
        end if
        q = min(1.0_dp, max(0.0_dp, q))
    end function gammq

    real(dp) function gser(a, x) result(p)
        real(dp), intent(in) :: a, x
        integer, parameter :: itmax = 10000
        real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
        real(dp) :: ap, del, s
        integer :: n

        ap = a
        del = 1.0_dp/a
        s = del
        do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            s = s + del
            if (abs(del) <= abs(s)*eps) exit
        end do
        p = s*exp(-x + a*log(x) - log_gamma(a))
    end function gser

    real(dp) function gcf(a, x) result(q)
        real(dp), intent(in) :: a, x
        integer, parameter :: itmax = 10000
        real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
        real(dp), parameter :: fpmin = tiny(1.0_dp)/epsilon(1.0_dp)
        real(dp) :: b, c, d, h, an, del
        integer :: i

        b = x + 1.0_dp - a
        c = 1.0_dp/fpmin
        d = 1.0_dp/max(abs(b),fpmin)
        if (b < 0.0_dp) d = -d
        h = d
        do i = 1, itmax
            an = -real(i,dp)*(real(i,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = sign(fpmin,d)
            c = b + an/c
            if (abs(c) < fpmin) c = sign(fpmin,c)
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
        end do
        q = exp(-x + a*log(x) - log_gamma(a))*h
    end function gcf

end module hyper2_tests
