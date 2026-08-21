! Computational translation of the R package bivpois.
! License: GPL-2.0-or-later.
module bivpois_fit
    use bivpois_kinds, only : dp
    use bivpois_math, only : mean_counts, normal_upper_tail, chisq1_upper_tail
    use bivpois_distribution, only : dbp_scalar
    implicit none
    private

    type, public :: bp_mle2_result
        real(dp) :: lambda(3) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        logical :: converged = .false.
    end type bp_mle2_result

    type, public :: bp_mle_result
        real(dp) :: lambda(3) = 0.0_dp
        real(dp) :: rho = 0.0_dp
        real(dp) :: ci(2, 2) = 0.0_dp
        real(dp) :: loglik(2) = -huge(1.0_dp)
        real(dp) :: pvalue(3) = 1.0_dp
        real(dp) :: var_observed = huge(1.0_dp)
        real(dp) :: var_asymptotic = huge(1.0_dp)
        logical :: converged = .false.
    end type bp_mle_result

    type, public :: bp_profile_result
        real(dp), allocatable :: lambda3(:)
        real(dp), allocatable :: loglik(:)
        real(dp) :: ci(2) = 0.0_dp
        real(dp) :: mle = 0.0_dp
        real(dp) :: max_loglik = -huge(1.0_dp)
    end type bp_profile_result

    public :: profile_loglik, bp_mle2, bp_mle, lambda3_profile

contains

    real(dp) function profile_loglik(x1, x2, lambda3) result(ll)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: lambda3
        real(dp) :: m1, m2, lambda(3)
        integer :: i

        if (size(x1) == 0 .or. size(x2) /= size(x1)) error stop "profile_loglik: invalid input"
        m1 = mean_counts(x1)
        m2 = mean_counts(x2)
        if (lambda3 < 0.0_dp .or. lambda3 > min(m1, m2)) then
            ll = -huge(1.0_dp)
            return
        end if
        lambda = [m1 - lambda3, m2 - lambda3, lambda3]
        ll = 0.0_dp
        do i = 1, size(x1)
            ll = ll + dbp_scalar(x1(i), x2(i), lambda, .true.)
        end do
    end function profile_loglik

    function bp_mle2(x1, x2, tol) result(res)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in), optional :: tol
        type(bp_mle2_result) :: res
        real(dp) :: m1, m2, upper, xtol, l3, ll

        if (size(x1) == 0 .or. size(x2) /= size(x1)) error stop "bp_mle2: invalid input"
        if (any(x1 < 0) .or. any(x2 < 0)) error stop "bp_mle2: counts must be nonnegative"
        m1 = mean_counts(x1)
        m2 = mean_counts(x2)
        xtol = 1.0e-6_dp
        if (present(tol)) xtol = tol

        ! The R code uses min(mean1,mean2)-0.05. Preserve that when possible,
        ! but fall back to the mathematically valid open boundary for sparse data.
        upper = min(m1, m2) - 0.05_dp
        if (upper <= 0.0_dp) upper = max(0.0_dp, min(m1, m2) * (1.0_dp - 1.0e-10_dp))
        if (upper <= 0.0_dp) then
            l3 = 0.0_dp
            ll = profile_loglik(x1, x2, l3)
        else
            call maximize_profile(x1, x2, 0.0_dp, upper, xtol, l3, ll)
        end if
        res%lambda = [m1 - l3, m2 - l3, l3]
        res%loglik = ll
        res%converged = .true.
    end function bp_mle2

    function bp_mle(x1, x2, tol) result(res)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in), optional :: tol
        type(bp_mle_result) :: res
        type(bp_mle2_result) :: fast
        real(dp) :: m1, m2, l3, l0, l1, stat, d2ll, upper
        real(dp) :: tau, d2, num, logp1, logp2, term
        real(dp) :: se1, se2, z1, z2
        integer :: r, s, ma

        if (present(tol)) then
            fast = bp_mle2(x1, x2, tol)
        else
            fast = bp_mle2(x1, x2)
        end if
        m1 = mean_counts(x1)
        m2 = mean_counts(x2)
        l3 = fast%lambda(3)
        l1 = fast%loglik
        l0 = profile_loglik(x1, x2, 0.0_dp)
        res%lambda = fast%lambda
        if (m1 > 0.0_dp .and. m2 > 0.0_dp) res%rho = l3 / sqrt(m1 * m2)
        res%loglik = [l0, l1]
        stat = max(0.0_dp, 2.0_dp * (l1 - l0))
        res%pvalue(1) = chisq1_upper_tail(stat)

        upper = min(m1, m2) - 0.05_dp
        if (upper <= 0.0_dp) upper = max(0.0_dp, min(m1, m2) * (1.0_dp - 1.0e-10_dp))
        d2ll = profile_second_derivative(x1, x2, l3, 0.0_dp, max(upper, l3))
        if (d2ll < 0.0_dp) res%var_observed = -1.0_dp / d2ll

        ! Kawamura asymptotic variance calculation, translated from bp.mle().
        ma = max(maxval(x1), maxval(x2)) + 20
        tau = 0.0_dp
        do s = 1, ma
            do r = 1, ma
                logp1 = dbp_scalar(r - 1, s - 1, res%lambda, .true.)
                logp2 = dbp_scalar(r, s, res%lambda, .true.)
                if (logp2 > -0.5_dp * huge(1.0_dp) .and. logp1 > -0.5_dp * huge(1.0_dp)) then
                    term = exp(2.0_dp * logp1 - logp2)
                    if (term < huge(1.0_dp)) tau = tau + term
                end if
            end do
        end do
        d2 = -m1 + l3 - m2 + l3 + (m1 * m2 - l3 * l3) * (tau - 1.0_dp)
        num = (m1 - l3) * (m2 - l3) + l3 * (m1 + m2 - 2.0_dp * l3) * &
              (l3 * (tau - 1.0_dp) - 1.0_dp)
        if (abs(d2) > tiny(1.0_dp)) then
            res%var_asymptotic = (num / d2) / real(size(x1), dp)
            if (res%var_asymptotic <= 0.0_dp) res%var_asymptotic = huge(1.0_dp)
        end if

        if (res%var_observed < huge(1.0_dp) / 2.0_dp) then
            se1 = sqrt(res%var_observed)
            z1 = l3 / se1
            res%pvalue(2) = normal_upper_tail(z1)
            res%ci(1, :) = [l3 - 1.959964_dp * se1, l3 + 1.959964_dp * se1]
        else
            res%pvalue(2) = 1.0_dp
            res%ci(1, :) = [-huge(1.0_dp), huge(1.0_dp)]
        end if
        if (res%var_asymptotic < huge(1.0_dp) / 2.0_dp) then
            se2 = sqrt(res%var_asymptotic)
            z2 = l3 / se2
            res%pvalue(3) = normal_upper_tail(z2)
            res%ci(2, :) = [l3 - 1.959964_dp * se2, l3 + 1.959964_dp * se2]
        else
            res%pvalue(3) = 1.0_dp
            res%ci(2, :) = [-huge(1.0_dp), huge(1.0_dp)]
        end if
        res%converged = fast%converged
    end function bp_mle

    subroutine maximize_profile(x1, x2, lower, upper, tol, xmax, fmax)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: lower, upper, tol
        real(dp), intent(out) :: xmax, fmax
        real(dp), parameter :: gr = 0.6180339887498948482_dp
        real(dp) :: a, b, c, d, fc, fd, fl, fu

        if (upper < lower) error stop "maximize_profile: invalid interval"
        if (upper <= lower) then
            xmax = lower
            fmax = profile_loglik(x1, x2, lower)
            return
        end if
        a = lower
        b = upper
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        fc = profile_loglik(x1, x2, c)
        fd = profile_loglik(x1, x2, d)
        do while (abs(b - a) > tol * (1.0_dp + abs(a) + abs(b)))
            if (fc > fd) then
                b = d
                d = c
                fd = fc
                c = b - gr * (b - a)
                fc = profile_loglik(x1, x2, c)
            else
                a = c
                c = d
                fc = fd
                d = a + gr * (b - a)
                fd = profile_loglik(x1, x2, d)
            end if
        end do
        xmax = 0.5_dp * (a + b)
        fmax = profile_loglik(x1, x2, xmax)
        fl = profile_loglik(x1, x2, lower)
        if (fl > fmax) then
            xmax = lower
            fmax = fl
        end if
        fu = profile_loglik(x1, x2, upper)
        if (fu > fmax) then
            xmax = upper
            fmax = fu
        end if
    end subroutine maximize_profile

    real(dp) function profile_second_derivative(x1, x2, x, lower, upper) result(d2)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: x, lower, upper
        real(dp) :: h, f0, fm, fp, f1, f2, f3

        h = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(x))
        h = min(h, 0.2_dp * max(upper - lower, h))
        f0 = profile_loglik(x1, x2, x)
        if (x - h >= lower .and. x + h <= upper) then
            fm = profile_loglik(x1, x2, x - h)
            fp = profile_loglik(x1, x2, x + h)
            d2 = (fp - 2.0_dp * f0 + fm) / (h * h)
        else if (x + 3.0_dp * h <= upper) then
            f1 = profile_loglik(x1, x2, x + h)
            f2 = profile_loglik(x1, x2, x + 2.0_dp * h)
            f3 = profile_loglik(x1, x2, x + 3.0_dp * h)
            d2 = (2.0_dp * f0 - 5.0_dp * f1 + 4.0_dp * f2 - f3) / (h * h)
        else if (x - 3.0_dp * h >= lower) then
            f1 = profile_loglik(x1, x2, x - h)
            f2 = profile_loglik(x1, x2, x - 2.0_dp * h)
            f3 = profile_loglik(x1, x2, x - 3.0_dp * h)
            d2 = (2.0_dp * f0 - 5.0_dp * f1 + 4.0_dp * f2 - f3) / (h * h)
        else
            d2 = 0.0_dp
        end if
    end function profile_second_derivative

    function lambda3_profile(x1, x2, step) result(prof)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in), optional :: step
        type(bp_profile_result) :: prof
        real(dp) :: h, m1, m2, upper, cutoff
        integer :: ngrid, i, first, last, imax

        if (size(x1) == 0 .or. size(x2) /= size(x1)) error stop "lambda3_profile: invalid input"
        h = 0.01_dp
        if (present(step)) h = step
        if (h <= 0.0_dp) error stop "lambda3_profile: step must be positive"
        m1 = mean_counts(x1)
        m2 = mean_counts(x2)
        upper = max(0.0_dp, min(m1, m2) - 0.1_dp)
        ngrid = floor(upper / h) + 1
        ngrid = max(ngrid, 1)
        allocate(prof%lambda3(ngrid), prof%loglik(ngrid))
        do i = 1, ngrid
            prof%lambda3(i) = real(i - 1, dp) * h
            prof%loglik(i) = profile_loglik(x1, x2, prof%lambda3(i))
        end do
        imax = maxloc(prof%loglik, dim=1)
        prof%mle = prof%lambda3(imax)
        prof%max_loglik = prof%loglik(imax)
        cutoff = prof%max_loglik - 1.920729_dp
        first = 0
        last = 0
        do i = 1, ngrid
            if (prof%loglik(i) >= cutoff) then
                if (first == 0) first = i
                last = i
            end if
        end do
        if (first > 0) then
            prof%ci = [prof%lambda3(first), prof%lambda3(last)]
        else
            prof%ci = [prof%mle, prof%mle]
        end if
    end function lambda3_profile

end module bivpois_fit
