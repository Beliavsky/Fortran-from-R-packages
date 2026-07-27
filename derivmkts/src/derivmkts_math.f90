! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_math
    use derivmkts_kinds, only: dp, pi
    implicit none
    private
    public :: norm_pdf, norm_cdf, bivar_norm_cdf, clamp, mean_value, sample_sd
    public :: covariance, variance, linear_interp, cholesky_lower

contains

    pure elemental real(dp) function clamp(x, lo, hi) result(y)
        real(dp), intent(in) :: x, lo, hi
        y = min(max(x, lo), hi)
    end function clamp

    pure elemental real(dp) function norm_pdf(x) result(y)
        real(dp), intent(in) :: x
        y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
    end function norm_pdf

    pure elemental real(dp) function norm_cdf(x) result(y)
        real(dp), intent(in) :: x
        y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
    end function norm_cdf

    real(dp) function bivar_norm_cdf(a, b, rho) result(p)
        real(dp), intent(in) :: a, b, rho
        real(dp) :: rr, lower
        rr = clamp(rho, -0.999999999999_dp, 0.999999999999_dp)
        if (a <= -10.0_dp .or. b <= -10.0_dp) then
            p = 0.0_dp
        else if (a >= 10.0_dp) then
            p = norm_cdf(b)
        else if (b >= 10.0_dp) then
            p = norm_cdf(a)
        else if (abs(rr) < 1.0e-14_dp) then
            p = norm_cdf(a)*norm_cdf(b)
        else
            lower = -10.0_dp
            p = adaptive_simpson(lower, a, 1.0e-12_dp, 20)
        end if
    contains
        real(dp) function integrand(x) result(y)
            real(dp), intent(in) :: x
            y = norm_pdf(x)*norm_cdf((b - rr*x)/sqrt(1.0_dp - rr*rr))
        end function integrand

        real(dp) function adaptive_simpson(x0, x1, tol, maxdepth) result(res)
            real(dp), intent(in) :: x0, x1, tol
            integer, intent(in) :: maxdepth
            real(dp) :: fa, fb, fc, whole
            if (x1 <= x0) then
                res = 0.0_dp
                return
            end if
            fa = integrand(x0)
            fb = integrand(x1)
            fc = integrand(0.5_dp*(x0+x1))
            whole = (x1-x0)*(fa + 4.0_dp*fc + fb)/6.0_dp
            res = recurse(x0, x1, fa, fb, fc, whole, tol, maxdepth)
        end function adaptive_simpson

        recursive real(dp) function recurse(x0, x1, fa, fb, fc, whole, tol, depth) result(res)
            real(dp), intent(in) :: x0, x1, fa, fb, fc, whole, tol
            integer, intent(in) :: depth
            real(dp) :: mid, lmid, rmid, fl, fr, left, right, delta
            mid = 0.5_dp*(x0+x1)
            lmid = 0.5_dp*(x0+mid)
            rmid = 0.5_dp*(mid+x1)
            fl = integrand(lmid)
            fr = integrand(rmid)
            left = (mid-x0)*(fa + 4.0_dp*fl + fc)/6.0_dp
            right = (x1-mid)*(fc + 4.0_dp*fr + fb)/6.0_dp
            delta = left + right - whole
            if (depth <= 0 .or. abs(delta) <= 15.0_dp*tol) then
                res = left + right + delta/15.0_dp
            else
                res = recurse(x0, mid, fa, fc, fl, left, 0.5_dp*tol, depth-1) + &
                    recurse(mid, x1, fc, fb, fr, right, 0.5_dp*tol, depth-1)
            end if
        end function recurse
    end function bivar_norm_cdf

    pure real(dp) function mean_value(x) result(m)
        real(dp), intent(in) :: x(:)
        if (size(x) == 0) then
            m = 0.0_dp
        else
            m = sum(x)/real(size(x), dp)
        end if
    end function mean_value

    pure real(dp) function variance(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: m
        if (size(x) <= 1) then
            v = 0.0_dp
        else
            m = mean_value(x)
            v = sum((x-m)**2)/real(size(x)-1, dp)
        end if
    end function variance

    pure real(dp) function sample_sd(x) result(s)
        real(dp), intent(in) :: x(:)
        s = sqrt(max(variance(x), 0.0_dp))
    end function sample_sd

    pure real(dp) function covariance(x, y) result(c)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: mx, my
        integer :: n
        n = min(size(x), size(y))
        if (n <= 1) then
            c = 0.0_dp
        else
            mx = sum(x(1:n))/real(n, dp)
            my = sum(y(1:n))/real(n, dp)
            c = sum((x(1:n)-mx)*(y(1:n)-my))/real(n-1, dp)
        end if
    end function covariance

    pure real(dp) function linear_interp(x, y, x0) result(y0)
        real(dp), intent(in) :: x(:), y(:), x0
        integer :: i, n
        n = min(size(x), size(y))
        if (n < 2) then
            y0 = merge(y(1), 0.0_dp, n == 1)
            return
        end if
        if (x0 <= x(1)) then
            i = 1
        else if (x0 >= x(n)) then
            i = n-1
        else
            do i = 1, n-1
                if (x0 <= x(i+1)) exit
            end do
        end if
        y0 = y(i) + (y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i))
    end function linear_interp

    subroutine cholesky_lower(a, l, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: l(:, :)
        logical, intent(out) :: ok
        integer :: i, j, k, n
        real(dp) :: s
        n = size(a, 1)
        l = 0.0_dp
        ok = size(a, 2) == n .and. size(l,1) == n .and. size(l,2) == n
        if (.not. ok) return
        do i = 1, n
            do j = 1, i
                s = a(i,j)
                do k = 1, j-1
                    s = s - l(i,k)*l(j,k)
                end do
                if (i == j) then
                    if (s <= 0.0_dp) then
                        ok = .false.
                        return
                    end if
                    l(i,j) = sqrt(s)
                else
                    l(i,j) = s/l(j,j)
                end if
            end do
        end do
    end subroutine cholesky_lower

end module derivmkts_math
