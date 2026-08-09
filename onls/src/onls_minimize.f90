! SPDX-License-Identifier: GPL-2.0-or-later
module onls_minimize
    use onls_kinds, only : dp
    implicit none
    private
    public :: scalar_fn_ctx, brent_minimize

    abstract interface
        function scalar_fn_ctx(x, ctx) result(f)
            import dp
            real(dp), intent(in) :: x
            class(*), intent(inout) :: ctx
            real(dp) :: f
        end function scalar_fn_ctx
    end interface
contains
    subroutine brent_minimize(fn, ctx, ax, bx, xmin, fmin, tol, maxiter)
        procedure(scalar_fn_ctx) :: fn
        class(*), intent(inout) :: ctx
        real(dp), intent(in) :: ax, bx, tol
        real(dp), intent(out) :: xmin, fmin
        integer, intent(in), optional :: maxiter
        real(dp), parameter :: cgold = 0.3819660112501051_dp
        real(dp), parameter :: zeps = 1.0e-15_dp
        real(dp) :: a, b, d, e, etemp, fu, fv, fw, fx
        real(dp) :: p, q, r, tol1, tol2, u, v, w, x, xm
        integer :: iter, itmax

        itmax = 200
        if (present(maxiter)) itmax = maxiter
        a = min(ax,bx)
        b = max(ax,bx)
        x = a + cgold * (b - a)
        w = x
        v = x
        fx = fn(x, ctx)
        fw = fx
        fv = fx
        d = 0.0_dp
        e = 0.0_dp
        do iter = 1, itmax
            xm = 0.5_dp * (a + b)
            tol1 = tol * abs(x) + zeps
            tol2 = 2.0_dp * tol1
            if (abs(x - xm) <= tol2 - 0.5_dp * (b - a)) exit
            if (abs(e) > tol1) then
                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = 2.0_dp * (q - r)
                if (q > 0.0_dp) p = -p
                q = abs(q)
                etemp = e
                e = d
                if (abs(p) >= abs(0.5_dp*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
                    if (x >= xm) then
                        e = a - x
                    else
                        e = b - x
                    end if
                    d = cgold * e
                else
                    d = p / q
                    u = x + d
                    if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
                end if
            else
                if (x >= xm) then
                    e = a - x
                else
                    e = b - x
                end if
                d = cgold * e
            end if
            if (abs(d) >= tol1) then
                u = x + d
            else
                u = x + sign(tol1, d)
            end if
            fu = fn(u, ctx)
            if (fu <= fx) then
                if (u >= x) then
                    a = x
                else
                    b = x
                end if
                v = w
                fv = fw
                w = x
                fw = fx
                x = u
                fx = fu
            else
                if (u < x) then
                    a = u
                else
                    b = u
                end if
                if (fu <= fw .or. abs(w-x) <= epsilon(1.0_dp)) then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                else if (fu <= fv .or. abs(v-x) <= epsilon(1.0_dp) .or. &
                         abs(v-w) <= epsilon(1.0_dp)) then
                    v = u
                    fv = fu
                end if
            end if
        end do
        xmin = x
        fmin = fx
    end subroutine brent_minimize
end module onls_minimize
