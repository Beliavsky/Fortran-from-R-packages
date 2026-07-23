! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_optimize
    use longmemo_kinds, only : dp
    implicit none
    private

    public :: minimize_scalar, nelder_mead

    abstract interface
        function scalar_objective(x, context) result(f)
            import dp
            real(dp), intent(in) :: x
            class(*), intent(in) :: context
            real(dp) :: f
        end function scalar_objective

        function vector_objective(x, context) result(f)
            import dp
            real(dp), intent(in) :: x(:)
            class(*), intent(in) :: context
            real(dp) :: f
        end function vector_objective
    end interface

contains

    subroutine minimize_scalar(fun, context, lower, upper, xmin, fmin, tol, max_iter)
        procedure(scalar_objective) :: fun
        class(*), intent(in) :: context
        real(dp), intent(in) :: lower, upper
        real(dp), intent(out) :: xmin, fmin
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        real(dp), parameter :: golden = 0.3819660112501051518_dp
        real(dp) :: a, b, x, w, v, fx, fw, fv, d, e, p, q, r
        real(dp) :: u, fu, xm, tol1, tol2, tolerance
        integer :: iter, niter

        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        niter = 500
        if (present(max_iter)) niter = max_iter

        a = min(lower, upper)
        b = max(lower, upper)
        x = a + golden*(b - a)
        w = x
        v = x
        fx = fun(x, context)
        fw = fx
        fv = fx
        d = 0.0_dp
        e = 0.0_dp

        do iter = 1, niter
            xm = 0.5_dp*(a + b)
            tol1 = tolerance*abs(x) + 1.0e-12_dp
            tol2 = 2.0_dp*tol1
            if (abs(x - xm) <= tol2 - 0.5_dp*(b - a)) exit

            if (abs(e) > tol1) then
                r = (x - w)*(fx - fv)
                q = (x - v)*(fx - fw)
                p = (x - v)*q - (x - w)*r
                q = 2.0_dp*(q - r)
                if (q > 0.0_dp) p = -p
                q = abs(q)
                r = e
                e = d
                if (abs(p) < abs(0.5_dp*q*r) .and. p > q*(a - x) .and. p < q*(b - x)) then
                    d = p/q
                    u = x + d
                    if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
                else
                    e = merge(a - x, b - x, x >= xm)
                    d = golden*e
                end if
            else
                e = merge(a - x, b - x, x >= xm)
                d = golden*e
            end if

            u = x + merge(d, sign(tol1, d), abs(d) >= tol1)
            fu = fun(u, context)
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
                if (fu <= fw .or. same_value(w, x)) then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                else if (fu <= fv .or. same_value(v, x) .or. same_value(v, w)) then
                    v = u
                    fv = fu
                end if
            end if
        end do

        xmin = x
        fmin = fx
    end subroutine minimize_scalar


    subroutine nelder_mead(fun, context, start, xmin, fmin, converged, step, tol, max_iter)
        procedure(vector_objective) :: fun
        class(*), intent(in) :: context
        real(dp), intent(in) :: start(:)
        real(dp), allocatable, intent(out) :: xmin(:)
        real(dp), intent(out) :: fmin
        logical, intent(out) :: converged
        real(dp), intent(in), optional :: step(:), tol
        integer, intent(in), optional :: max_iter
        real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)
        real(dp) :: alpha, gamma, rho, sigma, tolerance, spread, fspread
        integer :: n, j, best, worst, second_worst, iter, niter

        n = size(start)
        allocate(simplex(n, n + 1), f(n + 1), centroid(n), xr(n), xe(n), xc(n), xmin(n))
        simplex(:, 1) = start
        do j = 1, n
            simplex(:, j + 1) = start
            if (present(step)) then
                simplex(j, j + 1) = simplex(j, j + 1) + step(j)
            else
                simplex(j, j + 1) = simplex(j, j + 1) + merge(0.05_dp*abs(start(j)), 0.05_dp, abs(start(j)) > tiny(1.0_dp))
            end if
        end do
        do j = 1, n + 1
            f(j) = fun(simplex(:, j), context)
        end do

        alpha = 1.0_dp
        gamma = 2.0_dp
        rho = 0.5_dp
        sigma = 0.5_dp
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        niter = 4000
        if (present(max_iter)) niter = max_iter
        converged = .false.

        do iter = 1, niter
            call order_simplex(f, best, worst, second_worst)
            spread = 0.0_dp
            do j = 1, n + 1
                spread = max(spread, maxval(abs(simplex(:, j) - simplex(:, best))))
            end do
            fspread = maxval(abs(f - f(best)))
            if (spread <= tolerance*(1.0_dp + maxval(abs(simplex(:, best)))) .and. &
                fspread <= tolerance*(1.0_dp + abs(f(best)))) then
                converged = .true.
                exit
            end if

            centroid = 0.0_dp
            do j = 1, n + 1
                if (j /= worst) centroid = centroid + simplex(:, j)
            end do
            centroid = centroid/real(n, dp)

            xr = centroid + alpha*(centroid - simplex(:, worst))
            if (fun(xr, context) < f(best)) then
                xe = centroid + gamma*(xr - centroid)
                if (fun(xe, context) < fun(xr, context)) then
                    simplex(:, worst) = xe
                    f(worst) = fun(xe, context)
                else
                    simplex(:, worst) = xr
                    f(worst) = fun(xr, context)
                end if
            else if (fun(xr, context) < f(second_worst)) then
                simplex(:, worst) = xr
                f(worst) = fun(xr, context)
            else
                if (fun(xr, context) < f(worst)) then
                    xc = centroid + rho*(xr - centroid)
                else
                    xc = centroid + rho*(simplex(:, worst) - centroid)
                end if
                if (fun(xc, context) < min(f(worst), fun(xr, context))) then
                    simplex(:, worst) = xc
                    f(worst) = fun(xc, context)
                else
                    do j = 1, n + 1
                        if (j /= best) then
                            simplex(:, j) = simplex(:, best) + sigma*(simplex(:, j) - simplex(:, best))
                            f(j) = fun(simplex(:, j), context)
                        end if
                    end do
                end if
            end if
        end do

        call order_simplex(f, best, worst, second_worst)
        xmin = simplex(:, best)
        fmin = f(best)
    end subroutine nelder_mead


    subroutine order_simplex(f, best, worst, second_worst)
        real(dp), intent(in) :: f(:)
        integer, intent(out) :: best, worst, second_worst
        integer :: j

        best = 1
        worst = 1
        do j = 2, size(f)
            if (f(j) < f(best)) best = j
            if (f(j) > f(worst)) worst = j
        end do
        second_worst = best
        do j = 1, size(f)
            if (j /= worst) then
                if (second_worst == worst .or. f(j) > f(second_worst)) second_worst = j
            end if
        end do
    end subroutine order_simplex


    pure logical function same_value(a, b)
        real(dp), intent(in) :: a, b

        same_value = a <= b .and. a >= b
    end function same_value

end module longmemo_optimize
