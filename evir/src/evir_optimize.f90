! SPDX-License-Identifier: GPL-2.0-or-later
module evir_optimize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use evir_kinds, only : dp
    use evir_math, only : invert_matrix, safe_nan
    implicit none
    private

    abstract interface
        function objective_function(x, context) result(f)
            import dp
            real(dp), intent(in) :: x(:)
            class(*), intent(in) :: context
            real(dp) :: f
        end function objective_function
    end interface

    public :: objective_function, nelder_mead, numerical_hessian
    public :: minimize_scalar_bounded

contains

    subroutine nelder_mead(objective, context, x0, xbest, fbest, converged, &
        max_iter, x_tol, f_tol, initial_step)
        procedure(objective_function) :: objective
        class(*), intent(in) :: context
        real(dp), intent(in) :: x0(:)
        real(dp), intent(out) :: xbest(size(x0))
        real(dp), intent(out) :: fbest
        logical, intent(out) :: converged
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: x_tol, f_tol, initial_step

        integer :: n, miter, iter, i, j, ilo, ihi, inhi
        real(dp) :: xtol, ftol, step, alpha, gamma, rho, shrink
        real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)
        real(dp) :: fr, fe, fc, spread_x, spread_f

        n = size(x0)
        miter = 4000
        if (present(max_iter)) miter = max_iter
        xtol = 1.0e-9_dp
        if (present(x_tol)) xtol = x_tol
        ftol = 1.0e-10_dp
        if (present(f_tol)) ftol = f_tol
        step = 0.1_dp
        if (present(initial_step)) step = initial_step
        alpha = 1.0_dp
        gamma = 2.0_dp
        rho = 0.5_dp
        shrink = 0.5_dp

        allocate(simplex(n, n+1), f(n+1), centroid(n), xr(n), xe(n), xc(n))
        simplex(:, 1) = x0
        do j = 2, n+1
            simplex(:, j) = x0
            i = j - 1
            simplex(i, j) = x0(i) + step*max(1.0_dp, abs(x0(i)))
            if (abs(simplex(i, j)-x0(i)) <= epsilon(1.0_dp)) simplex(i, j) = x0(i) + step
        end do
        do j = 1, n+1
            f(j) = objective(simplex(:, j), context)
            if (.not. (f(j) < huge(1.0_dp))) f(j) = huge(1.0_dp)/100.0_dp
        end do

        converged = .false.
        do iter = 1, miter
            call order_vertices(f, ilo, ihi, inhi)
            spread_f = maxval(abs(f - f(ilo)))
            spread_x = 0.0_dp
            do j = 1, n+1
                spread_x = max(spread_x, maxval(abs(simplex(:, j)-simplex(:, ilo))))
            end do
            if (spread_f <= ftol*(1.0_dp+abs(f(ilo))) .and. &
                spread_x <= xtol*(1.0_dp+maxval(abs(simplex(:, ilo))))) then
                converged = .true.
                exit
            end if

            centroid = 0.0_dp
            do j = 1, n+1
                if (j /= ihi) centroid = centroid + simplex(:, j)
            end do
            centroid = centroid / real(n, dp)
            xr = centroid + alpha*(centroid-simplex(:, ihi))
            fr = objective(xr, context)

            if (fr < f(ilo)) then
                xe = centroid + gamma*(xr-centroid)
                fe = objective(xe, context)
                if (fe < fr) then
                    simplex(:, ihi) = xe
                    f(ihi) = fe
                else
                    simplex(:, ihi) = xr
                    f(ihi) = fr
                end if
            else if (fr < f(inhi)) then
                simplex(:, ihi) = xr
                f(ihi) = fr
            else
                if (fr < f(ihi)) then
                    xc = centroid + rho*(xr-centroid)
                else
                    xc = centroid + rho*(simplex(:, ihi)-centroid)
                end if
                fc = objective(xc, context)
                if (fc < min(fr, f(ihi))) then
                    simplex(:, ihi) = xc
                    f(ihi) = fc
                else
                    do j = 1, n+1
                        if (j == ilo) cycle
                        simplex(:, j) = simplex(:, ilo) + shrink*(simplex(:, j)-simplex(:, ilo))
                        f(j) = objective(simplex(:, j), context)
                    end do
                end if
            end if
        end do
        call order_vertices(f, ilo, ihi, inhi)
        xbest = simplex(:, ilo)
        fbest = f(ilo)
    contains
        subroutine order_vertices(values, low, high, next_high)
            real(dp), intent(in) :: values(:)
            integer, intent(out) :: low, high, next_high
            integer :: k
            low = 1
            high = 1
            do k = 2, size(values)
                if (values(k) < values(low)) low = k
                if (values(k) > values(high)) high = k
            end do
            if (high == 1) then
                next_high = 2
            else
                next_high = 1
            end if
            do k = 1, size(values)
                if (k /= high .and. values(k) > values(next_high)) next_high = k
            end do
        end subroutine order_vertices
    end subroutine nelder_mead

    subroutine numerical_hessian(objective, context, x, hessian, covariance, status)
        procedure(objective_function) :: objective
        class(*), intent(in) :: context
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: hessian(size(x), size(x))
        real(dp), intent(out) :: covariance(size(x), size(x))
        integer, intent(out) :: status

        integer :: n, i, j, inv_status
        logical :: bad
        real(dp) :: f0, fpp, fpm, fmp, fmm, fp, fm
        real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), h(:)

        n = size(x)
        allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), h(n))
        do i = 1, n
            h(i) = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(x(i)))
        end do
        f0 = objective(x, context)
        hessian = 0.0_dp
        bad = .not. valid_objective(f0)
        do i = 1, n
            xp = x
            xm = x
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            fp = objective(xp, context)
            fm = objective(xm, context)
            if (valid_objective(fp) .and. valid_objective(fm) .and. valid_objective(f0)) then
                hessian(i, i) = (fp - 2.0_dp*f0 + fm)/(h(i)*h(i))
            else
                bad = .true.
                hessian(i, i) = 0.0_dp
            end if
            do j = i+1, n
                xpp = x
                xpm = x
                xmp = x
                xmm = x
                xpp(i) = xpp(i)+h(i); xpp(j) = xpp(j)+h(j)
                xpm(i) = xpm(i)+h(i); xpm(j) = xpm(j)-h(j)
                xmp(i) = xmp(i)-h(i); xmp(j) = xmp(j)+h(j)
                xmm(i) = xmm(i)-h(i); xmm(j) = xmm(j)-h(j)
                fpp = objective(xpp, context)
                fpm = objective(xpm, context)
                fmp = objective(xmp, context)
                fmm = objective(xmm, context)
                if (valid_objective(fpp) .and. valid_objective(fpm) .and. &
                    valid_objective(fmp) .and. valid_objective(fmm)) then
                    hessian(i, j) = (fpp-fpm-fmp+fmm)/(4.0_dp*h(i)*h(j))
                    hessian(j, i) = hessian(i, j)
                else
                    bad = .true.
                    hessian(i, j) = 0.0_dp
                    hessian(j, i) = 0.0_dp
                end if
            end do
        end do
        if (bad) then
            status = 1
            covariance = safe_nan()
            return
        end if
        call invert_matrix(hessian, covariance, inv_status)
        if (inv_status /= 0) then
            status = 1
            covariance = safe_nan()
        else
            status = 0
        end if
    end subroutine numerical_hessian

    subroutine minimize_scalar_bounded(objective, context, lower, upper, xmin, fmin, converged, tol, max_iter)
        procedure(objective_function) :: objective
        class(*), intent(in) :: context
        real(dp), intent(in) :: lower, upper
        real(dp), intent(out) :: xmin, fmin
        logical, intent(out) :: converged
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter

        real(dp), parameter :: gr = 0.6180339887498948482_dp
        real(dp) :: a, b, c, d, fc, fd, tolerance
        integer :: iter, miter
        real(dp) :: xv(1)

        tolerance = 1.0e-9_dp
        if (present(tol)) tolerance = tol
        miter = 500
        if (present(max_iter)) miter = max_iter
        a = lower
        b = upper
        c = b - gr*(b-a)
        d = a + gr*(b-a)
        xv(1) = c; fc = objective(xv, context)
        xv(1) = d; fd = objective(xv, context)
        converged = .false.
        do iter = 1, miter
            if (abs(b-a) <= tolerance*(1.0_dp+abs(a)+abs(b))) then
                converged = .true.
                exit
            end if
            if (fc < fd) then
                b = d
                d = c
                fd = fc
                c = b - gr*(b-a)
                xv(1) = c; fc = objective(xv, context)
            else
                a = c
                c = d
                fc = fd
                d = a + gr*(b-a)
                xv(1) = d; fd = objective(xv, context)
            end if
        end do
        if (fc < fd) then
            xmin = c
            fmin = fc
        else
            xmin = d
            fmin = fd
        end if
    end subroutine minimize_scalar_bounded

    pure logical function valid_objective(f) result(ok)
        real(dp), intent(in) :: f
        ok = ieee_is_finite(f) .and. abs(f) < huge(1.0_dp)/1000.0_dp
    end function valid_objective

end module evir_optimize
