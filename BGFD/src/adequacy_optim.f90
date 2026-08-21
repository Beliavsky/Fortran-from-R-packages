! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_optim
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use adequacy_kinds, only: dp
    use adequacy_interfaces, only: objective_fn
    use adequacy_math, only: numeric_gradient, numeric_hessian, invert_matrix
    implicit none
    private

    type, public :: optimize_result
        real(dp), allocatable :: par(:)
        real(dp), allocatable :: hessian(:, :)
        real(dp), allocatable :: history(:)
        real(dp) :: value = huge(1.0_dp)
        integer :: iterations = 0
        integer :: convergence = 1
    end type optimize_result

    public :: pso_optimize, nelder_mead_optimize, bfgs_optimize
    public :: cg_optimize, sann_optimize

contains

    subroutine pso_optimize(func, data, lower, upper, result, swarm_size, tol_var, min_history, prop, max_iter)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), lower(:), upper(:)
        type(optimize_result), intent(out) :: result
        integer, intent(in), optional :: swarm_size, min_history, max_iter
        real(dp), intent(in), optional :: tol_var, prop
        integer :: s, npar, nhist, maxit, i, j, iter, nwin, best_idx
        real(dp) :: tol, frac, omega, phi_p, phi_g, r1, r2, span, hist_var
        real(dp), allocatable :: x(:, :), v(:, :), pbest(:, :), f(:), fp(:), hist(:)

        npar = size(lower)
        if (size(upper) /= npar) error stop 'pso_optimize: bound dimensions differ'
        s = 350
        if (present(swarm_size)) s = swarm_size
        nhist = 500
        if (present(min_history)) nhist = min_history
        maxit = 10000
        if (present(max_iter)) maxit = max_iter
        tol = 1.0e-4_dp
        if (present(tol_var)) tol = tol_var
        frac = 0.2_dp
        if (present(prop)) frac = prop

        allocate(x(s, npar), v(s, npar), pbest(s, npar), f(s), fp(s), hist(maxit))
        do i = 1, s
            do j = 1, npar
                call random_number(r1)
                x(i, j) = lower(j) + r1 * (upper(j) - lower(j))
                call random_number(r2)
                span = maxval(upper) - minval(lower)
                v(i, j) = (2.0_dp*r2 - 1.0_dp) * abs(span)
            end do
            f(i) = safe_eval(func, x(i, :), data)
        end do
        pbest = x
        fp = f
        best_idx = minloc(fp, dim=1)
        allocate(result%par(npar))
        result%par = pbest(best_idx, :)
        result%value = fp(best_idx)

        omega = 0.5_dp
        phi_p = 0.5_dp
        phi_g = 0.5_dp
        do iter = 1, maxit
            do i = 1, s
                do j = 1, npar
                    call random_number(r1)
                    call random_number(r2)
                    v(i, j) = omega*v(i, j) + phi_p*r1*(pbest(i, j)-x(i, j)) + &
                              phi_g*r2*(result%par(j)-x(i, j))
                    x(i, j) = x(i, j) + v(i, j)
                    if (x(i, j) < lower(j) .or. x(i, j) > upper(j)) then
                        call random_number(r1)
                        x(i, j) = lower(j) + r1 * (upper(j) - lower(j))
                        v(i, j) = 0.0_dp
                    end if
                end do
                f(i) = safe_eval(func, x(i, :), data)
                if (f(i) <= fp(i)) then
                    pbest(i, :) = x(i, :)
                    fp(i) = f(i)
                end if
                if (fp(i) < result%value) then
                    result%value = fp(i)
                    result%par = pbest(i, :)
                end if
            end do
            hist(iter) = result%value
            if (iter >= nhist) then
                nwin = max(2, ceiling(frac * real(iter, dp)))
                hist_var = variance_tail(hist(iter-nwin+1:iter))
                if (hist_var <= tol) then
                    result%convergence = 0
                    exit
                end if
            end if
        end do
        result%iterations = min(iter, maxit)
        if (iter > maxit) result%convergence = 1
        allocate(result%history(result%iterations))
        result%history = hist(1:result%iterations)
        allocate(result%hessian(npar, npar))
        call numeric_hessian(func, result%par, data, result%hessian)
    end subroutine pso_optimize

    subroutine nelder_mead_optimize(func, data, start, result, tol, max_iter)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), start(:)
        type(optimize_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        integer :: n, nv, i, iter, worst, second_worst, best, maxit
        real(dp) :: eps, alpha, gamma, rho, sigma, fr, fe, fc
        real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)

        n = size(start)
        nv = n + 1
        eps = 1.0e-8_dp
        if (present(tol)) eps = tol
        maxit = 2000
        if (present(max_iter)) maxit = max_iter
        alpha = 1.0_dp
        gamma = 2.0_dp
        rho = 0.5_dp
        sigma = 0.5_dp
        allocate(simplex(nv, n), f(nv), centroid(n), xr(n), xe(n), xc(n))
        simplex(1, :) = start
        do i = 1, n
            simplex(i+1, :) = start
            simplex(i+1, i) = simplex(i+1, i) + 0.05_dp*max(1.0_dp, abs(start(i)))
        end do
        do i = 1, nv
            f(i) = safe_eval(func, simplex(i, :), data)
        end do

        do iter = 1, maxit
            best = minloc(f, dim=1)
            worst = maxloc(f, dim=1)
            second_worst = second_largest_index(f, worst)
            if (maxval(abs(f - f(best))) <= eps*(1.0_dp + abs(f(best)))) exit
            centroid = 0.0_dp
            do i = 1, nv
                if (i /= worst) centroid = centroid + simplex(i, :)
            end do
            centroid = centroid / real(n, dp)
            xr = centroid + alpha*(centroid - simplex(worst, :))
            fr = safe_eval(func, xr, data)
            if (fr < f(best)) then
                xe = centroid + gamma*(xr - centroid)
                fe = safe_eval(func, xe, data)
                if (fe < fr) then
                    simplex(worst, :) = xe
                    f(worst) = fe
                else
                    simplex(worst, :) = xr
                    f(worst) = fr
                end if
            else if (fr < f(second_worst)) then
                simplex(worst, :) = xr
                f(worst) = fr
            else
                if (fr < f(worst)) then
                    xc = centroid + rho*(xr - centroid)
                else
                    xc = centroid + rho*(simplex(worst, :) - centroid)
                end if
                fc = safe_eval(func, xc, data)
                if (fc < min(fr, f(worst))) then
                    simplex(worst, :) = xc
                    f(worst) = fc
                else
                    best = minloc(f, dim=1)
                    do i = 1, nv
                        if (i == best) cycle
                        simplex(i, :) = simplex(best, :) + sigma*(simplex(i, :) - simplex(best, :))
                        f(i) = safe_eval(func, simplex(i, :), data)
                    end do
                end if
            end if
        end do
        best = minloc(f, dim=1)
        allocate(result%par(n), result%hessian(n, n))
        result%par = simplex(best, :)
        result%value = f(best)
        result%iterations = min(iter, maxit)
        result%convergence = merge(0, 1, iter <= maxit)
        call numeric_hessian(func, result%par, data, result%hessian)
    end subroutine nelder_mead_optimize

    subroutine bfgs_optimize(func, data, start, result, tol, max_iter)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), start(:)
        type(optimize_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        integer :: n, i, iter, maxit
        real(dp) :: eps, fcur, fnew, step, ys, rho
        real(dp), allocatable :: x(:), xn(:), g(:), gn(:), p(:), s(:), y(:), h(:, :), eye(:, :), tmp(:, :)

        n = size(start)
        eps = 1.0e-8_dp
        if (present(tol)) eps = tol
        maxit = 1000
        if (present(max_iter)) maxit = max_iter
        allocate(x(n), xn(n), g(n), gn(n), p(n), s(n), y(n), h(n,n), eye(n,n), tmp(n,n))
        x = start
        h = 0.0_dp
        eye = 0.0_dp
        do i = 1, n
            h(i, i) = 1.0_dp
            eye(i, i) = 1.0_dp
        end do
        fcur = safe_eval(func, x, data)
        call numeric_gradient(func, x, data, g)
        do iter = 1, maxit
            if (sqrt(sum(g*g)) <= eps) exit
            p = -matmul(h, g)
            step = backtracking_step(func, data, x, fcur, g, p)
            if (step <= 1.0e-14_dp) exit
            xn = x + step*p
            fnew = safe_eval(func, xn, data)
            call numeric_gradient(func, xn, data, gn)
            s = xn - x
            y = gn - g
            ys = dot_product(y, s)
            if (ys > 1.0e-14_dp) then
                rho = 1.0_dp / ys
                tmp = eye - rho*outer_product(s, y)
                h = matmul(tmp, matmul(h, transpose(tmp))) + rho*outer_product(s, s)
            else
                h = eye
            end if
            x = xn
            g = gn
            fcur = fnew
        end do
        allocate(result%par(n), result%hessian(n,n))
        result%par = x
        result%value = fcur
        result%iterations = min(iter, maxit)
        result%convergence = merge(0, 1, sqrt(sum(g*g)) <= max(eps, 1.0e-6_dp))
        call numeric_hessian(func, result%par, data, result%hessian)
    end subroutine bfgs_optimize

    subroutine cg_optimize(func, data, start, result, tol, max_iter)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), start(:)
        type(optimize_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        integer :: n, iter, maxit
        real(dp) :: eps, beta, step, fcur
        real(dp), allocatable :: x(:), xn(:), g(:), gn(:), d(:)

        n = size(start)
        eps = 1.0e-8_dp
        if (present(tol)) eps = tol
        maxit = 1000
        if (present(max_iter)) maxit = max_iter
        allocate(x(n), xn(n), g(n), gn(n), d(n))
        x = start
        fcur = safe_eval(func, x, data)
        call numeric_gradient(func, x, data, g)
        d = -g
        do iter = 1, maxit
            if (sqrt(sum(g*g)) <= eps) exit
            if (dot_product(g, d) >= 0.0_dp) d = -g
            step = backtracking_step(func, data, x, fcur, g, d)
            if (step <= 1.0e-14_dp) exit
            xn = x + step*d
            fcur = safe_eval(func, xn, data)
            call numeric_gradient(func, xn, data, gn)
            beta = max(0.0_dp, dot_product(gn, gn-g) / max(dot_product(g,g), tiny(1.0_dp)))
            d = -gn + beta*d
            x = xn
            g = gn
        end do
        allocate(result%par(n), result%hessian(n,n))
        result%par = x
        result%value = fcur
        result%iterations = min(iter, maxit)
        result%convergence = merge(0, 1, sqrt(sum(g*g)) <= max(eps, 1.0e-6_dp))
        call numeric_hessian(func, result%par, data, result%hessian)
    end subroutine cg_optimize

    subroutine sann_optimize(func, data, start, result, temp0, cooling, max_iter)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), start(:)
        type(optimize_result), intent(out) :: result
        real(dp), intent(in), optional :: temp0, cooling
        integer, intent(in), optional :: max_iter
        integer :: n, i, iter, maxit
        real(dp) :: temp, cool, fcur, fnew, fbest, u, z
        real(dp), allocatable :: x(:), xn(:), best(:)

        n = size(start)
        maxit = 10000
        if (present(max_iter)) maxit = max_iter
        temp = 10.0_dp
        if (present(temp0)) temp = temp0
        cool = 0.995_dp
        if (present(cooling)) cool = cooling
        allocate(x(n), xn(n), best(n))
        x = start
        fcur = safe_eval(func, x, data)
        best = x
        fbest = fcur
        do iter = 1, maxit
            xn = x
            do i = 1, n
                call random_number(u)
                call random_number(z)
                z = sqrt(-2.0_dp*log(max(u, tiny(1.0_dp)))) * cos(2.0_dp*acos(-1.0_dp)*z)
                xn(i) = x(i) + z*max(1.0_dp, abs(x(i)))*0.1_dp
            end do
            fnew = safe_eval(func, xn, data)
            call random_number(u)
            if (fnew <= fcur .or. u < exp(min(0.0_dp, (fcur-fnew)/max(temp, tiny(1.0_dp))))) then
                x = xn
                fcur = fnew
                if (fcur < fbest) then
                    best = x
                    fbest = fcur
                end if
            end if
            temp = max(temp*cool, 1.0e-10_dp)
        end do
        allocate(result%par(n), result%hessian(n,n))
        result%par = best
        result%value = fbest
        result%iterations = maxit
        result%convergence = 0
        call numeric_hessian(func, result%par, data, result%hessian)
    end subroutine sann_optimize

    function safe_eval(func, par, data) result(v)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: par(:), data(:)
        real(dp) :: v
        v = func(par, data)
        if (.not. ieee_is_finite(v)) v = huge(1.0_dp) / 1000.0_dp
    end function safe_eval

    pure function variance_tail(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: v, m
        if (size(x) < 2) then
            v = 0.0_dp
        else
            m = sum(x) / real(size(x), dp)
            v = sum((x-m)**2) / real(size(x)-1, dp)
        end if
    end function variance_tail

    integer function second_largest_index(x, largest) result(idx)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: largest
        integer :: i
        real(dp) :: val
        idx = merge(2, 1, largest == 1)
        val = x(idx)
        do i = 1, size(x)
            if (i == largest) cycle
            if (x(i) > val) then
                idx = i
                val = x(i)
            end if
        end do
    end function second_largest_index

    pure function outer_product(a, b) result(c)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: c(size(a), size(b))
        integer :: i, j
        do i = 1, size(a)
            do j = 1, size(b)
                c(i, j) = a(i)*b(j)
            end do
        end do
    end function outer_product

    function backtracking_step(func, data, x, fx, g, p) result(step)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: data(:), x(:), fx, g(:), p(:)
        real(dp) :: step, trial
        real(dp) :: xn(size(x))
        step = 1.0_dp
        do while (step > 1.0e-14_dp)
            xn = x + step*p
            trial = safe_eval(func, xn, data)
            if (trial <= fx + 1.0e-4_dp*step*dot_product(g, p)) return
            step = 0.5_dp*step
        end do
    end function backtracking_step

end module adequacy_optim
