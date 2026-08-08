! SPDX-License-Identifier: MIT
! Translated from the computational core of the R package trust 0.1-9.
module trust_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_negative_inf
    use trust_kinds, only : dp
    use trust_types, only : trust_options, trust_result, trust_history, trust_objective, &
        trust_step_newton, trust_step_easy_easy, trust_step_hard_easy, trust_step_hard_hard, &
        trust_ok, trust_err_bad_input, trust_err_initial_objective, trust_err_objective, &
        trust_err_eigensolver
    use trust_linalg, only : symmetric_eigen_jacobi, vector_norm
    implicit none
    private
    public :: trust_optimize, trust_step_name
contains

    function trust_step_name(step_type) result(name)
        integer, intent(in) :: step_type
        character(len=12) :: name
        select case (step_type)
        case (trust_step_newton)
            name = 'Newton'
        case (trust_step_easy_easy)
            name = 'easy-easy'
        case (trust_step_hard_easy)
            name = 'hard-easy'
        case (trust_step_hard_hard)
            name = 'hard-hard'
        case default
            name = 'unknown'
        end select
    end function trust_step_name

    subroutine trust_optimize(objfun, parinit, options, result)
        procedure(trust_objective) :: objfun
        real(dp), intent(in) :: parinit(:)
        type(trust_options), intent(in) :: options
        type(trust_result), intent(out) :: result

        real(dp), allocatable :: theta(:), theta_try(:), g(:), b(:, :), p(:), &
            evals(:), evecs(:, :), gq(:), scale(:), grad_raw(:), hess_raw(:, :)
        real(dp) :: f, ftry, value_raw, radius, preddiff, rho, stepnorm
        real(dp) :: neg_inf
        integer :: d, iiter, stat, eigstat, step_type
        logical :: accept, terminate, have_scale

        d = size(parinit)
        call initialize_result(result, d)
        if (d < 1 .or. .not. all_finite(parinit)) then
            call fail(result, trust_err_bad_input, 'parinit must be a nonempty finite vector')
            return
        end if
        if (.not. ieee_is_finite(options%rinit) .or. options%rinit <= 0.0_dp .or. &
            .not. ieee_is_finite(options%rmax) .or. options%rmax < options%rinit .or. &
            options%iterlim < 1 .or. options%fterm < 0.0_dp .or. options%mterm < 0.0_dp) then
            call fail(result, trust_err_bad_input, 'invalid trust-region controls')
            return
        end if

        allocate(theta(d), theta_try(d), g(d), b(d, d), p(d), evals(d), evecs(d, d), &
            gq(d), scale(d), grad_raw(d), hess_raw(d, d))
        theta = parinit
        radius = options%rinit
        accept = .true.
        terminate = .false.
        have_scale = allocated(options%parscale)
        if (have_scale) then
            if (size(options%parscale) /= d .or. any(options%parscale <= 0.0_dp) .or. &
                .not. all_finite(options%parscale) .or. &
                .not. all_finite(1.0_dp / options%parscale)) then
                call fail(result, trust_err_bad_input, 'parscale must be positive, finite, and conformable')
                return
            end if
            scale = options%parscale
        else
            scale = 1.0_dp
        end if

        call evaluate_objective(objfun, theta, value_raw, grad_raw, hess_raw, stat)
        if (stat /= 0) then
            result%argument = theta
            call fail(result, trust_err_initial_objective, 'objective failed at initial point')
            return
        end if
        if (.not. ieee_is_finite(value_raw)) then
            result%argument = theta
            call fail(result, trust_err_bad_input, 'initial point is not feasible')
            return
        end if
        if (.not. valid_derivatives(grad_raw, hess_raw, d)) then
            result%argument = theta
            call fail(result, trust_err_bad_input, 'objective returned invalid gradient or Hessian')
            return
        end if

        if (options%save_history) call allocate_history(result%history, d, options%iterlim)
        neg_inf = ieee_value(0.0_dp, ieee_negative_inf)

        do iiter = 1, options%iterlim
            if (options%save_history) then
                result%history%argument(:, iiter) = theta
                result%history%radius(iiter) = radius
                if (accept) then
                    result%history%value(iiter) = value_raw
                else
                    result%history%value(iiter) = merge(-f, f, .not. options%minimize)
                end if
            end if

            if (accept) then
                f = value_raw
                g = grad_raw / scale
                b = hess_raw
                b = b / spread(scale, 2, d)
                b = b / spread(scale, 1, d)
                if (.not. options%minimize) then
                    f = -f
                    g = -g
                    b = -b
                end if
                call symmetric_eigen_jacobi(b, evals, evecs, eigstat)
                if (eigstat /= 0) then
                    call fail(result, trust_err_eigensolver, 'symmetric eigensolver did not converge')
                    result%argument = theta
                    result%iterations = iiter - 1
                    return
                end if
                gq = matmul(transpose(evecs), g)
            end if

            call solve_trust_subproblem(evals, evecs, gq, radius, p, step_type)
            preddiff = dot_product(p, g + 0.5_dp * matmul(b, p))
            theta_try = theta + p / scale
            stepnorm = vector_norm(p)

            call evaluate_objective(objfun, theta_try, ftry, grad_raw, hess_raw, stat)
            if (stat /= 0) then
                result%argument = theta
                result%iterations = iiter
                result%converged = .false.
                call fail(result, trust_err_objective, 'objective failed at a trial point')
                if (options%save_history) then
                    result%history%n = iiter - 1
                    call trim_history(result%history)
                end if
                call finalize_at_point(objfun, theta, result)
                return
            end if
            if (ieee_is_nan(ftry) .or. &
                (.not. ieee_is_finite(ftry) .and. options%minimize .and. ftry < 0.0_dp) .or. &
                (.not. ieee_is_finite(ftry) .and. .not. options%minimize .and. ftry > 0.0_dp)) then
                call fail(result, trust_err_bad_input, 'objective returned invalid infinite or NaN value')
                result%argument = theta
                result%iterations = iiter
                return
            end if
            if (ieee_is_finite(ftry) .and. .not. valid_derivatives(grad_raw, hess_raw, d)) then
                call fail(result, trust_err_bad_input, 'objective returned invalid gradient or Hessian')
                result%argument = theta
                result%iterations = iiter
                return
            end if

            if (.not. options%minimize) ftry = -ftry
            if (abs(preddiff) > tiny(1.0_dp)) then
                rho = (ftry - f) / preddiff
            else if (abs(ftry - f) <= tiny(1.0_dp)) then
                rho = 1.0_dp
            else
                rho = neg_inf
            end if

            if (ieee_is_finite(ftry)) then
                terminate = abs(ftry - f) < options%fterm .or. abs(preddiff) < options%mterm
            else
                terminate = .false.
                rho = neg_inf
            end if

            if (terminate) then
                if (ftry < f) then
                    accept = .true.
                    theta = theta_try
                    value_raw = merge(-ftry, ftry, .not. options%minimize)
                end if
            else if (rho < 0.25_dp) then
                accept = .false.
                radius = radius / 4.0_dp
            else
                accept = .true.
                theta = theta_try
                value_raw = merge(-ftry, ftry, .not. options%minimize)
                if (rho > 0.75_dp .and. step_type /= trust_step_newton) &
                    radius = min(2.0_dp * radius, options%rmax)
            end if

            if (options%save_history) then
                result%history%argument_try(:, iiter) = theta_try
                result%history%value_try(iiter) = merge(-ftry, ftry, .not. options%minimize)
                result%history%accepted(iiter) = accept
                result%history%predicted_difference(iiter) = &
                    merge(-preddiff, preddiff, .not. options%minimize)
                result%history%step_norm(iiter) = stepnorm
                result%history%step_type(iiter) = step_type
                result%history%rho(iiter) = rho
            end if

            if (terminate) exit

            if (accept) then
                ! The trial derivatives already belong to the new accepted point.
                ! value_raw has been restored to the user's minimize/maximize sign above.
                if (.not. valid_derivatives(grad_raw, hess_raw, d)) then
                    call fail(result, trust_err_bad_input, 'objective returned invalid derivatives')
                    result%argument = theta
                    result%iterations = iiter
                    return
                end if
            end if
        end do

        result%iterations = min(iiter, options%iterlim)
        result%argument = theta
        result%converged = terminate
        result%status = trust_ok
        if (terminate) then
            result%message = 'converged'
        else
            result%message = 'iteration limit reached'
        end if
        call finalize_at_point(objfun, theta, result)
        if (options%save_history) then
            result%history%n = min(iiter, options%iterlim)
            call trim_history(result%history)
        end if
    end subroutine trust_optimize

    subroutine evaluate_objective(callback, x, value, gradient, hessian, status)
        procedure(trust_objective) :: callback
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value
        real(dp), intent(out) :: gradient(:)
        real(dp), intent(out) :: hessian(:, :)
        integer, intent(out) :: status
        call callback(x, value, gradient, hessian, status)
    end subroutine evaluate_objective

    subroutine finalize_at_point(callback, x, result)
        procedure(trust_objective) :: callback
        real(dp), intent(in) :: x(:)
        type(trust_result), intent(inout) :: result
        real(dp) :: v
        real(dp), allocatable :: g(:), h(:, :)
        integer :: stat, n
        n = size(x)
        allocate(g(n), h(n, n))
        call evaluate_objective(callback, x, v, g, h, stat)
        if (stat == 0 .and. ieee_is_finite(v) .and. valid_derivatives(g, h, n)) then
            result%value = v
            result%gradient = g
            result%hessian = h
        else if (result%status == trust_ok) then
            result%status = trust_err_objective
            result%message = 'objective failed at final point'
            result%converged = .false.
        end if
    end subroutine finalize_at_point

    subroutine solve_trust_subproblem(values, vectors, gq, radius, p, step_type)
        real(dp), intent(in) :: values(:), vectors(:, :), gq(:), radius
        real(dp), intent(out) :: p(:)
        integer, intent(out) :: step_type
        real(dp), allocatable :: beta(:), w(:)
        logical, allocatable :: imin(:)
        real(dp) :: lambda_min, c1, c2, c3, beta_dn, beta_up, root, utry
        real(dp) :: tol, denom
        integer :: n, first_min

        n = size(values)
        allocate(beta(n), w(n), imin(n))
        tol = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))

        if (all(values > 0.0_dp)) then
            w = -gq / values
            p = matmul(vectors, w)
            if (vector_norm(p) <= radius) then
                step_type = trust_step_newton
                return
            end if
        end if

        lambda_min = minval(values)
        beta = values - lambda_min
        imin = abs(beta) <= tol
        c1 = 0.0_dp
        where (.not. imin) w = gq / beta
        where (imin) w = 0.0_dp
        c1 = sum(pack(w * w, .not. imin))
        c2 = sum(pack(gq * gq, imin))
        c3 = sum(gq * gq)

        if (c2 > 0.0_dp .or. c1 > radius * radius) then
            if (c2 > 0.0_dp) then
                step_type = trust_step_easy_easy
            else
                step_type = trust_step_hard_easy
            end if
            beta_dn = sqrt(c2) / radius
            beta_up = sqrt(c3) / radius
            if (secular_value(beta, gq, beta_up, radius) <= 0.0_dp) then
                root = beta_up
            else if (secular_value(beta, gq, beta_dn, radius) >= 0.0_dp) then
                root = beta_dn
            else
                call secular_bisect(beta, gq, radius, beta_dn, beta_up, root)
            end if
            do n = 1, size(w)
                denom = beta(n) + root
                if (abs(denom) <= tiny(1.0_dp)) then
                    w(n) = 0.0_dp
                else
                    w(n) = gq(n) / denom
                end if
            end do
            p = -matmul(vectors, w)
        else
            step_type = trust_step_hard_hard
            w = 0.0_dp
            where (.not. imin) w = gq / beta
            p = -matmul(vectors, w)
            utry = sqrt(max(0.0_dp, radius * radius - dot_product(p, p)))
            if (utry > 0.0_dp) then
                first_min = find_first_true(imin)
                if (first_min > 0) p = p + utry * vectors(:, first_min)
            end if
        end if
    end subroutine solve_trust_subproblem

    real(dp) function secular_value(beta, gq, beep, radius) result(v)
        real(dp), intent(in) :: beta(:), gq(:), beep, radius
        real(dp) :: ss, den
        integer :: i
        if (beep <= tiny(1.0_dp)) then
            ss = 0.0_dp
            do i = 1, size(beta)
                if (abs(beta(i)) <= tiny(1.0_dp)) then
                    if (abs(gq(i)) > 0.0_dp) then
                        v = -1.0_dp / radius
                        return
                    end if
                else
                    ss = ss + (gq(i) / beta(i)) ** 2
                end if
            end do
        else
            ss = 0.0_dp
            do i = 1, size(beta)
                den = beta(i) + beep
                ss = ss + (gq(i) / den) ** 2
            end do
        end if
        if (ss <= tiny(1.0_dp)) then
            v = huge(1.0_dp)
        else
            v = sqrt(1.0_dp / ss) - 1.0_dp / radius
        end if
    end function secular_value

    subroutine secular_bisect(beta, gq, radius, lo_in, hi_in, root)
        real(dp), intent(in) :: beta(:), gq(:), radius, lo_in, hi_in
        real(dp), intent(out) :: root
        real(dp) :: lo, hi, mid, fm
        integer :: iter
        lo = lo_in
        hi = hi_in
        do iter = 1, 120
            mid = 0.5_dp * (lo + hi)
            fm = secular_value(beta, gq, mid, radius)
            if (abs(fm) <= 8.0_dp * epsilon(1.0_dp) / max(1.0_dp, radius)) exit
            if (fm > 0.0_dp) then
                hi = mid
            else
                lo = mid
            end if
            if (abs(hi - lo) <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        root = 0.5_dp * (lo + hi)
    end subroutine secular_bisect

    pure integer function find_first_true(mask) result(idx)
        logical, intent(in) :: mask(:)
        integer :: i
        idx = 0
        do i = 1, size(mask)
            if (mask(i)) then
                idx = i
                return
            end if
        end do
    end function find_first_true

    pure logical function all_finite(x) result(ok)
        real(dp), intent(in) :: x(:)
        integer :: i
        ok = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                ok = .false.
                return
            end if
        end do
    end function all_finite

    pure logical function valid_derivatives(g, h, d) result(ok)
        real(dp), intent(in) :: g(:), h(:, :)
        integer, intent(in) :: d
        integer :: i, j
        ok = size(g) == d .and. size(h, 1) == d .and. size(h, 2) == d
        if (.not. ok) return
        do i = 1, d
            if (.not. ieee_is_finite(g(i))) then
                ok = .false.
                return
            end if
            do j = 1, d
                if (.not. ieee_is_finite(h(i, j))) then
                    ok = .false.
                    return
                end if
            end do
        end do
    end function valid_derivatives

    subroutine initialize_result(result, d)
        type(trust_result), intent(out) :: result
        integer, intent(in) :: d
        allocate(result%argument(d), result%gradient(d), result%hessian(d, d))
        result%argument = 0.0_dp
        result%gradient = 0.0_dp
        result%hessian = 0.0_dp
    end subroutine initialize_result

    subroutine fail(result, code, message)
        type(trust_result), intent(inout) :: result
        integer, intent(in) :: code
        character(len=*), intent(in) :: message
        result%status = code
        result%message = message
        result%converged = .false.
    end subroutine fail

    subroutine allocate_history(hist, d, nmax)
        type(trust_history), intent(inout) :: hist
        integer, intent(in) :: d, nmax
        allocate(hist%argument(d, nmax), hist%argument_try(d, nmax), hist%step_type(nmax), &
            hist%accepted(nmax), hist%radius(nmax), hist%step_norm(nmax), hist%rho(nmax), &
            hist%value(nmax), hist%value_try(nmax), hist%predicted_difference(nmax))
        hist%argument = 0.0_dp
        hist%argument_try = 0.0_dp
        hist%step_type = 0
        hist%accepted = .false.
        hist%radius = 0.0_dp
        hist%step_norm = 0.0_dp
        hist%rho = 0.0_dp
        hist%value = 0.0_dp
        hist%value_try = 0.0_dp
        hist%predicted_difference = 0.0_dp
    end subroutine allocate_history

    subroutine trim_history(hist)
        type(trust_history), intent(inout) :: hist
        integer :: n
        n = hist%n
        if (.not. allocated(hist%radius)) return
        if (n < size(hist%radius)) then
            hist%argument = hist%argument(:, :n)
            hist%argument_try = hist%argument_try(:, :n)
            hist%step_type = hist%step_type(:n)
            hist%accepted = hist%accepted(:n)
            hist%radius = hist%radius(:n)
            hist%step_norm = hist%step_norm(:n)
            hist%rho = hist%rho(:n)
            hist%value = hist%value(:n)
            hist%value_try = hist%value_try(:n)
            hist%predicted_difference = hist%predicted_difference(:n)
        end if
    end subroutine trim_history

end module trust_core
