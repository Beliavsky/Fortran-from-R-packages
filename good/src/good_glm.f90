! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good_glm
    use good_kinds, only : dp
    use good_special, only : good_series_stats, normal_sf_two_sided, chi_square_sf, quiet_nan
    use good_linalg, only : invert_matrix, outer_product
    implicit none
    private

    type, public :: good_glm_fit
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: vcov(:, :)
        real(dp), allocatable :: hessian(:, :)
        real(dp), allocatable :: fitted_values(:)
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: iterations = 0
        logical :: converged = .false.
        character(len=8) :: link = 'log'
    end type good_glm_fit

    type, public :: good_prediction
        real(dp), allocatable :: fit(:)
        real(dp), allocatable :: se_fit(:)
    end type good_prediction

    type, public :: good_glm_summary
        real(dp), allocatable :: coefficient_table(:, :)
        real(dp), allocatable :: residuals(:)
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp), allocatable :: transformed(:)
        real(dp), allocatable :: restricted_loglik(:)
        real(dp), allocatable :: lrt(:)
        integer, allocatable :: df(:)
        real(dp), allocatable :: p_value(:)
        character(len=16), allocatable :: test_distribution(:)
    end type good_glm_summary

    public :: glm_good
    public :: predict_good
    public :: summary_good
    public :: good_loglik

contains

    subroutine glm_good(y, x, link, fit, start, max_iter, tol)
        integer, intent(in) :: y(:)
        real(dp), intent(in) :: x(:, :)
        character(len=*), intent(in), optional :: link
        type(good_glm_fit), intent(out) :: fit
        real(dp), intent(in), optional :: start(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol

        real(dp), allocatable :: theta(:), g(:), gnew(:), direction(:), trial(:)
        real(dp), allocatable :: hinv(:, :), ident(:, :), svec(:), yvec(:), tmp(:, :)
        real(dp), allocatable :: hess(:, :), gp(:), gm(:), tp(:), tm(:)
        real(dp) :: f, fnew, alpha, c1, dg, ys, rho, tolerance, step
        integer :: p, i, j, iter, itermax, status, invstat
        logical :: accepted
        character(len=8) :: link_name

        if (size(y) /= size(x, 1) .or. size(y) == 0) error stop 'glm_good: incompatible y and x'
        if (any(y < 0)) error stop 'glm_good: response values must be non-negative'

        p = size(x, 2) + 1
        link_name = 'log'
        if (present(link)) link_name = adjustl(link)
        if (.not. valid_link(link_name)) error stop 'glm_good: link must be log, logit, or identity'

        allocate(theta(p), g(p), gnew(p), direction(p), trial(p), svec(p), yvec(p))
        allocate(hinv(p, p), ident(p, p), tmp(p, p))
        ident = 0.0_dp
        do i = 1, p
            ident(i, i) = 1.0_dp
        end do

        if (present(start)) then
            if (size(start) /= p) error stop 'glm_good: start has wrong size'
            theta = start
        else
            theta = 0.0_dp
            theta(1) = -2.0_dp
            select case (trim(link_name))
            case ('log')
                if (all(abs(x(:, 1) - 1.0_dp) < 100.0_dp * epsilon(1.0_dp))) then
                    theta(2) = log(0.5_dp)
                else
                    error stop 'glm_good: log link without intercept requires feasible start'
                end if
            case ('identity')
                if (all(abs(x(:, 1) - 1.0_dp) < 100.0_dp * epsilon(1.0_dp))) then
                    theta(2) = 0.5_dp
                else
                    error stop 'glm_good: identity link without intercept requires feasible start'
                end if
            case ('logit')
                theta(2:) = 0.0_dp
            end select
        end if
        if (.not. feasible(theta, x, link_name)) error stop 'glm_good: initial value not feasible'

        itermax = 500
        if (present(max_iter)) itermax = max(1, max_iter)
        tolerance = 1.0e-7_dp
        if (present(tol)) tolerance = max(tol, 10.0_dp * epsilon(1.0_dp))

        call negloglik_grad(theta, y, x, link_name, f, g, status)
        if (status < 0) error stop 'glm_good: likelihood evaluation failed'
        hinv = ident
        c1 = 1.0e-4_dp
        fit%converged = .false.

        do iter = 1, itermax
            if (maxval(abs(g)) <= tolerance) then
                fit%converged = .true.
                exit
            end if

            direction = -matmul(hinv, g)
            dg = dot_product(g, direction)
            if (dg >= -epsilon(1.0_dp)) then
                direction = -g
                dg = -dot_product(g, g)
                hinv = ident
            end if

            alpha = 1.0_dp
            accepted = .false.
            do while (alpha >= 1.0e-12_dp)
                trial = theta + alpha * direction
                if (feasible(trial, x, link_name)) then
                    call negloglik_grad(trial, y, x, link_name, fnew, gnew, status)
                    if (status >= 0) then
                        if (fnew <= f + c1 * alpha * dg) then
                            accepted = .true.
                            exit
                        end if
                    end if
                end if
                alpha = 0.5_dp * alpha
            end do
            if (.not. accepted) exit

            svec = trial - theta
            yvec = gnew - g
            ys = dot_product(yvec, svec)
            if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(svec, svec) * dot_product(yvec, yvec)))) then
                rho = 1.0_dp / ys
                tmp = matmul(ident - rho * outer_product(svec, yvec), hinv)
                hinv = matmul(tmp, ident - rho * outer_product(yvec, svec)) + rho * outer_product(svec, svec)
            else
                hinv = ident
            end if

            if (maxval(abs(svec)) <= tolerance * (1.0_dp + maxval(abs(theta))) .and. &
                abs(fnew - f) <= tolerance * (1.0_dp + abs(f))) then
                theta = trial
                f = fnew
                g = gnew
                fit%converged = .true.
                exit
            end if

            theta = trial
            f = fnew
            g = gnew
        end do

        fit%iterations = min(iter, itermax)
        fit%link = link_name
        allocate(fit%coefficients(p), fit%fitted_values(size(y)), fit%hessian(p, p), fit%vcov(p, p))
        fit%coefficients = theta
        fit%loglik = -f
        do i = 1, size(y)
            fit%fitted_values(i) = fitted_mean(theta, x(i, :), link_name)
        end do

        allocate(hess(p, p), gp(p), gm(p), tp(p), tm(p))
        hess = 0.0_dp
        do j = 1, p
            step = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp) * max(1.0_dp, abs(theta(j)))
            tp = theta
            tm = theta
            tp(j) = tp(j) + step
            tm(j) = tm(j) - step
            if (feasible(tp, x, link_name) .and. feasible(tm, x, link_name)) then
                call negloglik_grad(tp, y, x, link_name, fnew, gp, status)
                call negloglik_grad(tm, y, x, link_name, fnew, gm, status)
                hess(:, j) = (gp - gm) / (2.0_dp * step)
            else if (feasible(tp, x, link_name)) then
                call negloglik_grad(tp, y, x, link_name, fnew, gp, status)
                hess(:, j) = (gp - g) / step
            else if (feasible(tm, x, link_name)) then
                call negloglik_grad(tm, y, x, link_name, fnew, gm, status)
                hess(:, j) = (g - gm) / step
            end if
        end do
        hess = 0.5_dp * (hess + transpose(hess))
        fit%hessian = hess
        call invert_matrix(hess, fit%vcov, invstat)
        if (invstat /= 0) then
            fit%vcov = quiet_nan()
            fit%converged = .false.
        end if
    end subroutine glm_good

    function good_loglik(y, x, coefficients, link, status) result(ll)
        integer, intent(in) :: y(:)
        real(dp), intent(in) :: x(:, :), coefficients(:)
        character(len=*), intent(in) :: link
        integer, intent(out), optional :: status
        real(dp) :: ll
        real(dp), allocatable :: g(:)
        real(dp) :: f
        integer :: istat

        allocate(g(size(coefficients)))
        call negloglik_grad(coefficients, y, x, link, f, g, istat)
        ll = -f
        if (present(status)) status = istat
    end function good_loglik

    subroutine predict_good(fit, xnew, prediction, with_se)
        type(good_glm_fit), intent(in) :: fit
        real(dp), intent(in) :: xnew(:, :)
        type(good_prediction), intent(out) :: prediction
        logical, intent(in), optional :: with_se

        real(dp) :: eta, z, log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        real(dp) :: dlogz, v
        real(dp), allocatable :: grad(:)
        integer :: i, status, p
        logical :: calc_se

        p = size(fit%coefficients)
        if (size(xnew, 2) /= p - 1) error stop 'predict_good: xnew has wrong number of columns'
        calc_se = .false.
        if (present(with_se)) calc_se = with_se
        allocate(prediction%fit(size(xnew, 1)))
        if (calc_se) allocate(prediction%se_fit(size(xnew, 1)))
        allocate(grad(p))

        do i = 1, size(xnew, 1)
            eta = dot_product(xnew(i, :), fit%coefficients(2:))
            z = inverse_link(eta, fit%link)
            call good_series_stats(z, fit%coefficients(1), log_norm, mean_n, var_n, mean_log_n, &
                                   cov_n_log_n, status)
            prediction%fit(i) = mean_n - 1.0_dp
            if (calc_se) then
                grad(1) = -cov_n_log_n
                dlogz = dlogz_deta(z, fit%link)
                grad(2:) = var_n * dlogz * xnew(i, :)
                v = dot_product(grad, matmul(fit%vcov, grad))
                prediction%se_fit(i) = sqrt(max(0.0_dp, v))
            end if
        end do
    end subroutine predict_good

    subroutine summary_good(fit, y, x, summary)
        type(good_glm_fit), intent(in) :: fit
        integer, intent(in) :: y(:)
        real(dp), intent(in) :: x(:, :)
        type(good_glm_summary), intent(out) :: summary

        type(good_glm_fit) :: null_fit
        real(dp), allocatable :: x0(:, :)
        real(dp) :: se, zstat, beta, zpar, sez, ll_logseries, ll_geom, pgeom
        integer :: i, p, n

        p = size(fit%coefficients)
        n = size(y)
        if (size(x, 1) /= n) error stop 'summary_good: x and y have incompatible rows'
        allocate(summary%coefficient_table(p, 4), summary%residuals(n))
        summary%residuals = real(y, dp) - fit%fitted_values
        summary%loglik = fit%loglik
        summary%aic = 2.0_dp * real(p, dp) - 2.0_dp * fit%loglik
        summary%bic = real(p, dp) * log(real(n, dp)) - 2.0_dp * fit%loglik
        do i = 1, p
            beta = fit%coefficients(i)
            se = sqrt(max(0.0_dp, fit%vcov(i, i)))
            if (se > 0.0_dp) then
                zstat = beta / se
            else
                zstat = quiet_nan()
            end if
            summary%coefficient_table(i, :) = [beta, se, zstat, normal_sf_two_sided(zstat)]
        end do

        if (p == 2) then
            allocate(summary%transformed(2))
            beta = fit%coefficients(2)
            se = sqrt(max(0.0_dp, fit%vcov(2, 2)))
            zpar = inverse_link(beta, fit%link)
            select case (trim(fit%link))
            case ('log')
                sez = zpar * se
            case ('logit')
                sez = zpar * (1.0_dp - zpar) * se
            case default
                sez = se
            end select
            summary%transformed = [zpar, sez]

            allocate(summary%restricted_loglik(2), summary%lrt(2), summary%df(2), summary%p_value(2))
            allocate(character(len=16) :: summary%test_distribution(2))
            ll_logseries = maximize_logseries(y)
            pgeom = real(n, dp) / (real(n, dp) + sum(real(y, dp)))
            pgeom = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, pgeom))
            ll_geom = real(n, dp) * log(pgeom) + sum(real(y, dp)) * log(1.0_dp - pgeom)
            summary%restricted_loglik = [ll_logseries, ll_geom]
            summary%lrt = max(0.0_dp, 2.0_dp * (fit%loglik - summary%restricted_loglik))
            summary%df = 1
            summary%p_value(1) = chi_square_sf(summary%lrt(1), 1)
            summary%p_value(2) = chi_square_sf(summary%lrt(2), 1)
            summary%test_distribution = [character(len=16) :: 'logarithmic', 'geometric']
        else
            allocate(x0(n, 1))
            x0 = 1.0_dp
            call glm_good(y, x0, 'log', null_fit)
            allocate(summary%restricted_loglik(1), summary%lrt(1), summary%df(1), summary%p_value(1))
            allocate(character(len=16) :: summary%test_distribution(1))
            summary%restricted_loglik(1) = null_fit%loglik
            summary%lrt(1) = max(0.0_dp, 2.0_dp * (fit%loglik - null_fit%loglik))
            summary%df(1) = p - 2
            summary%p_value(1) = chi_square_sf(summary%lrt(1), summary%df(1))
            summary%test_distribution(1) = 'good'
        end if
    end subroutine summary_good

    subroutine negloglik_grad(theta, y, x, link, f, grad, status)
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: y(:)
        real(dp), intent(in) :: x(:, :)
        character(len=*), intent(in) :: link
        real(dp), intent(out) :: f, grad(:)
        integer, intent(out) :: status

        real(dp) :: eta, z, log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        real(dp) :: nval, factor, ll
        integer :: i, istat

        if (.not. feasible(theta, x, link)) then
            f = huge(1.0_dp)
            grad = 0.0_dp
            status = -1
            return
        end if

        ll = 0.0_dp
        grad = 0.0_dp
        status = 0
        do i = 1, size(y)
            eta = dot_product(x(i, :), theta(2:))
            z = inverse_link(eta, link)
            call good_series_stats(z, theta(1), log_norm, mean_n, var_n, mean_log_n, &
                                   cov_n_log_n, istat)
            if (istat < 0) then
                f = huge(1.0_dp)
                grad = 0.0_dp
                status = -2
                return
            end if
            if (istat > 0) status = max(status, istat)
            nval = real(y(i) + 1, dp)
            ll = ll + nval * log(z) - theta(1) * log(nval) - log_norm
            grad(1) = grad(1) - (mean_log_n - log(nval))
            factor = (nval - mean_n) * dlogz_deta(z, link)
            grad(2:) = grad(2:) - factor * x(i, :)
        end do
        f = -ll
    end subroutine negloglik_grad

    pure logical function feasible(theta, x, link)
        real(dp), intent(in) :: theta(:), x(:, :)
        character(len=*), intent(in) :: link
        real(dp), allocatable :: eta(:)

        feasible = .false.
        if (size(theta) /= size(x, 2) + 1) return
        if (trim(link) == 'logit') then
            feasible = .true.
            return
        end if
        eta = matmul(x, theta(2:))
        select case (trim(link))
        case ('log')
            feasible = all(eta < -1.0e-12_dp)
        case ('identity')
            feasible = all(eta > 1.0e-12_dp .and. eta < 1.0_dp - 1.0e-12_dp)
        end select
    end function feasible

    pure logical function valid_link(link)
        character(len=*), intent(in) :: link
        valid_link = trim(link) == 'log' .or. trim(link) == 'logit' .or. trim(link) == 'identity'
    end function valid_link

    pure function inverse_link(eta, link) result(z)
        real(dp), intent(in) :: eta
        character(len=*), intent(in) :: link
        real(dp) :: z

        select case (trim(link))
        case ('log')
            z = exp(eta)
        case ('logit')
            if (eta >= 0.0_dp) then
                z = 1.0_dp / (1.0_dp + exp(-eta))
            else
                z = exp(eta) / (1.0_dp + exp(eta))
            end if
        case default
            z = eta
        end select
    end function inverse_link

    pure function dlogz_deta(z, link) result(v)
        real(dp), intent(in) :: z
        character(len=*), intent(in) :: link
        real(dp) :: v

        select case (trim(link))
        case ('log')
            v = 1.0_dp
        case ('logit')
            v = 1.0_dp - z
        case default
            v = 1.0_dp / z
        end select
    end function dlogz_deta

    pure function fitted_mean(theta, xrow, link) result(mu)
        real(dp), intent(in) :: theta(:), xrow(:)
        character(len=*), intent(in) :: link
        real(dp) :: mu
        real(dp) :: z, log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        integer :: status

        z = inverse_link(dot_product(xrow, theta(2:)), link)
        call good_series_stats(z, theta(1), log_norm, mean_n, var_n, mean_log_n, cov_n_log_n, status)
        mu = mean_n - 1.0_dp
    end function fitted_mean

    function maximize_logseries(y) result(llmax)
        integer, intent(in) :: y(:)
        real(dp) :: llmax
        real(dp) :: a, b, c, d, fc, fd, gr
        integer :: iter

        a = 1.0e-8_dp
        b = 1.0_dp - 1.0e-8_dp
        gr = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        fc = logseries_loglik(c, y)
        fd = logseries_loglik(d, y)
        do iter = 1, 200
            if (fc > fd) then
                b = d
                d = c
                fd = fc
                c = b - gr * (b - a)
                fc = logseries_loglik(c, y)
            else
                a = c
                c = d
                fc = fd
                d = a + gr * (b - a)
                fd = logseries_loglik(d, y)
            end if
            if (abs(b - a) < 1.0e-12_dp) exit
        end do
        llmax = max(fc, fd)
    end function maximize_logseries

    pure function logseries_loglik(z, y) result(ll)
        real(dp), intent(in) :: z
        integer, intent(in) :: y(:)
        real(dp) :: ll
        integer :: i

        ll = 0.0_dp
        do i = 1, size(y)
            ll = ll + real(y(i) + 1, dp) * log(z) - log(real(y(i) + 1, dp)) - log(-log(1.0_dp - z))
        end do
    end function logseries_loglik

end module good_glm
