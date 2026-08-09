! SPDX-License-Identifier: GPL-2.0-only
module nls2_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use nls2_kinds, only : dp
    use nls2_types, only : nls_control, nls_result, nls_model, nls_jacobian, plinear_basis, &
        nls2_ok, nls2_maxiter, nls2_singular, nls2_bad_input, nls2_model_error
    use nls2_linalg, only : least_squares, invert_spd, norm2_sq
    implicit none
    private
    public :: fit_nls, fit_plinear, evaluate_model, evaluate_plinear

contains

    subroutine evaluate_model(model, x, y, par, result, weights)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), y(:), par(:)
        type(nls_result), intent(out) :: result
        real(dp), intent(in), optional :: weights(:)
        real(dp), allocatable :: yhat(:), resid(:), sw(:)
        integer :: ierr, n

        n = size(y)
        call initialize_result(result, size(par), n)
        if (size(x,1) /= n) return
        allocate(yhat(n), resid(n), sw(n))
        call make_sqrt_weights(n, weights, sw, ierr)
        if (ierr /= 0) return
        call invoke_model(model, x, par, yhat, ierr)
        result%evaluations = 1
        if (ierr /= 0 .or. any(.not. is_finite_vec(yhat))) then
            result%status = nls2_model_error
            return
        end if
        resid = sw * (y - yhat)
        result%par = par
        result%fitted = yhat
        result%residuals = resid
        result%rss = norm2_sq(resid)
        result%status = nls2_ok
        result%converged = .true.
        result%iterations = 0
        result%fin_tol = 0.0_dp
        call finalize_stats(result, sw)
    end subroutine evaluate_model

    subroutine fit_nls(model, x, y, start, result, control, weights, jacobian, lower, upper)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), y(:), start(:)
        type(nls_result), intent(out) :: result
        type(nls_control), intent(in), optional :: control
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        procedure(nls_jacobian), optional :: jacobian
        type(nls_control) :: ctl
        real(dp), allocatable :: par(:), trial(:), yhat(:), ytrial(:), resid(:), rtrial(:)
        real(dp), allocatable :: jac(:,:), wjac(:,:), step(:), sw(:), rmat(:,:), jt_j(:,:), inv(:,:)
        real(dp) :: rss, rss_trial, fac, conv, proj_ss, rem_ss
        integer :: n, p, iter, ierr, evals
        logical :: rank_ok, inv_ok, accepted

        ctl = nls_control()
        if (present(control)) ctl = control
        n = size(y)
        p = size(start)
        call initialize_result(result, p, n)
        if (n <= p .or. size(x,1) /= n .or. p < 1 .or. ctl%maxiter < 1) return
        allocate(par(p), trial(p), yhat(n), ytrial(n), resid(n), rtrial(n), jac(n,p), wjac(n,p))
        allocate(step(p), sw(n), rmat(p,p), jt_j(p,p), inv(p,p))
        call make_sqrt_weights(n, weights, sw, ierr)
        if (ierr /= 0) return
        par = start
        call clamp_params(par, lower, upper, ierr)
        if (ierr /= 0) return
        call invoke_model(model, x, par, yhat, ierr)
        evals = 1
        if (ierr /= 0 .or. any(.not. is_finite_vec(yhat))) then
            result%status = nls2_model_error
            return
        end if
        resid = sw * (y - yhat)
        rss = norm2_sq(resid)
        result%status = nls2_maxiter

        do iter = 1, ctl%maxiter
            if (present(jacobian)) then
                call invoke_jacobian(jacobian, x, par, jac, ierr)
                if (ierr /= 0) then
                    result%status = nls2_model_error
                    exit
                end if
            else
                call numerical_jacobian(model, x, par, yhat, jac, ctl, ierr, evals, upper)
                if (ierr /= 0) then
                    result%status = nls2_model_error
                    exit
                end if
            end if
            wjac = jac
            wjac = spread(sw, 2, p) * wjac
            call least_squares(wjac, resid, step, rank_ok, rmat)
            if (.not. rank_ok) then
                result%status = nls2_singular
                exit
            end if
            proj_ss = norm2_sq(matmul(wjac, step))
            rem_ss = max(rss - proj_ss + ctl%scale_offset, tiny(1.0_dp))
            conv = sqrt(max(proj_ss, 0.0_dp) / rem_ss)
            result%fin_tol = conv
            if (conv <= ctl%tol) then
                result%status = nls2_ok
                result%converged = .true.
                result%iterations = iter - 1
                exit
            end if

            fac = 1.0_dp
            accepted = .false.
            do while (fac >= ctl%min_factor)
                trial = par + fac * step
                call clamp_params(trial, lower, upper, ierr)
                if (ierr /= 0) then
                    result%status = nls2_bad_input
                    exit
                end if
                call invoke_model(model, x, trial, ytrial, ierr)
                evals = evals + 1
                if (ierr == 0 .and. all(is_finite_vec(ytrial))) then
                    rtrial = sw * (y - ytrial)
                    rss_trial = norm2_sq(rtrial)
                    if (rss_trial <= rss) then
                        accepted = .true.
                        exit
                    end if
                end if
                fac = 0.5_dp * fac
            end do
            if (result%status == nls2_bad_input) exit
            if (.not. accepted) then
                result%status = nls2_singular
                exit
            end if

            par = trial
            yhat = ytrial
            resid = rtrial
            rss = rss_trial
            result%iterations = iter
            if (iter == ctl%maxiter) result%status = nls2_maxiter
        end do

        result%par = par
        result%fitted = yhat
        result%residuals = resid
        result%rss = rss
        result%evaluations = evals
        call finalize_stats(result, sw)
        if (result%status == nls2_ok) then
            if (present(jacobian)) then
                call invoke_jacobian(jacobian, x, par, jac, ierr)
            else
                call numerical_jacobian(model, x, par, yhat, jac, ctl, ierr, evals, upper)
            end if
            if (ierr == 0) then
                wjac = spread(sw, 2, p) * jac
                jt_j = matmul(transpose(wjac), wjac)
                call invert_spd(jt_j, inv, inv_ok)
                if (inv_ok) result%covariance = result%sigma**2 * inv
            end if
        end if
    end subroutine fit_nls

    subroutine evaluate_plinear(basis_fn, x, y, theta, result, n_linear, weights)
        procedure(plinear_basis) :: basis_fn
        real(dp), intent(in) :: x(:,:), y(:), theta(:)
        type(nls_result), intent(out) :: result
        integer, intent(in) :: n_linear
        real(dp), intent(in), optional :: weights(:)
        real(dp), allocatable :: basis(:,:), wb(:,:), beta(:), yhat(:), resid(:), sw(:)
        logical :: rank_ok
        integer :: n, ierr

        n = size(y)
        call initialize_result(result, size(theta), n, n_linear)
        if (size(x,1) /= n .or. n_linear < 1 .or. n < n_linear) return
        allocate(basis(n,n_linear), wb(n,n_linear), beta(n_linear), yhat(n), resid(n), sw(n))
        call make_sqrt_weights(n, weights, sw, ierr)
        if (ierr /= 0) return
        call invoke_basis(basis_fn, x, theta, basis, ierr)
        result%evaluations = 1
        if (ierr /= 0 .or. any(.not. is_finite_mat(basis))) then
            result%status = nls2_model_error
            return
        end if
        wb = spread(sw, 2, n_linear) * basis
        call least_squares(wb, sw * y, beta, rank_ok)
        if (.not. rank_ok) then
            result%status = nls2_singular
            return
        end if
        yhat = matmul(basis, beta)
        resid = sw * (y - yhat)
        result%par = theta
        result%linear_par = beta
        result%fitted = yhat
        result%residuals = resid
        result%rss = norm2_sq(resid)
        result%status = nls2_ok
        result%converged = .true.
        result%fin_tol = 0.0_dp
        call finalize_stats(result, sw)
    end subroutine evaluate_plinear

    subroutine fit_plinear(basis_fn, x, y, start, n_linear, result, control, weights, lower, upper)
        procedure(plinear_basis) :: basis_fn
        real(dp), intent(in) :: x(:,:), y(:), start(:)
        integer, intent(in) :: n_linear
        type(nls_result), intent(out) :: result
        type(nls_control), intent(in), optional :: control
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        type(nls_control) :: ctl
        type(nls_result) :: cur, tr
        real(dp), allocatable :: theta(:), trial(:), jac(:,:), step(:), rmat(:,:)
        real(dp) :: h, fac, conv, proj_ss, rem_ss
        integer :: p, n, j, iter, ierr, evals
        logical :: rank_ok, accepted

        ctl = nls_control()
        if (present(control)) ctl = control
        n = size(y)
        p = size(start)
        call initialize_result(result, p, n, n_linear)
        if (p < 1 .or. n_linear < 1 .or. n <= p + n_linear) return
        allocate(theta(p), trial(p), jac(n,p), step(p), rmat(p,p))
        theta = start
        call clamp_params(theta, lower, upper, ierr)
        if (ierr /= 0) return
        call evaluate_plinear(basis_fn, x, y, theta, cur, n_linear, weights)
        evals = 1
        if (cur%status /= nls2_ok) then
            result = cur
            return
        end if

        do iter = 1, ctl%maxiter
            do j = 1, p
                h = ctl%diff_step * max(1.0_dp, abs(theta(j)))
                if (abs(h) <= tiny(1.0_dp)) h = max(ctl%diff_step, sqrt(epsilon(1.0_dp)))
                trial = theta
                trial(j) = theta(j) + h
                if (present(upper)) then
                    if (size(upper) == p .and. trial(j) > upper(j)) trial(j) = theta(j) - h
                end if
                call evaluate_plinear(basis_fn, x, y, trial, tr, n_linear, weights)
                evals = evals + 1
                if (tr%status /= nls2_ok) then
                    result%status = nls2_model_error
                    result = cur
                    return
                end if
                jac(:,j) = (tr%fitted - cur%fitted) / (trial(j) - theta(j))
            end do
            ! Weighted residuals are stored; the finite-difference fitted values are not weighted.
            call apply_weights_to_jac(jac, weights)
            call least_squares(jac, cur%residuals, step, rank_ok, rmat)
            if (.not. rank_ok) then
                result = cur
                result%status = nls2_singular
                result%converged = .false.
                result%evaluations = evals
                return
            end if
            proj_ss = norm2_sq(matmul(jac, step))
            rem_ss = max(cur%rss - proj_ss + ctl%scale_offset, tiny(1.0_dp))
            conv = sqrt(max(proj_ss, 0.0_dp) / rem_ss)
            cur%fin_tol = conv
            if (conv <= ctl%tol) then
                cur%status = nls2_ok
                cur%converged = .true.
                cur%iterations = iter - 1
                exit
            end if
            fac = 1.0_dp
            accepted = .false.
            do while (fac >= ctl%min_factor)
                trial = theta + fac * step
                call clamp_params(trial, lower, upper, ierr)
                if (ierr /= 0) exit
                call evaluate_plinear(basis_fn, x, y, trial, tr, n_linear, weights)
                evals = evals + 1
                if (tr%status == nls2_ok .and. tr%rss <= cur%rss) then
                    accepted = .true.
                    exit
                end if
                fac = 0.5_dp * fac
            end do
            if (.not. accepted) then
                result = cur
                result%status = nls2_singular
                result%converged = .false.
                result%evaluations = evals
                return
            end if
            theta = trial
            cur = tr
            cur%iterations = iter
            if (iter == ctl%maxiter) then
                cur%status = nls2_maxiter
                cur%converged = .false.
            end if
        end do
        cur%evaluations = evals
        result = cur
    end subroutine fit_plinear

    subroutine numerical_jacobian(model, x, par, y0, jac, ctl, ierr, evals, upper)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), par(:), y0(:)
        real(dp), intent(out) :: jac(:,:)
        type(nls_control), intent(in) :: ctl
        integer, intent(out) :: ierr
        integer, intent(inout) :: evals
        real(dp), intent(in), optional :: upper(:)
        real(dp), allocatable :: pp(:), pm(:), yp(:), ym(:)
        real(dp) :: h, den
        integer :: j, n, p, ie

        n = size(y0)
        p = size(par)
        allocate(pp(p), pm(p), yp(n), ym(n))
        ierr = 0
        do j = 1, p
            h = ctl%diff_step * max(1.0_dp, abs(par(j)))
            if (abs(h) <= tiny(1.0_dp)) h = max(ctl%diff_step, sqrt(epsilon(1.0_dp)))
            pp = par
            pm = par
            if (ctl%central_diff) then
                pp(j) = par(j) + h
                pm(j) = par(j) - h
                call invoke_model(model, x, pp, yp, ie)
                evals = evals + 1
                if (ie /= 0) then
                    ierr = ie
                    return
                end if
                call invoke_model(model, x, pm, ym, ie)
                evals = evals + 1
                if (ie /= 0) then
                    ierr = ie
                    return
                end if
                jac(:,j) = (yp - ym) / (2.0_dp*h)
            else
                pp(j) = par(j) + h
                if (present(upper)) then
                    if (size(upper) == p .and. pp(j) > upper(j)) pp(j) = par(j) - h
                end if
                den = pp(j) - par(j)
                call invoke_model(model, x, pp, yp, ie)
                evals = evals + 1
                if (ie /= 0 .or. abs(den) <= tiny(1.0_dp)) then
                    ierr = max(ie, 1)
                    return
                end if
                jac(:,j) = (yp - y0) / den
            end if
        end do
    end subroutine numerical_jacobian

    subroutine apply_weights_to_jac(jac, weights)
        real(dp), intent(inout) :: jac(:,:)
        real(dp), intent(in), optional :: weights(:)
        integer :: i
        if (.not. present(weights)) return
        do i = 1, size(jac,1)
            jac(i,:) = sqrt(weights(i)) * jac(i,:)
        end do
    end subroutine apply_weights_to_jac

    subroutine make_sqrt_weights(n, weights, sw, ierr)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: weights(:)
        real(dp), intent(out) :: sw(n)
        integer, intent(out) :: ierr
        ierr = 0
        sw = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n .or. any(weights < 0.0_dp)) then
                ierr = 1
                return
            end if
            sw = sqrt(weights)
        end if
    end subroutine make_sqrt_weights

    subroutine clamp_params(par, lower, upper, ierr)
        real(dp), intent(inout) :: par(:)
        real(dp), intent(in), optional :: lower(:), upper(:)
        integer, intent(out) :: ierr
        ierr = 0
        if (present(lower)) then
            if (size(lower) /= size(par)) then
                ierr = 1
                return
            end if
            par = max(par, lower)
        end if
        if (present(upper)) then
            if (size(upper) /= size(par)) then
                ierr = 1
                return
            end if
            par = min(par, upper)
        end if
        if (present(lower) .and. present(upper)) then
            if (any(lower > upper)) ierr = 1
        end if
    end subroutine clamp_params

    subroutine initialize_result(result, p, n, q)
        type(nls_result), intent(out) :: result
        integer, intent(in) :: p, n
        integer, intent(in), optional :: q
        integer :: nq
        nq = 0
        if (present(q)) nq = q
        allocate(result%par(p), result%fitted(n), result%residuals(n), result%covariance(p,p))
        if (nq > 0) allocate(result%linear_par(nq))
        result%par = 0.0_dp
        if (nq > 0) result%linear_par = 0.0_dp
        result%fitted = 0.0_dp
        result%residuals = 0.0_dp
        result%covariance = 0.0_dp
        result%rss = huge(1.0_dp)
        result%sigma = huge(1.0_dp)
        result%status = nls2_bad_input
        result%converged = .false.
    end subroutine initialize_result

    subroutine finalize_stats(result, sw)
        type(nls_result), intent(inout) :: result
        real(dp), intent(in) :: sw(:)
        integer :: n_eff, p_total, df
        n_eff = count(sw > 0.0_dp)
        p_total = size(result%par)
        if (allocated(result%linear_par)) p_total = p_total + size(result%linear_par)
        df = n_eff - p_total
        if (df > 0 .and. result%rss >= 0.0_dp) result%sigma = sqrt(result%rss / real(df,dp))
    end subroutine finalize_stats

    subroutine invoke_model(model, x, par, yhat, ierr)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: yhat(:)
        integer, intent(out) :: ierr
        call model(x, par, yhat, ierr)
    end subroutine invoke_model

    subroutine invoke_jacobian(jacobian, x, par, jac, ierr)
        procedure(nls_jacobian) :: jacobian
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: jac(:,:)
        integer, intent(out) :: ierr
        call jacobian(x, par, jac, ierr)
    end subroutine invoke_jacobian

    subroutine invoke_basis(basis_fn, x, theta, basis, ierr)
        procedure(plinear_basis) :: basis_fn
        real(dp), intent(in) :: x(:,:), theta(:)
        real(dp), intent(out) :: basis(:,:)
        integer, intent(out) :: ierr
        call basis_fn(x, theta, basis, ierr)
    end subroutine invoke_basis

    elemental logical function is_finite_scalar(x) result(ok)
        real(dp), intent(in) :: x
        ok = ieee_is_finite(x)
    end function is_finite_scalar

    pure function is_finite_vec(x) result(ok)
        real(dp), intent(in) :: x(:)
        logical :: ok(size(x))
        integer :: i
        do i = 1, size(x)
            ok(i) = is_finite_scalar(x(i))
        end do
    end function is_finite_vec

    pure function is_finite_mat(x) result(ok)
        real(dp), intent(in) :: x(:,:)
        logical :: ok(size(x,1), size(x,2))
        integer :: i, j
        do j = 1, size(x,2)
            do i = 1, size(x,1)
                ok(i,j) = is_finite_scalar(x(i,j))
            end do
        end do
    end function is_finite_mat

end module nls2_core
