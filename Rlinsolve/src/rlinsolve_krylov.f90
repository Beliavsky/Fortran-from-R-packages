module rlinsolve_krylov
    use rlinsolve_kinds, only : dp
    use rlinsolve_types, only : solver_result, sparse_matrix
    use rlinsolve_utils, only : matrix_operator, preconditioner_type, dense_operator, &
        sparse_operator, sparse_to_dense, op_matvec, op_matvec_t, prepare_system, &
        residual_norm, vec_norm2, setup_preconditioner, apply_preconditioner, &
        apply_preconditioner_t, estimate_preconditioned_bounds, triangular_lower_solve, &
        triangular_upper_solve
    implicit none
    private

    public :: lsolve_cg_dense, lsolve_cg_sparse
    public :: lsolve_bicg_dense, lsolve_bicg_sparse
    public :: lsolve_bicgstab_dense, lsolve_bicgstab_sparse
    public :: lsolve_cgs_dense, lsolve_cgs_sparse
    public :: lsolve_cheby_dense, lsolve_cheby_sparse
    public :: lsolve_gmres_dense, lsolve_gmres_sparse
    public :: lsolve_qmr_dense, lsolve_qmr_sparse

contains

    subroutine initialize_x(op, b, xinit, tol, x, status)
        type(matrix_operator), intent(in) :: op !! Prepared square coefficient operator used to size and test the initial iterate.
        real(dp), intent(in) :: b(:) !! Prepared right-hand-side vector.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial iterate; absence uses a deterministic zero vector.
        real(dp), intent(in) :: tol !! Residual threshold below which upstream replaces the start with a random vector.
        real(dp), allocatable, intent(out) :: x(:) !! Selected initial iterate of length op%ncol.
        integer, intent(out) :: status !! Zero on success and negative when xinit has invalid length.
        integer :: i

        status = 0
        allocate(x(op%ncol))
        if (present(xinit)) then
            if (size(xinit) /= op%ncol) then
                status = -1
                x = 0.0_dp
                return
            end if
            x = xinit
        else
            x = 0.0_dp
        end if
        if (residual_norm(op, x, b) < tol) then
            do i = 1, size(x)
                x(i) = sin(real(i, dp))
            end do
        end if
    end subroutine initialize_x

    subroutine finish_result(x, errors, iter, maxiter, flag, tol, result)
        real(dp), intent(in) :: x(:) !! Final approximation vector.
        real(dp), intent(in) :: errors(:) !! Allocated convergence-history workspace.
        integer, intent(in) :: iter !! Upstream-compatible zero-based outer iteration index.
        integer, intent(in) :: maxiter !! Maximum iteration count used to size errors.
        integer, intent(in) :: flag !! Upstream-style status flag: zero success, one no convergence, negative breakdown.
        real(dp), intent(in) :: tol !! Relative residual tolerance used for convergence.
        type(solver_result), intent(out) :: result !! Packaged Krylov solution and diagnostics.
        integer :: nerr

        allocate(result%x(size(x)))
        result%x = x
        result%iter = iter
        result%flag = flag
        nerr = min(maxiter, max(1, iter + 1))
        allocate(result%errors(nerr))
        result%errors = errors(1:nerr)
        result%converged = flag == 0
        if (nerr > 0 .and. result%errors(nerr) <= tol) result%converged = .true.
    end subroutine finish_result

    subroutine prepare_krylov(input_op, input_b, require_symmetric, adjsym, xinit, tol, &
            preconditioner, op, b, x, pre, status)
        type(matrix_operator), intent(in) :: input_op !! Original dense or sparse coefficient operator.
        real(dp), intent(in) :: input_b(:) !! Original right-hand side with length input_op%nrow.
        logical, intent(in) :: require_symmetric !! Whether nonsymmetric square input is converted to normal equations.
        logical, intent(in) :: adjsym !! Whether the symmetry conversion is allowed.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial iterate for the prepared square system.
        real(dp), intent(in) :: tol !! Iterative stopping tolerance.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner matching the prepared system.
        type(matrix_operator), intent(out) :: op !! Prepared square coefficient operator.
        real(dp), allocatable, intent(out) :: b(:) !! Prepared right-hand-side vector.
        real(dp), allocatable, intent(out) :: x(:) !! Prepared initial iterate.
        type(preconditioner_type), intent(out) :: pre !! Reusable LU-based preconditioner.
        integer, intent(out) :: status !! Zero on success and negative for preprocessing errors.

        call prepare_system(input_op, input_b, require_symmetric, adjsym, op, b, status)
        if (status /= 0) return
        call initialize_x(op, b, xinit, tol, x, status)
        if (status /= 0) return
        call setup_preconditioner(op%ncol, preconditioner, pre, status)
    end subroutine prepare_krylov

    subroutine cg_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Preconditioned conjugate-gradient output.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M applied through M^{-1}; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: z(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: rho
        real(dp) :: rho_old
        real(dp) :: tol
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 10000
        if (present(maxiter)) maxit = maxiter
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .true., dosym, xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        allocate(r(op%nrow), z(op%nrow), p(op%nrow), q(op%nrow), errors(maxit))
        call op_matvec(op, x, q)
        r = b - q
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        rho_old = 1.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            call apply_preconditioner(pre, r, z, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            rho = dot_product(r, z)
            if (iter > 0) then
                if (abs(rho_old) <= tiny(1.0_dp)) then
                    flag = -1
                    exit
                end if
                beta = rho / rho_old
                p = z + beta * p
            else
                p = z
            end if
            call op_matvec(op, p, q)
            denom = dot_product(p, q)
            if (abs(denom) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            alpha = rho / denom
            x = x + alpha * p
            r = r - alpha * q
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err <= tol) then
                flag = 0
                exit
            end if
            rho_old = rho
            iter = iter + 1
        end do
        if (iter >= maxit .and. err > tol) flag = 1
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine cg_core

    subroutine bicg_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Preconditioned BiCG output with breakdown flag.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: pt(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: qt(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: rt(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: zt(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: rho
        real(dp) :: rho_old
        real(dp) :: tol
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 10000
        if (present(maxiter)) maxit = maxiter
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .false., .false., xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        allocate(r(op%nrow), rt(op%nrow), z(op%nrow), zt(op%nrow))
        allocate(p(op%nrow), pt(op%nrow), q(op%nrow), qt(op%nrow), errors(maxit))
        call op_matvec(op, x, q)
        r = b - q
        rt = r
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        rho = 0.0_dp
        rho_old = 1.0_dp
        p = 0.0_dp
        pt = 0.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            call apply_preconditioner(pre, r, z, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            call apply_preconditioner_t(pre, rt, zt, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            rho = dot_product(z, rt)
            if (abs(rho) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            if (iter > 0) then
                beta = rho / rho_old
                p = z + beta * p
                pt = zt + beta * pt
            else
                p = z
                pt = zt
            end if
            call op_matvec(op, p, q)
            call op_matvec_t(op, pt, qt)
            denom = dot_product(pt, q)
            if (abs(denom) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            alpha = rho / denom
            x = x + alpha * p
            r = r - alpha * q
            rt = rt - alpha * qt
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err < tol) then
                flag = 0
                exit
            end if
            rho_old = rho
            iter = iter + 1
        end do
        if (iter >= maxit .and. err > tol) flag = 1
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine bicg_core

    subroutine bicgstab_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Preconditioned BiCGSTAB output with breakdown flags.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: phat(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: rt(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: shat(:)
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: x(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: omega
        real(dp) :: rho
        real(dp) :: rho_old
        real(dp) :: tol
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .false., .false., xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        allocate(r(op%nrow), rt(op%nrow), p(op%nrow), phat(op%nrow))
        allocate(v(op%nrow), s(op%nrow), shat(op%nrow), t(op%nrow), errors(maxit))
        call op_matvec(op, x, v)
        r = b - v
        rt = r
        p = 0.0_dp
        v = 0.0_dp
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        rho_old = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            rho = dot_product(rt, r)
            if (abs(rho) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            if (iter > 0) then
                if (abs(omega) <= tiny(1.0_dp)) then
                    flag = -2
                    exit
                end if
                beta = (rho / rho_old) * (alpha / omega)
                p = r + beta * (p - omega * v)
            else
                p = r
            end if
            call apply_preconditioner(pre, p, phat, istat)
            if (istat /= 0) then
                flag = -3
                exit
            end if
            call op_matvec(op, phat, v)
            denom = dot_product(rt, v)
            if (abs(denom) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            alpha = rho / denom
            s = r - alpha * v
            if (vec_norm2(s) / bnrm2 < tol) then
                x = x + alpha * phat
                err = vec_norm2(s) / bnrm2
                if (iter < maxit - 1) errors(iter + 2) = err
                flag = 0
                exit
            end if
            call apply_preconditioner(pre, s, shat, istat)
            if (istat /= 0) then
                flag = -3
                exit
            end if
            call op_matvec(op, shat, t)
            denom = dot_product(t, t)
            if (abs(denom) <= tiny(1.0_dp)) then
                flag = -2
                exit
            end if
            omega = dot_product(t, s) / denom
            x = x + alpha * phat + omega * shat
            r = s - omega * t
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err < tol) then
                flag = 0
                exit
            end if
            if (abs(omega) <= tiny(1.0_dp)) then
                flag = -2
                exit
            end if
            rho_old = rho
            iter = iter + 1
        end do
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine bicgstab_core

    subroutine cgs_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Preconditioned conjugate-gradient-squared output.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: phat(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: rt(:)
        real(dp), allocatable :: u(:)
        real(dp), allocatable :: uhat(:)
        real(dp), allocatable :: vhat(:)
        real(dp), allocatable :: x(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: rho
        real(dp) :: rho_old
        real(dp) :: tol
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 10000
        if (present(maxiter)) maxit = maxiter
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .true., dosym, xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        allocate(r(op%nrow), rt(op%nrow), u(op%nrow), p(op%nrow), q(op%nrow))
        allocate(phat(op%nrow), vhat(op%nrow), uhat(op%nrow), errors(maxit))
        call op_matvec(op, x, vhat)
        r = b - vhat
        rt = r
        q = 0.0_dp
        p = 0.0_dp
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        rho_old = 1.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            rho = dot_product(rt, r)
            if (abs(rho) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            if (iter > 0) then
                beta = rho / rho_old
                u = r + beta * q
                p = u + beta * (q + beta * p)
            else
                u = r
                p = u
            end if
            call apply_preconditioner(pre, p, phat, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            call op_matvec(op, phat, vhat)
            denom = dot_product(rt, vhat)
            if (abs(denom) <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            alpha = rho / denom
            q = u - alpha * vhat
            call apply_preconditioner(pre, u + q, uhat, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            x = x + alpha * uhat
            call op_matvec(op, uhat, vhat)
            r = r - alpha * vhat
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err <= tol) then
                flag = 0
                exit
            end if
            rho_old = rho
            iter = iter + 1
        end do
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine cgs_core

    subroutine cheby_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner, &
            adjsym, eigmin, eigmax)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Chebyshev semi-iteration output.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        real(dp), intent(in), optional :: eigmin !! Optional smallest eigenvalue of M^{-1}A; estimated when absent.
        real(dp), intent(in), optional :: eigmax !! Optional largest eigenvalue of M^{-1}A; estimated when absent.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: ap(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: z(:)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: c
        real(dp) :: d
        real(dp) :: emax
        real(dp) :: emin
        real(dp) :: err
        real(dp) :: tol
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 10000
        if (present(maxiter)) maxit = maxiter
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .true., dosym, xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        if (present(eigmin) .and. present(eigmax)) then
            emin = eigmin
            emax = eigmax
        else
            call estimate_preconditioned_bounds(op, pre, emin, emax, istat)
            if (istat /= 0) then
                result%flag = -20
                return
            end if
        end if
        if (emin <= 0.0_dp .or. emax <= emin) then
            result%flag = -21
            return
        end if
        allocate(r(op%nrow), z(op%nrow), p(op%nrow), ap(op%nrow), errors(maxit))
        call op_matvec(op, x, ap)
        r = b - ap
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        c = (emax - emin) / 2.0_dp
        d = (emax + emin) / 2.0_dp
        alpha = 0.0_dp
        p = 0.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            call apply_preconditioner(pre, r, z, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            if (iter > 0) then
                beta = (c * alpha / 2.0_dp)**2
                alpha = 1.0_dp / (d - beta)
                p = z + beta * p
            else
                p = z
                alpha = 2.0_dp / d
            end if
            x = x + alpha * p
            call op_matvec(op, p, ap)
            r = r - alpha * ap
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err <= tol) then
                flag = 0
                exit
            end if
            iter = iter + 1
        end do
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine cheby_core

    pure subroutine givens(a, b, c, s)
        real(dp), intent(in) :: a !! First scalar in the two-vector to rotate.
        real(dp), intent(in) :: b !! Second scalar in the two-vector to rotate.
        real(dp), intent(out) :: c !! Cosine-like coefficient of the Givens rotation.
        real(dp), intent(out) :: s !! Sine-like coefficient of the Givens rotation.
        real(dp) :: r

        if (abs(b) <= tiny(1.0_dp)) then
            c = 1.0_dp
            s = 0.0_dp
        else
            r = hypot(a, b)
            c = a / r
            s = b / r
        end if
    end subroutine givens

    subroutine backsolve_upper(h, rhs, n, x, status)
        real(dp), intent(in) :: h(:, :) !! Upper-triangular Hessenberg leading block after Givens rotations.
        real(dp), intent(in) :: rhs(:) !! Right-hand side for the triangular least-squares problem.
        integer, intent(in) :: n !! Leading block order to solve.
        real(dp), intent(out) :: x(:) !! Solution vector; only x(1:n) is defined.
        integer, intent(out) :: status !! Zero on success and negative for a zero diagonal pivot.
        integer :: i

        status = 0
        x = 0.0_dp
        do i = n, 1, -1
            if (abs(h(i, i)) <= tiny(1.0_dp)) then
                status = -1
                return
            end if
            x(i) = rhs(i)
            if (i < n) x(i) = x(i) - dot_product(h(i, i + 1:n), x(i + 1:n))
            x(i) = x(i) / h(i, i)
        end do
    end subroutine backsolve_upper

    subroutine gmres_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner, restart)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Restarted left-preconditioned GMRES output.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of restart cycles; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        integer, intent(in), optional :: restart !! Krylov subspace dimension per cycle; default n-1 and must be at least 2.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: av(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: cs(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: h(:, :)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: svec(:)
        real(dp), allocatable :: sn(:)
        real(dp), allocatable :: v(:, :)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp) :: beta0
        real(dp) :: bnrm2
        real(dp) :: c
        real(dp) :: err
        real(dp) :: s
        real(dp) :: temp
        real(dp) :: tol
        integer :: flag
        integer :: i
        integer :: iter
        integer :: istat
        integer :: j
        integer :: k
        integer :: maxit
        integer :: m

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        m = max(2, input_op%ncol - 1)
        if (present(restart)) m = restart
        if (tol < 0.0_dp .or. maxit < 1 .or. m < 2) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .false., .false., xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        if (m >= op%ncol .and. op%ncol > 2) m = op%ncol - 1
        if (op%ncol <= 2) m = op%ncol
        allocate(v(op%nrow, m + 1), h(m + 1, m), cs(m), sn(m), svec(m + 1))
        allocate(w(op%nrow), av(op%nrow), r(op%nrow), y(m), errors(maxit))
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = residual_norm(op, x, b) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        iter = 0
        flag = 1
        do while (iter < maxit)
            call op_matvec(op, x, av)
            call apply_preconditioner(pre, b - av, r, istat)
            if (istat /= 0) then
                flag = -2
                exit
            end if
            beta0 = vec_norm2(r)
            if (beta0 <= tiny(1.0_dp)) then
                err = 0.0_dp
                flag = 0
                exit
            end if
            v = 0.0_dp
            h = 0.0_dp
            cs = 0.0_dp
            sn = 0.0_dp
            svec = 0.0_dp
            v(:, 1) = r / beta0
            svec(1) = beta0
            k = m
            do j = 1, m
                call op_matvec(op, v(:, j), av)
                call apply_preconditioner(pre, av, w, istat)
                if (istat /= 0) then
                    flag = -2
                    exit
                end if
                do i = 1, j
                    h(i, j) = dot_product(w, v(:, i))
                    w = w - h(i, j) * v(:, i)
                end do
                h(j + 1, j) = vec_norm2(w)
                if (h(j + 1, j) > tiny(1.0_dp)) v(:, j + 1) = w / h(j + 1, j)
                do i = 1, j - 1
                    temp = cs(i) * h(i, j) + sn(i) * h(i + 1, j)
                    h(i + 1, j) = -sn(i) * h(i, j) + cs(i) * h(i + 1, j)
                    h(i, j) = temp
                end do
                call givens(h(j, j), h(j + 1, j), c, s)
                cs(j) = c
                sn(j) = s
                h(j, j) = c * h(j, j) + s * h(j + 1, j)
                h(j + 1, j) = 0.0_dp
                temp = c * svec(j)
                svec(j + 1) = -s * svec(j)
                svec(j) = temp
                err = abs(svec(j + 1)) / bnrm2
                if (err <= tol) then
                    k = j
                    exit
                end if
            end do
            if (flag < 0) exit
            call backsolve_upper(h, svec, k, y, istat)
            if (istat /= 0) then
                flag = -3
                exit
            end if
            x = x + matmul(v(:, 1:k), y(1:k))
            err = residual_norm(op, x, b) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err <= tol) then
                flag = 0
                exit
            end if
            iter = iter + 1
        end do
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine gmres_core

    subroutine split_lu_no_pivot(a, l, u, status)
        real(dp), intent(in) :: a(:, :) !! Square preconditioner matrix to split as L*U without row pivoting.
        real(dp), intent(out) :: l(:, :) !! Unit lower-triangular factor.
        real(dp), intent(out) :: u(:, :) !! Upper-triangular factor.
        integer, intent(out) :: status !! Zero on success and negative when a required pivot is zero.
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        real(dp) :: accum

        n = size(a, 1)
        status = 0
        l = 0.0_dp
        u = 0.0_dp
        do i = 1, n
            l(i, i) = 1.0_dp
        end do
        do k = 1, n
            do j = k, n
                accum = 0.0_dp
                do i = 1, k - 1
                    accum = accum + l(k, i) * u(i, j)
                end do
                u(k, j) = a(k, j) - accum
            end do
            if (abs(u(k, k)) <= tiny(1.0_dp)) then
                status = -1
                return
            end if
            do i = k + 1, n
                accum = 0.0_dp
                do j = 1, k - 1
                    accum = accum + l(i, j) * u(j, k)
                end do
                l(i, k) = (a(i, k) - accum) / u(k, k)
            end do
        end do
    end subroutine split_lu_no_pivot

    subroutine qmr_core(input_op, input_b, result, xinit, reltol, maxiter, preconditioner)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse A, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand side corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! QMR output with upstream-style breakdown flags.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square M; default identity.
        type(matrix_operator) :: op
        type(preconditioner_type) :: pre
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: dvec(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: l(:, :)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: pt(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: svec(:)
        real(dp), allocatable :: u(:, :)
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: vt(:)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: wt(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: yt(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: zt(:)
        real(dp) :: beta
        real(dp) :: bnrm2
        real(dp) :: delta
        real(dp) :: ep
        real(dp) :: err
        real(dp) :: eta
        real(dp) :: gamma
        real(dp) :: gamma_old
        real(dp) :: rho
        real(dp) :: rho_old
        real(dp) :: theta
        real(dp) :: theta_old
        real(dp) :: tol
        real(dp) :: xi
        integer :: flag
        integer :: iter
        integer :: istat
        integer :: maxit
        integer :: n

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_krylov(input_op, input_b, .false., .false., xinit, tol, &
            preconditioner, op, b, x, pre, istat)
        if (istat /= 0) then
            result%flag = istat
            return
        end if
        n = op%nrow
        allocate(l(n, n), u(n, n))
        call split_lu_no_pivot(pre%mat, l, u, istat)
        if (istat /= 0) then
            result%flag = -7
            return
        end if
        allocate(r(n), vt(n), y(n), wt(n), z(n), v(n), w(n), yt(n), zt(n))
        allocate(p(n), q(n), pt(n), dvec(n), svec(n), errors(maxit))
        call op_matvec(op, x, pt)
        r = b - pt
        bnrm2 = vec_norm2(b)
        if (bnrm2 <= tiny(1.0_dp)) bnrm2 = 1.0_dp
        err = vec_norm2(r) / bnrm2
        errors = 0.0_dp
        errors(1) = err
        vt = r
        call triangular_lower_solve(l, vt, y, istat)
        if (istat /= 0) then
            result%flag = -7
            return
        end if
        rho = vec_norm2(y)
        wt = r
        call triangular_lower_solve(transpose(u), wt, z, istat)
        if (istat /= 0) then
            result%flag = -7
            return
        end if
        xi = vec_norm2(z)
        gamma = 1.0_dp
        eta = -1.0_dp
        theta = 0.0_dp
        p = 0.0_dp
        q = 0.0_dp
        dvec = 0.0_dp
        svec = 0.0_dp
        ep = 1.0_dp
        beta = 1.0_dp
        delta = 1.0_dp
        iter = 0
        flag = 1
        do while (iter < maxit)
            if (rho <= tiny(1.0_dp)) then
                flag = -1
                exit
            end if
            if (xi <= tiny(1.0_dp)) then
                flag = -6
                exit
            end if
            v = vt / rho
            y = y / rho
            w = wt / xi
            z = z / xi
            delta = dot_product(z, y)
            if (abs(delta) <= tiny(1.0_dp)) then
                flag = -4
                exit
            end if
            call triangular_upper_solve(u, y, yt, istat)
            if (istat /= 0) then
                flag = -7
                exit
            end if
            call triangular_upper_solve(transpose(l), z, zt, istat)
            if (istat /= 0) then
                flag = -7
                exit
            end if
            if (iter > 0) then
                p = yt - (xi * delta / ep) * p
                q = zt - (rho * delta / ep) * q
            else
                p = yt
                q = zt
            end if
            call op_matvec(op, p, pt)
            ep = dot_product(q, pt)
            if (abs(ep) <= tiny(1.0_dp)) then
                flag = -5
                exit
            end if
            beta = ep / delta
            if (abs(beta) <= tiny(1.0_dp)) then
                flag = -2
                exit
            end if
            vt = pt - beta * v
            call triangular_lower_solve(l, vt, y, istat)
            if (istat /= 0) then
                flag = -7
                exit
            end if
            rho_old = rho
            rho = vec_norm2(y)
            call op_matvec_t(op, q, zt)
            wt = zt - beta * w
            call triangular_lower_solve(transpose(u), wt, z, istat)
            if (istat /= 0) then
                flag = -7
                exit
            end if
            xi = vec_norm2(z)
            gamma_old = gamma
            theta_old = theta
            theta = rho / (gamma_old * beta)
            gamma = 1.0_dp / sqrt(1.0_dp + theta * theta)
            if (gamma <= tiny(1.0_dp)) then
                flag = -3
                exit
            end if
            eta = -eta * rho_old * gamma * gamma / (beta * gamma_old * gamma_old)
            if (iter > 0) then
                dvec = eta * p + (theta_old * gamma)**2 * dvec
                svec = eta * pt + (theta_old * gamma)**2 * svec
            else
                dvec = eta * p
                svec = eta * pt
            end if
            x = x + dvec
            r = r - svec
            err = vec_norm2(r) / bnrm2
            if (iter < maxit - 1) errors(iter + 2) = err
            if (err < tol) then
                flag = 0
                exit
            end if
            iter = iter + 1
        end do
        call finish_result(x, errors, iter, maxit, flag, tol, result)
    end subroutine qmr_core

    subroutine lsolve_cg_dense(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Conjugate-gradient result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.

        call cg_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner, adjsym)
    end subroutine lsolve_cg_dense

    subroutine lsolve_cg_sparse(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Conjugate-gradient result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call cg_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m, adjsym)
        else
            call cg_core(sparse_operator(a), b, result, xinit, reltol, maxiter, adjsym=adjsym)
        end if
    end subroutine lsolve_cg_sparse

    subroutine lsolve_bicg_dense(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! BiCG result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.

        call bicg_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner)
    end subroutine lsolve_bicg_dense

    subroutine lsolve_bicg_sparse(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! BiCG result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call bicg_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m)
        else
            call bicg_core(sparse_operator(a), b, result, xinit, reltol, maxiter)
        end if
    end subroutine lsolve_bicg_sparse

    subroutine lsolve_bicgstab_dense(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! BiCGSTAB result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.

        call bicgstab_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner)
    end subroutine lsolve_bicgstab_dense

    subroutine lsolve_bicgstab_sparse(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! BiCGSTAB result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call bicgstab_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m)
        else
            call bicgstab_core(sparse_operator(a), b, result, xinit, reltol, maxiter)
        end if
    end subroutine lsolve_bicgstab_sparse

    subroutine lsolve_cgs_dense(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! CGS result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.

        call cgs_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner, adjsym)
    end subroutine lsolve_cgs_dense

    subroutine lsolve_cgs_sparse(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! CGS result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call cgs_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m, adjsym)
        else
            call cgs_core(sparse_operator(a), b, result, xinit, reltol, maxiter, adjsym=adjsym)
        end if
    end subroutine lsolve_cgs_sparse

    subroutine lsolve_cheby_dense(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Chebyshev result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.
        real(dp), intent(in), optional :: eigmin !! Optional smallest spectral bound for M^{-1}A.
        real(dp), intent(in), optional :: eigmax !! Optional largest spectral bound for M^{-1}A.

        call cheby_core(dense_operator(a), b, result, xinit, reltol, maxiter, &
            preconditioner, adjsym, eigmin, eigmax)
    end subroutine lsolve_cheby_dense

    subroutine lsolve_cheby_sparse(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Chebyshev result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Convert nonsymmetric square systems to normal equations when true.
        real(dp), intent(in), optional :: eigmin !! Optional smallest spectral bound for M^{-1}A.
        real(dp), intent(in), optional :: eigmax !! Optional largest spectral bound for M^{-1}A.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call cheby_core(sparse_operator(a), b, result, xinit, reltol, maxiter, &
                m, adjsym, eigmin, eigmax)
        else
            call cheby_core(sparse_operator(a), b, result, xinit, reltol, maxiter, &
                adjsym=adjsym, eigmin=eigmin, eigmax=eigmax)
        end if
    end subroutine lsolve_cheby_sparse

    subroutine lsolve_gmres_dense(a, b, result, xinit, reltol, maxiter, preconditioner, restart)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Restarted GMRES result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum restart cycles; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.
        integer, intent(in), optional :: restart !! Restart dimension; default n-1.

        call gmres_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner, restart)
    end subroutine lsolve_gmres_dense

    subroutine lsolve_gmres_sparse(a, b, result, xinit, reltol, maxiter, preconditioner, restart)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! Restarted GMRES result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum restart cycles; default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        integer, intent(in), optional :: restart !! Restart dimension; default n-1.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call gmres_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m, restart)
        else
            call gmres_core(sparse_operator(a), b, result, xinit, reltol, maxiter, restart=restart)
        end if
    end subroutine lsolve_gmres_sparse

    subroutine lsolve_qmr_dense(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! QMR result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional preconditioner M; default identity.

        call qmr_core(dense_operator(a), b, result, xinit, reltol, maxiter, preconditioner)
    end subroutine lsolve_qmr_dense

    subroutine lsolve_qmr_sparse(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector.
        type(solver_result), intent(out) :: result !! QMR result.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector.
        real(dp), intent(in), optional :: reltol !! Relative tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iterations; default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional sparse preconditioner M; default identity.
        real(dp), allocatable :: m(:, :)

        if (present(preconditioner)) then
            m = sparse_to_dense(preconditioner)
            call qmr_core(sparse_operator(a), b, result, xinit, reltol, maxiter, m)
        else
            call qmr_core(sparse_operator(a), b, result, xinit, reltol, maxiter)
        end if
    end subroutine lsolve_qmr_sparse

end module rlinsolve_krylov
