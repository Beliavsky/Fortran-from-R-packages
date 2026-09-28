module rlinsolve_stationary
    use rlinsolve_kinds, only : dp
    use rlinsolve_types, only : solver_result, sparse_matrix
    use rlinsolve_utils, only : matrix_operator, dense_operator, sparse_operator, op_matvec, &
        op_to_dense, prepare_system, ensure_nonzero_diagonal, vec_norm2, residual_norm, &
        triangular_lower_solve, triangular_upper_solve
    implicit none
    private

    public :: lsolve_jacobi_dense, lsolve_jacobi_sparse
    public :: lsolve_gs_dense, lsolve_gs_sparse
    public :: lsolve_sor_dense, lsolve_sor_sparse
    public :: lsolve_ssor_dense, lsolve_ssor_sparse

contains

    subroutine initialize_x(op, b, xinit, tol, x, status)
        type(matrix_operator), intent(in) :: op !! Prepared square coefficient operator used to size and test the initial iterate.
        real(dp), intent(in) :: b(:) !! Prepared right-hand-side vector.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial iterate; absence uses a deterministic zero vector.
        real(dp), intent(in) :: tol !! Residual threshold triggering deterministic start replacement.
        real(dp), allocatable, intent(out) :: x(:) !! Initial iterate selected for the solver.
        integer, intent(out) :: status !! Zero on success and negative if xinit has the wrong length.
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

    pure subroutine get_diagonal(op, d)
        type(matrix_operator), intent(in) :: op !! Square operator whose diagonal entries are extracted.
        real(dp), intent(out) :: d(:) !! Diagonal vector with length min(op%nrow,op%ncol).
        integer :: i
        integer :: k

        d = 0.0_dp
        if (op%is_sparse) then
            if (.not. allocated(op%sparse%val)) return
            do k = 1, size(op%sparse%val)
                if (op%sparse%row(k) == op%sparse%col(k)) then
                    i = op%sparse%row(k)
                    if (i <= size(d)) d(i) = d(i) + op%sparse%val(k)
                end if
            end do
        else
            do i = 1, size(d)
                d(i) = op%dense(i, i)
            end do
        end if
    end subroutine get_diagonal

    subroutine finalize_stationary(x, errors, iter, maxiter, tol, result)
        real(dp), intent(in) :: x(:) !! Final approximation vector.
        real(dp), intent(in) :: errors(:) !! Full allocated convergence-history workspace.
        integer, intent(in) :: iter !! Upstream-compatible zero-based iteration index on exit.
        integer, intent(in) :: maxiter !! Maximum allowed iteration count used to size history.
        real(dp), intent(in) :: tol !! Relative residual tolerance used to determine convergence.
        type(solver_result), intent(out) :: result !! Packaged solution, history, iteration count, flag, and convergence indicator.
        integer :: nerr

        allocate(result%x(size(x)))
        result%x = x
        result%iter = iter
        nerr = min(maxiter, max(1, iter + 1))
        allocate(result%errors(nerr))
        result%errors = errors(1:nerr)
        result%converged = result%errors(nerr) < tol
        if (result%converged) then
            result%flag = 0
        else
            result%flag = 1
        end if
    end subroutine finalize_stationary

    subroutine jacobi_core(input_op, input_b, result, xinit, reltol, maxiter, weight, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse input coefficient operator, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand-side vector corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Jacobi solution, error history, and convergence information.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        real(dp), intent(in), optional :: weight !! Weighted-Jacobi relaxation in (0,1]; default 2/3.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square systems use normal equations; default true.
        type(matrix_operator) :: op
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: ax(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: tol
        real(dp) :: w
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: status
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        w = 2.0_dp / 3.0_dp
        if (present(weight)) w = weight
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1 .or. w <= 0.0_dp .or. w > 1.0_dp) then
            result%flag = -10
            return
        end if
        call prepare_system(input_op, input_b, .true., dosym, op, b, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        call ensure_nonzero_diagonal(op)
        call initialize_x(op, b, xinit, tol, x, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        allocate(d(op%nrow), ax(op%nrow), xnew(op%ncol), errors(maxit))
        call get_diagonal(op, d)
        denom = residual_norm(op, x, b)
        if (denom <= tiny(1.0_dp)) denom = 1.0_dp
        errors = 0.0_dp
        iter = 0
        do while (iter < maxit)
            call op_matvec(op, x, ax)
            xnew = x
            do i = 1, op%ncol
                if (abs(d(i)) <= tiny(1.0_dp)) cycle
                xnew(i) = x(i) + w * (b(i) - ax(i)) / d(i)
            end do
            err = residual_norm(op, xnew, b) / denom
            errors(iter + 1) = err
            x = xnew
            if (err < tol) exit
            iter = iter + 1
        end do
        call finalize_stationary(x, errors, iter, maxit, tol, result)
    end subroutine jacobi_core

    subroutine gs_core(input_op, input_b, result, xinit, reltol, maxiter, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse input coefficient operator, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand-side vector corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! Forward Gauss-Seidel solution and convergence history.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square systems use normal equations; default true.
        type(matrix_operator) :: op
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: sigma
        real(dp) :: tol
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: status
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1) then
            result%flag = -10
            return
        end if
        call prepare_system(input_op, input_b, .true., dosym, op, b, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        call ensure_nonzero_diagonal(op)
        a = op_to_dense(op)
        call initialize_x(op, b, xinit, tol, x, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        allocate(xnew(size(x)), errors(maxit))
        denom = residual_norm(op, x, b)
        if (denom <= tiny(1.0_dp)) denom = 1.0_dp
        errors = 0.0_dp
        iter = 0
        do while (iter < maxit)
            xnew = x
            do i = 1, size(x)
                sigma = 0.0_dp
                if (i > 1) sigma = sigma + dot_product(a(i, 1:i - 1), xnew(1:i - 1))
                if (i < size(x)) sigma = sigma + dot_product(a(i, i + 1:size(x)), x(i + 1:size(x)))
                xnew(i) = (b(i) - sigma) / a(i, i)
            end do
            err = residual_norm(op, xnew, b) / denom
            errors(iter + 1) = err
            x = xnew
            if (err < tol) exit
            iter = iter + 1
        end do
        call finalize_stationary(x, errors, iter, maxit, tol, result)
    end subroutine gs_core

    subroutine sor_core(input_op, input_b, result, xinit, reltol, maxiter, weight, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse input coefficient operator, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand-side vector corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! SOR solution and convergence history.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        real(dp), intent(in), optional :: weight !! SOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square systems use normal equations; default true.
        type(matrix_operator) :: op
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: sigma
        real(dp) :: tol
        real(dp) :: w
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: status
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        w = 1.0_dp
        if (present(weight)) w = weight
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1 .or. w <= 0.0_dp .or. w >= 2.0_dp) then
            result%flag = -10
            return
        end if
        call prepare_system(input_op, input_b, .true., dosym, op, b, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        call ensure_nonzero_diagonal(op)
        a = op_to_dense(op)
        call initialize_x(op, b, xinit, tol, x, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        allocate(xnew(size(x)), errors(maxit))
        denom = residual_norm(op, x, b)
        if (denom <= tiny(1.0_dp)) denom = 1.0_dp
        errors = 0.0_dp
        iter = 0
        do while (iter < maxit)
            xnew = x
            do i = 1, size(x)
                sigma = 0.0_dp
                if (i > 1) sigma = sigma + dot_product(a(i, 1:i - 1), xnew(1:i - 1))
                if (i < size(x)) sigma = sigma + dot_product(a(i, i + 1:size(x)), x(i + 1:size(x)))
                xnew(i) = (1.0_dp - w) * x(i) + w * (b(i) - sigma) / a(i, i)
            end do
            err = residual_norm(op, xnew, b) / denom
            errors(iter + 1) = err
            x = xnew
            if (err < tol) exit
            iter = iter + 1
        end do
        call finalize_stationary(x, errors, iter, maxit, tol, result)
    end subroutine sor_core

    subroutine ssor_core(input_op, input_b, result, xinit, reltol, maxiter, weight, adjsym)
        type(matrix_operator), intent(in) :: input_op !! Dense or sparse input coefficient operator, possibly overdetermined.
        real(dp), intent(in) :: input_b(:) !! Right-hand-side vector corresponding to input_op rows.
        type(solver_result), intent(out) :: result !! SSOR solution and convergence history.
        real(dp), intent(in), optional :: xinit(:) !! Optional starting vector of length input_op%ncol.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum iteration count; default 1000.
        real(dp), intent(in), optional :: weight !! SSOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square systems use normal equations; default true.
        type(matrix_operator) :: op
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: c1(:)
        real(dp), allocatable :: c2(:)
        real(dp), allocatable :: dmat(:, :)
        real(dp), allocatable :: errors(:)
        real(dp), allocatable :: lower(:, :)
        real(dp), allocatable :: lstar(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: s1(:, :)
        real(dp), allocatable :: s2(:, :)
        real(dp), allocatable :: upper(:, :)
        real(dp), allocatable :: ustar(:, :)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: xhalf(:)
        real(dp), allocatable :: xnew(:)
        real(dp) :: denom
        real(dp) :: err
        real(dp) :: tol
        real(dp) :: w
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: n
        integer :: status
        logical :: dosym

        tol = 1.0e-5_dp
        if (present(reltol)) tol = reltol
        maxit = 1000
        if (present(maxiter)) maxit = maxiter
        w = 1.0_dp
        if (present(weight)) w = weight
        dosym = .true.
        if (present(adjsym)) dosym = adjsym
        if (tol < 0.0_dp .or. maxit < 1 .or. w <= 0.0_dp .or. w >= 2.0_dp) then
            result%flag = -10
            return
        end if
        call prepare_system(input_op, input_b, .true., dosym, op, b, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        call ensure_nonzero_diagonal(op)
        a = op_to_dense(op)
        n = size(a, 1)
        allocate(dmat(n, n), lower(n, n), upper(n, n), lstar(n, n), ustar(n, n))
        allocate(s1(n, n), s2(n, n), c1(n), c2(n), rhs(n), xhalf(n), xnew(n))
        dmat = 0.0_dp
        lower = 0.0_dp
        upper = 0.0_dp
        do i = 1, n
            dmat(i, i) = a(i, i)
            if (i > 1) lower(i, 1:i - 1) = a(i, 1:i - 1)
            if (i < n) upper(i, i + 1:n) = a(i, i + 1:n)
        end do
        lstar = dmat + w * lower
        ustar = dmat + w * upper
        s1 = -w * lower + (1.0_dp - w) * dmat
        s2 = -w * upper + (1.0_dp - w) * dmat
        call triangular_lower_solve(lstar, b, c1, status)
        if (status /= 0) then
            result%flag = -11
            return
        end if
        rhs = matmul(dmat, c1)
        call triangular_upper_solve(ustar, rhs, c2, status)
        if (status /= 0) then
            result%flag = -12
            return
        end if
        c2 = w * (2.0_dp - w) * c2
        call initialize_x(op, b, xinit, tol, x, status)
        if (status /= 0) then
            result%flag = status
            return
        end if
        allocate(errors(maxit))
        errors = 0.0_dp
        denom = residual_norm(op, x, b)
        if (denom <= tiny(1.0_dp)) denom = 1.0_dp
        iter = 0
        do while (iter < maxit)
            rhs = matmul(s2, x)
            call triangular_lower_solve(lstar, rhs, xhalf, status)
            if (status /= 0) exit
            rhs = matmul(s1, xhalf)
            call triangular_upper_solve(ustar, rhs, xnew, status)
            if (status /= 0) exit
            xnew = xnew + c2
            err = residual_norm(op, xnew, b) / denom
            errors(iter + 1) = err
            x = xnew
            if (err < tol) exit
            iter = iter + 1
        end do
        if (status /= 0) then
            result%flag = -13
            return
        end if
        call finalize_stationary(x, errors, iter, maxit, tol, result)
    end subroutine ssor_core

    subroutine lsolve_jacobi_dense(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! Weighted-Jacobi solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! Jacobi relaxation weight in (0,1]; default 2/3.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call jacobi_core(dense_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_jacobi_dense

    subroutine lsolve_jacobi_sparse(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! Weighted-Jacobi solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! Jacobi relaxation weight in (0,1]; default 2/3.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call jacobi_core(sparse_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_jacobi_sparse

    subroutine lsolve_gs_dense(a, b, result, xinit, reltol, maxiter, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! Forward Gauss-Seidel solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call gs_core(dense_operator(a), b, result, xinit, reltol, maxiter, adjsym)
    end subroutine lsolve_gs_dense

    subroutine lsolve_gs_sparse(a, b, result, xinit, reltol, maxiter, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! Forward Gauss-Seidel solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call gs_core(sparse_operator(a), b, result, xinit, reltol, maxiter, adjsym)
    end subroutine lsolve_gs_sparse

    subroutine lsolve_sor_dense(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! SOR solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! SOR relaxation parameter w in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call sor_core(dense_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_sor_dense

    subroutine lsolve_sor_sparse(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! SOR solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! SOR relaxation parameter w in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call sor_core(sparse_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_sor_sparse

    subroutine lsolve_ssor_dense(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! SSOR solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! SSOR relaxation parameter w in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call ssor_core(dense_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_ssor_dense

    subroutine lsolve_ssor_sparse(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A with m rows and n columns.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length m.
        type(solver_result), intent(out) :: result !! SSOR solver output.
        real(dp), intent(in), optional :: xinit(:) !! Optional initial vector of length n.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance; default 1e-5.
        integer, intent(in), optional :: maxiter !! Maximum number of iterations; default 1000.
        real(dp), intent(in), optional :: weight !! SSOR relaxation parameter w in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.

        call ssor_core(sparse_operator(a), b, result, xinit, reltol, maxiter, weight, adjsym)
    end subroutine lsolve_ssor_sparse

end module rlinsolve_stationary
