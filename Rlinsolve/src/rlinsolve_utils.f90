module rlinsolve_utils
    use ieee_arithmetic, only : ieee_is_finite
    use rlinsolve_kinds, only : dp
    use rlinsolve_types, only : sparse_matrix
    implicit none
    private

    type, public :: matrix_operator
        integer :: nrow = 0
        integer :: ncol = 0
        logical :: is_sparse = .false.
        real(dp), allocatable :: dense(:, :)
        type(sparse_matrix) :: sparse
    end type matrix_operator

    type, public :: lu_factorization
        real(dp), allocatable :: lu(:, :)
        integer, allocatable :: piv(:)
        logical :: ok = .false.
    end type lu_factorization

    type, public :: preconditioner_type
        real(dp), allocatable :: mat(:, :)
        type(lu_factorization) :: fac
        logical :: identity = .true.
    end type preconditioner_type

    public :: dense_operator, sparse_operator, sparse_from_dense, sparse_to_dense
    public :: op_matvec, op_matvec_t, op_to_dense, vec_norm2, residual_norm
    public :: prepare_system, ensure_nonzero_diagonal, setup_preconditioner
    public :: apply_preconditioner, apply_preconditioner_t, lu_factor, lu_solve
    public :: lu_solve_t, triangular_lower_solve, triangular_upper_solve
    public :: estimate_preconditioned_bounds, is_diagonally_dominant

contains

    pure function dense_operator(a) result(op)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix to wrap as a linear operator.
        type(matrix_operator) :: op

        op%nrow = size(a, 1)
        op%ncol = size(a, 2)
        op%is_sparse = .false.
        allocate(op%dense(op%nrow, op%ncol))
        op%dense = a
    end function dense_operator

    pure function sparse_operator(a) result(op)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix to wrap as a linear operator.
        type(matrix_operator) :: op

        op%nrow = a%nrow
        op%ncol = a%ncol
        op%is_sparse = .true.
        op%sparse%nrow = a%nrow
        op%sparse%ncol = a%ncol
        if (allocated(a%row)) then
            allocate(op%sparse%row(size(a%row)))
            allocate(op%sparse%col(size(a%col)))
            allocate(op%sparse%val(size(a%val)))
            op%sparse%row = a%row
            op%sparse%col = a%col
            op%sparse%val = a%val
        end if
    end function sparse_operator

    pure function sparse_from_dense(a, threshold) result(sp)
        real(dp), intent(in) :: a(:, :) !! Dense matrix to convert to COO storage.
        real(dp), intent(in), optional :: threshold !! Entries with magnitude at or below this threshold are omitted; default zero.
        type(sparse_matrix) :: sp
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        tol = 0.0_dp
        if (present(threshold)) tol = max(0.0_dp, threshold)
        n = count(abs(a) > tol)
        sp%nrow = size(a, 1)
        sp%ncol = size(a, 2)
        allocate(sp%row(n), sp%col(n), sp%val(n))
        k = 0
        do j = 1, sp%ncol
            do i = 1, sp%nrow
                if (abs(a(i, j)) > tol) then
                    k = k + 1
                    sp%row(k) = i
                    sp%col(k) = j
                    sp%val(k) = a(i, j)
                end if
            end do
        end do
    end function sparse_from_dense

    pure function sparse_to_dense(sp) result(a)
        type(sparse_matrix), intent(in) :: sp !! Sparse COO matrix to expand into dense storage, summing duplicate entries.
        real(dp), allocatable :: a(:, :)
        integer :: k

        allocate(a(sp%nrow, sp%ncol))
        a = 0.0_dp
        if (.not. allocated(sp%val)) return
        do k = 1, size(sp%val)
            if (sp%row(k) >= 1 .and. sp%row(k) <= sp%nrow .and. &
                sp%col(k) >= 1 .and. sp%col(k) <= sp%ncol) then
                a(sp%row(k), sp%col(k)) = a(sp%row(k), sp%col(k)) + sp%val(k)
            end if
        end do
    end function sparse_to_dense

    pure function op_to_dense(op) result(a)
        type(matrix_operator), intent(in) :: op !! Linear operator whose full dense matrix representation is requested.
        real(dp), allocatable :: a(:, :)

        if (op%is_sparse) then
            a = sparse_to_dense(op%sparse)
        else
            allocate(a(op%nrow, op%ncol))
            a = op%dense
        end if
    end function op_to_dense

    pure subroutine op_matvec(op, x, y)
        type(matrix_operator), intent(in) :: op !! Matrix operator applied without transposition.
        real(dp), intent(in) :: x(:) !! Input vector of length op%ncol.
        real(dp), intent(out) :: y(:) !! Output vector of length op%nrow receiving op*x.
        integer :: k

        y = 0.0_dp
        if (op%is_sparse) then
            if (.not. allocated(op%sparse%val)) return
            do k = 1, size(op%sparse%val)
                y(op%sparse%row(k)) = y(op%sparse%row(k)) + &
                    op%sparse%val(k) * x(op%sparse%col(k))
            end do
        else
            y = matmul(op%dense, x)
        end if
    end subroutine op_matvec

    pure subroutine op_matvec_t(op, x, y)
        type(matrix_operator), intent(in) :: op !! Matrix operator applied as its transpose.
        real(dp), intent(in) :: x(:) !! Input vector of length op%nrow.
        real(dp), intent(out) :: y(:) !! Output vector of length op%ncol receiving transpose(op)*x.
        integer :: k

        y = 0.0_dp
        if (op%is_sparse) then
            if (.not. allocated(op%sparse%val)) return
            do k = 1, size(op%sparse%val)
                y(op%sparse%col(k)) = y(op%sparse%col(k)) + &
                    op%sparse%val(k) * x(op%sparse%row(k))
            end do
        else
            y = matmul(transpose(op%dense), x)
        end if
    end subroutine op_matvec_t

    pure real(dp) function vec_norm2(x) result(value)
        real(dp), intent(in) :: x(:) !! Real vector whose Euclidean norm is requested.

        value = sqrt(max(0.0_dp, dot_product(x, x)))
    end function vec_norm2

    pure real(dp) function residual_norm(op, x, b) result(value)
        type(matrix_operator), intent(in) :: op !! Coefficient operator in the residual b-A*x.
        real(dp), intent(in) :: x(:) !! Candidate solution vector of length op%ncol.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length op%nrow.
        real(dp), allocatable :: ax(:)

        allocate(ax(op%nrow))
        call op_matvec(op, x, ax)
        value = vec_norm2(b - ax)
    end function residual_norm

    pure logical function all_finite_matrix(a) result(ok)
        real(dp), intent(in) :: a(:, :) !! Matrix checked for finite IEEE values.
        integer :: i
        integer :: j

        ok = .true.
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (.not. ieee_is_finite(a(i, j))) then
                    ok = .false.
                    return
                end if
            end do
        end do
    end function all_finite_matrix

    pure logical function all_finite_vector(x) result(ok)
        real(dp), intent(in) :: x(:) !! Vector checked for finite IEEE values.
        integer :: i

        ok = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                ok = .false.
                return
            end if
        end do
    end function all_finite_vector

    subroutine prepare_system(input_op, input_b, require_symmetric, adjsym, op, b, status)
        type(matrix_operator), intent(in) :: input_op !! Original coefficient matrix, possibly rectangular.
        real(dp), intent(in) :: input_b(:) !! Original right-hand-side vector with length input_op%nrow.
        logical, intent(in) :: require_symmetric !! Whether a square nonsymmetric system should be converted to normal equations.
        logical, intent(in) :: adjsym !! Whether normal-equation conversion is permitted for a nonsymmetric square matrix.
        type(matrix_operator), intent(out) :: op !! Prepared square coefficient operator used by an iterative method.
        real(dp), allocatable, intent(out) :: b(:) !! Prepared right-hand-side vector matching op%nrow.
        integer, intent(out) :: status !! Zero on success; negative when preparation fails.
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: ata(:, :)
        real(dp), allocatable :: atb(:)
        real(dp) :: asym

        status = 0
        if (size(input_b) /= input_op%nrow) then
            status = -1
            return
        end if
        if (input_op%nrow < input_op%ncol) then
            status = -2
            return
        end if
        a = op_to_dense(input_op)
        if (.not. all_finite_matrix(a) .or. .not. all_finite_vector(input_b)) then
            status = -3
            return
        end if

        if (input_op%nrow > input_op%ncol) then
            allocate(ata(input_op%ncol, input_op%ncol))
            allocate(atb(input_op%ncol))
            ata = matmul(transpose(a), a)
            atb = matmul(transpose(a), input_b)
            op = dense_operator(ata)
            call move_alloc(atb, b)
            return
        end if

        if (require_symmetric) then
            asym = sqrt(sum((a - transpose(a))**2))
            if (asym > 1.0e-10_dp) then
                if (.not. adjsym) then
                    status = -4
                    return
                end if
                allocate(ata(input_op%ncol, input_op%ncol))
                allocate(atb(input_op%ncol))
                ata = matmul(transpose(a), a)
                atb = matmul(transpose(a), input_b)
                op = dense_operator(ata)
                call move_alloc(atb, b)
                return
            end if
        end if

        if (input_op%is_sparse) then
            op = sparse_operator(input_op%sparse)
        else
            op = dense_operator(input_op%dense)
        end if
        allocate(b(size(input_b)))
        b = input_b
    end subroutine prepare_system

    subroutine ensure_nonzero_diagonal(op)
        type(matrix_operator), intent(inout) :: op !! Square operator whose zero diagonal entries are perturbed.
        real(dp), allocatable :: a(:, :)
        integer :: i
        logical :: changed

        a = op_to_dense(op)
        changed = .false.
        do i = 1, min(size(a, 1), size(a, 2))
            if (abs(a(i, i)) <= tiny(1.0_dp)) then
                a(i, i) = 1.0e-5_dp
                changed = .true.
            end if
        end do
        if (changed) op = dense_operator(a)
    end subroutine ensure_nonzero_diagonal

    pure logical function is_diagonally_dominant(op, strict) result(ok)
        type(matrix_operator), intent(in) :: op !! Square matrix operator tested for row diagonal dominance.
        logical, intent(in), optional :: strict !! If true require strict dominance; otherwise weak dominance is accepted.
        real(dp), allocatable :: a(:, :)
        real(dp) :: offsum
        logical :: want_strict
        integer :: i

        want_strict = .false.
        if (present(strict)) want_strict = strict
        a = op_to_dense(op)
        ok = .true.
        do i = 1, min(size(a, 1), size(a, 2))
            offsum = sum(abs(a(i, :))) - abs(a(i, i))
            if (want_strict) then
                if (abs(a(i, i)) <= offsum) then
                    ok = .false.
                    return
                end if
            else
                if (abs(a(i, i)) < offsum) then
                    ok = .false.
                    return
                end if
            end if
        end do
    end function is_diagonally_dominant

    subroutine lu_factor(a, fac)
        real(dp), intent(in) :: a(:, :) !! Square dense matrix to factor with partial pivoting.
        type(lu_factorization), intent(out) :: fac !! LU factors and row-swap history for subsequent solves.
        integer :: n
        integer :: i
        integer :: k
        integer :: p
        real(dp) :: pivot_abs
        real(dp), allocatable :: tmp(:)

        n = size(a, 1)
        fac%ok = .false.
        if (size(a, 2) /= n) return
        allocate(fac%lu(n, n), fac%piv(n), tmp(n))
        fac%lu = a
        fac%piv = [(i, i = 1, n)]
        do k = 1, n - 1
            p = k
            pivot_abs = abs(fac%lu(k, k))
            do i = k + 1, n
                if (abs(fac%lu(i, k)) > pivot_abs) then
                    p = i
                    pivot_abs = abs(fac%lu(i, k))
                end if
            end do
            fac%piv(k) = p
            if (pivot_abs <= tiny(1.0_dp)) return
            if (p /= k) then
                tmp = fac%lu(k, :)
                fac%lu(k, :) = fac%lu(p, :)
                fac%lu(p, :) = tmp
            end if
            do i = k + 1, n
                fac%lu(i, k) = fac%lu(i, k) / fac%lu(k, k)
                fac%lu(i, k + 1:n) = fac%lu(i, k + 1:n) - &
                    fac%lu(i, k) * fac%lu(k, k + 1:n)
            end do
        end do
        fac%piv(n) = n
        if (n > 0 .and. abs(fac%lu(n, n)) <= tiny(1.0_dp)) return
        fac%ok = .true.
    end subroutine lu_factor

    subroutine lu_solve(fac, rhs, x, status)
        type(lu_factorization), intent(in) :: fac !! Previously computed partial-pivoted LU factorization.
        real(dp), intent(in) :: rhs(:) !! Right-hand side of A*x=rhs.
        real(dp), intent(out) :: x(:) !! Solution vector overwritten with the computed result.
        integer, intent(out) :: status !! Zero on success and negative if the factorization or dimensions are invalid.
        integer :: i
        integer :: k
        integer :: n
        real(dp) :: tmp

        status = 0
        n = size(rhs)
        if (.not. fac%ok .or. size(x) /= n .or. size(fac%lu, 1) /= n) then
            status = -1
            x = 0.0_dp
            return
        end if
        x = rhs
        do k = 1, n - 1
            if (fac%piv(k) /= k) then
                tmp = x(k)
                x(k) = x(fac%piv(k))
                x(fac%piv(k)) = tmp
            end if
        end do
        do i = 2, n
            x(i) = x(i) - dot_product(fac%lu(i, 1:i - 1), x(1:i - 1))
        end do
        do i = n, 1, -1
            if (i < n) x(i) = x(i) - dot_product(fac%lu(i, i + 1:n), x(i + 1:n))
            x(i) = x(i) / fac%lu(i, i)
        end do
    end subroutine lu_solve

    subroutine lu_solve_t(fac, rhs, x, status)
        type(lu_factorization), intent(in) :: fac !! Previously computed factorization of A used to solve transpose(A)*x=rhs.
        real(dp), intent(in) :: rhs(:) !! Right-hand side of the transposed system.
        real(dp), intent(out) :: x(:) !! Solution vector for transpose(A)*x=rhs.
        integer, intent(out) :: status !! Zero on success and negative if the factorization or dimensions are invalid.
        integer :: i
        integer :: k
        integer :: n
        real(dp) :: tmp

        status = 0
        n = size(rhs)
        if (.not. fac%ok .or. size(x) /= n .or. size(fac%lu, 1) /= n) then
            status = -1
            x = 0.0_dp
            return
        end if
        x = rhs
        do i = 1, n
            if (i > 1) x(i) = x(i) - dot_product(fac%lu(1:i - 1, i), x(1:i - 1))
            x(i) = x(i) / fac%lu(i, i)
        end do
        do i = n - 1, 1, -1
            x(i) = x(i) - dot_product(fac%lu(i + 1:n, i), x(i + 1:n))
        end do
        do k = n - 1, 1, -1
            if (fac%piv(k) /= k) then
                tmp = x(k)
                x(k) = x(fac%piv(k))
                x(fac%piv(k)) = tmp
            end if
        end do
    end subroutine lu_solve_t

    subroutine setup_preconditioner(n, matrix, pre, status)
        integer, intent(in) :: n !! Required square order of the preconditioner.
        real(dp), intent(in), optional :: matrix(:, :) !! Optional explicit preconditioner M; absence means identity.
        type(preconditioner_type), intent(out) :: pre !! Prepared preconditioner with reusable LU factors.
        integer, intent(out) :: status !! Zero on success and negative for shape or singularity errors.
        integer :: i

        status = 0
        if (.not. present(matrix)) then
            pre%identity = .true.
            allocate(pre%mat(n, n))
            pre%mat = 0.0_dp
            do i = 1, n
                pre%mat(i, i) = 1.0_dp
            end do
            return
        end if
        if (size(matrix, 1) /= n .or. size(matrix, 2) /= n) then
            status = -1
            return
        end if
        pre%identity = .false.
        allocate(pre%mat(n, n))
        pre%mat = matrix
        call lu_factor(pre%mat, pre%fac)
        if (.not. pre%fac%ok) status = -2
    end subroutine setup_preconditioner

    subroutine apply_preconditioner(pre, rhs, x, status)
        type(preconditioner_type), intent(in) :: pre !! Prepared M used to compute M^{-1}*rhs.
        real(dp), intent(in) :: rhs(:) !! Vector to which the inverse preconditioner is applied.
        real(dp), intent(out) :: x(:) !! Result M^{-1}*rhs.
        integer, intent(out) :: status !! Zero on success and negative if the preconditioner solve fails.

        if (pre%identity) then
            x = rhs
            status = 0
        else
            call lu_solve(pre%fac, rhs, x, status)
        end if
    end subroutine apply_preconditioner

    subroutine apply_preconditioner_t(pre, rhs, x, status)
        type(preconditioner_type), intent(in) :: pre !! Prepared M used to compute transpose(M)^{-1}*rhs.
        real(dp), intent(in) :: rhs(:) !! Vector to which the transposed inverse preconditioner is applied.
        real(dp), intent(out) :: x(:) !! Result transpose(M)^{-1}*rhs.
        integer, intent(out) :: status !! Zero on success and negative if the preconditioner solve fails.

        if (pre%identity) then
            x = rhs
            status = 0
        else
            call lu_solve_t(pre%fac, rhs, x, status)
        end if
    end subroutine apply_preconditioner_t

    subroutine triangular_lower_solve(a, rhs, x, status)
        real(dp), intent(in) :: a(:, :) !! Square lower-triangular matrix including its diagonal.
        real(dp), intent(in) :: rhs(:) !! Right-hand side of the triangular system.
        real(dp), intent(out) :: x(:) !! Solution of a*x=rhs.
        integer, intent(out) :: status !! Zero on success and negative if a diagonal pivot vanishes.
        integer :: i
        integer :: n

        n = size(rhs)
        status = 0
        x = rhs
        do i = 1, n
            if (abs(a(i, i)) <= tiny(1.0_dp)) then
                status = -1
                x = 0.0_dp
                return
            end if
            if (i > 1) x(i) = x(i) - dot_product(a(i, 1:i - 1), x(1:i - 1))
            x(i) = x(i) / a(i, i)
        end do
    end subroutine triangular_lower_solve

    subroutine triangular_upper_solve(a, rhs, x, status)
        real(dp), intent(in) :: a(:, :) !! Square upper-triangular matrix including its diagonal.
        real(dp), intent(in) :: rhs(:) !! Right-hand side of the triangular system.
        real(dp), intent(out) :: x(:) !! Solution of a*x=rhs.
        integer, intent(out) :: status !! Zero on success and negative if a diagonal pivot vanishes.
        integer :: i
        integer :: n

        n = size(rhs)
        status = 0
        x = rhs
        do i = n, 1, -1
            if (abs(a(i, i)) <= tiny(1.0_dp)) then
                status = -1
                x = 0.0_dp
                return
            end if
            if (i < n) x(i) = x(i) - dot_product(a(i, i + 1:n), x(i + 1:n))
            x(i) = x(i) / a(i, i)
        end do
    end subroutine triangular_upper_solve

    subroutine estimate_preconditioned_bounds(op, pre, eigmin, eigmax, status)
        type(matrix_operator), intent(in) :: op !! Square positive-definite coefficient operator A.
        type(preconditioner_type), intent(in) :: pre !! Positive-definite preconditioner M in M^{-1}A.
        real(dp), intent(out) :: eigmin !! Estimated smallest eigenvalue of M^{-1}A.
        real(dp), intent(out) :: eigmax !! Estimated largest eigenvalue of M^{-1}A.
        integer, intent(out) :: status !! Zero on success and negative if an internal solve fails.
        type(lu_factorization) :: afac
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: av(:)
        real(dp), allocatable :: mv(:)
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: w(:)
        real(dp) :: denom
        real(dp) :: lambda
        integer :: i
        integer :: istat
        integer :: n

        status = 0
        n = op%ncol
        allocate(v(n), w(n), av(n), mv(n))
        v = 1.0_dp / sqrt(real(max(1, n), dp))
        lambda = 0.0_dp
        do i = 1, 60
            call op_matvec(op, v, av)
            call apply_preconditioner(pre, av, w, istat)
            if (istat /= 0) then
                status = -1
                return
            end if
            denom = vec_norm2(w)
            if (denom <= tiny(1.0_dp)) then
                status = -2
                return
            end if
            v = w / denom
            lambda = dot_product(v, w)
        end do
        eigmax = abs(lambda)

        a = op_to_dense(op)
        call lu_factor(a, afac)
        if (.not. afac%ok) then
            status = -3
            return
        end if
        v = 1.0_dp / sqrt(real(max(1, n), dp))
        lambda = 0.0_dp
        do i = 1, 60
            mv = matmul(pre%mat, v)
            call lu_solve(afac, mv, w, istat)
            if (istat /= 0) then
                status = -4
                return
            end if
            denom = vec_norm2(w)
            if (denom <= tiny(1.0_dp)) then
                status = -5
                return
            end if
            v = w / denom
            lambda = dot_product(v, w)
        end do
        if (abs(lambda) <= tiny(1.0_dp)) then
            status = -6
            return
        end if
        eigmin = abs(1.0_dp / lambda)
        if (eigmin > eigmax) then
            denom = eigmin
            eigmin = eigmax
            eigmax = denom
        end if
    end subroutine estimate_preconditioned_bounds

end module rlinsolve_utils
