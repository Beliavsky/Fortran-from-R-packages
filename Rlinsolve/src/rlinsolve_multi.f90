module rlinsolve_multi
    use rlinsolve_kinds, only : dp
    use rlinsolve_types, only : solver_result, solver_matrix_result, sparse_matrix
    use rlinsolve_stationary, only : lsolve_jacobi_dense, lsolve_jacobi_sparse, &
        lsolve_gs_dense, lsolve_gs_sparse, lsolve_sor_dense, lsolve_sor_sparse, &
        lsolve_ssor_dense, lsolve_ssor_sparse
    use rlinsolve_krylov, only : lsolve_cg_dense, lsolve_cg_sparse, &
        lsolve_bicg_dense, lsolve_bicg_sparse, lsolve_bicgstab_dense, &
        lsolve_bicgstab_sparse, lsolve_cgs_dense, lsolve_cgs_sparse, &
        lsolve_cheby_dense, lsolve_cheby_sparse, lsolve_gmres_dense, &
        lsolve_gmres_sparse, lsolve_qmr_dense, lsolve_qmr_sparse
    implicit none
    private

    public :: lsolve_jacobi_dense_many, lsolve_jacobi_sparse_many
    public :: lsolve_gs_dense_many, lsolve_gs_sparse_many
    public :: lsolve_sor_dense_many, lsolve_sor_sparse_many
    public :: lsolve_ssor_dense_many, lsolve_ssor_sparse_many
    public :: lsolve_cg_dense_many, lsolve_cg_sparse_many
    public :: lsolve_bicg_dense_many, lsolve_bicg_sparse_many
    public :: lsolve_bicgstab_dense_many, lsolve_bicgstab_sparse_many
    public :: lsolve_cgs_dense_many, lsolve_cgs_sparse_many
    public :: lsolve_cheby_dense_many, lsolve_cheby_sparse_many
    public :: lsolve_gmres_dense_many, lsolve_gmres_sparse_many
    public :: lsolve_qmr_dense_many, lsolve_qmr_sparse_many

contains

    subroutine initialize_many(n, nrhs, result)
        integer, intent(in) :: n !! Number of unknowns in each solved system.
        integer, intent(in) :: nrhs !! Number of right-hand-side columns to solve.
        type(solver_matrix_result), intent(out) :: result !! Allocated multi-column result container.

        allocate(result%x(n, nrhs), result%iter(nrhs), result%flag(nrhs))
        allocate(result%converged(nrhs), result%column(nrhs))
        result%x = 0.0_dp
        result%iter = 0
        result%flag = 0
        result%converged = .false.
    end subroutine initialize_many

    subroutine store_one(j, one, result)
        integer, intent(in) :: j !! One-based right-hand-side column index being stored.
        type(solver_result), intent(in) :: one !! Single-column solver result copied into the aggregate output.
        type(solver_matrix_result), intent(inout) :: result !! Aggregate result receiving the column solution and diagnostics.

        result%x(:, j) = one%x
        result%iter(j) = one%iter
        result%flag(j) = one%flag
        result%converged(j) = one%converged
        result%column(j) = one
    end subroutine store_one

    subroutine lsolve_jacobi_dense_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! Weighted-Jacobi relaxation in (0,1]; default 2/3.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_jacobi_dense(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_jacobi_dense_many

    subroutine lsolve_jacobi_sparse_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! Weighted-Jacobi relaxation in (0,1]; default 2/3.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_jacobi_sparse(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_jacobi_sparse_many

    subroutine lsolve_gs_dense_many(a, b, result, xinit, reltol, maxiter, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_gs_dense(a, b(:, j), one, xinit, reltol, maxiter, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_gs_dense_many

    subroutine lsolve_gs_sparse_many(a, b, result, xinit, reltol, maxiter, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_gs_sparse(a, b(:, j), one, xinit, reltol, maxiter, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_gs_sparse_many

    subroutine lsolve_sor_dense_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! SOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_sor_dense(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_sor_dense_many

    subroutine lsolve_sor_sparse_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! SOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_sor_sparse(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_sor_sparse_many

    subroutine lsolve_ssor_dense_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! SSOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_ssor_dense(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_ssor_dense_many

    subroutine lsolve_ssor_sparse_many(a, b, result, xinit, reltol, maxiter, weight, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: weight !! SSOR relaxation parameter in (0,2); default 1.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_ssor_sparse(a, b(:, j), one, xinit, reltol, maxiter, weight, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_ssor_sparse_many

    subroutine lsolve_cg_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cg_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cg_dense_many

    subroutine lsolve_cg_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cg_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cg_sparse_many

    subroutine lsolve_bicg_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_bicg_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_bicg_dense_many

    subroutine lsolve_bicg_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_bicg_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_bicg_sparse_many

    subroutine lsolve_bicgstab_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_bicgstab_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_bicgstab_dense_many

    subroutine lsolve_bicgstab_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_bicgstab_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_bicgstab_sparse_many

    subroutine lsolve_cgs_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cgs_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cgs_dense_many

    subroutine lsolve_cgs_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cgs_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cgs_sparse_many

    subroutine lsolve_cheby_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        real(dp), intent(in), optional :: eigmin !! Optional smallest eigenvalue bound for M^{-1}A.
        real(dp), intent(in), optional :: eigmax !! Optional largest eigenvalue bound for M^{-1}A.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cheby_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cheby_dense_many

    subroutine lsolve_cheby_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 10000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        logical, intent(in), optional :: adjsym !! Whether nonsymmetric square A is converted to normal equations; default true.
        real(dp), intent(in), optional :: eigmin !! Optional smallest eigenvalue bound for M^{-1}A.
        real(dp), intent(in), optional :: eigmax !! Optional largest eigenvalue bound for M^{-1}A.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_cheby_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, adjsym, eigmin, eigmax)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_cheby_sparse_many

    subroutine lsolve_gmres_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner, restart)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        integer, intent(in), optional :: restart !! Restart dimension; default n-1.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_gmres_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, restart)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_gmres_dense_many

    subroutine lsolve_gmres_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner, restart)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        integer, intent(in), optional :: restart !! Restart dimension; default n-1.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_gmres_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner, restart)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_gmres_sparse_many

    subroutine lsolve_qmr_dense_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        real(dp), intent(in), optional :: preconditioner(:, :) !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(size(a, 2), size(b, 2), result)
        if (size(b, 1) /= size(a, 1)) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_qmr_dense(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_qmr_dense_many

    subroutine lsolve_qmr_sparse_many(a, b, result, xinit, reltol, maxiter, preconditioner)
        type(sparse_matrix), intent(in) :: a !! Sparse COO coefficient matrix A shared by all right-hand-side columns.
        real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides; each column is solved independently.
        type(solver_matrix_result), intent(out) :: result !! Multi-column solutions and per-column diagnostics.
        real(dp), intent(in), optional :: xinit(:) !! Optional common starting vector used for every right-hand-side column.
        real(dp), intent(in), optional :: reltol !! Relative residual tolerance passed to each single-column solve.
        integer, intent(in), optional :: maxiter !! Maximum iteration count passed to each column; upstream default 1000.
        type(sparse_matrix), intent(in), optional :: preconditioner !! Optional square preconditioner M; default identity.
        type(solver_result) :: one
        integer :: j

        call initialize_many(a%ncol, size(b, 2), result)
        if (size(b, 1) /= a%nrow) then
            result%flag = -1
            return
        end if
        do j = 1, size(b, 2)
            call lsolve_qmr_sparse(a, b(:, j), one, xinit, reltol, maxiter, preconditioner)
            call store_one(j, one, result)
        end do
    end subroutine lsolve_qmr_sparse_many

end module rlinsolve_multi
