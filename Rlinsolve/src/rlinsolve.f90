module rlinsolve
    use rlinsolve_kinds, only : dp
    use rlinsolve_types, only : sparse_matrix, solver_result, solver_matrix_result
    use rlinsolve_utils, only : sparse_from_dense, sparse_to_dense
    use rlinsolve_stationary, only : lsolve_jacobi_dense, lsolve_jacobi_sparse, &
        lsolve_gs_dense, lsolve_gs_sparse, lsolve_sor_dense, lsolve_sor_sparse, &
        lsolve_ssor_dense, lsolve_ssor_sparse
    use rlinsolve_krylov, only : lsolve_cg_dense, lsolve_cg_sparse, &
        lsolve_bicg_dense, lsolve_bicg_sparse, lsolve_bicgstab_dense, &
        lsolve_bicgstab_sparse, lsolve_cgs_dense, lsolve_cgs_sparse, &
        lsolve_cheby_dense, lsolve_cheby_sparse, lsolve_gmres_dense, &
        lsolve_gmres_sparse, lsolve_qmr_dense, lsolve_qmr_sparse
    use rlinsolve_multi, only : lsolve_jacobi_dense_many, lsolve_jacobi_sparse_many, &
        lsolve_gs_dense_many, lsolve_gs_sparse_many, lsolve_sor_dense_many, &
        lsolve_sor_sparse_many, lsolve_ssor_dense_many, lsolve_ssor_sparse_many, &
        lsolve_cg_dense_many, lsolve_cg_sparse_many, lsolve_bicg_dense_many, &
        lsolve_bicg_sparse_many, lsolve_bicgstab_dense_many, lsolve_bicgstab_sparse_many, &
        lsolve_cgs_dense_many, lsolve_cgs_sparse_many, lsolve_cheby_dense_many, &
        lsolve_cheby_sparse_many, lsolve_gmres_dense_many, lsolve_gmres_sparse_many, &
        lsolve_qmr_dense_many, lsolve_qmr_sparse_many
    implicit none
    private

    public :: dp, sparse_matrix, solver_result, solver_matrix_result
    public :: sparse_from_dense, sparse_to_dense
    public :: aux_fisch, aux_fisch_dense, aux_fisch_sparse
    public :: lsolve_jacobi, lsolve_gs, lsolve_sor, lsolve_ssor
    public :: lsolve_cg, lsolve_bicg, lsolve_bicgstab, lsolve_cgs
    public :: lsolve_cheby, lsolve_gmres, lsolve_qmr

    interface aux_fisch
        module procedure aux_fisch_dense
        module procedure aux_fisch_sparse
    end interface aux_fisch

    interface lsolve_jacobi
        procedure :: lsolve_jacobi_dense
        procedure :: lsolve_jacobi_sparse
        procedure :: lsolve_jacobi_dense_many
        procedure :: lsolve_jacobi_sparse_many
    end interface lsolve_jacobi

    interface lsolve_gs
        procedure :: lsolve_gs_dense
        procedure :: lsolve_gs_sparse
        procedure :: lsolve_gs_dense_many
        procedure :: lsolve_gs_sparse_many
    end interface lsolve_gs

    interface lsolve_sor
        procedure :: lsolve_sor_dense
        procedure :: lsolve_sor_sparse
        procedure :: lsolve_sor_dense_many
        procedure :: lsolve_sor_sparse_many
    end interface lsolve_sor

    interface lsolve_ssor
        procedure :: lsolve_ssor_dense
        procedure :: lsolve_ssor_sparse
        procedure :: lsolve_ssor_dense_many
        procedure :: lsolve_ssor_sparse_many
    end interface lsolve_ssor

    interface lsolve_cg
        procedure :: lsolve_cg_dense
        procedure :: lsolve_cg_sparse
        procedure :: lsolve_cg_dense_many
        procedure :: lsolve_cg_sparse_many
    end interface lsolve_cg

    interface lsolve_bicg
        procedure :: lsolve_bicg_dense
        procedure :: lsolve_bicg_sparse
        procedure :: lsolve_bicg_dense_many
        procedure :: lsolve_bicg_sparse_many
    end interface lsolve_bicg

    interface lsolve_bicgstab
        procedure :: lsolve_bicgstab_dense
        procedure :: lsolve_bicgstab_sparse
        procedure :: lsolve_bicgstab_dense_many
        procedure :: lsolve_bicgstab_sparse_many
    end interface lsolve_bicgstab

    interface lsolve_cgs
        procedure :: lsolve_cgs_dense
        procedure :: lsolve_cgs_sparse
        procedure :: lsolve_cgs_dense_many
        procedure :: lsolve_cgs_sparse_many
    end interface lsolve_cgs

    interface lsolve_cheby
        procedure :: lsolve_cheby_dense
        procedure :: lsolve_cheby_sparse
        procedure :: lsolve_cheby_dense_many
        procedure :: lsolve_cheby_sparse_many
    end interface lsolve_cheby

    interface lsolve_gmres
        procedure :: lsolve_gmres_dense
        procedure :: lsolve_gmres_sparse
        procedure :: lsolve_gmres_dense_many
        procedure :: lsolve_gmres_sparse_many
    end interface lsolve_gmres

    interface lsolve_qmr
        procedure :: lsolve_qmr_dense
        procedure :: lsolve_qmr_sparse
        procedure :: lsolve_qmr_dense_many
        procedure :: lsolve_qmr_sparse_many
    end interface lsolve_qmr

contains

    pure subroutine aux_fisch_dense(n, a)
        integer, intent(in) :: n !! Number of grid points in each spatial direction; valid values are integers greater than 2.
        real(dp), allocatable, intent(out) :: a(:, :) !! Dense n^2-by-n^2 five-point discrete Poisson matrix.
        integer :: col
        integer :: idx
        integer :: row
        integer :: total

        if (n <= 2) then
            allocate(a(0, 0))
            return
        end if
        total = n * n
        allocate(a(total, total))
        a = 0.0_dp
        do row = 1, n
            do col = 1, n
                idx = (row - 1) * n + col
                a(idx, idx) = 4.0_dp
                if (col > 1) a(idx, idx - 1) = -1.0_dp
                if (col < n) a(idx, idx + 1) = -1.0_dp
                if (row > 1) a(idx, idx - n) = -1.0_dp
                if (row < n) a(idx, idx + n) = -1.0_dp
            end do
        end do
    end subroutine aux_fisch_dense

    pure subroutine aux_fisch_sparse(n, a)
        integer, intent(in) :: n !! Number of grid points in each spatial direction; valid values are integers greater than 2.
        type(sparse_matrix), intent(out) :: a !! Sparse COO n^2-by-n^2 five-point discrete Poisson matrix.
        integer :: col
        integer :: idx
        integer :: k
        integer :: nnz
        integer :: row
        integer :: total

        a%nrow = 0
        a%ncol = 0
        if (n <= 2) then
            allocate(a%row(0), a%col(0), a%val(0))
            return
        end if
        total = n * n
        nnz = 5 * n * n - 4 * n
        a%nrow = total
        a%ncol = total
        allocate(a%row(nnz), a%col(nnz), a%val(nnz))
        k = 0
        do row = 1, n
            do col = 1, n
                idx = (row - 1) * n + col
                k = k + 1
                a%row(k) = idx
                a%col(k) = idx
                a%val(k) = 4.0_dp
                if (col > 1) then
                    k = k + 1
                    a%row(k) = idx
                    a%col(k) = idx - 1
                    a%val(k) = -1.0_dp
                end if
                if (col < n) then
                    k = k + 1
                    a%row(k) = idx
                    a%col(k) = idx + 1
                    a%val(k) = -1.0_dp
                end if
                if (row > 1) then
                    k = k + 1
                    a%row(k) = idx
                    a%col(k) = idx - n
                    a%val(k) = -1.0_dp
                end if
                if (row < n) then
                    k = k + 1
                    a%row(k) = idx
                    a%col(k) = idx + n
                    a%val(k) = -1.0_dp
                end if
            end do
        end do
    end subroutine aux_fisch_sparse

end module rlinsolve
