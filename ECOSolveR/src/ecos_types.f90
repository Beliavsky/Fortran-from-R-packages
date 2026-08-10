! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_types
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    integer, parameter, public :: ECOS_OPTIMAL = 0
    integer, parameter, public :: ECOS_PINF = 1
    integer, parameter, public :: ECOS_DINF = 2
    integer, parameter, public :: ECOS_INACC_OFFSET = 10
    integer, parameter, public :: ECOS_MAXIT = -1
    integer, parameter, public :: ECOS_NUMERICS = -2
    integer, parameter, public :: ECOS_OUTCONE = -3
    integer, parameter, public :: ECOS_SIGINT = -4
    integer, parameter, public :: ECOS_FATAL = -7

    integer, parameter, public :: MI_OPTIMAL_SOLN = 0
    integer, parameter, public :: MI_MAXITER_FEASIBLE_SOLN = 1
    integer, parameter, public :: MI_MAXITER_NO_SOLN = 2
    integer, parameter, public :: MI_INFEASIBLE = 3

    type, public :: ecos_dims
        integer :: l = 0
        integer, allocatable :: q(:)
        integer :: e = 0
    contains
        procedure :: cone_rows => dims_cone_rows
        procedure :: scalar_inequalities => dims_scalar_inequalities
    end type ecos_dims

    type, public :: ecos_settings
        real(dp) :: feastol = 1.0e-8_dp
        real(dp) :: reltol = 1.0e-8_dp
        real(dp) :: abstol = 1.0e-8_dp
        real(dp) :: feastol_inacc = 1.0e-4_dp
        real(dp) :: abstol_inacc = 5.0e-5_dp
        real(dp) :: reltol_inacc = 5.0e-5_dp
        integer :: maxit = 100
        logical :: verbose = .false.
        integer :: mi_max_iters = 1000
        real(dp) :: mi_abs_eps = 1.0e-6_dp
        real(dp) :: mi_rel_eps = 1.0e-6_dp
        real(dp) :: mi_int_tol = 1.0e-4_dp
        real(dp) :: regularization = 1.0e-9_dp
        logical :: sparse_kkt = .true.
        logical :: sparse_rcm = .false.
        logical :: sparse_amd = .true.
        logical :: equilibrate = .true.
        integer :: equilibration_iters = 3
        real(dp) :: equilibration_min = 1.0e-4_dp
        real(dp) :: equilibration_max = 1.0e4_dp
        logical :: dynamic_regularization = .true.
        integer :: max_regularization_updates = 6
        integer :: certificate_maxit = 80
        integer :: dense_diagnostic_limit = 256
        integer :: iterative_refinement = 3
        real(dp) :: refinement_tol = 1.0e-10_dp
    end type ecos_settings

    type, public :: ecos_csc_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: colptr(:)  ! 1-based pointer, size ncol+1
        integer, allocatable :: rowind(:)  ! 1-based row index
        real(dp), allocatable :: values(:)
    contains
        procedure :: to_dense => csc_to_dense
    end type ecos_csc_matrix

    type, public :: ecos_csr_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: rowptr(:)  ! 1-based pointer, size nrow+1
        integer, allocatable :: colind(:)  ! 1-based column index
        real(dp), allocatable :: values(:)
    contains
        procedure :: to_dense => csr_to_dense
    end type ecos_csr_matrix

    type, public :: ecos_problem
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: gmat(:,:)
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: amat(:,:)
        real(dp), allocatable :: b(:)
        type(ecos_csc_matrix) :: g_csc
        type(ecos_csr_matrix) :: g_csr
        type(ecos_csc_matrix) :: a_csc
        type(ecos_csr_matrix) :: a_csr
        logical :: sparse_storage = .false.
        type(ecos_dims) :: dims
        integer, allocatable :: bool_vars(:)
        integer, allocatable :: int_vars(:)
    contains
        procedure :: valid => problem_valid
        procedure :: nvar => problem_nvar
        procedure :: ncone => problem_ncone
        procedure :: neq => problem_neq
    end type ecos_problem

    type, public :: ecos_result
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: z(:)
        integer :: exitflag = ECOS_FATAL
        integer :: iter = 0
        integer :: mi_iter = -1
        integer :: numerr = 0
        character(len=96) :: infostring = 'Not solved'
        real(dp) :: pcost = huge(1.0_dp)
        real(dp) :: dcost = huge(1.0_dp)
        real(dp) :: pres = huge(1.0_dp)
        real(dp) :: dres = huge(1.0_dp)
        real(dp) :: pinf = 0.0_dp
        real(dp) :: dinf = 0.0_dp
        real(dp) :: pinfres = huge(1.0_dp)
        real(dp) :: dinfres = huge(1.0_dp)
        real(dp) :: gap = huge(1.0_dp)
        real(dp) :: relgap = huge(1.0_dp)
        real(dp) :: r0 = 0.0_dp
        logical :: sparse_backend_used = .false.
        integer :: kkt_nnz = 0
        integer :: ldl_nnz = 0
        integer :: iterative_refinements = 0
        integer :: symbolic_analyses = 0
        integer :: numeric_factorizations = 0
        integer :: regularization_updates = 0
        integer :: cached_symbolic_reuses = 0
        integer :: cached_warm_starts = 0
        integer :: bb_symbolic_reuses = 0
        real(dp) :: min_col_scale = 1.0_dp
        real(dp) :: max_col_scale = 1.0_dp
        real(dp) :: min_row_scale = 1.0_dp
        real(dp) :: max_row_scale = 1.0_dp
        real(dp) :: factor_fill_ratio = 0.0_dp
        real(dp) :: time_ordering = 0.0_dp
        real(dp) :: time_factorization = 0.0_dp
        real(dp) :: time_refinement = 0.0_dp
        logical :: primal_certificate_valid = .false.
        logical :: dual_certificate_valid = .false.
        real(dp), allocatable :: primal_certificate(:)
        real(dp), allocatable :: dual_certificate(:)
    end type ecos_result

    type, public :: ecos_sparse_cache
        logical :: symbolic_valid = .false.
        integer :: n = 0
        integer :: structure_hash = 0
        integer :: symbolic_nnz = 0
        integer, allocatable :: perm(:)
        integer, allocatable :: parent(:)
        integer, allocatable :: lp(:)
        real(dp), allocatable :: warm_x(:)
        real(dp), allocatable :: warm_y(:)
        real(dp), allocatable :: warm_sl(:)
        real(dp), allocatable :: warm_lam(:)
        logical :: warm_valid = .false.
    end type ecos_sparse_cache

    type, public :: ecos_workspace
        type(ecos_problem) :: problem
        type(ecos_settings) :: settings
        type(ecos_sparse_cache) :: sparse_cache
        logical :: initialized = .false.
    end type ecos_workspace

    public :: make_csc_matrix, csc_from_zero_based

contains

    integer function dims_cone_rows(this) result(m)
        class(ecos_dims), intent(in) :: this
        m = this%l + 3*this%e
        if (allocated(this%q)) m = m + sum(this%q)
    end function dims_cone_rows

    integer function dims_scalar_inequalities(this) result(m)
        class(ecos_dims), intent(in) :: this
        m = this%l + 3*this%e
        if (allocated(this%q)) m = m + size(this%q)
    end function dims_scalar_inequalities

    integer function problem_nvar(this) result(n)
        class(ecos_problem), intent(in) :: this
        if (allocated(this%c)) then
            n = size(this%c)
        else
            n = 0
        end if
    end function problem_nvar

    integer function problem_ncone(this) result(m)
        class(ecos_problem), intent(in) :: this
        if (allocated(this%h)) then
            m = size(this%h)
        else
            m = 0
        end if
    end function problem_ncone

    integer function problem_neq(this) result(p)
        class(ecos_problem), intent(in) :: this
        if (allocated(this%b)) then
            p = size(this%b)
        else
            p = 0
        end if
    end function problem_neq

    logical function problem_valid(this) result(ok)
        class(ecos_problem), intent(in) :: this
        integer :: n, m, p
        ok = .false.
        if (.not. allocated(this%c)) return
        n = size(this%c)
        if (n < 1) return
        m = this%dims%cone_rows()
        if (.not.allocated(this%h)) return
        if (size(this%h) /= m) return
        p = this%neq()
        if (this%sparse_storage) then
            if (m > 0) then
                if (this%g_csc%nrow /= m .or. this%g_csc%ncol /= n) return
                if (.not.allocated(this%g_csc%colptr)) return
                if (this%g_csr%nrow /= m .or. this%g_csr%ncol /= n) return
            end if
            if (p > 0) then
                if (this%a_csc%nrow /= p .or. this%a_csc%ncol /= n) return
                if (.not.allocated(this%a_csc%colptr)) return
                if (this%a_csr%nrow /= p .or. this%a_csr%ncol /= n) return
            end if
        else
            if (m > 0) then
                if (.not. allocated(this%gmat)) return
                if (size(this%gmat,1) /= m .or. size(this%gmat,2) /= n) return
            else if (allocated(this%gmat)) then
                if (size(this%gmat,2) /= n .or. size(this%gmat,1) /= 0) return
            end if
            if (p > 0) then
                if (.not. allocated(this%amat)) return
                if (size(this%amat,1) /= p .or. size(this%amat,2) /= n) return
            else if (allocated(this%amat)) then
                if (size(this%amat,2) /= n .or. size(this%amat,1) /= 0) return
            end if
        end if
        ok = .true.
    end function problem_valid

    subroutine csc_to_dense(this, a)
        class(ecos_csc_matrix), intent(in) :: this
        real(dp), allocatable, intent(out) :: a(:,:)
        integer :: j, k
        allocate(a(this%nrow, this%ncol))
        a = 0.0_dp
        do j = 1, this%ncol
            do k = this%colptr(j), this%colptr(j+1)-1
                a(this%rowind(k),j) = a(this%rowind(k),j) + this%values(k)
            end do
        end do
    end subroutine csc_to_dense

    subroutine csr_to_dense(this, a)
        class(ecos_csr_matrix), intent(in) :: this
        real(dp), allocatable, intent(out) :: a(:,:)
        integer :: i, k
        allocate(a(this%nrow, this%ncol))
        a = 0.0_dp
        do i = 1, this%nrow
            do k = this%rowptr(i), this%rowptr(i+1)-1
                a(i,this%colind(k)) = a(i,this%colind(k)) + this%values(k)
            end do
        end do
    end subroutine csr_to_dense

    subroutine make_csc_matrix(a, csc, tol)
        real(dp), intent(in) :: a(:,:)
        type(ecos_csc_matrix), intent(out) :: csc
        real(dp), intent(in), optional :: tol
        real(dp) :: threshold
        integer :: i, j, k, nnz
        threshold = 0.0_dp
        if (present(tol)) threshold = max(0.0_dp, tol)
        nnz = count(abs(a) > threshold)
        csc%nrow = size(a,1)
        csc%ncol = size(a,2)
        allocate(csc%colptr(csc%ncol+1), csc%rowind(nnz), csc%values(nnz))
        k = 1
        csc%colptr(1) = 1
        do j = 1, csc%ncol
            do i = 1, csc%nrow
                if (abs(a(i,j)) > threshold) then
                    csc%rowind(k) = i
                    csc%values(k) = a(i,j)
                    k = k + 1
                end if
            end do
            csc%colptr(j+1) = k
        end do
    end subroutine make_csc_matrix

    subroutine csc_from_zero_based(nrow, ncol, matbeg, matind, values, csc, ierr)
        integer, intent(in) :: nrow, ncol
        integer, intent(in) :: matbeg(:), matind(:)
        real(dp), intent(in) :: values(:)
        type(ecos_csc_matrix), intent(out) :: csc
        integer, intent(out) :: ierr
        ierr = 0
        if (size(matbeg) /= ncol+1 .or. size(matind) /= size(values)) then
            ierr = 1
            return
        end if
        csc%nrow = nrow
        csc%ncol = ncol
        allocate(csc%colptr(ncol+1), csc%rowind(size(matind)), csc%values(size(values)))
        csc%colptr = matbeg + 1
        csc%rowind = matind + 1
        csc%values = values
        if (any(csc%rowind < 1) .or. any(csc%rowind > nrow)) ierr = 2
    end subroutine csc_from_zero_based

end module ecos_types
