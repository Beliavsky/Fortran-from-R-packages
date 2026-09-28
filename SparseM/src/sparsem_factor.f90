module sparsem_factor
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_value
    use sparsem_kinds, only: dp
    use sparsem_types, only: csr_matrix, csc_matrix, csr_cholesky, csr_cholesky_stats, slm_result, slm_summary_result
    use sparsem_types, only: mslm_result, mslm_summary_result
    use sparsem_core, only: csr_to_dense, csr_from_dense, csr_from_csc, csr_transpose, csr_matmul, csr_matvec
    use sparsem_ng_order, only: minimum_degree_order
    use sparsem_ng_symbolic, only: ng_symbolic_analysis, ng_supernodal_structure
    use sparsem_ng_numeric, only: ng_numeric_factor, ng_forwardsolve, ng_backsolve, ng_solve
    implicit none
    private

    public :: csr_chol, chol_statistics, chol_as_csr, chol_determinant
    public :: chol_forwardsolve, chol_backsolve, chol_solve
    public :: chol_forwardsolve_matrix, chol_backsolve_matrix, chol_solve_matrix
    public :: chol_forwardsolve_csr, chol_backsolve_csr, chol_solve_csr
    public :: csr_solve, csr_solve_matrix, csr_solve_csr_rhs, csr_inverse, csr_determinant
    public :: slm, slm_fit, slm_wfit, slm_fit_csr
    public :: slm_coef, slm_fitted, slm_residuals, slm_deviance, slm_extract_aic
    public :: slm_summary, mslm_fit_csr, mslm_summary

contains

    function csr_chol(a, info, ordering, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(factor)
        type(csr_matrix), intent(in) :: a !! Square symmetric positive-definite CSR matrix to factor.
        integer, intent(out), optional :: info !! Zero on success; 17+k signals k tiny-pivot replacements; negative values report input errors.
        character(len=*), intent(in), optional :: ordering !! Ordering strategy: "minimum_degree" (default) or "natural".
        integer, intent(in), optional :: cache_kb !! Cache size in KiB used for Ng-Peyton supernode splitting; default 1024 and minimum 32.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold used by Ng-Peyton factorization; default 1e-30.
        real(dp), intent(in), optional :: large !! Replacement diagonal value for a tiny pivot; default 1e128 and must be positive.
        integer, intent(in), optional :: level !! Ng-Peyton loop-unrolling level: 1, 2, 4, or 8; default 8 as in SparseM.
        real(dp), intent(in), optional :: eps !! Maximum allowed absolute asymmetry; default machine epsilon, and zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional maximum compressed symbolic-index storage; exceeding it returns SparseM status 5.
        integer, intent(in), optional :: nnzlmax !! Optional maximum numerical factor entries; exceeding it returns SparseM status 4.
        integer, intent(in), optional :: tmpmax !! Optional maximum real block workspace; exceeding it returns SparseM status 8.
        type(csr_cholesky) :: factor
        type(csc_matrix) :: lower_csc
        real(dp), allocatable :: lnz(:)
        integer, allocatable :: split(:)
        integer, allocatable :: col_count(:)
        integer, allocatable :: etree_parent(:)
        integer, allocatable :: inv_perm(:)
        integer, allocatable :: lindx(:)
        integer, allocatable :: perm(:)
        integer, allocatable :: supernode(:)
        integer, allocatable :: xlindx(:)
        integer, allocatable :: xlnz(:)
        integer, allocatable :: xsuper(:)
        integer :: cache_size
        integer :: i
        integer :: nnzl
        integer :: nnzl_limit
        integer :: nofsub
        integer :: nsub
        integer :: nsub_limit
        integer :: nsuper
        integer :: tmp_limit
        integer :: numeric_level
        integer :: stat
        real(dp) :: asymmetry
        real(dp) :: large_value
        real(dp) :: symmetry_tol
        real(dp) :: tiny_value
        character(len=16) :: selected_ordering

        if (present(info)) info = 0
        if (a%nrow /= a%ncol) then
            call initialize_empty_cholesky(factor)
            factor%info = -1
            if (present(info)) info = factor%info
            return
        end if

        cache_size = 1024
        if (present(cache_kb)) cache_size = cache_kb
        if (cache_size < 32) then
            call initialize_empty_cholesky(factor)
            factor%info = -5
            if (present(info)) info = factor%info
            return
        end if

        numeric_level = 8
        if (present(level)) numeric_level = level
        if (all(numeric_level /= [1, 2, 4, 8])) then
            call initialize_empty_cholesky(factor)
            factor%info = -4
            if (present(info)) info = factor%info
            return
        end if

        tiny_value = 1.0e-30_dp
        if (present(tiny)) tiny_value = tiny
        large_value = 1.0e128_dp
        if (present(large)) large_value = large
        symmetry_tol = epsilon(1.0_dp)
        if (present(eps)) symmetry_tol = eps
        if (tiny_value < 0.0_dp .or. large_value <= 0.0_dp .or. symmetry_tol < 0.0_dp) then
            call initialize_empty_cholesky(factor)
            factor%info = -6
            if (present(info)) info = factor%info
            return
        end if

        nsub_limit = size(a%values)
        if (present(nsubmax)) nsub_limit = nsubmax
        nnzl_limit = max(4 * size(a%values), int(floor(0.2_dp * real(size(a%values), dp)**1.3_dp)))
        if (present(nnzlmax)) nnzl_limit = nnzlmax
        tmp_limit = 50 * a%nrow
        if (present(tmpmax)) tmp_limit = tmpmax
        if (nsub_limit < 0 .or. nnzl_limit < 0 .or. tmp_limit < 0) then
            call initialize_empty_cholesky(factor)
            factor%info = -6
            if (present(info)) info = factor%info
            return
        end if

        if (symmetry_tol > 0.0_dp) then
            asymmetry = csr_max_asymmetry(a)
            if (asymmetry > symmetry_tol) then
                call initialize_empty_cholesky(factor)
                factor%info = -2
                if (present(info)) info = factor%info
                return
            end if
        end if

        selected_ordering = "minimum_degree"
        if (present(ordering)) selected_ordering = adjustl(ordering)
        select case (trim(selected_ordering))
        case ("minimum_degree", "mmd", "ng-peyton", "ng_peyton")
            call minimum_degree_order(a, perm, inv_perm, nofsub, stat)
            if (stat /= 0) then
                call initialize_empty_cholesky(factor)
                factor%info = stat
                if (present(info)) info = stat
                return
            end if
            call ng_symbolic_analysis(a, perm, inv_perm, etree_parent, col_count, nsuper, supernode, xsuper, &
                nnzl, nsub, stat)
            if (stat /= 0) then
                call initialize_empty_cholesky(factor)
                factor%info = stat
                if (present(info)) info = stat
                return
            end if
            if (nnzl > nnzl_limit) then
                call initialize_empty_cholesky(factor)
                factor%info = 4
                if (present(info)) info = factor%info
                return
            end if
            if (nsub > nsub_limit) then
                call initialize_empty_cholesky(factor)
                factor%info = 5
                if (present(info)) info = factor%info
                return
            end if
            call ng_supernodal_structure(a, perm, inv_perm, col_count, nsuper, xsuper, supernode, nsub, &
                xlindx, lindx, xlnz, stat)
            if (stat /= 0) then
                call initialize_empty_cholesky(factor)
                factor%info = stat
                if (present(info)) info = stat
                return
            end if
            call ng_numeric_factor(a, perm, inv_perm, nsuper, xsuper, supernode, xlindx, lindx, xlnz, &
                lnz, split, factor%tmp_size, stat, cache_kb=cache_size, tiny=tiny_value, large=large_value, &
                level=numeric_level, tmpmax=tmp_limit)
            if (stat == -1) stat = 9
            if (stat == -2) stat = 10
            if (stat < 0 .or. (stat > 0 .and. stat < 17)) then
                call initialize_empty_cholesky(factor)
                factor%info = stat
                factor%numeric_info = stat
                if (present(info)) info = stat
                return
            end if
            call supernodal_factor_to_csc(a%nrow, nsuper, xsuper, xlindx, lindx, xlnz, lnz, lower_csc)
            factor%n = a%nrow
            factor%lower = csr_from_csc(lower_csc)
            factor%upper = csr_transpose(factor%lower)
            factor%ordering = "minimum_degree"
            write (factor%numeric_kernel, '(a,i0)') "ng_peyton_level", numeric_level
            factor%nofsub = nofsub
            factor%nnzl = nnzl
            factor%nsub = nsub
            factor%nsuper = nsuper
            factor%etree_parent = etree_parent
            factor%col_count = col_count
            factor%supernode = supernode
            factor%xsuper = xsuper
            factor%xlindx = xlindx
            factor%lindx = lindx
            factor%xlnz = xlnz
            factor%lnz = lnz
            factor%split = split
            factor%factor_log_determinant = 0.0_dp
            do i = 1, factor%n
                factor%factor_log_determinant = factor%factor_log_determinant + log(lnz(xlnz(i)))
            end do
            factor%factor_determinant = exp_or_positive_inf(factor%factor_log_determinant)
            factor%log_determinant = 2.0_dp * factor%factor_log_determinant
            factor%determinant = exp_or_positive_inf(factor%log_determinant)
            factor%info = 0
            factor%numeric_info = stat
            if (stat >= 17) factor%tiny_pivots = stat - 16
        case ("natural")
            allocate(perm(a%nrow), inv_perm(a%nrow))
            do i = 1, a%nrow
                perm(i) = i
                inv_perm(i) = i
            end do
            factor = csr_chol_natural(a, stat)
            factor%ordering = "natural"
            factor%numeric_kernel = "row_sparse"
            factor%nofsub = 0
            factor%nnzl = size(factor%lower%values)
            factor%nsub = 0
            factor%nsuper = 0
            allocate(factor%etree_parent(0), factor%col_count(0), factor%supernode(0), factor%xsuper(1))
            allocate(factor%xlindx(1), factor%lindx(0), factor%xlnz(1))
            allocate(factor%lnz(0), factor%split(0))
            factor%tmp_size = 0
            factor%xsuper(1) = 1
            factor%xlindx(1) = 1
            factor%xlnz(1) = 1
        case default
            call initialize_empty_cholesky(factor)
            factor%info = -3
            if (present(info)) info = factor%info
            return
        end select

        factor%perm = perm
        factor%inv_perm = inv_perm
        factor%cache_kb = cache_size
        factor%level = numeric_level
        factor%tiny = tiny_value
        factor%large = large_value
        factor%symmetry_tolerance = symmetry_tol
        if (present(info)) info = stat
    end function csr_chol

    pure function chol_statistics(factor) result(stats)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor whose Ng-Peyton storage statistics are requested.
        type(csr_cholesky_stats) :: stats
        integer :: column
        integer :: column_length
        integer :: jlen
        integer :: ncols
        integer :: super
        integer :: super_nonzeros

        stats%nsuper = factor%nsuper
        stats%nnzl = factor%nnzl
        stats%nsub = factor%nsub
        stats%tmp_size = factor%tmp_size
        stats%supernodal = factor%nsuper > 0 .and. allocated(factor%xsuper) .and. &
            allocated(factor%xlindx) .and. allocated(factor%xlnz)
        if (.not. stats%supernodal) return
        if (size(factor%xsuper) < factor%nsuper + 1 .or. size(factor%xlindx) < factor%nsuper + 1 .or. &
                size(factor%xlnz) < factor%n + 1) then
            stats%supernodal = .false.
            return
        end if

        do super = 1, factor%nsuper
            ncols = factor%xsuper(super + 1) - factor%xsuper(super)
            jlen = factor%xlindx(super + 1) - factor%xlindx(super)
            super_nonzeros = ((2 * jlen - ncols + 1) * ncols) / 2
            stats%max_supernode_columns = max(stats%max_supernode_columns, ncols)
            stats%max_supernode_nonzeros = max(stats%max_supernode_nonzeros, super_nonzeros)
        end do

        do column = 1, factor%n
            column_length = factor%xlnz(column + 1) - factor%xlnz(column)
            stats%factorization_flops = stats%factorization_flops + &
                real(column_length * column_length - 1, dp)
            stats%triangular_solve_flops = stats%triangular_solve_flops + &
                real(2 * column_length - 1, dp)
        end do
        stats%triangular_solve_flops = 2.0_dp * stats%triangular_solve_flops
    end function chol_statistics

    function chol_as_csr(factor, upper_tri, info) result(a)
        type(csr_cholesky), intent(in) :: factor !! Cholesky object whose triangular numerical factor is requested.
        logical, intent(in), optional :: upper_tri !! Return upper R by default; false returns the stored lower transpose.
        integer, intent(out), optional :: info !! Factor status; zero for a usable factor and fatal status otherwise.
        type(csr_matrix) :: a
        logical :: use_upper

        if (present(info)) info = factor%info
        if (factor%info /= 0) then
            a%nrow = 0
            a%ncol = 0
            allocate(a%values(0), a%col_ind(0), a%row_ptr(1))
            a%row_ptr(1) = 1
            return
        end if

        use_upper = .true.
        if (present(upper_tri)) use_upper = upper_tri
        if (use_upper) then
            a = factor%upper
        else
            a = factor%lower
        end if
    end function chol_as_csr

    function chol_determinant(factor, logarithm, sign, info) result(value)
        type(csr_cholesky), intent(in) :: factor !! Cholesky object whose triangular-factor determinant is requested.
        logical, intent(in), optional :: logarithm !! Return the log modulus when true/default; determinant otherwise.
        integer, intent(out), optional :: sign !! Sign of the triangular-factor determinant; one for a usable factor.
        integer, intent(out), optional :: info !! Fatal factor status or retained tiny-pivot warning status.
        real(dp) :: value
        logical :: use_log
        integer :: stat

        stat = factor%info
        if (stat == 0 .and. factor%numeric_info >= 17) stat = factor%numeric_info
        if (present(info)) info = stat
        if (factor%info /= 0) then
            value = 0.0_dp
            if (present(sign)) sign = 0
            return
        end if

        if (present(sign)) sign = 1
        use_log = .true.
        if (present(logarithm)) use_log = logarithm
        if (use_log) then
            value = factor%factor_log_determinant
        else
            value = factor%factor_determinant
        end if
    end function chol_determinant

    pure function exp_or_positive_inf(log_value) result(value)
        real(dp), intent(in) :: log_value !! Natural logarithm of a nonnegative quantity to exponentiate safely.
        real(dp) :: value

        if (log_value > log(huge(1.0_dp))) then
            value = ieee_value(0.0_dp, ieee_positive_inf)
        else
            value = exp(log_value)
        end if
    end function exp_or_positive_inf

    pure function csr_max_asymmetry(a) result(max_error)
        type(csr_matrix), intent(in) :: a !! Square CSR matrix whose A(i,j)-A(j,i) differences are inspected.
        real(dp) :: max_error
        real(dp) :: a_ij
        real(dp) :: a_ji
        integer :: column
        integer :: pos
        integer :: row

        max_error = 0.0_dp
        do row = 1, a%nrow
            do pos = a%row_ptr(row), a%row_ptr(row + 1) - 1
                column = a%col_ind(pos)
                if (column < 1 .or. column > a%ncol) then
                    max_error = huge(1.0_dp)
                    return
                end if
                a_ij = csr_entry_sum(a, row, column)
                a_ji = csr_entry_sum(a, column, row)
                max_error = max(max_error, abs(a_ij - a_ji))
            end do
        end do
    end function csr_max_asymmetry

    pure function csr_entry_sum(a, row, column) result(value)
        type(csr_matrix), intent(in) :: a !! CSR matrix searched for all entries at one coordinate.
        integer, intent(in) :: row !! One-based row index to search.
        integer, intent(in) :: column !! One-based column index whose duplicate entries are summed.
        real(dp) :: value
        integer :: pos

        value = 0.0_dp
        if (row < 1 .or. row > a%nrow) return
        do pos = a%row_ptr(row), a%row_ptr(row + 1) - 1
            if (a%col_ind(pos) == column) value = value + a%values(pos)
        end do
    end function csr_entry_sum

    function csr_chol_natural(a, info) result(factor)
        type(csr_matrix), intent(in) :: a !! Square symmetric positive-definite CSR matrix already arranged in factorization order.
        integer, intent(out) :: info !! Zero on success; positive index of the first non-positive pivot; -1 for nonsquare input.
        type(csr_cholesky) :: factor
        real(dp), allocatable :: current_row(:)
        real(dp), allocatable :: original_row(:)
        real(dp), allocatable :: diagonal(:)
        real(dp) :: pivot
        real(dp) :: s
        integer :: capacity
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz
        integer :: pos

        factor%n = a%nrow
        factor%info = 0
        factor%determinant = 0.0_dp
        factor%log_determinant = 0.0_dp
        factor%factor_determinant = 0.0_dp
        factor%factor_log_determinant = 0.0_dp
        if (a%nrow /= a%ncol) then
            call initialize_empty_cholesky(factor)
            factor%info = -1
            info = factor%info
            return
        end if

        factor%lower%nrow = a%nrow
        factor%lower%ncol = a%ncol
        allocate(factor%lower%row_ptr(a%nrow + 1))
        factor%lower%row_ptr(1) = 1
        capacity = max(1, max(a%nrow, size(a%values)))
        allocate(factor%lower%values(capacity))
        allocate(factor%lower%col_ind(capacity))
        allocate(current_row(a%nrow), source=0.0_dp)
        allocate(original_row(a%nrow), source=0.0_dp)
        allocate(diagonal(a%nrow), source=0.0_dp)
        nnz = 0

        do i = 1, a%nrow
            current_row = 0.0_dp
            original_row = 0.0_dp
            do pos = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col_ind(pos)
                if (j <= i) original_row(j) = a%values(pos)
            end do
            do k = 1, i - 1
                s = original_row(k)
                do pos = factor%lower%row_ptr(k), factor%lower%row_ptr(k + 1) - 1
                    j = factor%lower%col_ind(pos)
                    if (j >= k) exit
                    if (current_row(j) /= 0.0_dp) then
                        s = s - current_row(j) * factor%lower%values(pos)
                    end if
                end do
                if (s /= 0.0_dp) then
                    current_row(k) = s / diagonal(k)
                    call append_cholesky_entry(factor%lower, nnz, capacity, k, current_row(k))
                end if
            end do

            pivot = original_row(i) - dot_product(current_row(1:i - 1), current_row(1:i - 1))
            if (pivot <= 0.0_dp) then
                factor%info = i
                factor%lower%row_ptr(i + 1:) = nnz + 1
                call trim_csr_storage(factor%lower, nnz)
                factor%upper = csr_transpose(factor%lower)
                info = factor%info
                return
            end if
            diagonal(i) = sqrt(pivot)
            call append_cholesky_entry(factor%lower, nnz, capacity, i, diagonal(i))
            factor%lower%row_ptr(i + 1) = nnz + 1
        end do

        call trim_csr_storage(factor%lower, nnz)
        factor%upper = csr_transpose(factor%lower)
        if (factor%n == 0) then
            factor%factor_log_determinant = 0.0_dp
            factor%factor_determinant = 1.0_dp
            factor%log_determinant = 0.0_dp
            factor%determinant = 1.0_dp
        else
            factor%factor_log_determinant = sum(log(diagonal))
            factor%factor_determinant = exp_or_positive_inf(factor%factor_log_determinant)
            factor%log_determinant = 2.0_dp * factor%factor_log_determinant
            factor%determinant = exp_or_positive_inf(factor%log_determinant)
        end if
        info = 0
    end function csr_chol_natural

    function chol_forwardsolve(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor; stored triangular factors are in factor%perm ordering.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector in original matrix ordering, of length factor%n.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for RHS size mismatch.
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: ordered_b(:)
        real(dp), allocatable :: ordered_x(:)
        real(dp) :: diagonal
        real(dp) :: s
        integer :: column
        integer :: i
        integer :: pos

        if (present(info)) info = 0
        if (factor%info /= 0) then
            if (present(info)) info = factor%info
            allocate(x(0))
            return
        end if
        if (size(b) /= factor%n) then
            if (present(info)) info = -2
            allocate(x(0))
            return
        end if

        allocate(ordered_b(factor%n), ordered_x(factor%n), x(factor%n))
        if (allocated(factor%perm)) then
            ordered_b = b(factor%perm)
        else
            ordered_b = b
        end if
        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            ordered_x = ordered_b
            call ng_forwardsolve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, &
                factor%lnz, ordered_x)
        else
            do i = 1, factor%n
                s = ordered_b(i)
                diagonal = 0.0_dp
                do pos = factor%lower%row_ptr(i), factor%lower%row_ptr(i + 1) - 1
                    column = factor%lower%col_ind(pos)
                    if (column < i) then
                        s = s - factor%lower%values(pos) * ordered_x(column)
                    else if (column == i) then
                        diagonal = factor%lower%values(pos)
                    end if
                end do
                ordered_x(i) = s / diagonal
            end do
        end if

        if (allocated(factor%perm)) then
            do i = 1, factor%n
                x(factor%perm(i)) = ordered_x(i)
            end do
        else
            x = ordered_x
        end if
    end function chol_forwardsolve

    function chol_backsolve(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor; solves the transposed stored factor with permutation handling.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector in original matrix ordering, of length factor%n.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for RHS size mismatch.
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: ordered_b(:)
        real(dp), allocatable :: ordered_x(:)
        real(dp) :: diagonal
        real(dp) :: s
        integer :: column
        integer :: i
        integer :: pos

        if (present(info)) info = 0
        if (factor%info /= 0) then
            if (present(info)) info = factor%info
            allocate(x(0))
            return
        end if
        if (size(b) /= factor%n) then
            if (present(info)) info = -2
            allocate(x(0))
            return
        end if

        allocate(ordered_b(factor%n), ordered_x(factor%n), x(factor%n))
        if (allocated(factor%perm)) then
            ordered_b = b(factor%perm)
        else
            ordered_b = b
        end if
        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            ordered_x = ordered_b
            call ng_backsolve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, &
                factor%lnz, ordered_x)
        else
            do i = factor%n, 1, -1
                s = ordered_b(i)
                diagonal = 0.0_dp
                do pos = factor%upper%row_ptr(i), factor%upper%row_ptr(i + 1) - 1
                    column = factor%upper%col_ind(pos)
                    if (column == i) then
                        diagonal = factor%upper%values(pos)
                    else if (column > i) then
                        s = s - factor%upper%values(pos) * ordered_x(column)
                    end if
                end do
                ordered_x(i) = s / diagonal
            end do
        end if

        if (allocated(factor%perm)) then
            do i = 1, factor%n
                x(factor%perm(i)) = ordered_x(i)
            end do
        else
            x = ordered_x
        end if
    end function chol_backsolve

    function chol_solve(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor for A in the system A x = b.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector of length factor%n.
        integer, intent(out), optional :: info !! Zero on success; negative for size issues; otherwise factorization status.
        real(dp), allocatable :: x(:)
        real(dp) :: y(factor%n)
        integer :: stat

        if (factor%info /= 0) then
            if (present(info)) info = factor%info
            allocate(x(0))
            return
        end if
        if (size(b) /= factor%n) then
            if (present(info)) info = -2
            allocate(x(0))
            return
        end if
        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            allocate(x(factor%n))
            if (allocated(factor%perm)) then
                y = b(factor%perm)
            else
                y = b
            end if
            call ng_solve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, factor%lnz, y)
            if (allocated(factor%perm)) then
                do stat = 1, factor%n
                    x(factor%perm(stat)) = y(stat)
                end do
            else
                x = y
            end if
            stat = 0
        else
            y = chol_forwardsolve(factor, b, stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                allocate(x(0))
                return
            end if
            x = chol_backsolve(factor, y, stat)
        end if
        if (present(info)) info = stat
    end function chol_solve

    function chol_forwardsolve_matrix(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor containing lower-triangular L with A = L L-transpose.
        real(dp), intent(in) :: b(:,:) !! Dense right-hand-side matrix with factor%n rows and one or more RHS columns.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for RHS row mismatch.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: col(:)
        integer :: i
        integer :: j
        integer :: stat

        if (present(info)) info = 0
        if (factor%info /= 0) then
            allocate(x(0, 0))
            if (present(info)) info = factor%info
            return
        end if
        if (size(b, 1) /= factor%n) then
            allocate(x(0, 0))
            if (present(info)) info = -2
            return
        end if

        allocate(x(factor%n, size(b, 2)))
        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            allocate(work(factor%n))
            do j = 1, size(b, 2)
                if (allocated(factor%perm)) then
                    work = b(factor%perm, j)
                else
                    work = b(:, j)
                end if
                call ng_forwardsolve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, &
                    factor%lnz, work)
                if (allocated(factor%perm)) then
                    do i = 1, factor%n
                        x(factor%perm(i), j) = work(i)
                    end do
                else
                    x(:, j) = work
                end if
            end do
        else
            do j = 1, size(b, 2)
                col = chol_forwardsolve(factor, b(:, j), stat)
                if (stat /= 0) then
                    if (present(info)) info = stat
                    return
                end if
                x(:, j) = col
            end do
        end if
    end function chol_forwardsolve_matrix

    function chol_backsolve_matrix(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor containing lower-triangular L; solves L-transpose X = B.
        real(dp), intent(in) :: b(:,:) !! Dense right-hand-side matrix with factor%n rows and one or more RHS columns.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for RHS row mismatch.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: col(:)
        integer :: i
        integer :: j
        integer :: stat

        if (present(info)) info = 0
        if (factor%info /= 0) then
            allocate(x(0, 0))
            if (present(info)) info = factor%info
            return
        end if
        if (size(b, 1) /= factor%n) then
            allocate(x(0, 0))
            if (present(info)) info = -2
            return
        end if

        allocate(x(factor%n, size(b, 2)))
        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            allocate(work(factor%n))
            do j = 1, size(b, 2)
                if (allocated(factor%perm)) then
                    work = b(factor%perm, j)
                else
                    work = b(:, j)
                end if
                call ng_backsolve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, &
                    factor%lnz, work)
                if (allocated(factor%perm)) then
                    do i = 1, factor%n
                        x(factor%perm(i), j) = work(i)
                    end do
                else
                    x(:, j) = work
                end if
            end do
        else
            do j = 1, size(b, 2)
                col = chol_backsolve(factor, b(:, j), stat)
                if (stat /= 0) then
                    if (present(info)) info = stat
                    return
                end if
                x(:, j) = col
            end do
        end if
    end function chol_backsolve_matrix

    function chol_solve_matrix(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor for A in the matrix system A X = B.
        real(dp), intent(in) :: b(:,:) !! Dense right-hand-side matrix with factor%n rows.
        integer, intent(out), optional :: info !! Zero on success; negative for size issues; otherwise factorization status.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: y(:,:)
        integer :: i
        integer :: j
        integer :: stat

        if (present(info)) info = 0
        if (factor%info /= 0) then
            allocate(x(0, 0))
            if (present(info)) info = factor%info
            return
        end if
        if (size(b, 1) /= factor%n) then
            allocate(x(0, 0))
            if (present(info)) info = -2
            return
        end if

        if (allocated(factor%lnz) .and. size(factor%lnz) > 0 .and. factor%nsuper > 0) then
            allocate(x(factor%n, size(b, 2)), work(factor%n))
            do j = 1, size(b, 2)
                if (allocated(factor%perm)) then
                    work = b(factor%perm, j)
                else
                    work = b(:, j)
                end if
                call ng_solve(factor%nsuper, factor%xsuper, factor%xlindx, factor%lindx, factor%xlnz, &
                    factor%lnz, work)
                if (allocated(factor%perm)) then
                    do i = 1, factor%n
                        x(factor%perm(i), j) = work(i)
                    end do
                else
                    x(:, j) = work
                end if
            end do
        else
            y = chol_forwardsolve_matrix(factor, b, stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                allocate(x(0, 0))
                return
            end if
            x = chol_backsolve_matrix(factor, y, stat)
            if (present(info)) info = stat
        end if
    end function chol_solve_matrix

    function chol_forwardsolve_csr(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor used for the lower-triangular solve.
        type(csr_matrix), intent(in) :: b !! Sparse right-hand-side matrix with one equation vector per column.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for row mismatch.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: dense_b(:,:)

        dense_b = csr_to_dense(b)
        x = chol_forwardsolve_matrix(factor, dense_b, info)
    end function chol_forwardsolve_csr

    function chol_backsolve_csr(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor used for the upper-triangular solve.
        type(csr_matrix), intent(in) :: b !! Sparse right-hand-side matrix with one equation vector per column.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for row mismatch.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: dense_b(:,:)

        dense_b = csr_to_dense(b)
        x = chol_backsolve_matrix(factor, dense_b, info)
    end function chol_backsolve_csr

    function chol_solve_csr(factor, b, info) result(x)
        type(csr_cholesky), intent(in) :: factor !! Cholesky factor used for the two-pass linear solve.
        type(csr_matrix), intent(in) :: b !! Sparse right-hand-side matrix with one equation vector per column.
        integer, intent(out), optional :: info !! Zero on success, factor status if invalid, or -2 for row mismatch.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: dense_b(:,:)

        dense_b = csr_to_dense(b)
        x = chol_solve_matrix(factor, dense_b, info)
    end function chol_solve_csr

    function csr_solve(a, b, info, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(x)
        type(csr_matrix), intent(in) :: a !! Symmetric positive-definite CSR coefficient matrix.
        real(dp), intent(in) :: b(:) !! Dense right-hand-side vector.
        integer, intent(out), optional :: info !! Zero on success, 17+k for tiny-pivot replacement, or a fatal status.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB forwarded to Ng-Peyton Cholesky.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold forwarded to Cholesky.
        real(dp), intent(in), optional :: large !! Positive tiny-pivot replacement value forwarded to Cholesky.
        integer, intent(in), optional :: level !! Ng-Peyton numerical level 1, 2, 4, or 8 forwarded to Cholesky.
        real(dp), intent(in), optional :: eps !! Symmetry tolerance forwarded to Cholesky; zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional compressed-symbolic-index ceiling forwarded to Cholesky.
        integer, intent(in), optional :: nnzlmax !! Optional numerical-factor-entry ceiling forwarded to Cholesky.
        integer, intent(in), optional :: tmpmax !! Optional real-workspace ceiling forwarded to Cholesky.
        real(dp), allocatable :: x(:)
        type(csr_cholesky) :: factor
        integer :: chol_stat
        integer :: solve_stat

        factor = csr_chol(a, chol_stat, cache_kb=cache_kb, tiny=tiny, large=large, level=level, eps=eps, &
            nsubmax=nsubmax, nnzlmax=nnzlmax, tmpmax=tmpmax)
        if (chol_stat < 0 .or. (chol_stat > 0 .and. chol_stat < 17)) then
            if (present(info)) info = chol_stat
            allocate(x(0))
            return
        end if
        x = chol_solve(factor, b, solve_stat)
        if (present(info)) then
            if (solve_stat /= 0) then
                info = solve_stat
            else
                info = chol_stat
            end if
        end if
    end function csr_solve

    function csr_solve_matrix(a, b, info, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(x)
        type(csr_matrix), intent(in) :: a !! Symmetric positive-definite CSR coefficient matrix.
        real(dp), intent(in) :: b(:,:) !! Dense right-hand-side matrix with one equation vector per column.
        integer, intent(out), optional :: info !! Zero on success, 17+k for tiny-pivot replacement, or a fatal status.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB forwarded to Ng-Peyton Cholesky.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold forwarded to Cholesky.
        real(dp), intent(in), optional :: large !! Positive tiny-pivot replacement value forwarded to Cholesky.
        integer, intent(in), optional :: level !! Ng-Peyton numerical level 1, 2, 4, or 8 forwarded to Cholesky.
        real(dp), intent(in), optional :: eps !! Symmetry tolerance forwarded to Cholesky; zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional compressed-symbolic-index ceiling forwarded to Cholesky.
        integer, intent(in), optional :: nnzlmax !! Optional numerical-factor-entry ceiling forwarded to Cholesky.
        integer, intent(in), optional :: tmpmax !! Optional real-workspace ceiling forwarded to Cholesky.
        real(dp), allocatable :: x(:,:)
        type(csr_cholesky) :: factor
        integer :: chol_stat
        integer :: solve_stat

        factor = csr_chol(a, chol_stat, cache_kb=cache_kb, tiny=tiny, large=large, level=level, eps=eps, &
            nsubmax=nsubmax, nnzlmax=nnzlmax, tmpmax=tmpmax)
        if (chol_stat < 0 .or. (chol_stat > 0 .and. chol_stat < 17)) then
            if (present(info)) info = chol_stat
            allocate(x(0, 0))
            return
        end if
        x = chol_solve_matrix(factor, b, solve_stat)
        if (present(info)) then
            if (solve_stat /= 0) then
                info = solve_stat
            else
                info = chol_stat
            end if
        end if
    end function csr_solve_matrix

    function csr_solve_csr_rhs(a, b, info, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(x)
        type(csr_matrix), intent(in) :: a !! Symmetric positive-definite CSR coefficient matrix.
        type(csr_matrix), intent(in) :: b !! Sparse right-hand-side matrix with one equation vector per column.
        integer, intent(out), optional :: info !! Zero on success, 17+k for tiny-pivot replacement, or a fatal status.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB forwarded to Ng-Peyton Cholesky.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold forwarded to Cholesky.
        real(dp), intent(in), optional :: large !! Positive tiny-pivot replacement value forwarded to Cholesky.
        integer, intent(in), optional :: level !! Ng-Peyton numerical level 1, 2, 4, or 8 forwarded to Cholesky.
        real(dp), intent(in), optional :: eps !! Symmetry tolerance forwarded to Cholesky; zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional compressed-symbolic-index ceiling forwarded to Cholesky.
        integer, intent(in), optional :: nnzlmax !! Optional numerical-factor-entry ceiling forwarded to Cholesky.
        integer, intent(in), optional :: tmpmax !! Optional real-workspace ceiling forwarded to Cholesky.
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: dense_b(:,:)

        dense_b = csr_to_dense(b)
        x = csr_solve_matrix(a, dense_b, info, cache_kb=cache_kb, tiny=tiny, large=large, level=level, eps=eps, &
            nsubmax=nsubmax, nnzlmax=nnzlmax, tmpmax=tmpmax)
    end function csr_solve_csr_rhs

    function csr_inverse(a, info, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(inverse)
        type(csr_matrix), intent(in) :: a !! Symmetric positive-definite CSR matrix whose inverse is requested.
        integer, intent(out), optional :: info !! Zero on success, 17+k for tiny-pivot replacement, or a fatal status.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB forwarded to Ng-Peyton Cholesky.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold forwarded to Cholesky.
        real(dp), intent(in), optional :: large !! Positive tiny-pivot replacement value forwarded to Cholesky.
        integer, intent(in), optional :: level !! Ng-Peyton numerical level 1, 2, 4, or 8 forwarded to Cholesky.
        real(dp), intent(in), optional :: eps !! Symmetry tolerance forwarded to Cholesky; zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional compressed-symbolic-index ceiling forwarded to Cholesky.
        integer, intent(in), optional :: nnzlmax !! Optional numerical-factor-entry ceiling forwarded to Cholesky.
        integer, intent(in), optional :: tmpmax !! Optional real-workspace ceiling forwarded to Cholesky.
        type(csr_matrix) :: inverse
        type(csr_cholesky) :: factor
        real(dp), allocatable :: identity(:,:)
        real(dp), allocatable :: dense_inverse(:,:)
        integer :: chol_stat
        integer :: i
        integer :: solve_stat

        factor = csr_chol(a, chol_stat, cache_kb=cache_kb, tiny=tiny, large=large, level=level, eps=eps, &
            nsubmax=nsubmax, nnzlmax=nnzlmax, tmpmax=tmpmax)
        if (chol_stat < 0 .or. (chol_stat > 0 .and. chol_stat < 17)) then
            if (present(info)) info = chol_stat
            inverse%nrow = 0
            inverse%ncol = 0
            allocate(inverse%values(0), inverse%col_ind(0), inverse%row_ptr(1))
            inverse%row_ptr(1) = 1
            return
        end if

        allocate(identity(a%nrow, a%nrow), source=0.0_dp)
        do i = 1, a%nrow
            identity(i, i) = 1.0_dp
        end do
        dense_inverse = chol_solve_matrix(factor, identity, solve_stat)
        if (solve_stat /= 0) then
            if (present(info)) info = solve_stat
            inverse%nrow = 0
            inverse%ncol = 0
            allocate(inverse%values(0), inverse%col_ind(0), inverse%row_ptr(1))
            inverse%row_ptr(1) = 1
            return
        end if
        inverse = csr_from_dense(dense_inverse)
        if (present(info)) info = chol_stat
    end function csr_inverse

    function csr_determinant(a, logarithm, sign, info, cache_kb, tiny, large, level, eps, nsubmax, nnzlmax, tmpmax) result(value)
        type(csr_matrix), intent(in) :: a !! Symmetric positive-definite CSR matrix.
        logical, intent(in), optional :: logarithm !! Return log determinant when true; determinant otherwise.
        integer, intent(out), optional :: sign !! Determinant sign; one on successful positive-definite factorization.
        integer, intent(out), optional :: info !! Zero on success, 17+k for tiny-pivot replacement, or a fatal status.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB forwarded to Ng-Peyton Cholesky.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold forwarded to Cholesky.
        real(dp), intent(in), optional :: large !! Positive tiny-pivot replacement value forwarded to Cholesky.
        integer, intent(in), optional :: level !! Ng-Peyton numerical level 1, 2, 4, or 8 forwarded to Cholesky.
        real(dp), intent(in), optional :: eps !! Symmetry tolerance forwarded to Cholesky; zero disables the check.
        integer, intent(in), optional :: nsubmax !! Optional compressed-symbolic-index ceiling forwarded to Cholesky.
        integer, intent(in), optional :: nnzlmax !! Optional numerical-factor-entry ceiling forwarded to Cholesky.
        integer, intent(in), optional :: tmpmax !! Optional real-workspace ceiling forwarded to Cholesky.
        real(dp) :: value
        type(csr_cholesky) :: factor
        integer :: stat

        factor = csr_chol(a, stat, cache_kb=cache_kb, tiny=tiny, large=large, level=level, eps=eps, &
            nsubmax=nsubmax, nnzlmax=nnzlmax, tmpmax=tmpmax)
        if (present(info)) info = stat
        if (stat < 0 .or. (stat > 0 .and. stat < 17)) then
            value = 0.0_dp
            if (present(sign)) sign = 0
            return
        end if
        if (present(sign)) sign = 1
        if (present(logarithm)) then
            if (logarithm) then
                value = factor%log_determinant
            else
                value = factor%determinant
            end if
        else
            value = factor%log_determinant
        end if
    end function csr_determinant

    function slm(x, y, info) result(fit)
        type(csr_matrix), intent(in) :: x !! Sparse model matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
        integer, intent(out), optional :: info !! Zero on success; nonzero for dimension or factorization failure.
        type(slm_result) :: fit
        fit = slm_fit_csr(x, y, info)
    end function slm

    function slm_fit(x, y, info) result(fit)
        type(csr_matrix), intent(in) :: x !! Sparse model matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
        integer, intent(out), optional :: info !! Zero on success; nonzero for dimension or factorization failure.
        type(slm_result) :: fit
        fit = slm_fit_csr(x, y, info)
    end function slm_fit

    function slm_wfit(x, y, weights, info) result(fit)
        type(csr_matrix), intent(in) :: x !! Sparse model matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
        real(dp), intent(in) :: weights(:) !! Nonnegative observation weights of length x%nrow.
        integer, intent(out), optional :: info !! Zero on success; -1 dimension mismatch; -2 negative weights; otherwise factorization status.
        type(slm_result) :: fit
        real(dp), allocatable :: xd(:,:)
        real(dp), allocatable :: yw(:)
        real(dp) :: w
        integer :: i
        type(csr_matrix) :: xw
        integer :: stat

        if (present(info)) info = 0
        if (size(y) /= x%nrow .or. size(weights) /= x%nrow) then
            if (present(info)) info = -1
            fit%df_residual = -1
            return
        end if
        if (any(weights < 0.0_dp)) then
            if (present(info)) info = -2
            fit%df_residual = -1
            return
        end if
        xd = csr_to_dense(x)
        allocate(yw(size(y)))
        do i = 1, x%nrow
            w = sqrt(weights(i))
            xd(i, :) = w * xd(i, :)
            yw(i) = w * y(i)
        end do
        xw = csr_from_dense(xd)
        fit = slm_fit_csr(xw, yw, stat)
        fit%weighted = .true.
        if (present(info)) info = stat
    end function slm_wfit

    function slm_fit_csr(x, y, info) result(fit)
        type(csr_matrix), intent(in) :: x !! Sparse model matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
        integer, intent(out), optional :: info !! Zero on success; -1 for dimension mismatch; otherwise Cholesky status.
        type(slm_result) :: fit
        type(csr_matrix) :: xt
        type(csr_matrix) :: xtx
        real(dp), allocatable :: xty(:)
        integer :: stat

        if (present(info)) info = 0
        if (size(y) /= x%nrow) then
            if (present(info)) info = -1
            fit%df_residual = -1
            return
        end if
        xt = csr_transpose(x)
        xtx = csr_matmul(xt, x, stat)
        if (stat /= 0) then
            if (present(info)) info = stat
            fit%df_residual = -1
            return
        end if
        xty = csr_matvec(xt, y, stat)
        fit%chol = csr_chol(xtx, stat)
        if (stat /= 0) then
            if (present(info)) info = stat
            fit%df_residual = x%nrow - x%ncol
            return
        end if
        fit%coefficients = chol_solve(fit%chol, xty, stat)
        fit%fitted = csr_matvec(x, fit%coefficients, stat)
        allocate(fit%residuals(size(y)))
        fit%residuals = y - fit%fitted
        fit%df_residual = x%nrow - x%ncol
        fit%weighted = .false.
        if (present(info)) info = stat
    end function slm_fit_csr

    pure function slm_coef(fit) result(values)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit whose coefficients are requested.
        real(dp), allocatable :: values(:)
        values = fit%coefficients
    end function slm_coef

    pure function slm_fitted(fit) result(values)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit whose fitted values are requested.
        real(dp), allocatable :: values(:)
        values = fit%fitted
    end function slm_fitted

    pure function slm_residuals(fit) result(values)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit whose residuals are requested.
        real(dp), allocatable :: values(:)
        values = fit%residuals
    end function slm_residuals

    pure function slm_deviance(fit) result(value)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit whose residual sum of squares is requested.
        real(dp) :: value
        value = sum(fit%residuals * fit%residuals)
    end function slm_deviance

    pure function slm_extract_aic(fit, scale, k) result(values)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit used to compute effective degrees of freedom and AIC criterion.
        real(dp), intent(in), optional :: scale !! Positive known residual scale; zero or absent uses the log-RSS form.
        real(dp), intent(in), optional :: k !! Penalty multiplier; default is 2.
        real(dp) :: values(2)
        real(dp) :: sc
        real(dp) :: penalty
        real(dp) :: rss
        real(dp) :: dev
        integer :: n
        integer :: edf

        sc = 0.0_dp
        penalty = 2.0_dp
        if (present(scale)) sc = scale
        if (present(k)) penalty = k
        n = size(fit%residuals)
        edf = n - fit%df_residual
        rss = slm_deviance(fit)
        if (sc > 0.0_dp) then
            dev = rss / sc - real(n, dp)
        else
            dev = real(n, dp) * log(rss / real(n, dp))
        end if
        values(1) = real(edf, dp)
        values(2) = dev + penalty * real(edf, dp)
    end function slm_extract_aic

    function slm_summary(fit, y, intercept, info) result(summary)
        type(slm_result), intent(in) :: fit !! Sparse linear-model fit to summarize.
        real(dp), intent(in) :: y(:) !! Original unweighted response, used to validate observation count and model sums of squares.
        logical, intent(in), optional :: intercept !! True when the fitted model contains an intercept; default is true.
        integer, intent(out), optional :: info !! Zero on success; -1 for inconsistent sizes; positive for covariance solve failure.
        type(slm_summary_result) :: summary
        real(dp), allocatable :: e(:)
        real(dp), allocatable :: col(:)
        real(dp) :: mss
        real(dp) :: fitted_mean
        real(dp) :: resvar
        real(dp) :: denom
        logical :: has_intercept
        integer :: i
        integer :: j
        integer :: p
        integer :: df_int
        integer :: stat

        if (present(info)) info = 0
        if (size(y) /= size(fit%residuals) .or. fit%df_residual <= 0) then
            if (present(info)) info = -1
            return
        end if
        p = size(fit%coefficients)
        summary%nobs = size(y)
        summary%df_model = p
        summary%df_residual = fit%df_residual
        summary%rss = slm_deviance(fit)
        resvar = summary%rss / real(summary%df_residual, dp)
        summary%sigma = sqrt(resvar)
        allocate(summary%covariance(p, p), source=0.0_dp)
        allocate(e(p), source=0.0_dp)
        do i = 1, p
            e = 0.0_dp
            e(i) = 1.0_dp
            col = chol_solve(fit%chol, e, stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                return
            end if
            summary%covariance(:, i) = col * resvar
        end do
        allocate(summary%std_error(p), summary%t_value(p), summary%p_value(p))
        allocate(summary%correlation(p, p), source=0.0_dp)
        do i = 1, p
            summary%std_error(i) = sqrt(max(0.0_dp, summary%covariance(i, i)))
            if (summary%std_error(i) > 0.0_dp) then
                summary%t_value(i) = fit%coefficients(i) / summary%std_error(i)
                summary%p_value(i) = student_t_two_sided_p(summary%t_value(i), summary%df_residual)
            else
                summary%t_value(i) = 0.0_dp
                if (fit%coefficients(i) == 0.0_dp) then
                    summary%p_value(i) = 1.0_dp
                else
                    summary%p_value(i) = 0.0_dp
                end if
            end if
        end do
        do j = 1, p
            do i = 1, p
                denom = summary%std_error(i) * summary%std_error(j)
                if (denom > 0.0_dp) summary%correlation(i, j) = summary%covariance(i, j) / denom
            end do
        end do

        has_intercept = .true.
        if (present(intercept)) has_intercept = intercept
        if (has_intercept) then
            fitted_mean = sum(fit%fitted) / real(size(fit%fitted), dp)
            mss = sum((fit%fitted - fitted_mean) * (fit%fitted - fitted_mean))
            df_int = 1
        else
            mss = sum(fit%fitted * fit%fitted)
            df_int = 0
        end if
        if (mss + summary%rss > 0.0_dp) then
            summary%r_squared = mss / (mss + summary%rss)
            summary%adj_r_squared = 1.0_dp - (1.0_dp - summary%r_squared) * &
                real(size(y) - df_int, dp) / real(summary%df_residual, dp)
        end if
        if (p > df_int .and. resvar > 0.0_dp) then
            summary%has_f_statistic = .true.
            summary%f_num_df = p - df_int
            summary%f_den_df = summary%df_residual
            summary%f_statistic = (mss / real(summary%f_num_df, dp)) / resvar
            summary%f_p_value = f_upper_tail_p(summary%f_statistic, summary%f_num_df, summary%f_den_df)
        end if
    end function slm_summary


    function mslm_fit_csr(x, y, info) result(fit)
        type(csr_matrix), intent(in) :: x !! Sparse model matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:,:) !! Dense response matrix with one observation per row and one response per column.
        integer, intent(out), optional :: info !! Zero on success; -1 for dimension mismatch; otherwise sparse-product or Cholesky status.
        type(mslm_result) :: fit
        type(csr_matrix) :: xt
        type(csr_matrix) :: xtx
        real(dp), allocatable :: xty(:)
        real(dp), allocatable :: coef(:)
        real(dp), allocatable :: fitted(:)
        integer :: stat
        integer :: j
        integer :: nresp

        if (present(info)) info = 0
        if (size(y, 1) /= x%nrow) then
            if (present(info)) info = -1
            fit%df_residual = -1
            return
        end if
        nresp = size(y, 2)
        xt = csr_transpose(x)
        xtx = csr_matmul(xt, x, stat)
        if (stat /= 0) then
            if (present(info)) info = stat
            fit%df_residual = -1
            return
        end if
        fit%chol = csr_chol(xtx, stat)
        if (stat /= 0) then
            if (present(info)) info = stat
            fit%df_residual = x%nrow - x%ncol
            return
        end if
        allocate(fit%coefficients(x%ncol, nresp))
        allocate(fit%fitted(x%nrow, nresp))
        allocate(fit%residuals(x%nrow, nresp))
        do j = 1, nresp
            xty = csr_matvec(xt, y(:, j), stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                return
            end if
            coef = chol_solve(fit%chol, xty, stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                return
            end if
            fitted = csr_matvec(x, coef, stat)
            if (stat /= 0) then
                if (present(info)) info = stat
                return
            end if
            fit%coefficients(:, j) = coef
            fit%fitted(:, j) = fitted
            fit%residuals(:, j) = y(:, j) - fitted
        end do
        fit%df_residual = x%nrow - x%ncol
        fit%weighted = .false.
        if (present(info)) info = 0
    end function mslm_fit_csr

    function mslm_summary(fit, y, intercept, info) result(summary)
        type(mslm_result), intent(in) :: fit !! Multi-response sparse linear-model fit to summarize response by response.
        real(dp), intent(in) :: y(:,:) !! Original unweighted response matrix used for per-response sums of squares.
        logical, intent(in), optional :: intercept !! True when the fitted model contains an intercept; default is true.
        integer, intent(out), optional :: info !! Zero on success; -1 for inconsistent dimensions; otherwise a per-response summary status.
        type(mslm_summary_result) :: summary
        type(slm_result) :: one
        integer :: stat
        integer :: j
        integer :: nresp

        if (present(info)) info = 0
        if (.not. allocated(fit%coefficients) .or. .not. allocated(fit%fitted) .or. .not. allocated(fit%residuals)) then
            if (present(info)) info = -1
            return
        end if
        nresp = size(fit%coefficients, 2)
        if (size(y, 1) /= size(fit%residuals, 1) .or. size(y, 2) /= nresp) then
            if (present(info)) info = -1
            return
        end if
        if (size(fit%fitted, 2) /= nresp .or. size(fit%residuals, 2) /= nresp) then
            if (present(info)) info = -1
            return
        end if
        allocate(summary%responses(nresp))
        do j = 1, nresp
            one%coefficients = fit%coefficients(:, j)
            one%fitted = fit%fitted(:, j)
            one%residuals = fit%residuals(:, j)
            one%chol = fit%chol
            one%df_residual = fit%df_residual
            one%weighted = fit%weighted
            if (present(intercept)) then
                summary%responses(j) = slm_summary(one, y(:, j), intercept, stat)
            else
                summary%responses(j) = slm_summary(one, y(:, j), info=stat)
            end if
            if (stat /= 0) then
                if (present(info)) info = stat
                return
            end if
        end do
        if (present(info)) info = 0
    end function mslm_summary


    subroutine append_cholesky_entry(lower, nnz, capacity, column, value)
        type(csr_matrix), intent(inout) :: lower !! Lower-triangular CSR factor whose current row is being assembled.
        integer, intent(inout) :: nnz !! Number of factor entries already stored; incremented on return.
        integer, intent(inout) :: capacity !! Allocated capacity of the factor value/index arrays; grown when necessary.
        integer, intent(in) :: column !! One-based column index for the new factor entry.
        real(dp), intent(in) :: value !! Numerical value of the new lower-triangular factor entry.
        real(dp), allocatable :: new_values(:)
        integer, allocatable :: new_columns(:)
        integer :: new_capacity

        if (nnz == capacity) then
            new_capacity = max(capacity + 1, 2 * capacity)
            allocate(new_values(new_capacity))
            allocate(new_columns(new_capacity))
            if (nnz > 0) then
                new_values(1:nnz) = lower%values(1:nnz)
                new_columns(1:nnz) = lower%col_ind(1:nnz)
            end if
            call move_alloc(new_values, lower%values)
            call move_alloc(new_columns, lower%col_ind)
            capacity = new_capacity
        end if
        nnz = nnz + 1
        lower%values(nnz) = value
        lower%col_ind(nnz) = column
    end subroutine append_cholesky_entry

    subroutine trim_csr_storage(matrix, nnz)
        type(csr_matrix), intent(inout) :: matrix !! CSR matrix whose overallocated value/index arrays are reduced to the used length.
        integer, intent(in) :: nnz !! Number of valid entries to retain from the start of the CSR value/index arrays.
        real(dp), allocatable :: values(:)
        integer, allocatable :: columns(:)

        allocate(values(nnz))
        allocate(columns(nnz))
        if (nnz > 0) then
            values = matrix%values(1:nnz)
            columns = matrix%col_ind(1:nnz)
        end if
        call move_alloc(values, matrix%values)
        call move_alloc(columns, matrix%col_ind)
    end subroutine trim_csr_storage

    subroutine supernodal_factor_to_csc(n, nsuper, xsuper, xlindx, lindx, xlnz, lnz, lower)
        integer, intent(in) :: n !! Number of rows and columns in the lower-triangular factor.
        integer, intent(in) :: nsuper !! Number of fundamental supernodes represented by the structure arrays.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers with active extent nsuper+1.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode row-index list.
        integer, intent(in) :: lindx(:) !! Compressed row indices for the supernodal factor.
        integer, intent(in) :: xlnz(:) !! Numerical column pointers for the lower factor.
        real(dp), intent(in) :: lnz(:) !! Numerical Cholesky values in Ng-Peyton supernodal storage.
        type(csc_matrix), intent(out) :: lower !! Equivalent lower-triangular factor in ordinary CSC storage.
        integer :: column
        integer :: first_column
        integer :: offset
        integer :: row_count
        integer :: super

        lower%nrow = n
        lower%ncol = n
        allocate(lower%values(size(lnz)), lower%row_ind(size(lnz)), lower%col_ptr(n + 1))
        lower%values = lnz
        lower%col_ptr = xlnz
        do super = 1, nsuper
            first_column = xsuper(super)
            do column = first_column, xsuper(super + 1) - 1
                offset = column - first_column
                row_count = xlnz(column + 1) - xlnz(column)
                lower%row_ind(xlnz(column):xlnz(column + 1) - 1) = &
                    lindx(xlindx(super) + offset:xlindx(super) + offset + row_count - 1)
            end do
        end do
    end subroutine supernodal_factor_to_csc

    subroutine initialize_empty_cholesky(factor)
        type(csr_cholesky), intent(inout) :: factor !! Cholesky object to initialize with empty sparse lower/upper matrices.

        factor%n = 0
        factor%ordering = "natural"
        factor%numeric_kernel = "row_sparse"
        factor%nofsub = 0
        factor%nnzl = 0
        factor%nsub = 0
        factor%nsuper = 0
        factor%tmp_size = 0
        factor%info = 0
        factor%determinant = 0.0_dp
        factor%log_determinant = 0.0_dp
        factor%factor_determinant = 0.0_dp
        factor%factor_log_determinant = 0.0_dp
        if (allocated(factor%perm)) deallocate(factor%perm)
        if (allocated(factor%inv_perm)) deallocate(factor%inv_perm)
        if (allocated(factor%etree_parent)) deallocate(factor%etree_parent)
        if (allocated(factor%col_count)) deallocate(factor%col_count)
        if (allocated(factor%supernode)) deallocate(factor%supernode)
        if (allocated(factor%xsuper)) deallocate(factor%xsuper)
        if (allocated(factor%xlindx)) deallocate(factor%xlindx)
        if (allocated(factor%lindx)) deallocate(factor%lindx)
        if (allocated(factor%xlnz)) deallocate(factor%xlnz)
        if (allocated(factor%lnz)) deallocate(factor%lnz)
        if (allocated(factor%split)) deallocate(factor%split)
        allocate(factor%perm(0), factor%inv_perm(0))
        allocate(factor%etree_parent(0), factor%col_count(0), factor%supernode(0), factor%xsuper(1))
        allocate(factor%xlindx(1), factor%lindx(0), factor%xlnz(1))
        allocate(factor%lnz(0), factor%split(0))
        factor%xsuper(1) = 1
        factor%xlindx(1) = 1
        factor%xlnz(1) = 1
        factor%lower%nrow = 0
        factor%lower%ncol = 0
        factor%upper%nrow = 0
        factor%upper%ncol = 0
        allocate(factor%lower%values(0), factor%lower%col_ind(0), factor%lower%row_ptr(1))
        allocate(factor%upper%values(0), factor%upper%col_ind(0), factor%upper%row_ptr(1))
        factor%lower%row_ptr(1) = 1
        factor%upper%row_ptr(1) = 1
    end subroutine initialize_empty_cholesky

    pure real(dp) function student_t_two_sided_p(t_value, degrees_freedom) result(probability)
        real(dp), intent(in) :: t_value !! Student-t statistic whose two-sided tail probability is requested.
        integer, intent(in) :: degrees_freedom !! Positive residual degrees of freedom of the Student-t distribution.
        real(dp) :: x

        if (degrees_freedom <= 0) then
            probability = 1.0_dp
            return
        end if
        x = real(degrees_freedom, dp) / (real(degrees_freedom, dp) + t_value * t_value)
        probability = regularized_beta(x, 0.5_dp * real(degrees_freedom, dp), 0.5_dp)
        probability = min(1.0_dp, max(0.0_dp, probability))
    end function student_t_two_sided_p

    pure real(dp) function f_upper_tail_p(f_value, numerator_df, denominator_df) result(probability)
        real(dp), intent(in) :: f_value !! Nonnegative F statistic whose upper-tail probability is requested.
        integer, intent(in) :: numerator_df !! Positive numerator degrees of freedom of the F distribution.
        integer, intent(in) :: denominator_df !! Positive denominator degrees of freedom of the F distribution.
        real(dp) :: x

        if (f_value < 0.0_dp .or. numerator_df <= 0 .or. denominator_df <= 0) then
            probability = 1.0_dp
            return
        end if
        x = real(denominator_df, dp) / &
            (real(denominator_df, dp) + real(numerator_df, dp) * f_value)
        probability = regularized_beta(x, 0.5_dp * real(denominator_df, dp), &
            0.5_dp * real(numerator_df, dp))
        probability = min(1.0_dp, max(0.0_dp, probability))
    end function f_upper_tail_p

    pure real(dp) function regularized_beta(x, a, b) result(value)
        real(dp), intent(in) :: x !! Beta-distribution argument in the closed interval from zero to one.
        real(dp), intent(in) :: a !! Positive first beta shape parameter.
        real(dp), intent(in) :: b !! Positive second beta shape parameter.
        real(dp) :: prefactor

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (x >= 1.0_dp) then
            value = 1.0_dp
            return
        end if
        if (a <= 0.0_dp .or. b <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        prefactor = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
            a * log(x) + b * log(1.0_dp - x))
        if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
            value = prefactor * beta_continued_fraction(a, b, x) / a
        else
            value = 1.0_dp - prefactor * beta_continued_fraction(b, a, 1.0_dp - x) / b
        end if
    end function regularized_beta

    pure real(dp) function beta_continued_fraction(a, b, x) result(value)
        real(dp), intent(in) :: a !! Positive first beta shape parameter.
        real(dp), intent(in) :: b !! Positive second beta shape parameter.
        real(dp), intent(in) :: x !! Beta-distribution argument strictly between zero and one.
        integer, parameter :: max_iterations = 300
        real(dp), parameter :: tolerance = 10.0_dp * epsilon(1.0_dp)
        real(dp), parameter :: tiny_value = tiny(1.0_dp) / epsilon(1.0_dp)
        real(dp) :: qab
        real(dp) :: qap
        real(dp) :: qam
        real(dp) :: c
        real(dp) :: d
        real(dp) :: h
        real(dp) :: aa
        real(dp) :: delta
        integer :: m
        integer :: m2

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab * x / qap
        if (abs(d) < tiny_value) d = tiny_value
        d = 1.0_dp / d
        h = d
        do m = 1, max_iterations
            m2 = 2 * m
            aa = real(m, dp) * (b - real(m, dp)) * x / &
                ((qam + real(m2, dp)) * (a + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < tiny_value) d = tiny_value
            c = 1.0_dp + aa / c
            if (abs(c) < tiny_value) c = tiny_value
            d = 1.0_dp / d
            h = h * d * c

            aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
                ((a + real(m2, dp)) * (qap + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < tiny_value) d = tiny_value
            c = 1.0_dp + aa / c
            if (abs(c) < tiny_value) c = tiny_value
            d = 1.0_dp / d
            delta = d * c
            h = h * delta
            if (abs(delta - 1.0_dp) <= tolerance) exit
        end do
        value = h
    end function beta_continued_fraction

end module sparsem_factor
