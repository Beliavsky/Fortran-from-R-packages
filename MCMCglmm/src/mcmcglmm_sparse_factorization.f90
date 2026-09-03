! SPDX-License-Identifier: GPL-2.0-or-later
! Independent sparse factorization for translated MCMCglmm engines; see NOTICE.md.
module mcmcglmm_sparse_factorization
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_normal, rng_state
    use mcmcglmm_sparse, only : mcmcglmm_sparse_matrix, sparse_from_coo, sparse_validate
    implicit none
    private

    type :: sparse_factor_row
        integer, allocatable :: column(:)
    end type sparse_factor_row

    type, public :: sparse_cholesky_analysis
        !! Reusable lower-factor sparsity pattern for one symmetric matrix graph.
        integer :: order = 0
        integer, allocatable :: row_pointer(:)
        integer, allocatable :: column(:)
    end type sparse_cholesky_analysis

    type, public :: sparse_precision_cache
        !! Reusable ordering and symbolic Cholesky analysis for repeated precision draws.
        logical :: initialized = .false.
        logical :: reordered = .false.
        integer :: analysis_count = 0
        integer, allocatable :: permutation(:)
        type(sparse_cholesky_analysis) :: analysis
    end type sparse_precision_cache

    public :: sample_mvn_sparse_precision
    public :: sparse_cholesky_factor
    public :: sparse_cholesky_analyze
    public :: sparse_cholesky_factor_analyzed
    public :: sparse_cholesky_solve
    public :: sparse_cholesky_transpose_solve
    public :: sparse_reverse_cuthill_mckee
    public :: sparse_symmetric_permute

contains

    pure subroutine sparse_cholesky_factor(matrix, factor, info, symmetry_tolerance)
        !! Compute a lower CSR Cholesky factor using a natural-order left-looking algorithm.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Symmetric positive-definite matrix in canonical CSR.
        type(mcmcglmm_sparse_matrix), intent(out) :: factor !! Lower triangular factor L satisfying A=L*transpose(L).
        integer, intent(out) :: info !! Zero on success; nonzero for invalid, nonsymmetric, or non-SPD input.
        real(dp), optional, intent(in) :: symmetry_tolerance !! Nonnegative absolute symmetry tolerance.
        type(sparse_cholesky_analysis) :: analysis

        call sparse_cholesky_analyze(matrix, analysis, info)
        if (info /= 0) return
        call sparse_cholesky_factor_analyzed(matrix, analysis, factor, info, symmetry_tolerance)
    end subroutine sparse_cholesky_factor

    pure subroutine sparse_cholesky_analyze(matrix, analysis, info)
        !! Compute and store the natural-order symbolic Cholesky fill pattern.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Square CSR matrix whose graph is analyzed.
        type(sparse_cholesky_analysis), intent(out) :: analysis !! Reusable lower-factor sparsity pattern.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid matrix storage or shape.
        type(sparse_factor_row), allocatable :: row(:)
        logical, allocatable :: pattern(:)
        integer :: count_in_row
        integer :: i
        integer :: j
        integer :: k
        integer :: total_nonzero

        info = 0
        call sparse_validate(matrix, info)
        if (info /= 0 .or. matrix%nrow /= matrix%ncol) then
            info = 1
            return
        end if
        allocate(row(matrix%nrow), pattern(matrix%nrow))
        do i = 1, matrix%nrow
            pattern = .false.
            pattern(i) = .true.
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                j = matrix%column(k)
                if (j <= i) pattern(j) = .true.
            end do
            do j = 1, i - 1
                if (pattern(j)) cycle
                do k = 1, size(row(j)%column) - 1
                    if (pattern(row(j)%column(k))) then
                        pattern(j) = .true.
                        exit
                    end if
                end do
            end do
            count_in_row = count(pattern(:i))
            allocate(row(i)%column(count_in_row))
            count_in_row = 0
            do j = 1, i
                if (.not. pattern(j)) cycle
                count_in_row = count_in_row + 1
                row(i)%column(count_in_row) = j
            end do
        end do
        analysis%order = matrix%nrow
        total_nonzero = 0
        do i = 1, matrix%nrow
            total_nonzero = total_nonzero + size(row(i)%column)
        end do
        allocate(analysis%row_pointer(matrix%nrow + 1), analysis%column(total_nonzero))
        analysis%row_pointer(1) = 1
        k = 1
        do i = 1, matrix%nrow
            count_in_row = size(row(i)%column)
            analysis%column(k:k + count_in_row - 1) = row(i)%column
            k = k + count_in_row
            analysis%row_pointer(i + 1) = k
        end do
    end subroutine sparse_cholesky_analyze

    pure subroutine sparse_cholesky_factor_analyzed(matrix, analysis, factor, info, symmetry_tolerance)
        !! Numerically factor a matrix using a compatible reusable symbolic Cholesky pattern.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Symmetric positive-definite matrix.
        type(sparse_cholesky_analysis), intent(in) :: analysis !! Compatible lower-factor sparsity pattern.
        type(mcmcglmm_sparse_matrix), intent(out) :: factor !! Numerical lower CSR Cholesky factor.
        integer, intent(out) :: info !! Zero on success; five when the matrix graph exceeds the analysis.
        real(dp), optional, intent(in) :: symmetry_tolerance !! Nonnegative absolute symmetry tolerance.
        integer :: diagonal_position
        integer :: i
        integer :: j
        integer :: k
        integer :: left_position
        integer :: right_position
        real(dp) :: diagonal
        real(dp) :: entry
        real(dp) :: tolerance

        info = 0
        tolerance = 100.0_dp * epsilon(1.0_dp)
        if (present(symmetry_tolerance)) tolerance = symmetry_tolerance
        call sparse_validate(matrix, info)
        if (info /= 0 .or. matrix%nrow /= matrix%ncol .or. tolerance < 0.0_dp .or. &
            analysis%order /= matrix%nrow .or. .not. allocated(analysis%row_pointer) .or. &
            .not. allocated(analysis%column)) then
            info = 1
            return
        end if
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                j = matrix%column(k)
                if (abs(matrix%value(k) - sparse_entry(matrix, j, i)) > tolerance) then
                    info = 2
                    return
                end if
                if (j <= i .and. .not. analysis_contains(analysis, i, j)) then
                    info = 5
                    return
                end if
            end do
        end do
        factor%nrow = matrix%nrow
        factor%ncol = matrix%ncol
        factor%row_pointer = analysis%row_pointer
        factor%column = analysis%column
        allocate(factor%value(size(analysis%column)), source=0.0_dp)
        do i = 1, matrix%nrow
            diagonal_position = factor%row_pointer(i + 1) - 1
            do k = factor%row_pointer(i), diagonal_position - 1
                j = factor%column(k)
                entry = sparse_entry(matrix, i, j)
                left_position = factor%row_pointer(i)
                right_position = factor%row_pointer(j)
                do while (left_position < k .and. right_position < factor%row_pointer(j + 1) - 1)
                    if (factor%column(left_position) == factor%column(right_position)) then
                        entry = entry - factor%value(left_position) * factor%value(right_position)
                        left_position = left_position + 1
                        right_position = right_position + 1
                    else if (factor%column(left_position) < factor%column(right_position)) then
                        left_position = left_position + 1
                    else
                        right_position = right_position + 1
                    end if
                end do
                diagonal = factor%value(factor%row_pointer(j + 1) - 1)
                if (.not. ieee_is_finite(diagonal) .or. diagonal <= 0.0_dp) then
                    info = 3
                    return
                end if
                factor%value(k) = entry / diagonal
            end do
            entry = sparse_entry(matrix, i, i)
            if (diagonal_position > factor%row_pointer(i)) then
                entry = entry - sum(factor%value(factor%row_pointer(i):diagonal_position - 1)**2)
            end if
            if (.not. ieee_is_finite(entry) .or. entry <= 0.0_dp) then
                info = 3
                return
            end if
            factor%value(diagonal_position) = sqrt(entry)
        end do
    end subroutine sparse_cholesky_factor_analyzed

    pure subroutine sparse_cholesky_solve(factor, rhs, solution, info)
        !! Solve L*transpose(L)*x=rhs from a canonical lower CSR Cholesky factor.
        type(mcmcglmm_sparse_matrix), intent(in) :: factor !! Lower CSR Cholesky factor.
        real(dp), intent(in) :: rhs(:) !! Right-hand-side vector.
        real(dp), allocatable, intent(out) :: solution(:) !! Allocated solution vector.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid factor or shape.
        real(dp), allocatable :: intermediate(:)

        call sparse_forward_solve(factor, rhs, intermediate, info)
        if (info /= 0) then
            allocate(solution(0))
            return
        end if
        call sparse_cholesky_transpose_solve(factor, intermediate, solution, info)
    end subroutine sparse_cholesky_solve

    pure subroutine sparse_cholesky_transpose_solve(factor, rhs, solution, info)
        !! Solve transpose(L)*x=rhs from a canonical lower CSR Cholesky factor.
        type(mcmcglmm_sparse_matrix), intent(in) :: factor !! Lower CSR Cholesky factor.
        real(dp), intent(in) :: rhs(:) !! Right-hand-side vector.
        real(dp), allocatable, intent(out) :: solution(:) !! Allocated solution vector.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid factor or shape.
        integer :: i
        integer :: k
        real(dp) :: diagonal

        call validate_lower_factor(factor, info)
        if (info /= 0 .or. size(rhs) /= factor%nrow) then
            allocate(solution(0))
            if (info == 0) info = 4
            return
        end if
        solution = rhs
        do i = factor%nrow, 1, -1
            diagonal = factor%value(factor%row_pointer(i + 1) - 1)
            solution(i) = solution(i) / diagonal
            do k = factor%row_pointer(i), factor%row_pointer(i + 1) - 2
                solution(factor%column(k)) = solution(factor%column(k)) - factor%value(k) * solution(i)
            end do
        end do
    end subroutine sparse_cholesky_transpose_solve

    pure subroutine sample_mvn_sparse_precision(state, rhs, precision, sample, mean_value, info, use_rcm, cache)
        !! Sample N(P^{-1}rhs,P^{-1}) using a sparse Cholesky factor of precision P.
        type(rng_state), intent(inout) :: state !! Generator state consumed by standard-normal innovations.
        real(dp), intent(in) :: rhs(:) !! Precision-weighted mean vector.
        type(mcmcglmm_sparse_matrix), intent(in) :: precision !! Symmetric positive-definite CSR precision matrix.
        real(dp), allocatable, intent(out) :: sample(:) !! Allocated Gaussian draw.
        real(dp), allocatable, intent(out) :: mean_value(:) !! Allocated precision solve.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage or factorization failure.
        logical, optional, intent(in) :: use_rcm !! Apply deterministic reverse Cuthill-McKee ordering when true.
        type(sparse_precision_cache), optional, intent(inout) :: cache !! Optional reusable ordering and analysis.
        type(mcmcglmm_sparse_matrix) :: factor
        type(mcmcglmm_sparse_matrix) :: ordered_precision
        real(dp), allocatable :: innovation(:)
        real(dp), allocatable :: ordered_mean(:)
        real(dp), allocatable :: ordered_rhs(:)
        real(dp), allocatable :: ordered_sample(:)
        real(dp), allocatable :: normal_draw(:)
        integer, allocatable :: permutation(:)
        integer :: i
        logical :: reorder

        call sparse_validate(precision, info)
        if (info /= 0 .or. precision%nrow /= precision%ncol .or. size(rhs) /= precision%nrow) then
            allocate(sample(0), mean_value(0))
            info = 1
            return
        end if
        reorder = .false.
        if (present(use_rcm)) reorder = use_rcm
        if (present(cache)) then
            call prepare_ordered_factor(precision, reorder, cache, permutation, ordered_precision, factor, info)
        else
            if (reorder) then
                call sparse_reverse_cuthill_mckee(precision, permutation, info)
                if (info == 0) call sparse_symmetric_permute(precision, permutation, ordered_precision, info)
            else
                ordered_precision = precision
            end if
            if (info == 0) call sparse_cholesky_factor(ordered_precision, factor, info)
        end if
        if (info /= 0) then
            allocate(sample(0), mean_value(0))
            return
        end if
        allocate(ordered_rhs(size(rhs)))
        if (reorder) then
            ordered_rhs = rhs(permutation)
        else
            ordered_rhs = rhs
        end if
        call sparse_cholesky_solve(factor, ordered_rhs, ordered_mean, info)
        if (info /= 0) then
            allocate(sample(0))
            return
        end if
        allocate(normal_draw(size(rhs)))
        do i = 1, size(rhs)
            call rng_normal(state, normal_draw(i))
        end do
        call sparse_cholesky_transpose_solve(factor, normal_draw, innovation, info)
        if (info /= 0) then
            allocate(sample(0))
            return
        end if
        ordered_sample = ordered_mean + innovation
        if (reorder) then
            allocate(sample(size(rhs)), mean_value(size(rhs)))
            do i = 1, size(rhs)
                sample(permutation(i)) = ordered_sample(i)
                mean_value(permutation(i)) = ordered_mean(i)
            end do
        else
            sample = ordered_sample
            mean_value = ordered_mean
        end if
    end subroutine sample_mvn_sparse_precision

    pure subroutine prepare_ordered_factor(precision, reorder, cache, permutation, ordered_precision, factor, info)
        !! Reuse or safely rebuild ordering and symbolic analysis for one numerical precision matrix.
        type(mcmcglmm_sparse_matrix), intent(in) :: precision !! Current symmetric positive-definite precision.
        logical, intent(in) :: reorder !! Whether reverse Cuthill-McKee ordering is requested.
        type(sparse_precision_cache), intent(inout) :: cache !! Persistent ordering and symbolic analysis cache.
        integer, allocatable, intent(out) :: permutation(:) !! New-to-old ordering, empty for natural order.
        type(mcmcglmm_sparse_matrix), intent(out) :: ordered_precision !! Current precision after cached ordering.
        type(mcmcglmm_sparse_matrix), intent(out) :: factor !! Numerical lower Cholesky factor.
        integer, intent(out) :: info !! Zero on success; numerical factorization status otherwise.
        logical :: cache_usable

        cache_usable = cache%initialized .and. cache%reordered .eqv. reorder
        if (cache_usable) then
            if (cache%analysis%order /= precision%nrow) cache_usable = .false.
            if (reorder) then
                if (.not. allocated(cache%permutation)) then
                    cache_usable = .false.
                else if (size(cache%permutation) /= precision%nrow) then
                    cache_usable = .false.
                end if
            end if
        end if
        if (cache_usable) then
            if (reorder) then
                permutation = cache%permutation
                call sparse_symmetric_permute(precision, permutation, ordered_precision, info)
            else
                allocate(permutation(0))
                ordered_precision = precision
                info = 0
            end if
            if (info == 0) then
                call sparse_cholesky_factor_analyzed(ordered_precision, cache%analysis, factor, info)
                if (info /= 5) return
            end if
        end if

        if (reorder) then
            call sparse_reverse_cuthill_mckee(precision, permutation, info)
            if (info /= 0) return
            call sparse_symmetric_permute(precision, permutation, ordered_precision, info)
            if (info /= 0) return
        else
            if (allocated(permutation)) deallocate(permutation)
            allocate(permutation(0))
            ordered_precision = precision
        end if
        call sparse_cholesky_analyze(ordered_precision, cache%analysis, info)
        if (info /= 0) return
        cache%analysis_count = cache%analysis_count + 1
        call sparse_cholesky_factor_analyzed(ordered_precision, cache%analysis, factor, info)
        if (info /= 0) return
        cache%initialized = .true.
        cache%reordered = reorder
        if (reorder) then
            cache%permutation = permutation
        else if (allocated(cache%permutation)) then
            deallocate(cache%permutation)
        end if
    end subroutine prepare_ordered_factor

    pure subroutine sparse_reverse_cuthill_mckee(matrix, permutation, info)
        !! Compute a deterministic reverse Cuthill-McKee ordering of a symmetric CSR graph.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Symmetric square CSR matrix.
        integer, allocatable, intent(out) :: permutation(:) !! New-to-old index mapping.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid or nonsymmetric input.
        integer, allocatable :: degree(:)
        integer, allocatable :: neighbor(:)
        integer, allocatable :: order(:)
        integer, allocatable :: queue(:)
        logical, allocatable :: visited(:)
        integer :: candidate
        integer :: component_end
        integer :: component_start
        integer :: current
        integer :: i
        integer :: j
        integer :: k
        integer :: neighbor_count
        integer :: order_count
        integer :: queue_end
        integer :: queue_start
        integer :: start

        call sparse_validate(matrix, info)
        if (info /= 0 .or. matrix%nrow /= matrix%ncol) then
            allocate(permutation(0))
            info = 1
            return
        end if
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                j = matrix%column(k)
                if (abs(matrix%value(k) - sparse_entry(matrix, j, i)) > &
                    100.0_dp * epsilon(1.0_dp)) then
                    allocate(permutation(0))
                    info = 2
                    return
                end if
            end do
        end do
        allocate(degree(matrix%nrow), neighbor(matrix%nrow), order(matrix%nrow), queue(matrix%nrow), source=0)
        allocate(visited(matrix%nrow), source=.false.)
        visited = .false.
        do i = 1, matrix%nrow
            degree(i) = count(matrix%column(matrix%row_pointer(i):matrix%row_pointer(i + 1) - 1) /= i)
        end do
        order_count = 0
        do while (order_count < matrix%nrow)
            start = 0
            do i = 1, matrix%nrow
                if (visited(i)) cycle
                if (start == 0) then
                    start = i
                else if (degree(i) < degree(start)) then
                    start = i
                end if
            end do
            component_start = order_count + 1
            queue_start = 1
            queue_end = 1
            queue(1) = start
            visited(start) = .true.
            do while (queue_start <= queue_end)
                current = queue(queue_start)
                queue_start = queue_start + 1
                order_count = order_count + 1
                order(order_count) = current
                neighbor_count = 0
                do k = matrix%row_pointer(current), matrix%row_pointer(current + 1) - 1
                    candidate = matrix%column(k)
                    if (candidate == current .or. visited(candidate)) cycle
                    visited(candidate) = .true.
                    neighbor_count = neighbor_count + 1
                    neighbor(neighbor_count) = candidate
                end do
                call sort_nodes_by_degree(neighbor(:neighbor_count), degree)
                do i = 1, neighbor_count
                    queue_end = queue_end + 1
                    queue(queue_end) = neighbor(i)
                end do
            end do
            component_end = order_count
            do i = 0, (component_end - component_start) / 2
                j = order(component_start + i)
                order(component_start + i) = order(component_end - i)
                order(component_end - i) = j
            end do
        end do
        permutation = order
    end subroutine sparse_reverse_cuthill_mckee

    pure subroutine sparse_symmetric_permute(matrix, permutation, permuted, info)
        !! Apply new-to-old symmetric permutation B(i,j)=A(permutation(i),permutation(j)).
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Square CSR matrix to permute.
        integer, intent(in) :: permutation(:) !! New-to-old permutation of one through matrix order.
        type(mcmcglmm_sparse_matrix), intent(out) :: permuted !! Canonical symmetrically permuted CSR matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid matrix or permutation.
        integer, allocatable :: coordinate_column(:)
        integer, allocatable :: coordinate_row(:)
        real(dp), allocatable :: coordinate_value(:)
        integer, allocatable :: inverse_permutation(:)
        logical, allocatable :: seen(:)
        integer :: i
        integer :: k
        integer :: old_column
        integer :: old_row

        call sparse_validate(matrix, info)
        if (info /= 0 .or. matrix%nrow /= matrix%ncol .or. size(permutation) /= matrix%nrow) then
            info = 1
            return
        end if
        allocate(inverse_permutation(matrix%nrow), source=0)
        allocate(seen(matrix%nrow), source=.false.)
        do i = 1, matrix%nrow
            if (permutation(i) < 1 .or. permutation(i) > matrix%nrow) then
                info = 2
                return
            end if
            if (seen(permutation(i))) then
                info = 2
                return
            end if
            seen(permutation(i)) = .true.
            inverse_permutation(permutation(i)) = i
        end do
        allocate(coordinate_row(size(matrix%value)), coordinate_column(size(matrix%value)), &
            coordinate_value(size(matrix%value)))
        k = 0
        do old_row = 1, matrix%nrow
            do i = matrix%row_pointer(old_row), matrix%row_pointer(old_row + 1) - 1
                old_column = matrix%column(i)
                k = k + 1
                coordinate_row(k) = inverse_permutation(old_row)
                coordinate_column(k) = inverse_permutation(old_column)
                coordinate_value(k) = matrix%value(i)
            end do
        end do
        call sparse_from_coo(matrix%nrow, matrix%ncol, coordinate_row, coordinate_column, coordinate_value, &
            permuted, info)
    end subroutine sparse_symmetric_permute

    pure subroutine sort_nodes_by_degree(nodes, degree)
        !! Sort a small node list by increasing graph degree and then node number.
        integer, intent(inout) :: nodes(:) !! Node identifiers to sort in place.
        integer, intent(in) :: degree(:) !! Graph degree indexed by node identifier.
        integer :: i
        integer :: j
        integer :: value

        do i = 2, size(nodes)
            value = nodes(i)
            j = i - 1
            do while (j >= 1)
                if (degree(nodes(j)) < degree(value)) exit
                if (degree(nodes(j)) == degree(value) .and. nodes(j) < value) exit
                nodes(j + 1) = nodes(j)
                j = j - 1
            end do
            nodes(j + 1) = value
        end do
    end subroutine sort_nodes_by_degree

    pure subroutine sparse_forward_solve(factor, rhs, solution, info)
        !! Solve L*x=rhs from a canonical lower CSR Cholesky factor.
        type(mcmcglmm_sparse_matrix), intent(in) :: factor !! Lower CSR Cholesky factor.
        real(dp), intent(in) :: rhs(:) !! Right-hand-side vector.
        real(dp), allocatable, intent(out) :: solution(:) !! Allocated solution vector.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid factor or shape.
        integer :: i
        integer :: k
        real(dp) :: accumulated
        real(dp) :: diagonal

        call validate_lower_factor(factor, info)
        if (info /= 0 .or. size(rhs) /= factor%nrow) then
            allocate(solution(0))
            if (info == 0) info = 4
            return
        end if
        allocate(solution(factor%nrow))
        do i = 1, factor%nrow
            accumulated = rhs(i)
            do k = factor%row_pointer(i), factor%row_pointer(i + 1) - 2
                accumulated = accumulated - factor%value(k) * solution(factor%column(k))
            end do
            diagonal = factor%value(factor%row_pointer(i + 1) - 1)
            solution(i) = accumulated / diagonal
        end do
    end subroutine sparse_forward_solve

    pure subroutine validate_lower_factor(factor, info)
        !! Validate canonical CSR storage with one positive diagonal as the last entry of every row.
        type(mcmcglmm_sparse_matrix), intent(in) :: factor !! Candidate lower triangular Cholesky factor.
        integer, intent(out) :: info !! Zero when the factor structure is valid.
        integer :: i

        call sparse_validate(factor, info)
        if (info /= 0 .or. factor%nrow /= factor%ncol) then
            info = 1
            return
        end if
        do i = 1, factor%nrow
            if (factor%row_pointer(i) == factor%row_pointer(i + 1) .or. &
                factor%column(factor%row_pointer(i + 1) - 1) /= i .or. &
                factor%value(factor%row_pointer(i + 1) - 1) <= 0.0_dp) then
                info = 1
                return
            end if
        end do
    end subroutine validate_lower_factor

    pure logical function analysis_contains(analysis, row, column) result(found)
        !! Report whether one lower-triangular coordinate is present in a symbolic pattern.
        type(sparse_cholesky_analysis), intent(in) :: analysis !! Symbolic lower-factor pattern.
        integer, intent(in) :: row !! One-based factor row.
        integer, intent(in) :: column !! One-based factor column.
        integer :: k

        found = .false.
        do k = analysis%row_pointer(row), analysis%row_pointer(row + 1) - 1
            if (analysis%column(k) == column) then
                found = .true.
                return
            end if
            if (analysis%column(k) > column) return
        end do
    end function analysis_contains

    pure real(dp) function sparse_entry(matrix, row, column) result(value)
        !! Return one CSR entry or zero when the coordinate is structurally absent.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Canonical CSR matrix.
        integer, intent(in) :: row !! One-based row coordinate.
        integer, intent(in) :: column !! One-based column coordinate.
        integer :: k

        value = 0.0_dp
        do k = matrix%row_pointer(row), matrix%row_pointer(row + 1) - 1
            if (matrix%column(k) == column) then
                value = matrix%value(k)
                return
            end if
            if (matrix%column(k) > column) return
        end do
    end function sparse_entry

end module mcmcglmm_sparse_factorization
