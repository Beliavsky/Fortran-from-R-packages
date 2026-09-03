! SPDX-License-Identifier: GPL-2.0-or-later
! Sparse numerical foundation for translated MCMCglmm engines; see NOTICE.md.
module mcmcglmm_sparse
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use r_kinds, only : dp
    implicit none
    private

    type, public :: mcmcglmm_sparse_matrix
        !! One-based compressed sparse row matrix with sorted, unique column indices per row.
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: row_pointer(:)
        integer, allocatable :: column(:)
        real(dp), allocatable :: value(:)
    end type mcmcglmm_sparse_matrix

    type :: sparse_product_row
        integer, allocatable :: column(:)
        real(dp), allocatable :: value(:)
    end type sparse_product_row

    public :: sparse_crossproduct
    public :: sparse_from_coo
    public :: sparse_from_dense
    public :: sparse_is_initialized
    public :: sparse_matmul_matrix
    public :: sparse_matmul_vector
    public :: sparse_stacked_crossproduct
    public :: sparse_to_dense
    public :: sparse_transpose_matmul_matrix
    public :: sparse_transpose_matmul_sparse
    public :: sparse_validate

contains

    pure logical function sparse_is_initialized(matrix) result(initialized)
        !! Report whether a sparse matrix has positive dimensions and allocated CSR storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix being queried.

        initialized = matrix%nrow > 0 .and. matrix%ncol > 0 .and. allocated(matrix%row_pointer) .and. &
            allocated(matrix%column) .and. allocated(matrix%value)
    end function sparse_is_initialized

    pure subroutine sparse_validate(matrix, info)
        !! Validate CSR dimensions, pointers, sorted columns, and finite values.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix to validate.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage.
        integer :: i
        integer :: k

        info = 0
        if (.not. sparse_is_initialized(matrix)) then
            info = 1
            return
        end if
        if (size(matrix%row_pointer) /= matrix%nrow + 1 .or. &
            size(matrix%column) /= size(matrix%value)) then
            info = 1
            return
        end if
        if (matrix%row_pointer(1) /= 1 .or. &
            matrix%row_pointer(matrix%nrow + 1) /= size(matrix%value) + 1) then
            info = 2
            return
        end if
        if (any(matrix%row_pointer(2:) < matrix%row_pointer(:matrix%nrow)) .or. &
            any(matrix%column < 1) .or. any(matrix%column > matrix%ncol) .or. &
            any(.not. ieee_is_finite(matrix%value))) then
            info = 3
            return
        end if
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i) + 1, matrix%row_pointer(i + 1) - 1
                if (matrix%column(k) <= matrix%column(k - 1)) then
                    info = 4
                    return
                end if
            end do
        end do
    end subroutine sparse_validate

    pure subroutine sparse_from_dense(dense, matrix, tolerance, info)
        !! Compress a dense matrix, dropping entries whose magnitudes do not exceed tolerance.
        real(dp), intent(in) :: dense(:, :) !! Dense input matrix.
        type(mcmcglmm_sparse_matrix), intent(out) :: matrix !! Canonical CSR result.
        real(dp), optional, intent(in) :: tolerance !! Nonnegative drop tolerance; defaults to zero.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or tolerance.
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz
        real(dp) :: threshold

        info = 0
        threshold = 0.0_dp
        if (present(tolerance)) threshold = tolerance
        if (size(dense, 1) < 1 .or. size(dense, 2) < 1 .or. threshold < 0.0_dp .or. &
            any(.not. ieee_is_finite(dense))) then
            info = 1
            return
        end if
        matrix%nrow = size(dense, 1)
        matrix%ncol = size(dense, 2)
        nnz = count(abs(dense) > threshold)
        allocate(matrix%row_pointer(matrix%nrow + 1), matrix%column(nnz), matrix%value(nnz))
        k = 1
        matrix%row_pointer(1) = 1
        do i = 1, matrix%nrow
            do j = 1, matrix%ncol
                if (abs(dense(i, j)) <= threshold) cycle
                matrix%column(k) = j
                matrix%value(k) = dense(i, j)
                k = k + 1
            end do
            matrix%row_pointer(i + 1) = k
        end do
    end subroutine sparse_from_dense

    pure subroutine sparse_from_coo(nrow, ncol, row, column, value, matrix, info, tolerance)
        !! Build canonical CSR storage from one-based coordinate triplets, summing duplicate entries.
        integer, intent(in) :: nrow !! Number of matrix rows.
        integer, intent(in) :: ncol !! Number of matrix columns.
        integer, intent(in) :: row(:) !! One-based row index for each coordinate entry.
        integer, intent(in) :: column(:) !! One-based column index for each coordinate entry.
        real(dp), intent(in) :: value(:) !! Coordinate values; duplicates are summed.
        type(mcmcglmm_sparse_matrix), intent(out) :: matrix !! Canonical CSR result.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid triplets or tolerance.
        real(dp), optional, intent(in) :: tolerance !! Nonnegative post-summation drop tolerance.
        integer, allocatable :: count_by_row(:)
        integer, allocatable :: cursor(:)
        integer, allocatable :: raw_column(:)
        integer, allocatable :: raw_pointer(:)
        real(dp), allocatable :: raw_value(:)
        integer :: c
        integer :: i
        integer :: j
        integer :: k
        integer :: keep
        integer :: next_column
        integer :: start
        integer :: stop
        integer :: swap_column
        real(dp) :: accumulated
        real(dp) :: swap_value
        real(dp) :: threshold

        info = 0
        threshold = 0.0_dp
        if (present(tolerance)) threshold = tolerance
        if (nrow < 1 .or. ncol < 1 .or. size(row) /= size(column) .or. size(row) /= size(value) .or. &
            threshold < 0.0_dp .or. any(row < 1) .or. any(row > nrow) .or. any(column < 1) .or. &
            any(column > ncol) .or. any(.not. ieee_is_finite(value))) then
            info = 1
            return
        end if

        allocate(count_by_row(nrow), source=0)
        do k = 1, size(row)
            if (abs(value(k)) > threshold) count_by_row(row(k)) = count_by_row(row(k)) + 1
        end do
        allocate(raw_pointer(nrow + 1))
        raw_pointer(1) = 1
        do i = 1, nrow
            raw_pointer(i + 1) = raw_pointer(i) + count_by_row(i)
        end do
        allocate(raw_column(raw_pointer(nrow + 1) - 1), raw_value(raw_pointer(nrow + 1) - 1))
        allocate(cursor(nrow))
        cursor = raw_pointer(:nrow)
        do k = 1, size(row)
            if (abs(value(k)) <= threshold) cycle
            i = row(k)
            raw_column(cursor(i)) = column(k)
            raw_value(cursor(i)) = value(k)
            cursor(i) = cursor(i) + 1
        end do

        do i = 1, nrow
            start = raw_pointer(i)
            stop = raw_pointer(i + 1) - 1
            do k = start + 1, stop
                swap_column = raw_column(k)
                swap_value = raw_value(k)
                j = k - 1
                do while (j >= start)
                    if (raw_column(j) <= swap_column) exit
                    raw_column(j + 1) = raw_column(j)
                    raw_value(j + 1) = raw_value(j)
                    j = j - 1
                end do
                raw_column(j + 1) = swap_column
                raw_value(j + 1) = swap_value
            end do
        end do

        count_by_row = 0
        do i = 1, nrow
            k = raw_pointer(i)
            stop = raw_pointer(i + 1) - 1
            do while (k <= stop)
                c = raw_column(k)
                accumulated = 0.0_dp
                do while (k <= stop)
                    if (raw_column(k) /= c) exit
                    accumulated = accumulated + raw_value(k)
                    k = k + 1
                end do
                if (abs(accumulated) > threshold) count_by_row(i) = count_by_row(i) + 1
            end do
        end do
        matrix%nrow = nrow
        matrix%ncol = ncol
        allocate(matrix%row_pointer(nrow + 1))
        matrix%row_pointer(1) = 1
        do i = 1, nrow
            matrix%row_pointer(i + 1) = matrix%row_pointer(i) + count_by_row(i)
        end do
        allocate(matrix%column(matrix%row_pointer(nrow + 1) - 1))
        allocate(matrix%value(matrix%row_pointer(nrow + 1) - 1))
        keep = 1
        do i = 1, nrow
            k = raw_pointer(i)
            stop = raw_pointer(i + 1) - 1
            do while (k <= stop)
                next_column = raw_column(k)
                accumulated = 0.0_dp
                do while (k <= stop)
                    if (raw_column(k) /= next_column) exit
                    accumulated = accumulated + raw_value(k)
                    k = k + 1
                end do
                if (abs(accumulated) <= threshold) cycle
                matrix%column(keep) = next_column
                matrix%value(keep) = accumulated
                keep = keep + 1
            end do
        end do
    end subroutine sparse_from_coo

    pure subroutine sparse_to_dense(matrix, dense, info)
        !! Expand a canonical CSR matrix into a dense array.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse input matrix.
        real(dp), allocatable, intent(out) :: dense(:, :) !! Allocated dense result.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid CSR storage.
        integer :: i
        integer :: k

        call sparse_validate(matrix, info)
        if (info /= 0) then
            allocate(dense(0, 0))
            return
        end if
        allocate(dense(matrix%nrow, matrix%ncol), source=0.0_dp)
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                dense(i, matrix%column(k)) = matrix%value(k)
            end do
        end do
    end subroutine sparse_to_dense

    pure subroutine sparse_matmul_vector(matrix, x, y, info)
        !! Compute y = A*x directly from CSR storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix A.
        real(dp), intent(in) :: x(:) !! Dense vector of length A%ncol.
        real(dp), allocatable, intent(out) :: y(:) !! Allocated product of length A%nrow.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage or shape.
        integer :: i
        integer :: k

        call sparse_validate(matrix, info)
        if (info /= 0 .or. size(x) /= matrix%ncol) then
            allocate(y(0))
            if (info == 0) info = 5
            return
        end if
        allocate(y(matrix%nrow), source=0.0_dp)
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                y(i) = y(i) + matrix%value(k) * x(matrix%column(k))
            end do
        end do
    end subroutine sparse_matmul_vector

    pure subroutine sparse_matmul_matrix(matrix, x, y, info)
        !! Compute Y = A*X directly from CSR storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix A.
        real(dp), intent(in) :: x(:, :) !! Dense matrix with A%ncol rows.
        real(dp), allocatable, intent(out) :: y(:, :) !! Allocated product with A%nrow rows.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage or shape.
        integer :: i
        integer :: k

        call sparse_validate(matrix, info)
        if (info /= 0 .or. size(x, 1) /= matrix%ncol) then
            allocate(y(0, 0))
            if (info == 0) info = 5
            return
        end if
        allocate(y(matrix%nrow, size(x, 2)), source=0.0_dp)
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                y(i, :) = y(i, :) + matrix%value(k) * x(matrix%column(k), :)
            end do
        end do
    end subroutine sparse_matmul_matrix

    pure subroutine sparse_transpose_matmul_matrix(matrix, x, y, info)
        !! Compute Y = transpose(A)*X directly from CSR storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix A.
        real(dp), intent(in) :: x(:, :) !! Dense matrix with A%nrow rows.
        real(dp), allocatable, intent(out) :: y(:, :) !! Allocated product with A%ncol rows.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage or shape.
        integer :: i
        integer :: k

        call sparse_validate(matrix, info)
        if (info /= 0 .or. size(x, 1) /= matrix%nrow) then
            allocate(y(0, 0))
            if (info == 0) info = 5
            return
        end if
        allocate(y(matrix%ncol, size(x, 2)), source=0.0_dp)
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                y(matrix%column(k), :) = y(matrix%column(k), :) + matrix%value(k) * x(i, :)
            end do
        end do
    end subroutine sparse_transpose_matmul_matrix

    pure subroutine sparse_transpose_matmul_sparse(left, right, product, info)
        !! Compute transpose(left)*right as a dense matrix directly from two CSR inputs.
        type(mcmcglmm_sparse_matrix), intent(in) :: left !! Left sparse matrix sharing rows with right.
        type(mcmcglmm_sparse_matrix), intent(in) :: right !! Right sparse matrix sharing rows with left.
        real(dp), allocatable, intent(out) :: product(:, :) !! Allocated dense crossproduct.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid storage or unequal row counts.
        integer :: i
        integer :: j
        integer :: k

        call sparse_validate(left, info)
        if (info /= 0) then
            allocate(product(0, 0))
            return
        end if
        call sparse_validate(right, info)
        if (info /= 0 .or. left%nrow /= right%nrow) then
            allocate(product(0, 0))
            if (info == 0) info = 5
            return
        end if
        allocate(product(left%ncol, right%ncol), source=0.0_dp)
        do i = 1, left%nrow
            do j = left%row_pointer(i), left%row_pointer(i + 1) - 1
                do k = right%row_pointer(i), right%row_pointer(i + 1) - 1
                    product(left%column(j), right%column(k)) = &
                        product(left%column(j), right%column(k)) + left%value(j) * right%value(k)
                end do
            end do
        end do
    end subroutine sparse_transpose_matmul_sparse

    pure subroutine sparse_crossproduct(matrix, product, info)
        !! Compute the dense symmetric crossproduct transpose(A)*A without expanding A.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse matrix A.
        real(dp), allocatable, intent(out) :: product(:, :) !! Allocated A-transpose-A matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid CSR storage.
        integer :: i
        integer :: j
        integer :: k
        integer :: column_j
        integer :: column_k
        real(dp) :: contribution

        call sparse_validate(matrix, info)
        if (info /= 0) then
            allocate(product(0, 0))
            return
        end if
        allocate(product(matrix%ncol, matrix%ncol), source=0.0_dp)
        do i = 1, matrix%nrow
            do j = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                column_j = matrix%column(j)
                do k = j, matrix%row_pointer(i + 1) - 1
                    column_k = matrix%column(k)
                    contribution = matrix%value(j) * matrix%value(k)
                    product(column_j, column_k) = product(column_j, column_k) + contribution
                    if (column_j /= column_k) then
                        product(column_k, column_j) = product(column_k, column_j) + contribution
                    end if
                end do
            end do
        end do
    end subroutine sparse_crossproduct

    pure subroutine sparse_stacked_crossproduct(left, right, product, info)
        !! Form transpose([left right])*[left right] as canonical CSR without a dense design or crossproduct.
        type(mcmcglmm_sparse_matrix), intent(in) :: left !! First CSR column block.
        type(mcmcglmm_sparse_matrix), intent(in) :: right !! Second CSR column block with matching rows.
        type(mcmcglmm_sparse_matrix), intent(out) :: product !! Canonical CSR stacked crossproduct.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or unequal row counts.
        real(dp), allocatable :: accumulator(:)
        integer, allocatable :: column_count(:)
        integer, allocatable :: column_cursor(:)
        integer, allocatable :: column_pointer(:)
        integer, allocatable :: observation(:)
        real(dp), allocatable :: observation_value(:)
        type(sparse_product_row), allocatable :: row(:)
        logical, allocatable :: touched_marker(:)
        integer, allocatable :: touched_column(:)
        integer :: column
        integer :: i
        integer :: k
        integer :: keep
        integer :: number_columns
        integer :: number_nonzero
        integer :: position
        integer :: source_row
        integer :: total_nonzero
        integer :: touched_count
        real(dp) :: source_value

        call sparse_validate(left, info)
        if (info /= 0) return
        call sparse_validate(right, info)
        if (info /= 0 .or. left%nrow /= right%nrow) then
            if (info == 0) info = 5
            return
        end if
        number_columns = left%ncol + right%ncol
        number_nonzero = size(left%value) + size(right%value)
        allocate(column_count(number_columns), source=0)
        do i = 1, left%nrow
            do k = left%row_pointer(i), left%row_pointer(i + 1) - 1
                column_count(left%column(k)) = column_count(left%column(k)) + 1
            end do
            do k = right%row_pointer(i), right%row_pointer(i + 1) - 1
                column = left%ncol + right%column(k)
                column_count(column) = column_count(column) + 1
            end do
        end do
        allocate(column_pointer(number_columns + 1))
        column_pointer(1) = 1
        do column = 1, number_columns
            column_pointer(column + 1) = column_pointer(column) + column_count(column)
        end do
        allocate(column_cursor(number_columns), observation(number_nonzero), observation_value(number_nonzero))
        column_cursor = column_pointer(:number_columns)
        do i = 1, left%nrow
            do k = left%row_pointer(i), left%row_pointer(i + 1) - 1
                column = left%column(k)
                position = column_cursor(column)
                observation(position) = i
                observation_value(position) = left%value(k)
                column_cursor(column) = position + 1
            end do
            do k = right%row_pointer(i), right%row_pointer(i + 1) - 1
                column = left%ncol + right%column(k)
                position = column_cursor(column)
                observation(position) = i
                observation_value(position) = right%value(k)
                column_cursor(column) = position + 1
            end do
        end do
        allocate(accumulator(number_columns), source=0.0_dp)
        allocate(touched_marker(number_columns), source=.false.)
        allocate(touched_column(number_columns), row(number_columns))
        do column = 1, number_columns
            touched_count = 0
            do position = column_pointer(column), column_pointer(column + 1) - 1
                source_row = observation(position)
                source_value = observation_value(position)
                do k = left%row_pointer(source_row), left%row_pointer(source_row + 1) - 1
                    call accumulate_sparse_product(left%column(k), source_value * left%value(k), accumulator, &
                        touched_marker, touched_column, touched_count)
                end do
                do k = right%row_pointer(source_row), right%row_pointer(source_row + 1) - 1
                    call accumulate_sparse_product(left%ncol + right%column(k), source_value * right%value(k), &
                        accumulator, touched_marker, touched_column, touched_count)
                end do
            end do
            call sort_integer_ascending(touched_column(:touched_count))
            keep = 0
            do k = 1, touched_count
                if (accumulator(touched_column(k)) /= 0.0_dp) keep = keep + 1
            end do
            allocate(row(column)%column(keep), row(column)%value(keep))
            keep = 0
            do k = 1, touched_count
                i = touched_column(k)
                if (accumulator(i) /= 0.0_dp) then
                    keep = keep + 1
                    row(column)%column(keep) = i
                    row(column)%value(keep) = accumulator(i)
                end if
                accumulator(i) = 0.0_dp
                touched_marker(i) = .false.
            end do
        end do
        product%nrow = number_columns
        product%ncol = number_columns
        total_nonzero = 0
        do i = 1, number_columns
            total_nonzero = total_nonzero + size(row(i)%column)
        end do
        allocate(product%row_pointer(number_columns + 1), product%column(total_nonzero), &
            product%value(total_nonzero))
        product%row_pointer(1) = 1
        position = 1
        do i = 1, number_columns
            keep = size(row(i)%column)
            product%column(position:position + keep - 1) = row(i)%column
            product%value(position:position + keep - 1) = row(i)%value
            position = position + keep
            product%row_pointer(i + 1) = position
        end do
        info = 0
    end subroutine sparse_stacked_crossproduct

    pure subroutine accumulate_sparse_product(column, contribution, accumulator, marker, touched, touched_count)
        !! Accumulate one crossproduct contribution while recording newly touched output columns.
        integer, intent(in) :: column !! One-based output column.
        real(dp), intent(in) :: contribution !! Value added to the output coordinate.
        real(dp), intent(inout) :: accumulator(:) !! Dense one-row numerical workspace.
        logical, intent(inout) :: marker(:) !! Flags for columns currently present in touched.
        integer, intent(inout) :: touched(:) !! Workspace of columns touched in the current output row.
        integer, intent(inout) :: touched_count !! Number of populated touched entries.

        if (.not. marker(column)) then
            touched_count = touched_count + 1
            touched(touched_count) = column
            marker(column) = .true.
        end if
        accumulator(column) = accumulator(column) + contribution
    end subroutine accumulate_sparse_product

    pure subroutine sort_integer_ascending(values)
        !! Sort an integer vector in ascending order using Shell insertion gaps.
        integer, intent(inout) :: values(:) !! Values sorted in place.
        integer :: gap
        integer :: i
        integer :: j
        integer :: saved

        gap = size(values) / 2
        do while (gap > 0)
            do i = gap + 1, size(values)
                saved = values(i)
                j = i
                do while (j > gap)
                    if (values(j - gap) <= saved) exit
                    values(j) = values(j - gap)
                    j = j - gap
                end do
                values(j) = saved
            end do
            gap = gap / 2
        end do
    end subroutine sort_integer_ascending

end module mcmcglmm_sparse
