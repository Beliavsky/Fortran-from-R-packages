module e1071_sparse
    use e1071_kinds, only: dp
    implicit none
    private

    type, public :: matrix_csr
        real(dp), allocatable :: values(:)
        integer, allocatable :: col_index(:)
        integer, allocatable :: row_pointer(:)
        integer :: nrow = 0
        integer :: ncol = 0
    end type matrix_csr

    public :: dense_to_csr, csr_to_dense
    public :: read_matrix_csr, write_matrix_csr

contains

    subroutine dense_to_csr(x, matrix, tolerance)
        real(dp), intent(in) :: x(:, :) !! Dense row-by-column matrix converted to compressed sparse-row storage.
        type(matrix_csr), intent(out) :: matrix !! CSR result with one-based column indices and row pointers.
        real(dp), intent(in), optional :: tolerance !! Entries with absolute value at most this threshold are omitted; default zero.
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: n

        tol = 0.0_dp
        if (present(tolerance)) tol = max(0.0_dp, tolerance)
        n = count(abs(x) > tol)
        allocate(matrix%values(n), matrix%col_index(n), matrix%row_pointer(size(x, 1) + 1))
        matrix%nrow = size(x, 1)
        matrix%ncol = size(x, 2)
        n = 0
        matrix%row_pointer(1) = 1
        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                if (abs(x(i, j)) <= tol) cycle
                n = n + 1
                matrix%values(n) = x(i, j)
                matrix%col_index(n) = j
            end do
            matrix%row_pointer(i + 1) = n + 1
        end do
    end subroutine dense_to_csr

    subroutine csr_to_dense(matrix, x)
        type(matrix_csr), intent(in) :: matrix !! CSR matrix with one-based column indices and row pointers.
        real(dp), allocatable, intent(out) :: x(:, :) !! Dense row-by-column matrix reconstructed from the sparse entries.
        integer :: i
        integer :: k

        call validate_csr(matrix)
        allocate(x(matrix%nrow, matrix%ncol))
        x = 0.0_dp
        do i = 1, matrix%nrow
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                x(i, matrix%col_index(k)) = matrix%values(k)
            end do
        end do
    end subroutine csr_to_dense

    subroutine read_matrix_csr(filename, matrix, y, has_y, ncol)
        character(len=*), intent(in) :: filename !! LIBSVM-style text file containing `index:value` sparse rows and optional labels.
        type(matrix_csr), intent(out) :: matrix !! Parsed CSR matrix using one-based feature indices exactly as in the text file.
        real(dp), allocatable, intent(out), optional :: y(:) !! Optional numeric leading label from each row; empty when labels
        !! are absent.
        logical, intent(out), optional :: has_y !! True when the first token of every nonempty row is a response rather than
        !! `index:value`.
        integer, intent(in), optional :: ncol !! Optional minimum output column count; columns implied by the file are never
        !! discarded.
        character(len=16384) :: line
        character(len=512) :: token
        integer :: unit
        integer :: ios
        integer :: nrow
        integer :: nnz
        integer :: maxcol
        integer :: pos
        integer :: colon
        integer :: row
        integer :: k
        integer :: col
        real(dp) :: value
        real(dp) :: label_value
        logical :: row_has_y
        logical :: file_has_y

        nrow = 0
        nnz = 0
        maxcol = 0
        file_has_y = .false.
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "read_matrix_csr: unable to open input file"
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            nrow = nrow + 1
            pos = 1
            call next_token(line, pos, token)
            row_has_y = index(trim(token), ':') == 0
            if (nrow == 1) then
                file_has_y = row_has_y
            else if (row_has_y .neqv. file_has_y) then
                error stop "read_matrix_csr: inconsistent presence of response labels"
            end if
            if (.not. row_has_y) pos = 1
            do
                call next_token(line, pos, token)
                if (len_trim(token) == 0) exit
                colon = index(token, ':')
                if (colon <= 1) error stop "read_matrix_csr: malformed sparse token"
                read(token(:colon - 1), *, iostat=ios) col
                if (ios /= 0 .or. col < 1) error stop "read_matrix_csr: invalid column index"
                maxcol = max(maxcol, col)
                nnz = nnz + 1
            end do
        end do
        close(unit)
        if (present(ncol)) maxcol = max(maxcol, ncol)
        allocate(matrix%values(nnz), matrix%col_index(nnz), matrix%row_pointer(nrow + 1))
        if (present(y)) then
            if (file_has_y) then
                allocate(y(nrow))
            else
                allocate(y(0))
            end if
        end if
        matrix%nrow = nrow
        matrix%ncol = maxcol
        matrix%row_pointer(1) = 1
        open(newunit=unit, file=filename, status="old", action="read")
        row = 0
        k = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            row = row + 1
            pos = 1
            if (file_has_y) then
                call next_token(line, pos, token)
                read(token, *, iostat=ios) label_value
                if (ios /= 0) error stop "read_matrix_csr: nonnumeric response label"
                if (present(y)) y(row) = label_value
            end if
            do
                call next_token(line, pos, token)
                if (len_trim(token) == 0) exit
                colon = index(token, ':')
                read(token(:colon - 1), *, iostat=ios) col
                if (ios /= 0) error stop "read_matrix_csr: invalid column index"
                read(token(colon + 1:), *, iostat=ios) value
                if (ios /= 0) error stop "read_matrix_csr: invalid sparse value"
                k = k + 1
                matrix%col_index(k) = col
                matrix%values(k) = value
            end do
            matrix%row_pointer(row + 1) = k + 1
        end do
        close(unit)
        if (present(has_y)) has_y = file_has_y
    end subroutine read_matrix_csr

    subroutine write_matrix_csr(matrix, filename, y)
        type(matrix_csr), intent(in) :: matrix !! CSR matrix written in LIBSVM-compatible one-based `index:value` syntax.
        character(len=*), intent(in) :: filename !! Output text-file path, created or replaced by this routine.
        real(dp), intent(in), optional :: y(:) !! Optional numeric response label written before each sparse row.
        integer :: unit
        integer :: i
        integer :: k
        logical :: need_last_zero

        call validate_csr(matrix)
        if (present(y)) then
            if (size(y) /= matrix%nrow) error stop "write_matrix_csr: response length mismatch"
        end if
        need_last_zero = .false.
        if (size(matrix%col_index) == 0) then
            need_last_zero = matrix%ncol > 0
        else
            need_last_zero = maxval(matrix%col_index) < matrix%ncol
        end if
        open(newunit=unit, file=filename, status="replace", action="write")
        do i = 1, matrix%nrow
            if (present(y)) write(unit, '(g0,1x)', advance='no') y(i)
            do k = matrix%row_pointer(i), matrix%row_pointer(i + 1) - 1
                write(unit, '(i0,a,g0,1x)', advance='no') matrix%col_index(k), ':', matrix%values(k)
            end do
            if (need_last_zero .and. i == 1) then
                write(unit, '(i0,a,g0,1x)', advance='no') matrix%ncol, ':', 0.0_dp
                need_last_zero = .false.
            end if
            write(unit, '()')
        end do
        close(unit)
    end subroutine write_matrix_csr

    subroutine validate_csr(matrix)
        type(matrix_csr), intent(in) :: matrix !! CSR matrix whose dimensions, pointers, and column indices are checked.

        if (matrix%nrow < 0 .or. matrix%ncol < 0) error stop "validate_csr: negative dimension"
        if (.not. allocated(matrix%values)) error stop "validate_csr: values are not allocated"
        if (.not. allocated(matrix%col_index)) error stop "validate_csr: col_index is not allocated"
        if (.not. allocated(matrix%row_pointer)) error stop "validate_csr: row_pointer is not allocated"
        if (size(matrix%values) /= size(matrix%col_index)) error stop "validate_csr: value/index length mismatch"
        if (size(matrix%row_pointer) /= matrix%nrow + 1) error stop "validate_csr: row_pointer length mismatch"
        if (matrix%row_pointer(1) /= 1) error stop "validate_csr: first row pointer must be one"
        if (matrix%row_pointer(matrix%nrow + 1) /= size(matrix%values) + 1) then
            error stop "validate_csr: terminal row pointer is inconsistent"
        end if
        if (size(matrix%col_index) > 0) then
            if (any(matrix%col_index < 1) .or. any(matrix%col_index > matrix%ncol)) then
                error stop "validate_csr: column index out of range"
            end if
        end if
    end subroutine validate_csr

    subroutine next_token(line, position, token)
        character(len=*), intent(in) :: line !! Input line scanned left-to-right for whitespace-separated tokens.
        integer, intent(inout) :: position !! One-based scan position, advanced past the returned token.
        character(len=*), intent(out) :: token !! Next token, or an empty string when the line is exhausted.
        integer :: first
        integer :: last
        integer :: n

        token = ''
        n = len_trim(line)
        do while (position <= n)
            if (line(position:position) /= ' ' .and. line(position:position) /= char(9)) exit
            position = position + 1
        end do
        if (position > n) return
        first = position
        do while (position <= n)
            if (line(position:position) == ' ' .or. line(position:position) == char(9)) exit
            position = position + 1
        end do
        last = position - 1
        token = line(first:last)
    end subroutine next_token

end module e1071_sparse
