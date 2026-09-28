module sparsem_core
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use sparsem_kinds, only: dp
    use sparsem_types, only: csr_matrix, csc_matrix, coo_matrix, ssr_matrix, ssc_matrix
    implicit none
    private

    public :: csr_from_dense, csc_from_dense, coo_from_dense, ssr_from_dense, ssc_from_dense
    public :: csr_to_dense, csc_to_dense, coo_to_dense, ssr_to_dense, ssc_to_dense
    public :: csr_from_csc, csr_from_coo, csr_from_ssr, csr_from_ssc
    public :: csc_from_csr, coo_from_csr, ssr_from_csr, ssc_from_csr
    public :: csr_transpose, csr_shape, csr_diff_rows, csr_diag, csr_set_diag
    public :: csr_norm, csr_rbind, csr_cbind, csr_submatrix, csr_assign_submatrix
    public :: coo_submatrix, coo_assign_submatrix
    public :: csr_add, csr_subtract, csr_hadamard, csr_divide, csr_power
    public :: csr_intdiv, csr_modulo, csr_relation, csr_complement
    public :: csr_matmul, csr_matvec, csr_vecmat, csr_kronecker

contains

    pure function csr_from_dense(a, eps) result(out)
        real(dp), intent(in) :: a(:,:) !! Dense source matrix stored in ordinary Fortran column-major form.
        real(dp), intent(in), optional :: eps !! Entries with absolute value below this threshold are omitted; default is zero.
        type(csr_matrix) :: out
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz

        tol = 0.0_dp
        if (present(eps)) tol = max(0.0_dp, eps)

        out%nrow = size(a, 1)
        out%ncol = size(a, 2)
        nnz = 0
        do i = 1, out%nrow
            do j = 1, out%ncol
                if (keep_value(a(i, j), tol)) nnz = nnz + 1
            end do
        end do

        allocate(out%values(nnz), out%col_ind(nnz), out%row_ptr(out%nrow + 1))
        k = 1
        out%row_ptr(1) = 1
        do i = 1, out%nrow
            do j = 1, out%ncol
                if (keep_value(a(i, j), tol)) then
                    out%values(k) = a(i, j)
                    out%col_ind(k) = j
                    k = k + 1
                end if
            end do
            out%row_ptr(i + 1) = k
        end do
    end function csr_from_dense

    pure function csc_from_dense(a, eps) result(out)
        real(dp), intent(in) :: a(:,:) !! Dense source matrix.
        real(dp), intent(in), optional :: eps !! Magnitude threshold used to drop finite near-zero entries; default is zero.
        type(csc_matrix) :: out
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz

        tol = 0.0_dp
        if (present(eps)) tol = max(0.0_dp, eps)
        out%nrow = size(a, 1)
        out%ncol = size(a, 2)
        nnz = 0
        do j = 1, out%ncol
            do i = 1, out%nrow
                if (keep_value(a(i, j), tol)) nnz = nnz + 1
            end do
        end do
        allocate(out%values(nnz), out%row_ind(nnz), out%col_ptr(out%ncol + 1))
        k = 1
        out%col_ptr(1) = 1
        do j = 1, out%ncol
            do i = 1, out%nrow
                if (keep_value(a(i, j), tol)) then
                    out%values(k) = a(i, j)
                    out%row_ind(k) = i
                    k = k + 1
                end if
            end do
            out%col_ptr(j + 1) = k
        end do
    end function csc_from_dense

    pure function coo_from_dense(a, eps) result(out)
        real(dp), intent(in) :: a(:,:) !! Dense source matrix.
        real(dp), intent(in), optional :: eps !! Magnitude threshold used to drop finite near-zero entries; default is zero.
        type(coo_matrix) :: out
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz

        tol = 0.0_dp
        if (present(eps)) tol = max(0.0_dp, eps)
        out%nrow = size(a, 1)
        out%ncol = size(a, 2)
        nnz = 0
        do i = 1, out%nrow
            do j = 1, out%ncol
                if (keep_value(a(i, j), tol)) nnz = nnz + 1
            end do
        end do
        allocate(out%values(nnz), out%row_ind(nnz), out%col_ind(nnz))
        k = 1
        do i = 1, out%nrow
            do j = 1, out%ncol
                if (keep_value(a(i, j), tol)) then
                    out%values(k) = a(i, j)
                    out%row_ind(k) = i
                    out%col_ind(k) = j
                    k = k + 1
                end if
            end do
        end do
    end function coo_from_dense

    function ssr_from_dense(a, eps, info) result(out)
        real(dp), intent(in) :: a(:,:) !! Dense matrix expected to be square and symmetric.
        real(dp), intent(in), optional :: eps !! Symmetry and storage tolerance; default is 100 machine epsilons.
        integer, intent(out), optional :: info !! Zero on success, one if nonsquare, two if nonsymmetric.
        type(ssr_matrix) :: out
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz

        tol = 100.0_dp * epsilon(1.0_dp)
        if (present(eps)) tol = max(0.0_dp, eps)
        if (present(info)) info = 0
        if (size(a, 1) /= size(a, 2)) then
            if (present(info)) info = 1
            out%n = 0
            allocate(out%values(0), out%col_ind(0), out%row_ptr(1))
            out%row_ptr = 1
            return
        end if
        if (maxval(abs(a - transpose(a))) > tol) then
            if (present(info)) info = 2
            out%n = 0
            allocate(out%values(0), out%col_ind(0), out%row_ptr(1))
            out%row_ptr = 1
            return
        end if
        out%n = size(a, 1)
        nnz = 0
        do i = 1, out%n
            do j = 1, i
                if (keep_value(a(i, j), tol)) nnz = nnz + 1
            end do
        end do
        allocate(out%values(nnz), out%col_ind(nnz), out%row_ptr(out%n + 1))
        k = 1
        out%row_ptr(1) = 1
        do i = 1, out%n
            do j = 1, i
                if (keep_value(a(i, j), tol)) then
                    out%values(k) = a(i, j)
                    out%col_ind(k) = j
                    k = k + 1
                end if
            end do
            out%row_ptr(i + 1) = k
        end do
    end function ssr_from_dense

    function ssc_from_dense(a, eps, info) result(out)
        real(dp), intent(in) :: a(:,:) !! Dense matrix expected to be square and symmetric.
        real(dp), intent(in), optional :: eps !! Symmetry and storage tolerance; default is 100 machine epsilons.
        integer, intent(out), optional :: info !! Zero on success, one if nonsquare, two if nonsymmetric.
        type(ssc_matrix) :: out
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: nnz

        tol = 100.0_dp * epsilon(1.0_dp)
        if (present(eps)) tol = max(0.0_dp, eps)
        if (present(info)) info = 0
        if (size(a, 1) /= size(a, 2)) then
            if (present(info)) info = 1
            out%n = 0
            allocate(out%values(0), out%row_ind(0), out%col_ptr(1))
            out%col_ptr = 1
            return
        end if
        if (maxval(abs(a - transpose(a))) > tol) then
            if (present(info)) info = 2
            out%n = 0
            allocate(out%values(0), out%row_ind(0), out%col_ptr(1))
            out%col_ptr = 1
            return
        end if
        out%n = size(a, 1)
        nnz = 0
        do j = 1, out%n
            do i = j, out%n
                if (keep_value(a(i, j), tol)) nnz = nnz + 1
            end do
        end do
        allocate(out%values(nnz), out%row_ind(nnz), out%col_ptr(out%n + 1))
        k = 1
        out%col_ptr(1) = 1
        do j = 1, out%n
            do i = j, out%n
                if (keep_value(a(i, j), tol)) then
                    out%values(k) = a(i, j)
                    out%row_ind(k) = i
                    k = k + 1
                end if
            end do
            out%col_ptr(j + 1) = k
        end do
    end function ssc_from_dense

    pure function csr_to_dense(a) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix to expand into dense storage.
        real(dp), allocatable :: out(:,:)
        integer :: i
        integer :: k

        allocate(out(a%nrow, a%ncol), source=0.0_dp)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                out(i, a%col_ind(k)) = a%values(k)
            end do
        end do
    end function csr_to_dense

    pure function csc_to_dense(a) result(out)
        type(csc_matrix), intent(in) :: a !! CSC matrix to expand into dense storage.
        real(dp), allocatable :: out(:,:)
        integer :: j
        integer :: k

        allocate(out(a%nrow, a%ncol), source=0.0_dp)
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                out(a%row_ind(k), j) = a%values(k)
            end do
        end do
    end function csc_to_dense

    pure function coo_to_dense(a) result(out)
        type(coo_matrix), intent(in) :: a !! COO matrix to expand into dense storage; duplicate entries are summed.
        real(dp), allocatable :: out(:,:)
        integer :: k

        allocate(out(a%nrow, a%ncol), source=0.0_dp)
        do k = 1, size(a%values)
            out(a%row_ind(k), a%col_ind(k)) = out(a%row_ind(k), a%col_ind(k)) + a%values(k)
        end do
    end function coo_to_dense

    pure function ssr_to_dense(a) result(out)
        type(ssr_matrix), intent(in) :: a !! Symmetric sparse-row lower-triangle storage to expand.
        real(dp), allocatable :: out(:,:)
        integer :: i
        integer :: j
        integer :: k

        allocate(out(a%n, a%n), source=0.0_dp)
        do i = 1, a%n
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col_ind(k)
                out(i, j) = a%values(k)
                out(j, i) = a%values(k)
            end do
        end do
    end function ssr_to_dense

    pure function ssc_to_dense(a) result(out)
        type(ssc_matrix), intent(in) :: a !! Symmetric sparse-column lower-triangle storage to expand.
        real(dp), allocatable :: out(:,:)
        integer :: i
        integer :: j
        integer :: k

        allocate(out(a%n, a%n), source=0.0_dp)
        do j = 1, a%n
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                i = a%row_ind(k)
                out(i, j) = a%values(k)
                out(j, i) = a%values(k)
            end do
        end do
    end function ssc_to_dense

    pure function csr_from_csc(a) result(out)
        type(csc_matrix), intent(in) :: a !! CSC source matrix.
        type(csr_matrix) :: out
        out = csr_from_dense(csc_to_dense(a))
    end function csr_from_csc

    pure function csr_from_coo(a) result(out)
        type(coo_matrix), intent(in) :: a !! COO source matrix; duplicate coordinates are summed.
        type(csr_matrix) :: out
        out = csr_from_dense(coo_to_dense(a))
    end function csr_from_coo

    pure function csr_from_ssr(a) result(out)
        type(ssr_matrix), intent(in) :: a !! Symmetric sparse-row source matrix.
        type(csr_matrix) :: out
        out = csr_from_dense(ssr_to_dense(a))
    end function csr_from_ssr

    pure function csr_from_ssc(a) result(out)
        type(ssc_matrix), intent(in) :: a !! Symmetric sparse-column source matrix.
        type(csr_matrix) :: out
        out = csr_from_dense(ssc_to_dense(a))
    end function csr_from_ssc

    pure function csc_from_csr(a) result(out)
        type(csr_matrix), intent(in) :: a !! CSR source matrix.
        type(csc_matrix) :: out
        out = csc_from_dense(csr_to_dense(a))
    end function csc_from_csr

    pure function coo_from_csr(a) result(out)
        type(csr_matrix), intent(in) :: a !! CSR source matrix.
        type(coo_matrix) :: out
        integer :: i
        integer :: k
        integer :: p

        out%nrow = a%nrow
        out%ncol = a%ncol
        allocate(out%values(size(a%values)), out%row_ind(size(a%values)), out%col_ind(size(a%values)))
        p = 1
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                out%values(p) = a%values(k)
                out%row_ind(p) = i
                out%col_ind(p) = a%col_ind(k)
                p = p + 1
            end do
        end do
    end function coo_from_csr

    function ssr_from_csr(a, eps, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix expected to represent a square symmetric matrix.
        real(dp), intent(in), optional :: eps !! Symmetry/storage tolerance forwarded to the conversion.
        integer, intent(out), optional :: info !! Conversion status: zero success, one nonsquare, two nonsymmetric.
        type(ssr_matrix) :: out
        real(dp), allocatable :: d(:,:)
        d = csr_to_dense(a)
        if (present(eps)) then
            if (present(info)) then
                out = ssr_from_dense(d, eps, info)
            else
                out = ssr_from_dense(d, eps)
            end if
        else
            if (present(info)) then
                out = ssr_from_dense(d, info=info)
            else
                out = ssr_from_dense(d)
            end if
        end if
    end function ssr_from_csr

    function ssc_from_csr(a, eps, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix expected to represent a square symmetric matrix.
        real(dp), intent(in), optional :: eps !! Symmetry/storage tolerance forwarded to the conversion.
        integer, intent(out), optional :: info !! Conversion status: zero success, one nonsquare, two nonsymmetric.
        type(ssc_matrix) :: out
        real(dp), allocatable :: d(:,:)
        d = csr_to_dense(a)
        if (present(eps)) then
            if (present(info)) then
                out = ssc_from_dense(d, eps, info)
            else
                out = ssc_from_dense(d, eps)
            end if
        else
            if (present(info)) then
                out = ssc_from_dense(d, info=info)
            else
                out = ssc_from_dense(d)
            end if
        end if
    end function ssc_from_csr

    pure function csr_transpose(a) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix to transpose.
        type(csr_matrix) :: out
        integer, allocatable :: counts(:)
        integer, allocatable :: next(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: p

        out%nrow = a%ncol
        out%ncol = a%nrow
        allocate(out%values(size(a%values)), out%col_ind(size(a%values)), out%row_ptr(out%nrow + 1))
        allocate(counts(out%nrow), source=0)
        do k = 1, size(a%values)
            counts(a%col_ind(k)) = counts(a%col_ind(k)) + 1
        end do
        out%row_ptr(1) = 1
        do i = 1, out%nrow
            out%row_ptr(i + 1) = out%row_ptr(i) + counts(i)
        end do
        allocate(next(out%nrow))
        next = out%row_ptr(1:out%nrow)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col_ind(k)
                p = next(j)
                out%values(p) = a%values(k)
                out%col_ind(p) = i
                next(j) = p + 1
            end do
        end do
    end function csr_transpose

    pure subroutine csr_shape(a, nrow, ncol)
        type(csr_matrix), intent(in) :: a !! CSR matrix whose dimensions are requested.
        integer, intent(out) :: nrow !! Number of matrix rows.
        integer, intent(out) :: ncol !! Number of matrix columns.
        nrow = a%nrow
        ncol = a%ncol
    end subroutine csr_shape

    function csr_diff_rows(a, lag, differences, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix whose rows are differenced.
        integer, intent(in), optional :: lag !! Positive row lag; default is one.
        integer, intent(in), optional :: differences !! Positive number of repeated differences; default is one.
        integer, intent(out), optional :: info !! Zero on success and one for invalid lag/difference settings.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)
        real(dp), allocatable :: tmp(:,:)
        integer :: l
        integer :: nd
        integer :: k

        l = 1
        nd = 1
        if (present(lag)) l = lag
        if (present(differences)) nd = differences
        if (present(info)) info = 0
        if (l < 1 .or. nd < 1 .or. l * nd >= a%nrow) then
            if (present(info)) info = 1
            out%nrow = 0
            out%ncol = a%ncol
            allocate(out%values(0), out%col_ind(0), out%row_ptr(1))
            out%row_ptr = 1
            return
        end if
        d = csr_to_dense(a)
        do k = 1, nd
            allocate(tmp(size(d, 1) - l, size(d, 2)))
            tmp = d(l + 1:size(d, 1), :) - d(1:size(d, 1) - l, :)
            call move_alloc(tmp, d)
        end do
        out = csr_from_dense(d)
    end function csr_diff_rows

    pure function csr_diag(a) result(d)
        type(csr_matrix), intent(in) :: a !! CSR matrix from which to extract the main diagonal.
        real(dp), allocatable :: d(:)
        integer :: i
        integer :: k
        integer :: n

        n = min(a%nrow, a%ncol)
        allocate(d(n), source=0.0_dp)
        do i = 1, n
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col_ind(k) == i) then
                    d(i) = a%values(k)
                    exit
                end if
            end do
        end do
    end function csr_diag

    function csr_set_diag(a, d, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix whose main diagonal is replaced.
        real(dp), intent(in) :: d(:) !! Replacement diagonal; length one broadcasts, otherwise length must equal min dimensions.
        integer, intent(out), optional :: info !! Zero on success and one if the replacement length is incompatible.
        type(csr_matrix) :: out
        real(dp), allocatable :: dense(:,:)
        integer :: i
        integer :: n

        if (present(info)) info = 0
        n = min(a%nrow, a%ncol)
        if (size(d) /= 1 .and. size(d) /= n) then
            if (present(info)) info = 1
            out = a
            return
        end if
        dense = csr_to_dense(a)
        do i = 1, n
            if (size(d) == 1) then
                dense(i, i) = d(1)
            else
                dense(i, i) = d(i)
            end if
        end do
        out = csr_from_dense(dense)
    end function csr_set_diag

    function csr_norm(a, norm_type, info) result(value)
        type(csr_matrix), intent(in) :: a !! CSR matrix whose sparse-entry norm is requested.
        character(len=*), intent(in), optional :: norm_type !! SparseM norm code: sup/M, HS/F, or l1; default is sup.
        integer, intent(out), optional :: info !! Zero on success and one for an unsupported norm code.
        real(dp) :: value
        character(len=:), allocatable :: kind

        if (present(info)) info = 0
        kind = 'sup'
        if (present(norm_type)) kind = trim(norm_type)
        select case (kind)
        case ('sup', 'M')
            if (size(a%values) == 0) then
                value = 0.0_dp
            else
                value = maxval(abs(a%values))
            end if
        case ('HS', 'F')
            value = sqrt(sum(a%values * a%values))
        case ('l1')
            value = sum(abs(a%values))
        case default
            value = 0.0_dp
            if (present(info)) info = 1
        end select
    end function csr_norm

    function csr_rbind(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Upper matrix in the row binding operation.
        type(csr_matrix), intent(in) :: b !! Lower matrix; column count must match the first matrix.
        integer, intent(out), optional :: info !! Zero on success and one when column counts differ.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (a%ncol /= b%ncol) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        allocate(d(a%nrow + b%nrow, a%ncol))
        d(1:a%nrow, :) = csr_to_dense(a)
        d(a%nrow + 1:, :) = csr_to_dense(b)
        out = csr_from_dense(d)
    end function csr_rbind

    function csr_cbind(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Left matrix in the column binding operation.
        type(csr_matrix), intent(in) :: b !! Right matrix; row count must match the first matrix.
        integer, intent(out), optional :: info !! Zero on success and one when row counts differ.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (a%nrow /= b%nrow) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        allocate(d(a%nrow, a%ncol + b%ncol))
        d(:, 1:a%ncol) = csr_to_dense(a)
        d(:, a%ncol + 1:) = csr_to_dense(b)
        out = csr_from_dense(d)
    end function csr_cbind

    function csr_submatrix(a, rows, cols, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR source matrix.
        integer, intent(in) :: rows(:) !! One-based row indices; duplicates are retained in the requested order.
        integer, intent(in) :: cols(:) !! One-based column indices; duplicates are retained in the requested order.
        integer, intent(out), optional :: info !! Zero on success and one when any index is out of bounds.
        type(csr_matrix) :: out
        real(dp), allocatable :: dense(:,:)
        real(dp), allocatable :: sub(:,:)
        integer :: i
        integer :: j

        if (present(info)) info = 0
        if (any(rows < 1) .or. any(rows > a%nrow) .or. any(cols < 1) .or. any(cols > a%ncol)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        dense = csr_to_dense(a)
        allocate(sub(size(rows), size(cols)))
        do i = 1, size(rows)
            do j = 1, size(cols)
                sub(i, j) = dense(rows(i), cols(j))
            end do
        end do
        out = csr_from_dense(sub)
    end function csr_submatrix

    function csr_assign_submatrix(a, rows, cols, value, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix to copy and modify.
        integer, intent(in) :: rows(:) !! One-based target row indices.
        integer, intent(in) :: cols(:) !! One-based target column indices.
        type(csr_matrix), intent(in) :: value !! Replacement matrix matching the selected row and column counts.
        integer, intent(out), optional :: info !! Zero on success, one for bounds errors, two for replacement-shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: dense(:,:)
        real(dp), allocatable :: replacement(:,:)
        integer :: i
        integer :: j

        if (present(info)) info = 0
        if (any(rows < 1) .or. any(rows > a%nrow) .or. any(cols < 1) .or. any(cols > a%ncol)) then
            if (present(info)) info = 1
            out = a
            return
        end if
        if (value%nrow /= size(rows) .or. value%ncol /= size(cols)) then
            if (present(info)) info = 2
            out = a
            return
        end if
        dense = csr_to_dense(a)
        replacement = csr_to_dense(value)
        do i = 1, size(rows)
            do j = 1, size(cols)
                dense(rows(i), cols(j)) = replacement(i, j)
            end do
        end do
        out = csr_from_dense(dense)
    end function csr_assign_submatrix

    function coo_submatrix(a, rows, cols, info) result(out)
        type(coo_matrix), intent(in) :: a !! COO source matrix.
        integer, intent(in) :: rows(:) !! One-based row indices to retain.
        integer, intent(in) :: cols(:) !! One-based column indices to retain.
        integer, intent(out), optional :: info !! Zero on success and one when any index is out of bounds.
        type(coo_matrix) :: out
        type(csr_matrix) :: tmp
        tmp = csr_submatrix(csr_from_coo(a), rows, cols, info)
        out = coo_from_csr(tmp)
    end function coo_submatrix

    function coo_assign_submatrix(a, rows, cols, value, info) result(out)
        type(coo_matrix), intent(in) :: a !! COO matrix to copy and modify.
        integer, intent(in) :: rows(:) !! One-based target row indices.
        integer, intent(in) :: cols(:) !! One-based target column indices.
        type(coo_matrix), intent(in) :: value !! COO replacement matrix matching the selected shape.
        integer, intent(out), optional :: info !! Zero on success, one for bounds errors, two for replacement-shape mismatch.
        type(coo_matrix) :: out
        type(csr_matrix) :: tmp
        tmp = csr_assign_submatrix(csr_from_coo(a), rows, cols, csr_from_coo(value), info)
        out = coo_from_csr(tmp)
    end function coo_assign_submatrix

    function csr_add(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! First CSR addend.
        type(csr_matrix), intent(in) :: b !! Second CSR addend with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        d = csr_to_dense(a) + csr_to_dense(b)
        out = csr_from_dense(d)
    end function csr_add

    function csr_subtract(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! CSR minuend.
        type(csr_matrix), intent(in) :: b !! CSR subtrahend with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        d = csr_to_dense(a) - csr_to_dense(b)
        out = csr_from_dense(d)
    end function csr_subtract

    function csr_hadamard(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! First CSR matrix in the elementwise product.
        type(csr_matrix), intent(in) :: b !! Second CSR matrix with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        d = csr_to_dense(a) * csr_to_dense(b)
        out = csr_from_dense(d)
    end function csr_hadamard

    function csr_divide(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Numerator CSR matrix.
        type(csr_matrix), intent(in) :: b !! Denominator CSR matrix with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: da(:,:)
        real(dp), allocatable :: db(:,:)
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        da = csr_to_dense(a)
        db = csr_to_dense(b)
        allocate(d(a%nrow, a%ncol))
        d = da / db
        out = csr_from_dense(d)
    end function csr_divide

    function csr_power(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Base CSR matrix.
        type(csr_matrix), intent(in) :: b !! Exponent CSR matrix with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: da(:,:)
        real(dp), allocatable :: db(:,:)
        real(dp), allocatable :: d(:,:)

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        da = csr_to_dense(a)
        db = csr_to_dense(b)
        allocate(d(a%nrow, a%ncol))
        d = da ** db
        out = csr_from_dense(d)
    end function csr_power

    function csr_intdiv(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Dividend CSR matrix.
        type(csr_matrix), intent(in) :: b !! Divisor CSR matrix with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: da(:,:)
        real(dp), allocatable :: db(:,:)
        real(dp), allocatable :: d(:,:)
        integer :: i
        integer :: j

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        da = csr_to_dense(a)
        db = csr_to_dense(b)
        allocate(d(a%nrow, a%ncol))
        do i = 1, a%nrow
            do j = 1, a%ncol
                d(i, j) = floor(da(i, j) / db(i, j))
            end do
        end do
        out = csr_from_dense(d)
    end function csr_intdiv

    function csr_modulo(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Dividend CSR matrix.
        type(csr_matrix), intent(in) :: b !! Divisor CSR matrix with the same shape.
        integer, intent(out), optional :: info !! Zero on success and one for shape mismatch.
        type(csr_matrix) :: out
        real(dp), allocatable :: da(:,:)
        real(dp), allocatable :: db(:,:)
        real(dp), allocatable :: d(:,:)
        integer :: i
        integer :: j

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        da = csr_to_dense(a)
        db = csr_to_dense(b)
        allocate(d(a%nrow, a%ncol))
        do i = 1, a%nrow
            do j = 1, a%ncol
                d(i, j) = modulo(da(i, j), db(i, j))
            end do
        end do
        out = csr_from_dense(d)
    end function csr_modulo

    function csr_relation(a, b, relation, info) result(out)
        type(csr_matrix), intent(in) :: a !! Left CSR operand.
        type(csr_matrix), intent(in) :: b !! Right CSR operand with the same shape.
        character(len=*), intent(in) :: relation !! Relation code: gt, ge, lt, le, eq, or ne.
        integer, intent(out), optional :: info !! Zero on success, one for shape mismatch, two for an unknown relation code.
        type(csr_matrix) :: out
        real(dp), allocatable :: da(:,:)
        real(dp), allocatable :: db(:,:)
        real(dp), allocatable :: d(:,:)
        integer :: i
        integer :: j

        if (present(info)) info = 0
        if (.not. same_shape(a, b)) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        da = csr_to_dense(a)
        db = csr_to_dense(b)
        allocate(d(a%nrow, a%ncol), source=0.0_dp)
        select case (trim(relation))
        case ('gt')
            where (da > db) d = 1.0_dp
        case ('ge')
            where (da >= db) d = 1.0_dp
        case ('lt')
            where (da < db) d = 1.0_dp
        case ('le')
            where (da <= db) d = 1.0_dp
        case ('eq')
            where (da == db) d = 1.0_dp
        case ('ne')
            where (da /= db) d = 1.0_dp
        case default
            if (present(info)) info = 2
            d = 0.0_dp
        end select
        do i = 1, a%nrow
            do j = 1, a%ncol
                if (ieee_is_nan(da(i, j)) .or. ieee_is_nan(db(i, j))) d(i, j) = 0.0_dp
            end do
        end do
        out = csr_from_dense(d)
    end function csr_relation

    pure function csr_complement(a) result(out)
        type(csr_matrix), intent(in) :: a !! CSR matrix interpreted as logical, where zero is false and nonzero is true.
        type(csr_matrix) :: out
        real(dp), allocatable :: d(:,:)
        real(dp), allocatable :: x(:,:)

        x = csr_to_dense(a)
        allocate(d(a%nrow, a%ncol), source=0.0_dp)
        where (x == 0.0_dp) d = 1.0_dp
        out = csr_from_dense(d)
    end function csr_complement

    function csr_matmul(a, b, info) result(out)
        type(csr_matrix), intent(in) :: a !! Left CSR matrix.
        type(csr_matrix), intent(in) :: b !! Right CSR matrix with row count equal to the left column count.
        integer, intent(out), optional :: info !! Zero on success and one for nonconformable dimensions.
        type(csr_matrix) :: out
        real(dp), allocatable :: dense(:,:)
        real(dp), allocatable :: acc(:)
        integer :: i
        integer :: ka
        integer :: kb
        integer :: row_b

        if (present(info)) info = 0
        if (a%ncol /= b%nrow) then
            if (present(info)) info = 1
            out = empty_csr(0, 0)
            return
        end if
        allocate(dense(a%nrow, b%ncol), source=0.0_dp)
        allocate(acc(b%ncol))
        do i = 1, a%nrow
            acc = 0.0_dp
            do ka = a%row_ptr(i), a%row_ptr(i + 1) - 1
                row_b = a%col_ind(ka)
                do kb = b%row_ptr(row_b), b%row_ptr(row_b + 1) - 1
                    acc(b%col_ind(kb)) = acc(b%col_ind(kb)) + a%values(ka) * b%values(kb)
                end do
            end do
            dense(i, :) = acc
        end do
        out = csr_from_dense(dense)
    end function csr_matmul

    function csr_matvec(a, x, info) result(y)
        type(csr_matrix), intent(in) :: a !! CSR matrix in a matrix-vector product.
        real(dp), intent(in) :: x(:) !! Dense vector of length equal to the matrix column count.
        integer, intent(out), optional :: info !! Zero on success and one for a nonconformable vector length.
        real(dp), allocatable :: y(:)
        integer :: i
        integer :: k

        if (present(info)) info = 0
        if (size(x) /= a%ncol) then
            if (present(info)) info = 1
            allocate(y(0))
            return
        end if
        allocate(y(a%nrow), source=0.0_dp)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                y(i) = y(i) + a%values(k) * x(a%col_ind(k))
            end do
        end do
    end function csr_matvec

    function csr_vecmat(x, a, info) result(y)
        real(dp), intent(in) :: x(:) !! Dense row vector of length equal to the matrix row count.
        type(csr_matrix), intent(in) :: a !! CSR matrix in a vector-matrix product.
        integer, intent(out), optional :: info !! Zero on success and one for a nonconformable vector length.
        real(dp), allocatable :: y(:)
        integer :: i
        integer :: k

        if (present(info)) info = 0
        if (size(x) /= a%nrow) then
            if (present(info)) info = 1
            allocate(y(0))
            return
        end if
        allocate(y(a%ncol), source=0.0_dp)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                y(a%col_ind(k)) = y(a%col_ind(k)) + x(i) * a%values(k)
            end do
        end do
    end function csr_vecmat

    pure function csr_kronecker(a, b) result(out)
        type(csr_matrix), intent(in) :: a !! First CSR Kronecker factor.
        type(csr_matrix), intent(in) :: b !! Second CSR Kronecker factor.
        type(csr_matrix) :: out
        type(coo_matrix) :: coo
        integer :: ia
        integer :: ib
        integer :: ka
        integer :: kb
        integer :: p

        coo%nrow = a%nrow * b%nrow
        coo%ncol = a%ncol * b%ncol
        allocate(coo%values(size(a%values) * size(b%values)))
        allocate(coo%row_ind(size(coo%values)), coo%col_ind(size(coo%values)))
        p = 1
        do ia = 1, a%nrow
            do ka = a%row_ptr(ia), a%row_ptr(ia + 1) - 1
                do ib = 1, b%nrow
                    do kb = b%row_ptr(ib), b%row_ptr(ib + 1) - 1
                        coo%values(p) = a%values(ka) * b%values(kb)
                        coo%row_ind(p) = (ia - 1) * b%nrow + ib
                        coo%col_ind(p) = (a%col_ind(ka) - 1) * b%ncol + b%col_ind(kb)
                        p = p + 1
                    end do
                end do
            end do
        end do
        out = csr_from_coo(coo)
    end function csr_kronecker

    pure logical function keep_value(x, eps)
        real(dp), intent(in) :: x !! Candidate matrix value, including possible IEEE NaN or infinity.
        real(dp), intent(in) :: eps !! Nonnegative finite-value magnitude threshold.
        keep_value = abs(x) >= eps
        if (eps == 0.0_dp .and. x == 0.0_dp) keep_value = .false.
        if (ieee_is_nan(x)) keep_value = .true.
    end function keep_value

    pure logical function same_shape(a, b)
        type(csr_matrix), intent(in) :: a !! First CSR matrix in a shape comparison.
        type(csr_matrix), intent(in) :: b !! Second CSR matrix in a shape comparison.
        same_shape = a%nrow == b%nrow .and. a%ncol == b%ncol
    end function same_shape

    pure function empty_csr(nrow, ncol) result(out)
        integer, intent(in) :: nrow !! Desired nonnegative row count for an empty CSR matrix.
        integer, intent(in) :: ncol !! Desired nonnegative column count for an empty CSR matrix.
        type(csr_matrix) :: out
        out%nrow = max(0, nrow)
        out%ncol = max(0, ncol)
        allocate(out%values(0), out%col_ind(0), out%row_ptr(out%nrow + 1))
        out%row_ptr = 1
    end function empty_csr

end module sparsem_core
