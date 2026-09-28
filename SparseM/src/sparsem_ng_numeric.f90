! Supernodal numerical Cholesky routines for SparseM.
!
! These routines are modern free-form translations of the Ng-Peyton block
! factorization path distributed with SparseM. Algorithm/source credits include
! Esmond G. Ng, Barry W. Peyton, Guodong Zhang, Joseph W. H. Liu, and the Oak
! Ridge National Laboratory Mathematical Sciences Section. See NOTICE.md.
!
module sparsem_ng_numeric
    use sparsem_kinds, only: dp
    use sparsem_types, only: csr_matrix
    implicit none
    private

    public :: ng_numeric_factor, ng_forwardsolve, ng_backsolve, ng_solve

contains

    subroutine ng_numeric_factor(a, perm, inv_perm, nsuper, xsuper, supernode, xlindx, lindx, xlnz, &
            lnz, split, tmp_size, info, cache_kb, tiny, large, level, tmpmax)
        type(csr_matrix), intent(in) :: a !! Symmetric matrix in original CSR numbering whose lower triangle is factored.
        integer, intent(in) :: perm(:) !! New-to-old symbolic ordering used for the factor columns.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse ordering corresponding to perm.
        integer, intent(in) :: nsuper !! Number of fundamental supernodes in the symbolic factorization.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers with active extent nsuper+1.
        integer, intent(in) :: supernode(:) !! Supernode membership number for every ordered column.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode.
        integer, intent(in) :: lindx(:) !! Compressed row-index structure for the supernodal factor.
        integer, intent(in) :: xlnz(:) !! Column pointers into numerical lower-factor storage.
        real(dp), allocatable, intent(out) :: lnz(:) !! Numerical supernodal Cholesky storage on return.
        integer, allocatable, intent(out) :: split(:) !! Cache-block split sizes indexed from supernode first columns.
        integer, intent(out) :: tmp_size !! Floating-point workspace size computed for block factorization.
        integer, intent(out) :: info !! Zero on success; -1 structural error; 16+k when k tiny pivots are replaced.
        integer, intent(in), optional :: cache_kb !! Cache size in KiB used to split supernodes; default is 1024.
        real(dp), intent(in), optional :: tiny !! Relative tiny-pivot threshold; default is 1e-30.
        real(dp), intent(in), optional :: large !! Replacement value for tiny pivots; default is 1e128.
        integer, intent(in), optional :: level !! Numerical loop-unrolling level; valid values are 1, 2, 4, or 8, default 8.
        integer, intent(in), optional :: tmpmax !! Maximum temporary real workspace; status 8 avoids allocation beyond this limit.
        real(dp), allocatable :: temp(:)
        integer, allocatable :: index_map(:)
        integer, allocatable :: length(:)
        integer, allocatable :: link(:)
        integer, allocatable :: offset(:)
        integer, allocatable :: relative_index(:)
        integer :: cache_size
        integer :: n
        integer :: numeric_size
        integer :: unroll_level
        real(dp) :: large_value
        real(dp) :: tiny_value

        info = 0
        tmp_size = 0
        n = a%nrow
        if (a%nrow /= a%ncol .or. size(perm) /= n .or. size(inv_perm) /= n .or. &
                size(supernode) /= n .or. size(xlnz) /= n + 1 .or. size(xsuper) < nsuper + 1 .or. &
                size(xlindx) < nsuper + 1) then
            allocate(lnz(0), split(0))
            info = -1
            return
        end if
        if (n == 0) then
            allocate(lnz(0), split(0))
            return
        end if
        if (xlnz(1) /= 1 .or. xlnz(n + 1) < 1) then
            allocate(lnz(0), split(0))
            info = -1
            return
        end if

        numeric_size = xlnz(n + 1) - 1
        allocate(lnz(numeric_size), split(n), offset(n))
        lnz = 0.0_dp
        split = 0
        offset = 0
        call input_numeric_values(a, perm, inv_perm, nsuper, xsuper, xlindx, lindx, xlnz, lnz, offset)

        cache_size = 1024
        if (present(cache_kb)) cache_size = cache_kb
        tiny_value = 1.0e-30_dp
        if (present(tiny)) tiny_value = tiny
        large_value = 1.0e128_dp
        if (present(large)) large_value = large
        unroll_level = 8
        if (present(level)) unroll_level = level
        if (all(unroll_level /= [1, 2, 4, 8])) then
            info = -4
            return
        end if

        call factor_workspace_size(nsuper, xsuper, supernode, xlindx, lindx, tmp_size)
        if (present(tmpmax)) then
            if (tmp_size > tmpmax) then
                info = 8
                return
            end if
        end if
        call split_supernodes(n, nsuper, xsuper, xlindx, cache_size, split)
        allocate(link(max(1, nsuper)), length(max(1, nsuper)))
        allocate(index_map(n), relative_index(n), temp(max(1, tmp_size)))
        call block_factor(nsuper, xsuper, supernode, split, xlindx, lindx, xlnz, lnz, link, length, &
            index_map, relative_index, tmp_size, temp, info, tiny_value, large_value, unroll_level)
    end subroutine ng_numeric_factor

    subroutine input_numeric_values(a, perm, inv_perm, nsuper, xsuper, xlindx, lindx, xlnz, lnz, offset)
        type(csr_matrix), intent(in) :: a !! Symmetric input matrix in original CSR numbering.
        integer, intent(in) :: perm(:) !! New-to-old factor ordering used to choose each original matrix row.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse ordering used to map column indices.
        integer, intent(in) :: nsuper !! Number of supernodes represented by xsuper and xlindx.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers.
        integer, intent(in) :: xlindx(:) !! Pointers to each supernode row-index list in lindx.
        integer, intent(in) :: lindx(:) !! Compressed row-index lists for all supernodes.
        integer, intent(in) :: xlnz(:) !! Column pointers into lnz.
        real(dp), intent(inout) :: lnz(:) !! Numerical storage initialized from the matrix lower triangle.
        integer, intent(inout) :: offset(:) !! Workspace mapping ordered row indices to bottom-relative offsets.
        integer :: i
        integer :: ii
        integer :: j
        integer :: jlen
        integer :: jsuper
        integer :: last
        integer :: old_j

        do jsuper = 1, nsuper
            jlen = xlindx(jsuper + 1) - xlindx(jsuper)
            do ii = xlindx(jsuper), xlindx(jsuper + 1) - 1
                i = lindx(ii)
                jlen = jlen - 1
                offset(i) = jlen
            end do

            do j = xsuper(jsuper), xsuper(jsuper + 1) - 1
                lnz(xlnz(j):xlnz(j + 1) - 1) = 0.0_dp
                old_j = perm(j)
                last = xlnz(j + 1) - 1
                do ii = a%row_ptr(old_j), a%row_ptr(old_j + 1) - 1
                    i = inv_perm(a%col_ind(ii))
                    if (i >= j) lnz(last - offset(i)) = a%values(ii)
                end do
            end do
        end do
    end subroutine input_numeric_values

    subroutine factor_workspace_size(nsuper, xsuper, supernode, xlindx, lindx, tmp_size)
        integer, intent(in) :: nsuper !! Number of fundamental supernodes.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers.
        integer, intent(in) :: supernode(:) !! Supernode membership for every factor row index.
        integer, intent(in) :: xlindx(:) !! Pointers into the compressed row-index vector.
        integer, intent(in) :: lindx(:) !! Compressed row indices for the supernodes.
        integer, intent(out) :: tmp_size !! Required temporary real workspace size for block updates.
        integer :: bound
        integer :: child_length
        integer :: current_super
        integer :: i
        integer :: i_begin
        integer :: i_end
        integer :: k_super
        integer :: length_active
        integer :: ncols
        integer :: next_super
        integer :: tsize
        integer :: width

        tmp_size = 0
        do k_super = nsuper, 1, -1
            ncols = xsuper(k_super + 1) - xsuper(k_super)
            i_begin = xlindx(k_super) + ncols
            i_end = xlindx(k_super + 1) - 1
            length_active = i_end - i_begin + 1
            bound = length_active * (length_active + 1) / 2
            if (bound > tmp_size) then
                current_super = supernode(lindx(i_begin))
                child_length = xlindx(current_super + 1) - xlindx(current_super)
                width = 0
                do i = i_begin, i_end
                    next_super = supernode(lindx(i))
                    if (next_super == current_super) then
                        width = width + 1
                        if (i == i_end) then
                            if (child_length > length_active) then
                                tsize = length_active * width - (width - 1) * width / 2
                                tmp_size = max(tsize, tmp_size)
                            end if
                        end if
                    else
                        if (child_length > length_active) then
                            tsize = length_active * width - (width - 1) * width / 2
                            tmp_size = max(tsize, tmp_size)
                        end if
                        length_active = length_active - width
                        bound = length_active * (length_active + 1) / 2
                        if (bound <= tmp_size) exit
                        width = 1
                        current_super = next_super
                        child_length = xlindx(current_super + 1) - xlindx(current_super)
                    end if
                end do
            end if
        end do
    end subroutine factor_workspace_size

    subroutine split_supernodes(neqns, nsuper, xsuper, xlindx, cache_kb, split)
        integer, intent(in) :: neqns !! Number of factor columns represented by split.
        integer, intent(in) :: nsuper !! Number of supernodes to partition into cache-sized blocks.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers.
        integer, intent(in) :: xlindx(:) !! Row-index pointers whose differences give supernode heights.
        integer, intent(in) :: cache_kb !! Nominal cache size in KiB; nonpositive means effectively unlimited.
        integer, intent(out) :: split(:) !! Block widths stored consecutively from each supernode first column.
        integer :: cache_words
        integer :: current_column
        integer :: first_column
        integer :: height
        integer :: k_column
        integer :: k_super
        integer :: last_column
        integer :: ncols
        integer :: next_block
        integer :: used

        if (cache_kb <= 0) then
            cache_words = 2000000000
        else
            cache_words = int(real(cache_kb, dp) * 1024.0_dp / 8.0_dp * 0.9_dp)
        end if
        do k_column = 1, neqns
            split(k_column) = 0
        end do

        do k_super = 1, nsuper
            height = xlindx(k_super + 1) - xlindx(k_super)
            first_column = xsuper(k_super)
            last_column = xsuper(k_super + 1) - 1
            next_block = first_column
            current_column = first_column - 1
            do
                current_column = current_column + 1
                if (current_column < last_column) then
                    current_column = current_column + 1
                    ncols = 2
                    used = 4 * height - 1
                    height = height - 2
                else
                    ncols = 1
                    used = 3 * height
                    height = height - 1
                end if
                do while (used + height < cache_words .and. current_column < last_column)
                    current_column = current_column + 1
                    ncols = ncols + 1
                    used = used + height
                    height = height - 1
                end do
                split(next_block) = ncols
                next_block = next_block + 1
                if (current_column >= last_column) exit
            end do
        end do
    end subroutine split_supernodes

    subroutine block_factor(nsuper, xsuper, supernode, split, xlindx, lindx, xlnz, lnz, link, length, &
            index_map, relative_index, tmp_size, temp, info, tiny, large, level)
        integer, intent(in) :: nsuper !! Number of supernodes in the factorization.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers.
        integer, intent(in) :: supernode(:) !! Supernode membership of every factor row index.
        integer, intent(in) :: split(:) !! Cache-block widths for every supernode.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode.
        integer, intent(in) :: lindx(:) !! Compressed row indices for all supernodes.
        integer, intent(in) :: xlnz(:) !! Column pointers into numerical factor storage.
        real(dp), intent(inout) :: lnz(:) !! Numerical matrix storage overwritten by its Cholesky factor.
        integer, intent(inout) :: link(:) !! Workspace linking supernodes that update a target supernode.
        integer, intent(inout) :: length(:) !! Workspace holding active lengths of updating supernodes.
        integer, intent(inout) :: index_map(:) !! Workspace mapping global row indices to target-relative offsets.
        integer, intent(inout) :: relative_index(:) !! Workspace of gathered target-relative row offsets.
        integer, intent(in) :: tmp_size !! Available length of temporary real workspace.
        real(dp), intent(inout) :: temp(:) !! Temporary real workspace used for sparse block updates.
        integer, intent(out) :: info !! Zero on success; -2 for insufficient temp; 16+k for k tiny pivots replaced.
        real(dp), intent(in) :: tiny !! Relative tiny-pivot threshold applied to the maximum input diagonal.
        real(dp), intent(in) :: large !! Replacement matrix diagonal value used for tiny pivots.
        integer, intent(in) :: level !! Numerical loop-unrolling level: 1, 2, 4, or 8.
        real(dp) :: max_diagonal
        integer :: first_j_column
        integer :: first_k_column
        integer :: i
        integer :: index_difference
        integer :: i_length
        integer :: i_pointer
        integer :: j_length
        integer :: j_lower_pointer
        integer :: j_super
        integer :: j_index_pointer
        integer :: k_first
        integer :: k_last
        integer :: k_length
        integer :: k_lower_pointer
        integer :: k_super
        integer :: k_index_pointer
        integer :: last_j_column
        integer :: n_columns_update
        integer :: n_j_columns
        integer :: n_k_columns
        integer :: next_k_super
        integer :: next_column
        integer :: next_super
        integer :: store
        integer :: tiny_count

        info = 0
        tiny_count = 0
        next_column = -1
        link(1:nsuper) = 0
        if (tmp_size > 0) temp(1:tmp_size) = 0.0_dp

        max_diagonal = 0.0_dp
        do i = 1, xsuper(nsuper + 1) - 1
            first_j_column = xlnz(i)
            max_diagonal = max(max_diagonal, lnz(first_j_column))
        end do

        do j_super = 1, nsuper
            first_j_column = xsuper(j_super)
            n_j_columns = xsuper(j_super + 1) - first_j_column
            last_j_column = first_j_column + n_j_columns - 1
            j_length = xlnz(first_j_column + 1) - xlnz(first_j_column)
            j_index_pointer = xlindx(j_super)

            call load_index_map(j_length, lindx(j_index_pointer:), index_map)
            k_super = link(j_super)
            do while (k_super > 0)
                next_k_super = link(k_super)
                first_k_column = xsuper(k_super)
                n_k_columns = xsuper(k_super + 1) - first_k_column
                k_length = length(k_super)
                k_index_pointer = xlindx(k_super + 1) - k_length

                if (k_length /= j_length) then
                    do i = 0, k_length - 1
                        next_column = lindx(k_index_pointer + i)
                        if (next_column > last_j_column) exit
                    end do
                    n_columns_update = i

                    if (n_k_columns == 1) then
                        k_lower_pointer = xlnz(first_k_column + 1) - k_length
                        call multiply_indexed(k_length, n_columns_update, lindx(k_index_pointer:), &
                            lnz(k_lower_pointer:), xlnz, lnz, index_map)
                    else
                        k_first = lindx(k_index_pointer)
                        k_last = lindx(k_index_pointer + k_length - 1)
                        index_difference = index_map(k_first) - index_map(k_last)
                        if (index_difference < k_length) then
                            i_pointer = xlnz(k_first)
                            i_length = xlnz(k_first + 1) - i_pointer
                            call matrix_multiply(k_length, n_k_columns, n_columns_update, split(first_k_column:), &
                                xlnz(first_k_column:), lnz, lnz(i_pointer:), i_length, level)
                        else
                            store = k_length * n_columns_update - n_columns_update * (n_columns_update - 1) / 2
                            if (store > tmp_size) then
                                info = -2
                                return
                            end if
                            call matrix_multiply(k_length, n_k_columns, n_columns_update, split(first_k_column:), &
                                xlnz(first_k_column:), lnz, temp, k_length, level)
                            call gather_indices(k_length, lindx(k_index_pointer:), index_map, relative_index)
                            call assemble_update(k_length, n_columns_update, temp, relative_index, &
                                xlnz(first_j_column:), lnz, j_length)
                        end if
                    end if
                else
                    j_lower_pointer = xlnz(first_j_column)
                    call matrix_multiply(k_length, n_k_columns, n_j_columns, split(first_k_column:), &
                        xlnz(first_k_column:), lnz, lnz(j_lower_pointer:), j_length, level)
                    n_columns_update = n_j_columns
                    if (k_length > n_j_columns) next_column = lindx(j_index_pointer + n_j_columns)
                end if

                if (k_length > n_columns_update) then
                    next_super = supernode(next_column)
                    link(k_super) = link(next_super)
                    link(next_super) = k_super
                    length(k_super) = k_length - n_columns_update
                else
                    length(k_super) = 0
                end if
                k_super = next_k_super
            end do

            call cholesky_supernode(j_length, n_j_columns, split(first_j_column:), xlnz(first_j_column:), &
                lnz, max_diagonal, tiny_count, info, tiny, large, level)
            if (info /= 0) then
                info = -1
                return
            end if

            if (j_length > n_j_columns) then
                next_column = lindx(j_index_pointer + n_j_columns)
                next_super = supernode(next_column)
                link(j_super) = link(next_super)
                link(next_super) = j_super
                length(j_super) = j_length - n_j_columns
            else
                length(j_super) = 0
            end if
        end do

        if (tiny_count /= 0) info = 16 + tiny_count
    end subroutine block_factor

    subroutine load_index_map(jlen, lindx, index_map)
        integer, intent(in) :: jlen !! Number of row indices in the target supernode structure.
        integer, intent(in) :: lindx(:) !! Target supernode row-index list, beginning at its first entry.
        integer, intent(inout) :: index_map(:) !! Global-row map receiving bottom-relative offsets.
        integer :: current_length
        integer :: j
        integer :: row

        current_length = jlen
        do j = 1, jlen
            row = lindx(j)
            current_length = current_length - 1
            index_map(row) = current_length
        end do
    end subroutine load_index_map

    subroutine gather_indices(klen, lindx, index_map, relative_index)
        integer, intent(in) :: klen !! Number of active indices to gather.
        integer, intent(in) :: lindx(:) !! Global row indices of the active updating structure.
        integer, intent(in) :: index_map(:) !! Map from global row index to target-relative offset.
        integer, intent(out) :: relative_index(:) !! Gathered relative offsets for the first klen indices.
        integer :: i

        do i = 1, klen
            relative_index(i) = index_map(lindx(i))
        end do
    end subroutine gather_indices

    subroutine assemble_update(m, q, y, relative_index, xlnz, lnz, lda)
        integer, intent(in) :: m !! Number of rows represented in the temporary update block.
        integer, intent(in) :: q !! Number of target columns updated from the temporary block.
        real(dp), intent(inout) :: y(:) !! Packed update block, zeroed as its entries are assembled.
        integer, intent(in) :: relative_index(:) !! Bottom-relative target positions for update rows.
        integer, intent(in) :: xlnz(:) !! Column pointers beginning at the target first column.
        real(dp), intent(inout) :: lnz(:) !! Global numerical factor storage receiving the update.
        integer, intent(in) :: lda !! Length of the target first column.
        integer :: i_column
        integer :: i_lnz
        integer :: i_row
        integer :: i_y
        integer :: l_bottom
        integer :: y_column
        integer :: y_offset

        y_offset = 0
        do i_column = 1, q
            y_column = lda - relative_index(i_column)
            l_bottom = xlnz(y_column + 1) - 1
            do i_row = i_column, m
                i_lnz = l_bottom - relative_index(i_row)
                i_y = y_offset + i_row
                lnz(i_lnz) = lnz(i_lnz) + y(i_y)
                y(i_y) = 0.0_dp
            end do
            y_offset = i_y - i_column
        end do
    end subroutine assemble_update

    subroutine multiply_indexed(m, q, row_index, x, xlnz, y, relative_index)
        integer, intent(in) :: m !! Length of the active updating column.
        integer, intent(in) :: q !! Number of target columns updated by the one-column supernode.
        integer, intent(in) :: row_index(:) !! Active global row indices of the updating column.
        real(dp), intent(in) :: x(:) !! Active numerical values of the updating column.
        integer, intent(in) :: xlnz(:) !! Global column pointers for the target numerical storage.
        real(dp), intent(inout) :: y(:) !! Global numerical storage receiving the indexed outer-product update.
        integer, intent(in) :: relative_index(:) !! Global-row map to bottom-relative target offsets.
        real(dp) :: multiplier
        integer :: column
        integer :: i
        integer :: index
        integer :: k
        integer :: y_last

        do k = 1, q
            column = row_index(k)
            y_last = xlnz(column + 1) - 1
            multiplier = -x(k)
            do i = k, m
                index = row_index(i)
                index = y_last - relative_index(index)
                y(index) = y(index) + multiplier * x(i)
            end do
        end do
    end subroutine multiply_indexed

    subroutine matrix_multiply(m, n, q, split, xpnt, x, y, ldy, level)
        integer, intent(in) :: m !! Number of active rows in the updating block.
        integer, intent(in) :: n !! Number of updating columns in the source supernode.
        integer, intent(in) :: q !! Number of target columns updated.
        integer, intent(in) :: split(:) !! Cache-block widths beginning at the source supernode first column.
        integer, intent(in) :: xpnt(:) !! Absolute numerical column pointers beginning at the source first column.
        real(dp), intent(in) :: x(:) !! Global numerical supernodal storage containing the source block.
        real(dp), intent(inout) :: y(:) !! Target block receiving the negative product update.
        integer, intent(in) :: ldy !! Physical length of the first target column.
        integer, intent(in) :: level !! Numerical loop-unrolling level: 1, 2, 4, or 8.
        integer :: block
        integer :: first_column
        integer :: nn

        block = 1
        first_column = 1
        do while (first_column <= n)
            nn = split(block)
            if (level == 1) then
                call matrix_multiply_level1(m, nn, q, xpnt(first_column:), x, y, ldy)
            else
                call matrix_multiply_paired(m, nn, q, xpnt(first_column:), x, y, ldy, level)
            end if
            first_column = first_column + nn
            block = block + 1
        end do
    end subroutine matrix_multiply

    subroutine matrix_multiply_level1(m, n, q, xpnt, x, y, ldy)
        integer, intent(in) :: m !! Number of active rows in source and target blocks.
        integer, intent(in) :: n !! Number of source columns in this multiplication block.
        integer, intent(in) :: q !! Number of target columns to update.
        integer, intent(in) :: xpnt(:) !! Absolute pointers to source numerical columns.
        real(dp), intent(in) :: x(:) !! Global source numerical storage.
        real(dp), intent(inout) :: y(:) !! Packed target block receiving the update.
        integer, intent(in) :: ldy !! Physical length of the first packed target column.
        real(dp) :: multiplier
        integer :: i1
        integer :: i_y
        integer :: i_y_last
        integer :: i_y_start
        integer :: i_y_stop
        integer :: length_y
        integer :: mm
        integer :: x_column
        integer :: y_column

        mm = m
        i_y_last = 0
        length_y = ldy
        do y_column = 1, q
            i_y_start = i_y_last + 1
            i_y_stop = i_y_start + mm - 1
            i_y_last = i_y_last + length_y
            do x_column = 1, n
                i1 = xpnt(x_column + 1) - mm
                multiplier = -x(i1)
                do i_y = i_y_start, i_y_stop
                    y(i_y) = y(i_y) + multiplier * x(i1)
                    i1 = i1 + 1
                end do
            end do
            mm = mm - 1
            length_y = length_y - 1
        end do
    end subroutine matrix_multiply_level1


    subroutine matrix_multiply_paired(m, n, q, xpnt, x, y, ldy, level)
        integer, intent(in) :: m !! Number of active rows in source and target blocks.
        integer, intent(in) :: n !! Number of source columns in this multiplication block.
        integer, intent(in) :: q !! Number of packed target columns to update.
        integer, intent(in) :: xpnt(:) !! Absolute pointers to source numerical columns.
        real(dp), intent(in) :: x(:) !! Global source numerical storage.
        real(dp), intent(inout) :: y(:) !! Packed triangular target block receiving the negative update.
        integer, intent(in) :: ldy !! Physical length of the first packed target column.
        integer, intent(in) :: level !! Paired-target unrolling level: 2, 4, or 8.

        select case (level)
        case (2)
            call matrix_multiply_paired_level2(m, n, q, xpnt, x, y, ldy)
        case (4)
            call matrix_multiply_paired_level4(m, n, q, xpnt, x, y, ldy)
        case (8)
            call matrix_multiply_paired_level8(m, n, q, xpnt, x, y, ldy)
        end select
    end subroutine matrix_multiply_paired

    subroutine matrix_multiply_paired_level2(m, n, q, xpnt, x, y, ldy)
        integer, intent(in) :: m !! Number of active rows in source and target blocks.
        integer, intent(in) :: n !! Number of source columns in this multiplication block.
        integer, intent(in) :: q !! Number of packed target columns to update.
        integer, intent(in) :: xpnt(:) !! Absolute pointers to source numerical columns.
        real(dp), intent(in) :: x(:) !! Global source numerical storage shared by both target columns.
        real(dp), intent(inout) :: y(:) !! Packed triangular target block receiving the negative update.
        integer, intent(in) :: ldy !! Physical length of the first packed target column.
        integer :: i
        integer :: mm
        integer :: source_column
        integer :: target_column
        integer :: target_first
        integer :: target_second
        integer :: p1
        integer :: p2
        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: b1
        real(dp) :: b2
        real(dp) :: v1
        real(dp) :: v2

        if (q <= 0 .or. m <= 0 .or. n <= 0) return
        mm = m
        target_first = 1
        target_column = 1
        do while (target_column + 1 <= q .and. mm >= 2)
            target_second = target_first + ldy - target_column + 1
            source_column = 1
            do while (source_column + 1 <= n)
                p1 = xpnt(source_column + 1) - mm
                p2 = xpnt(source_column + 2) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                a2 = -x(p2)
                b2 = -x(p2 + 1)
                y(target_first) = ((y(target_first) + a1 * x(p1)) + a2 * x(p2))
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    p2 = p2 + 1
                    v2 = x(p2)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1 + a2 * v2
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1 + b2 * v2
                end do
                source_column = source_column + 2
            end do
            do while (source_column <= n)
                p1 = xpnt(source_column + 1) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                y(target_first) = y(target_first) + a1 * x(p1)
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1
                end do
                source_column = source_column + 1
            end do
            target_first = target_second + mm - 1
            target_column = target_column + 2
            mm = mm - 2
        end do
        if (target_column <= q .and. mm > 0) then
            call matrix_vector_update_unrolled(mm, n, y(target_first:), xpnt, x, 2)
        end if
    end subroutine matrix_multiply_paired_level2

    subroutine matrix_multiply_paired_level4(m, n, q, xpnt, x, y, ldy)
        integer, intent(in) :: m !! Number of active rows in source and target blocks.
        integer, intent(in) :: n !! Number of source columns in this multiplication block.
        integer, intent(in) :: q !! Number of packed target columns to update.
        integer, intent(in) :: xpnt(:) !! Absolute pointers to source numerical columns.
        real(dp), intent(in) :: x(:) !! Global source numerical storage shared by both target columns.
        real(dp), intent(inout) :: y(:) !! Packed triangular target block receiving the negative update.
        integer, intent(in) :: ldy !! Physical length of the first packed target column.
        integer :: i
        integer :: mm
        integer :: source_column
        integer :: target_column
        integer :: target_first
        integer :: target_second
        integer :: p1
        integer :: p2
        integer :: p3
        integer :: p4
        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: a3
        real(dp) :: a4
        real(dp) :: b1
        real(dp) :: b2
        real(dp) :: b3
        real(dp) :: b4
        real(dp) :: v1
        real(dp) :: v2
        real(dp) :: v3
        real(dp) :: v4

        if (q <= 0 .or. m <= 0 .or. n <= 0) return
        mm = m
        target_first = 1
        target_column = 1
        do while (target_column + 1 <= q .and. mm >= 2)
            target_second = target_first + ldy - target_column + 1
            source_column = 1
            do while (source_column + 3 <= n)
                p1 = xpnt(source_column + 1) - mm
                p2 = xpnt(source_column + 2) - mm
                p3 = xpnt(source_column + 3) - mm
                p4 = xpnt(source_column + 4) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                a2 = -x(p2)
                b2 = -x(p2 + 1)
                a3 = -x(p3)
                b3 = -x(p3 + 1)
                a4 = -x(p4)
                b4 = -x(p4 + 1)
                y(target_first) = ((((y(target_first) + a1 * x(p1)) + a2 * x(p2)) + a3 * x(p3)) + a4 * x(p4))
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    p2 = p2 + 1
                    v2 = x(p2)
                    p3 = p3 + 1
                    v3 = x(p3)
                    p4 = p4 + 1
                    v4 = x(p4)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1 + a2 * v2 + a3 * v3 + a4 * v4
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1 + b2 * v2 + b3 * v3 + b4 * v4
                end do
                source_column = source_column + 4
            end do
            do while (source_column <= n)
                p1 = xpnt(source_column + 1) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                y(target_first) = y(target_first) + a1 * x(p1)
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1
                end do
                source_column = source_column + 1
            end do
            target_first = target_second + mm - 1
            target_column = target_column + 2
            mm = mm - 2
        end do
        if (target_column <= q .and. mm > 0) then
            call matrix_vector_update_unrolled(mm, n, y(target_first:), xpnt, x, 4)
        end if
    end subroutine matrix_multiply_paired_level4

    subroutine matrix_multiply_paired_level8(m, n, q, xpnt, x, y, ldy)
        integer, intent(in) :: m !! Number of active rows in source and target blocks.
        integer, intent(in) :: n !! Number of source columns in this multiplication block.
        integer, intent(in) :: q !! Number of packed target columns to update.
        integer, intent(in) :: xpnt(:) !! Absolute pointers to source numerical columns.
        real(dp), intent(in) :: x(:) !! Global source numerical storage shared by both target columns.
        real(dp), intent(inout) :: y(:) !! Packed triangular target block receiving the negative update.
        integer, intent(in) :: ldy !! Physical length of the first packed target column.
        integer :: i
        integer :: mm
        integer :: source_column
        integer :: target_column
        integer :: target_first
        integer :: target_second
        integer :: p1
        integer :: p2
        integer :: p3
        integer :: p4
        integer :: p5
        integer :: p6
        integer :: p7
        integer :: p8
        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: a3
        real(dp) :: a4
        real(dp) :: a5
        real(dp) :: a6
        real(dp) :: a7
        real(dp) :: a8
        real(dp) :: b1
        real(dp) :: b2
        real(dp) :: b3
        real(dp) :: b4
        real(dp) :: b5
        real(dp) :: b6
        real(dp) :: b7
        real(dp) :: b8
        real(dp) :: v1
        real(dp) :: v2
        real(dp) :: v3
        real(dp) :: v4
        real(dp) :: v5
        real(dp) :: v6
        real(dp) :: v7
        real(dp) :: v8

        if (q <= 0 .or. m <= 0 .or. n <= 0) return
        mm = m
        target_first = 1
        target_column = 1
        do while (target_column + 1 <= q .and. mm >= 2)
            target_second = target_first + ldy - target_column + 1
            source_column = 1
            do while (source_column + 7 <= n)
                p1 = xpnt(source_column + 1) - mm
                p2 = xpnt(source_column + 2) - mm
                p3 = xpnt(source_column + 3) - mm
                p4 = xpnt(source_column + 4) - mm
                p5 = xpnt(source_column + 5) - mm
                p6 = xpnt(source_column + 6) - mm
                p7 = xpnt(source_column + 7) - mm
                p8 = xpnt(source_column + 8) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                a2 = -x(p2)
                b2 = -x(p2 + 1)
                a3 = -x(p3)
                b3 = -x(p3 + 1)
                a4 = -x(p4)
                b4 = -x(p4 + 1)
                a5 = -x(p5)
                b5 = -x(p5 + 1)
                a6 = -x(p6)
                b6 = -x(p6 + 1)
                a7 = -x(p7)
                b7 = -x(p7 + 1)
                a8 = -x(p8)
                b8 = -x(p8 + 1)
                y(target_first) = y(target_first) + a1 * x(p1) + a2 * x(p2) + a3 * x(p3) + a4 * x(p4) &
                    + a5 * x(p5) + a6 * x(p6) + a7 * x(p7) + a8 * x(p8)
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    p2 = p2 + 1
                    v2 = x(p2)
                    p3 = p3 + 1
                    v3 = x(p3)
                    p4 = p4 + 1
                    v4 = x(p4)
                    p5 = p5 + 1
                    v5 = x(p5)
                    p6 = p6 + 1
                    v6 = x(p6)
                    p7 = p7 + 1
                    v7 = x(p7)
                    p8 = p8 + 1
                    v8 = x(p8)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1 + a2 * v2 + a3 * v3 + a4 * v4 &
                        + a5 * v5 + a6 * v6 + a7 * v7 + a8 * v8
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1 + b2 * v2 + b3 * v3 + b4 * v4 &
                        + b5 * v5 + b6 * v6 + b7 * v7 + b8 * v8
                end do
                source_column = source_column + 8
            end do
            do while (source_column <= n)
                p1 = xpnt(source_column + 1) - mm
                a1 = -x(p1)
                b1 = -x(p1 + 1)
                y(target_first) = y(target_first) + a1 * x(p1)
                do i = 2, mm
                    p1 = p1 + 1
                    v1 = x(p1)
                    y(target_first + i - 1) = y(target_first + i - 1) + a1 * v1
                    y(target_second + i - 2) = y(target_second + i - 2) + b1 * v1
                end do
                source_column = source_column + 1
            end do
            target_first = target_second + mm - 1
            target_column = target_column + 2
            mm = mm - 2
        end do
        if (target_column <= q .and. mm > 0) then
            call matrix_vector_update_unrolled(mm, n, y(target_first:), xpnt, x, 8)
        end if
    end subroutine matrix_multiply_paired_level8

    subroutine cholesky_supernode(m, n, split, xpnt, x, max_diagonal, tiny_count, info, tiny, large, level)
        integer, intent(in) :: m !! Number of rows in the first column of the current supernode.
        integer, intent(in) :: n !! Number of columns in the current supernode.
        integer, intent(in) :: split(:) !! Cache-block widths beginning at this supernode first column.
        integer, intent(in) :: xpnt(:) !! Absolute numerical pointers beginning at this supernode first column.
        real(dp), intent(inout) :: x(:) !! Global numerical storage containing the current supernode.
        real(dp), intent(in) :: max_diagonal !! Maximum input diagonal used to scale the tiny-pivot threshold.
        integer, intent(inout) :: tiny_count !! Running count of tiny pivots replaced by the large value.
        integer, intent(out) :: info !! Zero unless partial Cholesky reports a failure.
        real(dp), intent(in) :: tiny !! Relative tiny-pivot threshold.
        real(dp), intent(in) :: large !! Replacement matrix diagonal value used for tiny pivots.
        integer, intent(in) :: level !! Numerical loop-unrolling level: 1, 2, 4, or 8.
        integer :: first_column
        integer :: j_block
        integer :: mm
        integer :: nn
        integer :: next_column
        integer :: q

        info = 0
        j_block = 0
        first_column = 1
        mm = m
        do while (first_column <= n)
            j_block = j_block + 1
            nn = split(j_block)
            call partial_cholesky(mm, nn, xpnt(first_column:), x, max_diagonal, tiny_count, info, tiny, large, level)
            if (info /= 0) return
            next_column = first_column + nn
            q = n - next_column + 1
            mm = mm - nn
            if (q > 0) then
                if (level == 1) then
                    call matrix_multiply_level1(mm, nn, q, xpnt(first_column:), x, x(xpnt(next_column):), mm)
                else
                    call matrix_multiply_paired(mm, nn, q, xpnt(first_column:), x, x(xpnt(next_column):), mm, level)
                end if
            end if
            first_column = next_column
        end do
    end subroutine cholesky_supernode

    subroutine partial_cholesky(m, n, xpnt, x, max_diagonal, tiny_count, info, tiny, large, level)
        integer, intent(in) :: m !! Number of active rows in the first column of the block.
        integer, intent(in) :: n !! Number of leading block columns to factor.
        integer, intent(in) :: xpnt(:) !! Absolute numerical column pointers for the block columns.
        real(dp), intent(inout) :: x(:) !! Global numerical supernodal storage containing the block.
        real(dp), intent(in) :: max_diagonal !! Maximum input diagonal used for tiny-pivot detection.
        integer, intent(inout) :: tiny_count !! Running number of tiny pivots replaced by the large value.
        integer, intent(out) :: info !! Zero on completion; reserved for numerical failure compatibility.
        real(dp), intent(in) :: tiny !! Relative tiny-pivot threshold.
        real(dp), intent(in) :: large !! Replacement matrix diagonal value used before square rooting a tiny pivot.
        integer, intent(in) :: level !! Numerical loop-unrolling level: 1, 2, 4, or 8.
        real(dp) :: diagonal
        integer :: j_column
        integer :: j_pointer
        integer :: mm

        info = 0
        mm = m
        j_pointer = xpnt(1)
        do j_column = 1, n
            if (j_column > 1) then
                if (level == 1) then
                    call matrix_vector_update(mm, j_column - 1, x(j_pointer:), xpnt, x)
                else
                    call matrix_vector_update_unrolled(mm, j_column - 1, x(j_pointer:), xpnt, x, level)
                end if
            end if
            diagonal = x(j_pointer)
            if (diagonal <= tiny * max_diagonal) then
                diagonal = large
                tiny_count = tiny_count + 1
            end if
            if (diagonal < 0.0_dp) then
                info = 1
                return
            end if
            diagonal = sqrt(diagonal)
            x(j_pointer) = diagonal
            diagonal = 1.0_dp / diagonal
            mm = mm - 1
            j_pointer = j_pointer + 1
            if (mm > 0) call scale_vector(mm, diagonal, x(j_pointer:))
            j_pointer = j_pointer + mm
        end do
    end subroutine partial_cholesky

    subroutine matrix_vector_update(m, n, y, apnt, a)
        integer, intent(in) :: m !! Number of active rows in the target vector.
        integer, intent(in) :: n !! Number of already-factorized columns contributing the update.
        real(dp), intent(inout) :: y(:) !! Target vector receiving the negative matrix-vector update.
        integer, intent(in) :: apnt(:) !! Absolute numerical column pointers for the source columns.
        real(dp), intent(in) :: a(:) !! Global numerical supernodal storage containing source coefficients.
        real(dp) :: multiplier
        integer :: i
        integer :: ii
        integer :: j

        do j = 1, n
            ii = apnt(j + 1) - m
            multiplier = -a(ii)
            do i = 1, m
                y(i) = y(i) + multiplier * a(ii)
                ii = ii + 1
            end do
        end do
    end subroutine matrix_vector_update


    subroutine matrix_vector_update_unrolled(m, n, y, apnt, a, level)
        integer, intent(in) :: m !! Number of active rows in the target vector.
        integer, intent(in) :: n !! Number of already-factorized columns contributing the update.
        real(dp), intent(inout) :: y(:) !! Target vector receiving the negative matrix-vector update.
        integer, intent(in) :: apnt(:) !! Absolute numerical column pointers for the source columns.
        real(dp), intent(in) :: a(:) !! Global numerical supernodal storage containing source coefficients.
        integer, intent(in) :: level !! Explicit source-column unrolling level: 2, 4, or 8.
        integer :: i
        integer :: p1
        integer :: p2
        integer :: p3
        integer :: p4
        integer :: p5
        integer :: p6
        integer :: p7
        integer :: p8
        integer :: j
        integer :: remainder
        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: a3
        real(dp) :: a4
        real(dp) :: a5
        real(dp) :: a6
        real(dp) :: a7
        real(dp) :: a8

        remainder = mod(n, level)
        do j = 1, remainder
            p1 = apnt(j + 1) - m
            a1 = -a(p1)
            do i = 1, m
                y(i) = y(i) + a1 * a(p1)
                p1 = p1 + 1
            end do
        end do

        j = remainder + 1
        select case (level)
        case (2)
            do while (j + 1 <= n)
                p1 = apnt(j + 1) - m
                p2 = apnt(j + 2) - m
                a1 = -a(p1)
                a2 = -a(p2)
                do i = 1, m
                    y(i) = (y(i) + a1 * a(p1)) + a2 * a(p2)
                    p1 = p1 + 1
                    p2 = p2 + 1
                end do
                j = j + 2
            end do
        case (4)
            do while (j + 3 <= n)
                p1 = apnt(j + 1) - m
                p2 = apnt(j + 2) - m
                p3 = apnt(j + 3) - m
                p4 = apnt(j + 4) - m
                a1 = -a(p1)
                a2 = -a(p2)
                a3 = -a(p3)
                a4 = -a(p4)
                do i = 1, m
                    y(i) = (((y(i) + a1 * a(p1)) + a2 * a(p2)) + a3 * a(p3)) + a4 * a(p4)
                    p1 = p1 + 1
                    p2 = p2 + 1
                    p3 = p3 + 1
                    p4 = p4 + 1
                end do
                j = j + 4
            end do
        case (8)
            do while (j + 7 <= n)
                p1 = apnt(j + 1) - m
                p2 = apnt(j + 2) - m
                p3 = apnt(j + 3) - m
                p4 = apnt(j + 4) - m
                p5 = apnt(j + 5) - m
                p6 = apnt(j + 6) - m
                p7 = apnt(j + 7) - m
                p8 = apnt(j + 8) - m
                a1 = -a(p1)
                a2 = -a(p2)
                a3 = -a(p3)
                a4 = -a(p4)
                a5 = -a(p5)
                a6 = -a(p6)
                a7 = -a(p7)
                a8 = -a(p8)
                do i = 1, m
                    y(i) = (((((((y(i) + a1 * a(p1)) + a2 * a(p2)) + a3 * a(p3)) + a4 * a(p4)) + &
                        a5 * a(p5)) + a6 * a(p6)) + a7 * a(p7)) + a8 * a(p8)
                    p1 = p1 + 1
                    p2 = p2 + 1
                    p3 = p3 + 1
                    p4 = p4 + 1
                    p5 = p5 + 1
                    p6 = p6 + 1
                    p7 = p7 + 1
                    p8 = p8 + 1
                end do
                j = j + 8
            end do
        end select
    end subroutine matrix_vector_update_unrolled

    subroutine scale_vector(n, multiplier, x)
        integer, intent(in) :: n !! Number of leading entries of x to scale.
        real(dp), intent(in) :: multiplier !! Scalar multiplier applied to every selected vector entry.
        real(dp), intent(inout) :: x(:) !! Vector whose first n values are multiplied in place.
        integer :: i

        do i = 1, n
            x(i) = multiplier * x(i)
        end do
    end subroutine scale_vector


    subroutine ng_forwardsolve(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
        integer, intent(in) :: nsuper !! Number of supernodes in the stored Cholesky factor.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers with active extent nsuper+1.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode row-index structure.
        integer, intent(in) :: lindx(:) !! Compressed row indices for the supernodal factor.
        integer, intent(in) :: xlnz(:) !! Numerical column pointers for the supernodal factor.
        real(dp), intent(in) :: lnz(:) !! Numerical lower Cholesky factor in Ng-Peyton storage.
        real(dp), intent(inout) :: rhs(:) !! Ordered right-hand side overwritten by solution of L x = rhs.
        real(dp) :: value
        integer :: first_column
        integer :: i
        integer :: index_pointer
        integer :: ix
        integer :: ix_start
        integer :: ix_stop
        integer :: j_column
        integer :: j_pointer
        integer :: j_super
        integer :: last_column

        if (nsuper <= 0) return
        first_column = xsuper(1)
        do j_super = 1, nsuper
            last_column = xsuper(j_super + 1) - 1
            ix_start = xlnz(first_column)
            j_pointer = xlindx(j_super)
            do j_column = first_column, last_column
                ix_stop = xlnz(j_column + 1) - 1
                if (rhs(j_column) /= 0.0_dp) then
                    value = rhs(j_column) / lnz(ix_start)
                    rhs(j_column) = value
                    index_pointer = j_pointer + 1
                    do ix = ix_start + 1, ix_stop
                        i = lindx(index_pointer)
                        rhs(i) = rhs(i) - value * lnz(ix)
                        index_pointer = index_pointer + 1
                    end do
                end if
                ix_start = ix_stop + 1
                j_pointer = j_pointer + 1
            end do
            first_column = last_column + 1
        end do
    end subroutine ng_forwardsolve

    subroutine ng_backsolve(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
        integer, intent(in) :: nsuper !! Number of supernodes in the stored Cholesky factor.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers with active extent nsuper+1.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode row-index structure.
        integer, intent(in) :: lindx(:) !! Compressed row indices for the supernodal factor.
        integer, intent(in) :: xlnz(:) !! Numerical column pointers for the supernodal factor.
        real(dp), intent(in) :: lnz(:) !! Numerical lower Cholesky factor in Ng-Peyton storage.
        real(dp), intent(inout) :: rhs(:) !! Ordered right-hand side overwritten by solution of L-transpose x = rhs.
        real(dp) :: value
        integer :: first_column
        integer :: i
        integer :: index_pointer
        integer :: ix
        integer :: ix_start
        integer :: ix_stop
        integer :: j_column
        integer :: j_pointer
        integer :: j_super
        integer :: last_column

        if (nsuper <= 0) return
        last_column = xsuper(nsuper + 1) - 1
        do j_super = nsuper, 1, -1
            first_column = xsuper(j_super)
            ix_stop = xlnz(last_column + 1) - 1
            j_pointer = xlindx(j_super) + (last_column - first_column)
            do j_column = last_column, first_column, -1
                ix_start = xlnz(j_column)
                index_pointer = j_pointer + 1
                value = rhs(j_column)
                do ix = ix_start + 1, ix_stop
                    i = lindx(index_pointer)
                    if (rhs(i) /= 0.0_dp) value = value - lnz(ix) * rhs(i)
                    index_pointer = index_pointer + 1
                end do
                if (value /= 0.0_dp) then
                    rhs(j_column) = value / lnz(ix_start)
                else
                    rhs(j_column) = 0.0_dp
                end if
                ix_stop = ix_start - 1
                j_pointer = j_pointer - 1
            end do
            last_column = first_column - 1
        end do
    end subroutine ng_backsolve

    subroutine ng_solve(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
        integer, intent(in) :: nsuper !! Number of supernodes in the stored Cholesky factor.
        integer, intent(in) :: xsuper(:) !! Supernode first-column pointers with active extent nsuper+1.
        integer, intent(in) :: xlindx(:) !! Pointers into lindx for each supernode row-index structure.
        integer, intent(in) :: lindx(:) !! Compressed row indices for the supernodal factor.
        integer, intent(in) :: xlnz(:) !! Numerical column pointers for the supernodal factor.
        real(dp), intent(in) :: lnz(:) !! Numerical lower Cholesky factor in Ng-Peyton storage.
        real(dp), intent(inout) :: rhs(:) !! Ordered right-hand side overwritten by solution of L L-transpose x = rhs.

        call ng_forwardsolve(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
        call ng_backsolve(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
    end subroutine ng_solve

end module sparsem_ng_numeric
