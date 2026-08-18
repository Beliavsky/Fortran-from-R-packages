! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_sets
    use partitions_kinds, only : i8
    use partitions_counts, only : set_partition_count, multiset_permutation_count
    use partitions_enumeration, only : block_parts
    use partitions_permutations, only : multiset_permutations
    implicit none
    private

    public :: set_partitions
    public :: restricted_set_partitions
    public :: multiset_sequences
    public :: multinomial_permutations
    public :: all_binomial
    public :: generalized_riffles
    public :: riffles

contains

    function set_partitions(block_sizes) result(out)
        integer, intent(in) :: block_sizes(:)
        integer, allocatable :: out(:,:)
        integer, allocatable :: shape_sizes(:), shape(:,:), digits(:), labels(:)
        integer(i8) :: nc8
        integer :: n, nb, max_size, nc, ndig, col, i
        integer, parameter :: sentinel_margin = 1000

        if (size(block_sizes) == 0) then
            allocate(out(0,1))
            return
        end if
        if (any(block_sizes <= 0)) error stop "set_partitions: block sizes must be positive"
        shape_sizes = block_sizes
        call sort_descending(shape_sizes)
        n = sum(shape_sizes)
        nb = size(shape_sizes)
        nc8 = set_partition_count(shape_sizes)
        call count_to_default(nc8, nc, "set_partitions")
        allocate(out(n, nc))
        out = 0
        max_size = maxval(shape_sizes)
        allocate(shape(max_size, nb), digits(n), labels(n))
        shape = huge(0) - sentinel_margin
        do i = 1, nb
            call reset_block(shape(:,i), shape_sizes(i))
        end do
        digits = [(i, i = 1, n)]
        ndig = n
        col = 0
        call recurse(1, 1)
        if (col /= nc) error stop "set_partitions: enumeration count mismatch"

    contains

        recursive subroutine recurse(t, e)
            integer, intent(in) :: t, e
            integer :: idx, chosen

            if (t > nb) then
                call save_partition()
                return
            end if

            if (e > shape_sizes(t)) then
                if (t > 1) then
                    if (shape_sizes(t - 1) == shape_sizes(t)) then
                        if (.not. tuple_less(shape(:,t - 1), shape(:,t), shape_sizes(t))) return
                    end if
                end if
                call recurse(t + 1, 1)
                return
            end if

            do idx = 1, ndig
                call get_digit(idx, chosen)
                if (e == 1) then
                    shape(e,t) = chosen
                    call recurse(t, e + 1)
                    shape(e,t) = huge(0) - shape_sizes(t) + e
                else if (chosen > shape(e - 1,t)) then
                    shape(e,t) = chosen
                    call recurse(t, e + 1)
                    shape(e,t) = huge(0) - shape_sizes(t) + e
                end if
                call put_digit(idx, chosen)
            end do
        end subroutine recurse

        subroutine get_digit(idx, chosen)
            integer, intent(in) :: idx
            integer, intent(out) :: chosen

            chosen = digits(idx)
            digits(idx) = digits(ndig)
            ndig = ndig - 1
        end subroutine get_digit

        subroutine put_digit(idx, chosen)
            integer, intent(in) :: idx, chosen

            if (idx == ndig + 1) then
                ndig = ndig + 1
                digits(ndig) = chosen
            else
                ndig = ndig + 1
                digits(ndig) = digits(idx)
                digits(idx) = chosen
            end if
        end subroutine put_digit

        subroutine save_partition()
            integer :: b, e1

            labels = 0
            do b = 1, nb
                do e1 = 1, shape_sizes(b)
                    labels(shape(e1,b)) = b
                end do
            end do
            col = col + 1
            if (col > nc) error stop "set_partitions: too many columns"
            out(:,col) = labels
        end subroutine save_partition

    end function set_partitions

    function restricted_set_partitions(block_sizes) result(out)
        integer, intent(in) :: block_sizes(:)
        integer, allocatable :: out(:,:), labels(:,:)
        integer :: n, nc, col, b, pos, i

        labels = set_partitions(block_sizes)
        n = size(labels,1)
        nc = size(labels,2)
        allocate(out(n,nc))
        do col = 1, nc
            pos = 0
            do b = 1, maxval(labels(:,col))
                do i = 1, n
                    if (labels(i,col) == b) then
                        pos = pos + 1
                        out(pos,col) = i
                    end if
                end do
            end do
        end do
    end function restricted_set_partitions

    function multiset_sequences(v, n_select) result(out)
        integer, intent(in) :: v(:)
        integer, intent(in), optional :: n_select
        integer, allocatable :: out(:,:)
        integer, allocatable :: sorted(:), values(:), counts(:), allocs(:,:), expanded(:), perms(:,:)
        integer(i8) :: total8
        integer :: n, nu, i, j, k, nsel, total, col, pcol, pos

        n = size(v)
        nsel = n
        if (present(n_select)) nsel = n_select
        if (nsel < 0 .or. nsel > n) error stop "multiset_sequences: invalid selection size"
        if (nsel == 0) then
            allocate(out(0,1))
            return
        end if

        sorted = v
        call sort_ascending(sorted)
        nu = 1
        do i = 2, n
            if (sorted(i) /= sorted(i - 1)) nu = nu + 1
        end do
        allocate(values(nu), counts(nu))
        values = 0
        counts = 0
        k = 1
        values(1) = sorted(1)
        counts(1) = 1
        do i = 2, n
            if (sorted(i) == values(k)) then
                counts(k) = counts(k) + 1
            else
                k = k + 1
                values(k) = sorted(i)
                counts(k) = 1
            end if
        end do

        if (nsel == n) then
            out = multiset_permutations(sorted)
            return
        end if
        if (nsel == 1) then
            allocate(out(1,nu))
            out(1,:) = values
            return
        end if

        allocs = block_parts(counts, nsel)
        total8 = 0_i8
        do j = 1, size(allocs,2)
            allocate(expanded(nsel))
            pos = 0
            do i = 1, nu
                do k = 1, allocs(i,j)
                    pos = pos + 1
                    expanded(pos) = values(i)
                end do
            end do
            total8 = total8 + multiset_permutation_count(expanded)
            deallocate(expanded)
        end do
        call count_to_default(total8, total, "multiset_sequences")
        allocate(out(nsel,total))
        col = 0
        do j = 1, size(allocs,2)
            allocate(expanded(nsel))
            pos = 0
            do i = 1, nu
                do k = 1, allocs(i,j)
                    pos = pos + 1
                    expanded(pos) = values(i)
                end do
            end do
            perms = multiset_permutations(expanded)
            do pcol = 1, size(perms,2)
                col = col + 1
                out(:,col) = perms(:,pcol)
            end do
            deallocate(expanded)
        end do
        if (col /= total) error stop "multiset_sequences: enumeration count mismatch"
    end function multiset_sequences

    function multinomial_permutations(group_sizes) result(out)
        integer, intent(in) :: group_sizes(:)
        integer, allocatable :: out(:,:), labels(:), seqs(:,:)
        integer :: n, g, i, j, col, pos, offset

        if (any(group_sizes < 0)) error stop "multinomial_permutations: negative group size"
        n = sum(group_sizes)
        if (n == 0) then
            allocate(out(0,1))
            return
        end if
        allocate(labels(n))
        pos = 0
        do g = 1, size(group_sizes)
            do i = 1, group_sizes(g)
                pos = pos + 1
                labels(pos) = g
            end do
        end do
        seqs = multiset_permutations(labels)
        allocate(out(n,size(seqs,2)))
        do col = 1, size(seqs,2)
            pos = 0
            offset = 0
            do g = 1, size(group_sizes)
                do j = 1, n
                    if (seqs(j,col) == g) then
                        pos = pos + 1
                        out(pos,col) = j
                    end if
                end do
                offset = offset + group_sizes(g)
            end do
        end do
    end function multinomial_permutations

    function all_binomial(n, k) result(out)
        integer, intent(in) :: n, k
        integer, allocatable :: out(:,:), full(:,:)

        if (n < 0 .or. k < 0 .or. k > n) error stop "all_binomial: invalid n or k"
        full = multinomial_permutations([k, n - k])
        allocate(out(k,size(full,2)))
        if (k > 0) out = full(1:k,:)
    end function all_binomial

    function generalized_riffles(group_sizes) result(out)
        integer, intent(in) :: group_sizes(:)
        integer, allocatable :: out(:,:), labels(:), seqs(:,:)
        integer :: n, g, i, col, pos, offset, rank

        if (any(group_sizes < 0)) error stop "generalized_riffles: negative group size"
        n = sum(group_sizes)
        if (n == 0) then
            allocate(out(0,1))
            return
        end if
        allocate(labels(n))
        pos = 0
        do g = 1, size(group_sizes)
            do i = 1, group_sizes(g)
                pos = pos + 1
                labels(pos) = g
            end do
        end do
        seqs = multiset_permutations(labels)
        allocate(out(n,size(seqs,2)))
        do col = 1, size(seqs,2)
            out(:,col) = seqs(:,col)
            offset = 0
            do g = 1, size(group_sizes)
                rank = 0
                do i = 1, n
                    if (seqs(i,col) == g) then
                        rank = rank + 1
                        out(i,col) = offset + rank
                    end if
                end do
                offset = offset + group_sizes(g)
            end do
        end do
    end function generalized_riffles

    function riffles(p, q) result(out)
        integer, intent(in) :: p
        integer, intent(in), optional :: q
        integer, allocatable :: out(:,:)
        integer :: qq

        qq = p
        if (present(q)) qq = q
        out = generalized_riffles([p,qq])
    end function riffles

    pure logical function tuple_less(a, b, n) result(value)
        integer, intent(in) :: a(:), b(:)
        integer, intent(in) :: n
        integer :: i

        value = .false.
        do i = 1, n
            if (a(i) /= b(i)) then
                value = a(i) < b(i)
                return
            end if
        end do
    end function tuple_less

    subroutine reset_block(block, n)
        integer, intent(inout) :: block(:)
        integer, intent(in) :: n
        integer :: i

        do i = 1, size(block)
            block(i) = huge(0) - n + i
        end do
    end subroutine reset_block

    pure subroutine sort_descending(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) >= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_descending

    pure subroutine sort_ascending(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_ascending

    subroutine count_to_default(n8, n, where)
        integer(i8), intent(in) :: n8
        integer, intent(out) :: n
        character(*), intent(in) :: where

        if (n8 < 0_i8 .or. n8 > int(huge(n), i8)) then
            error stop trim(where) // ": result has too many columns"
        end if
        n = int(n8)
    end subroutine count_to_default

end module partitions_sets
