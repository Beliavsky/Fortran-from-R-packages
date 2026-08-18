! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_permutations
    use partitions_kinds, only : i8
    use partitions_counts, only : factorial_i8, multiset_permutation_count
    implicit none
    private

    public :: next_permutation
    public :: permutations
    public :: plain_permutations
    public :: multiset_permutations

contains

    subroutine next_permutation(a, done)
        integer, intent(inout) :: a(:)
        logical, intent(out), optional :: done
        integer :: j, l, k, tmp
        logical :: finished

        finished = .false.
        if (size(a) <= 1) then
            finished = .true.
            if (present(done)) done = finished
            return
        end if

        j = size(a) - 1
        do while (j >= 1)
            if (a(j) < a(j + 1)) exit
            j = j - 1
        end do
        if (j < 1) then
            finished = .true.
            if (present(done)) done = finished
            return
        end if

        l = size(a)
        do while (a(j) >= a(l))
            l = l - 1
        end do
        tmp = a(l)
        a(l) = a(j)
        a(j) = tmp

        k = j + 1
        l = size(a)
        do while (k < l)
            tmp = a(l)
            a(l) = a(k)
            a(k) = tmp
            k = k + 1
            l = l - 1
        end do
        if (present(done)) done = finished
    end subroutine next_permutation

    function permutations(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        integer(i8) :: nc8
        integer :: nc, j, i
        logical :: done

        if (n < 0) error stop "permutations: n must be nonnegative"
        nc8 = factorial_i8(n)
        call count_to_default(nc8, nc, "permutations")
        allocate(out(n, nc))
        if (n == 0) return
        out(:,1) = [(i, i = 1, n)]
        do j = 2, nc
            out(:,j) = out(:,j - 1)
            call next_permutation(out(:,j), done)
            if (done) error stop "permutations: premature termination"
        end do
    end function permutations

    function multiset_permutations(v) result(out)
        integer, intent(in) :: v(:)
        integer, allocatable :: out(:,:)
        integer(i8) :: nc8
        integer :: nc, j
        logical :: done

        nc8 = multiset_permutation_count(v)
        call count_to_default(nc8, nc, "multiset_permutations")
        allocate(out(size(v), nc))
        if (size(v) == 0) return
        out(:,1) = v
        call sort_ascending(out(:,1))
        do j = 2, nc
            out(:,j) = out(:,j - 1)
            call next_permutation(out(:,j), done)
            if (done) error stop "multiset_permutations: premature termination"
        end do
    end function multiset_permutations

    function plain_permutations(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        integer, allocatable :: c(:), o(:)
        integer(i8) :: nc8
        integer :: nc, i, jj, j, q, s, i1, i2, tmp

        if (n < 0) error stop "plain_permutations: n must be nonnegative"
        nc8 = factorial_i8(n)
        call count_to_default(nc8, nc, "plain_permutations")
        allocate(out(n, nc))
        if (n == 0) return
        allocate(c(n), o(n))
        c = 0
        o = 1
        out(:,1) = [(jj, jj = 1, n)]

        do i = 2, nc
            out(:,i) = out(:,i - 1)
            j = n
            s = 0
            do
                q = c(j) + o(j)
                if (q >= 0) then
                    if (q /= j) then
                        i1 = j - c(j) + s
                        i2 = j - q + s
                        tmp = out(i1, i)
                        out(i1, i) = out(i2, i)
                        out(i2, i) = tmp
                        c(j) = q
                        exit
                    end if
                    s = s + 1
                end if
                o(j) = -o(j)
                j = j - 1
                if (j < 1) error stop "plain_permutations: algorithm failure"
            end do
        end do
    end function plain_permutations

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

end module partitions_permutations
