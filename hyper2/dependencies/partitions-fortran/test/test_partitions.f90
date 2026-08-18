! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_partitions
    use partitions
    implicit none
    integer, allocatable :: a(:,:), b(:,:), c(:,:), y(:), ar(:)
    integer :: n, j, k

    do n = 1, 20
        a = parts(n)
        call check(size(a,2) == int(partition_count(n)), "parts count")
        do j = 1, size(a,2)
            call check(sum(a(:,j)) == n, "parts sums")
            call check(all(a(:,j) >= 0), "parts nonnegative")
            do k = 1, n - 1
                call check(a(k,j) >= a(k+1,j), "parts decreasing")
            end do
        end do
        call check(all(a(:,size(a,2)) == 1), "parts final column")

        b = distinct_parts(n)
        call check(size(b,2) == int(distinct_partition_count(n)), "distinct count")
        do j = 1, size(b,2)
            call check(sum(b(:,j)) == n, "distinct sums")
            do k = 1, size(b,1) - 1
                if (b(k+1,j) > 0) call check(b(k,j) > b(k+1,j), "distinct positive entries")
            end do
        end do
    end do

    c = restricted_parts(12, 5, .false.)
    call check(size(c,2) == 13, "restricted count")
    do j = 1, size(c,2)
        call check(sum(c(:,j)) == 12, "restricted sum")
        call check(all(c(:,j) > 0), "restricted positive")
    end do

    c = restricted_parts(7, 4, .true.)
    do j = 1, size(c,2)
        call check(sum(c(:,j)) == 7, "restricted zero sum")
        call check(all(c(:,j) >= 0), "restricted zero nonnegative")
    end do

    y = conjugate([7,7,5,4,4,2,1])
    call check(all(y == [7,6,5,5,3,2,2]), "conjugate Andrews example")
    ar = [4,7,1,4,2,7,5]
    call check(all(conjugate(ar, .false.) == y), "conjugate unsorted")
    call check(durfee([7,7,5,4,4,2,1]) == 4, "durfee")
    call check(durfee(ar, .false.) == 4, "durfee unsorted")

    print *, "test_partitions: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_partitions
