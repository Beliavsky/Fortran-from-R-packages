! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_successors
    use partitions
    implicit none
    integer, allocatable :: allp(:,:), cur(:), allb(:,:), allc(:,:), alld(:,:), allr(:,:)
    integer :: j

    allp = parts(8)
    cur = first_part(8)
    call check(all(cur == allp(:,1)), "first part")
    do j = 2, size(allp,2)
        call next_part(cur)
        call check(all(cur == allp(:,j)), "next part")
    end do
    call check(is_last_part(cur), "last part")

    alld = distinct_parts(12)
    cur = first_distinct_part(12)
    do j = 2, size(alld,2)
        call next_distinct_part(cur)
        call check(all(cur == alld(:,j)), "next distinct")
    end do
    call check(is_last_distinct_part(cur), "last distinct")

    allr = restricted_parts(11, 3, .false.)
    cur = first_restricted_part(11, 3, .false.)
    do j = 2, size(allr,2)
        call next_restricted_part(cur)
        call check(all(cur == allr(:,j)), "next restricted")
    end do
    call check(is_last_restricted_part(cur), "last restricted")

    allb = block_parts([1,2,3,4], 4, .false.)
    cur = first_block_part([1,2,3,4], 4, .false.)
    do j = 2, size(allb,2)
        call next_block_part(cur, [1,2,3,4], 4, .false.)
        call check(all(cur == allb(:,j)), "next block exact")
    end do
    call check(is_last_block_part(cur, [1,2,3,4], 4, .false.), "last block exact")

    allb = block_parts([1,2,3,4], 4, .true.)
    cur = first_block_part([1,2,3,4], 4, .true.)
    do j = 2, size(allb,2)
        call next_block_part(cur, [1,2,3,4], 4, .true.)
        call check(all(cur == allb(:,j)), "next block fewer")
    end do
    call check(is_last_block_part(cur, [1,2,3,4], 4, .true.), "last block fewer")

    allc = compositions(7, 4, .true.)
    cur = first_composition(7, 4, .true.)
    do j = 2, size(allc,2)
        cur = next_composition(cur, .true., .true.)
        call check(all(cur == allc(:,j)), "next restricted composition zero")
    end do
    call check(is_last_composition(cur, .true., .true.), "last restricted composition zero")

    allc = compositions(7, 4, .false.)
    cur = first_composition(7, 4, .false.)
    do j = 2, size(allc,2)
        cur = next_composition(cur, .true., .false.)
        call check(all(cur == allc(:,j)), "next restricted composition positive")
    end do

    print *, "test_successors: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_successors
