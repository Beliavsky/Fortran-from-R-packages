! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_permutations
    use partitions
    implicit none
    integer, allocatable :: a(:,:), b(:,:), c(:,:)
    integer :: j

    a = permutations(5)
    call check(size(a,2) == 120, "permutation count")
    call check(all(a(:,1) == [1,2,3,4,5]), "permutation first")
    call check(all(a(:,120) == [5,4,3,2,1]), "permutation last")
    do j = 1, size(a,2)
        call check(sum(a(:,j)) == 15, "permutation sum")
    end do

    b = multiset_permutations([1,2,2,3])
    call check(size(b,2) == 12, "multiset permutation count")
    call check(all(b(:,1) == [1,2,2,3]), "multiset first")
    call check(all(b(:,12) == [3,2,2,1]), "multiset last")

    c = plain_permutations(5)
    call check(size(c,2) == 120, "plain count")
    do j = 2, size(c,2)
        call check(count(c(:,j) /= c(:,j-1)) == 2, "plain adjacent swap")
    end do

    print *, "test_permutations: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_permutations
