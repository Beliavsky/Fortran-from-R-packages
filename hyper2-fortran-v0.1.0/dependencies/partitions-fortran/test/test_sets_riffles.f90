! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_sets_riffles
    use partitions
    implicit none
    integer, allocatable :: a(:,:), rr(:,:), mm(:,:), multi(:,:)
    integer, parameter :: expected(4,6) = reshape([ &
        1,2,3,1, 1,1,2,3, 1,2,1,3, 2,1,3,1, 2,1,1,3, 2,3,1,1], [4,6])
    integer :: i, j

    a = set_partitions([2,1,1])
    call check(all(shape(a) == [4,6]), "setparts shape")
    call check(all(a == expected), "setparts exact upstream order")

    a = set_partitions([3,2])
    call check(size(a,2) == 10, "setparts 3,2 count")
    do j = 1, size(a,2)
        call check(count(a(:,j) == 1) == 3, "setparts block 1 size")
        call check(count(a(:,j) == 2) == 2, "setparts block 2 size")
    end do

    rr = restricted_set_partitions([2,1,1])
    call check(all(shape(rr) == [4,6]), "restricted setparts shape")
    do j = 1, size(rr,2)
        call check(all([(count(rr(:,j) == i) == 1, i=1,4)]), "restricted setparts permutation")
    end do

    multi = multiset_sequences([1,1,2,3], 3)
    call check(size(multi,1) == 3, "multiset selection rows")
    call check(size(multi,2) == 12, "multiset selection columns")

    mm = multinomial_permutations([2,2])
    call check(all(shape(mm) == [4,6]), "multinomial shape")
    a = all_binomial(4,2)
    call check(all(shape(a) == [2,6]), "allbinom shape")

    rr = generalized_riffles([2,2])
    call check(all(shape(rr) == [4,6]), "riffle shape")
    call check(all(rr(:,1) == [1,2,3,4]), "riffle first")
    rr = riffles(2,2)
    call check(all(shape(rr) == [4,6]), "riffles wrapper")

    print *, "test_sets_riffles: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_sets_riffles
