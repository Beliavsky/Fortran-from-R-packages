! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_counts
    use partitions
    implicit none
    integer(i8), allocatable :: pn(:), qn(:)
    integer :: i

    call check(partition_count(100) == 190569292_i8, "P(100)")
    call check(distinct_partition_count(100) == 444793_i8, "Q(100)")
    call check(restricted_partition_count(5, 12, .false.) == 13_i8, "R(5,12)")
    call check(block_partition_count([1,1,2,2,3,3,4,4], 5) == 474_i8, "S example")
    call check(set_partition_count([20,1]) == 21_i8, "set count avoids factorial overflow")
    call check(multiset_permutation_count([integer :: (1, i=1,20), 2]) == 21_i8, &
               "multiset count avoids factorial overflow")

    allocate(pn(0:10), qn(0:10))
    pn = partition_numbers(10)
    qn = distinct_partition_numbers(10)
    call check(all(pn == [1_i8,1_i8,2_i8,3_i8,5_i8,7_i8,11_i8,15_i8,22_i8,30_i8,42_i8]), &
               "partition number sequence")
    call check(all(qn == [1_i8,1_i8,1_i8,2_i8,2_i8,3_i8,4_i8,5_i8,6_i8,8_i8,10_i8]), &
               "distinct partition sequence")

    print *, "test_counts: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_counts
