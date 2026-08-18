! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_compositions
    use partitions
    implicit none
    integer, allocatable :: a(:,:), bits(:), comp(:), next(:), padded(:)
    integer(i8) :: code
    integer :: n, j

    do n = 1, 10
        a = compositions(n)
        call check(size(a,2) == 2**(n-1), "unrestricted composition count")
        do j = 1, size(a,2)
            call check(sum(a(:,j)) == n, "unrestricted composition sum")
        end do
    end do

    a = compositions(15, 4, .true.)
    call check(all([(sum(a(:,j)) == 15, j=1,size(a,2))]), "zero composition sums")
    a = compositions(15, 4, .false.)
    call check(all([(sum(a(:,j)) == 15, j=1,size(a,2))]), "positive composition sums")
    call check(all(a > 0), "positive composition entries")

    do code = 0_i8, 63_i8
        bits = to_binary(code, 6)
        call check(to_decimal(bits) == code, "binary round trip")
    end do

    comp = [3,2,1]
    bits = composition_to_binary(comp)
    call check(all(bits == [0,0,1,0,1]), "composition to binary")
    call check(all(binary_to_composition(bits) == comp), "binary to composition")

    padded = first_composition(6)
    a = compositions(6)
    do j = 1, size(a,2)
        call check(all(padded == a(:,j)), "unrestricted next order")
        if (j < size(a,2)) then
            next = next_composition(padded, .false.)
            padded = 0
            padded(1:size(next)) = next
        end if
    end do

    print *, "test_compositions: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_compositions
