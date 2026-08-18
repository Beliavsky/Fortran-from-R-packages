! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program test_compat
    use partitions
    implicit none
    integer, allocatable :: a(:,:), b(:), c(:,:)

    call check(p(100) == 190569292_i8, "p alias")
    call check(q(100) == 444793_i8, "q alias")
    call check(r(5,12) == 13_i8, "r alias")
    call check(s([1,1,2,2,3,3,4,4],5) == 474_i8, "s alias")
    a = diffparts(8)
    call check(all(shape(a) == [3,6]), "diffparts alias")
    b = firstpart(6)
    b = nextpart(b)
    call check(all(b == [5,1,0,0,0,0]), "nextpart alias")
    c = setparts([2,1,1])
    call check(all(shape(c) == [4,6]), "setparts alias")
    a = perms(4)
    call check(all(shape(a) == [4,24]), "perms alias")
    b = bintocomp(comptobin([3,2,1]))
    call check(all(b == [3,2,1]), "composition aliases")
    c = riffle(2,2)
    call check(all(shape(c) == [4,6]), "riffle alias")

    print *, "test_compat: PASS"

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", trim(label)
            error stop 1
        end if
    end subroutine check
end program test_compat
