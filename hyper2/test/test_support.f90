! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
    use hyper2_kinds, only : dp
    implicit none
    private
    public :: check, check_close, check_vec
contains
    subroutine check(ok, msg, failures)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: failures
        if (.not. ok) then
            failures = failures + 1
            write(*,'(a)') 'FAIL: '//trim(msg)
        end if
    end subroutine check

    subroutine check_close(x, y, tol, msg, failures)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: failures
        call check(abs(x-y) <= tol*max(1.0_dp,abs(y)), msg, failures)
        if (abs(x-y) > tol*max(1.0_dp,abs(y))) then
            write(*,'(a,es24.16,a,es24.16)') '  got=', x, ' expected=', y
        end if
    end subroutine check_close

    subroutine check_vec(x, y, tol, msg, failures)
        real(dp), intent(in) :: x(:), y(:), tol
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: failures
        call check(size(x)==size(y), msg//' size', failures)
        if (size(x)==size(y)) then
            call check(maxval(abs(x-y)) <= tol*max(1.0_dp,maxval(abs(y))), msg, failures)
        end if
    end subroutine check_vec
end module test_support
