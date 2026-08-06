! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module test_support
    use matlab_kinds, only : dp
    implicit none
    private

    public :: assert_true
    public :: assert_close
    public :: assert_all_close
    public :: assert_int_equal

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(*, '(a)') 'FAIL: ' // trim(message)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        call assert_true(abs(actual - expected) <= tolerance, message)
    end subroutine assert_close

    subroutine assert_all_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(..), expected(..), tolerance
        character(len=*), intent(in) :: message
        logical :: ok

        ok = .false.
        select rank (actual)
        rank (1)
            select rank (expected)
            rank (1)
                ok = size(actual) == size(expected)
                if (ok) ok = all(abs(actual - expected) <= tolerance)
            end select
        rank (2)
            select rank (expected)
            rank (2)
                ok = all(shape(actual) == shape(expected))
                if (ok) ok = all(abs(actual - expected) <= tolerance)
            end select
        rank (3)
            select rank (expected)
            rank (3)
                ok = all(shape(actual) == shape(expected))
                if (ok) ok = all(abs(actual - expected) <= tolerance)
            end select
        end select
        call assert_true(ok, message)
    end subroutine assert_all_close

    subroutine assert_int_equal(actual, expected, message)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: message

        call assert_true(actual == expected, message)
    end subroutine assert_int_equal
end module test_support
