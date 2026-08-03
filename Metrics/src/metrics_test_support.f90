module test_support
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use metrics_kinds, only : dp
    implicit none
    private
    public :: check_close, check_array_close, check_true, finish_tests

    integer :: failures = 0

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected
        real(dp), intent(in), optional :: tolerance
        real(dp) :: tol

        tol = 1.0e-12_dp
        if (present(tolerance)) tol = tolerance
        if (ieee_is_nan(expected)) then
            if (.not. ieee_is_nan(actual)) call fail(name, actual, expected)
        else if (.not. ieee_is_finite(expected)) then
            if (actual /= expected) call fail(name, actual, expected)
        else if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tol) then
            call fail(name, actual, expected)
        end if
    end subroutine check_close

    subroutine check_array_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:)
        real(dp), intent(in), optional :: tolerance
        integer :: i

        if (size(actual) /= size(expected)) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//trim(name)//' size mismatch'
            return
        end if
        do i = 1, size(actual)
            call check_close(trim(name)//'['//integer_string(i)//']', actual(i), expected(i), tolerance)
        end do
    end subroutine check_array_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: '//trim(name)
        end if
    end subroutine check_true

    subroutine finish_tests(name)
        character(len=*), intent(in) :: name
        if (failures > 0) then
            write(*, '(a, i0)') trim(name)//': failures = ', failures
            error stop 1
        end if
        write(*, '(a)') trim(name)//': PASS'
    end subroutine finish_tests

    subroutine fail(name, actual, expected)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected
        failures = failures + 1
        write(*, '(a, 2(1x, es24.16))') 'FAIL: '//trim(name)//' actual expected:', actual, expected
    end subroutine fail

    function integer_string(i) result(text)
        integer, intent(in) :: i
        character(len=32) :: text
        write(text, '(i0)') i
    end function integer_string

end module test_support
