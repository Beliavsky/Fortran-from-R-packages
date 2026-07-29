! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
    use evir, only : dp
    implicit none
    private
    public :: check, check_close, finish_tests
    integer :: failures = 0
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) then
            failures=failures+1
            print '(a)', 'FAIL: '//trim(message)
        end if
    end subroutine check
    subroutine check_close(actual,expected,tolerance,message)
        real(dp),intent(in)::actual,expected,tolerance
        character(len=*),intent(in)::message
        call check(abs(actual-expected)<=tolerance*(1.0_dp+abs(expected)),message)
    end subroutine check_close
    subroutine finish_tests()
        if(failures>0) then
            print '(a,i0)', 'test failures: ',failures
            error stop 1
        end if
        print '(a)', 'all tests passed'
    end subroutine finish_tests
end module test_support
