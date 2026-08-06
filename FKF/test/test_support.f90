! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
  use fkf_kinds, only : dp
  implicit none
  private
  public :: assert_close, assert_true, finish_test

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
      write(*, '(a,1x,es24.16,1x,a,1x,es24.16)') trim(label)//': got', actual, 'expected', expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a)') 'assertion failed: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine finish_test(name)
    character(len=*), intent(in) :: name
    write(*, '(a)') 'PASS: '//trim(name)
  end subroutine finish_test

end module test_support
