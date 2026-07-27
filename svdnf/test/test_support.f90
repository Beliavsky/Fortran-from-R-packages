! SPDX-License-Identifier: GPL-3.0-only
module test_support
  use svdnf_kinds, only : dp
  implicit none
  private
  public :: assert_true, assert_close, assert_vector_close
contains
  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'assertion failed: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,message)
    real(dp), intent(in) :: actual,expected,tolerance
    character(len=*), intent(in) :: message
    if (abs(actual-expected)>tolerance*(1.0_dp+abs(expected))) then
      write(*,'(a,3es24.15)') 'mismatch: ',actual,expected,abs(actual-expected)
      write(*,'(a)') trim(message)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual,expected,tolerance,message)
    real(dp), intent(in) :: actual(:),expected(:),tolerance
    character(len=*), intent(in) :: message
    integer :: i
    call assert_true(size(actual)==size(expected),message//' size')
    do i=1,size(actual)
      call assert_close(actual(i),expected(i),tolerance,message)
    end do
  end subroutine assert_vector_close
end module test_support
