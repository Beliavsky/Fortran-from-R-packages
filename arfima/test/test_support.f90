module test_support
  use arfima_kinds, only : dp
  implicit none
  private
  public :: assert_true, assert_close, assert_vector_close
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::message
    if(abs(actual-expected)>tolerance) then
      write(*,'(a,2es24.14,a,es12.4)') 'FAIL: '//trim(message)//' actual/expected=',actual,expected,' tol=',tolerance
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_vector_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual(:),expected(:),tolerance
    character(len=*),intent(in)::message
    if(size(actual)/=size(expected) .or. any(abs(actual-expected)>tolerance)) then
      write(*,'(a)') 'FAIL: '//trim(message)
      write(*,'(*(es14.6,1x))') actual
      write(*,'(*(es14.6,1x))') expected
      error stop 1
    end if
  end subroutine assert_vector_close
end module test_support
