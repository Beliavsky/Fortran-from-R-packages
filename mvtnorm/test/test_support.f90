! SPDX-License-Identifier: GPL-2.0-only
module test_support
  use mvtnorm_kinds, only : dp
  implicit none
contains
  subroutine assert_close(actual,expected,tol,label)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::label
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      write(*,'(a,3es25.16)') trim(label)//' mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(.not.condition) then
      write(*,'(a)') trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_matrix_close(a,b,tol,label)
    real(dp),intent(in)::a(:,:),b(:,:),tol
    character(len=*),intent(in)::label
    if(any(shape(a)/=shape(b))) then
      write(*,'(a)') trim(label)//' shape mismatch'; error stop 1
    end if
    if(maxval(abs(a-b))>tol*max(1.0_dp,maxval(abs(b)))) then
      write(*,'(a,es25.16)') trim(label)//' mismatch: ',maxval(abs(a-b)); error stop 1
    end if
  end subroutine assert_matrix_close
end module test_support
