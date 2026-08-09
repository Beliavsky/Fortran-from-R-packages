! SPDX-License-Identifier: BSD-2-Clause
module smoof_nk
  use smoof_kinds, only : dp
  implicit none
  private
  type, public :: nk_landscape
    integer :: n = 0
    integer, allocatable :: k(:)
    integer, allocatable :: link_start(:)
    integer, allocatable :: links(:)
    integer, allocatable :: value_start(:)
    real(dp), allocatable :: values(:)
  contains
    procedure :: evaluate => nk_evaluate
  end type nk_landscape
  public :: nk_evaluate_raw
contains
  pure real(dp) function nk_evaluate(self,x) result(f)
    class(nk_landscape),intent(in)::self
    integer,intent(in)::x(:)
    f=nk_evaluate_raw(x,self%k,self%link_start,self%links,self%value_start,self%values)
  end function nk_evaluate

  pure real(dp) function nk_evaluate_raw(x,k,link_start,links,value_start,values) result(f)
    integer,intent(in)::x(:),k(:),link_start(:),links(:),value_start(:)
    real(dp),intent(in)::values(:)
    integer::i,j,kk,off,idx
    f=0.0_dp
    do i=1,size(x)
      kk=k(i); off=x(i)*2**kk
      do j=1,kk
        off=off+x(links(link_start(i)+j-1))*2**(kk-j)
      end do
      idx=value_start(i)+off
      f=f+values(idx)
    end do
    f=f/real(size(x),dp)
  end function nk_evaluate_raw
end module smoof_nk
