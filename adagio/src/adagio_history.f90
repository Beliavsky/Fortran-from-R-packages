! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_history
  use adagio_kinds, only : dp
  implicit none
  private
  public :: history_buffer

  type :: history_buffer
     integer :: nvars=0
     integer :: ncalls=0
     logical :: store_inputs=.false.
     real(dp),allocatable :: input(:,:)
     real(dp),allocatable :: values(:)
   contains
     procedure :: reset => history_reset
     procedure :: record => history_record
  end type

contains

  subroutine history_reset(self)
    class(history_buffer),intent(inout)::self
    self%ncalls=0
    if(allocated(self%input))deallocate(self%input)
    if(allocated(self%values))deallocate(self%values)
  end subroutine

  subroutine history_record(self,x,value)
    class(history_buffer),intent(inout)::self
    real(dp),intent(in)::x(:),value
    real(dp),allocatable::ni(:,:),nv(:)
    integer::old
    old=self%ncalls;self%ncalls=old+1
    if(.not.allocated(self%values))then
       allocate(self%values(1));self%values(1)=value
    else
       allocate(nv(old+1));nv(1:old)=self%values;nv(old+1)=value;call move_alloc(nv,self%values)
    end if
    if(self%store_inputs)then
       if(self%nvars==0)self%nvars=size(x)
       if(.not.allocated(self%input))then
          allocate(self%input(1,self%nvars));self%input(1,:)=x
       else
          allocate(ni(old+1,self%nvars));ni(1:old,:)=self%input;ni(old+1,:)=x;call move_alloc(ni,self%input)
       end if
    end if
  end subroutine
end module adagio_history
