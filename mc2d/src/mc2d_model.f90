! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_model
  use mc2d_kinds, only : dp
  use mc2d_node, only : mc
  use mc2d_random, only : seed_random
  implicit none
  private
  public :: mcmodel, mcmodelcut, mccut_result, evalmcmod, evalmccut, evalmccut_reduce

  abstract interface
    subroutine model_callback(nsv,nsu,result)
      import mc
      integer,intent(in)::nsv,nsu
      type(mc),intent(out)::result
    end subroutine model_callback

    subroutine cut_setup_callback(nsv,nsu)
      integer,intent(in)::nsv,nsu
    end subroutine cut_setup_callback

    subroutine cut_column_callback(nsv,index,result)
      import mc
      integer,intent(in)::nsv,index
      type(mc),intent(out)::result
    end subroutine cut_column_callback

    subroutine cut_reduce_callback(column,statistic)
      import mc,dp
      type(mc),intent(in)::column
      real(dp),allocatable,intent(out)::statistic(:)
    end subroutine cut_reduce_callback
  end interface

  type :: mcmodel
    procedure(model_callback),pointer,nopass :: evaluate=>null()
  end type mcmodel

  type :: mcmodelcut
    procedure(model_callback),pointer,nopass :: evaluate=>null()
    procedure(cut_setup_callback),pointer,nopass :: setup=>null()
    procedure(cut_column_callback),pointer,nopass :: evaluate_column=>null()
  end type mcmodelcut

  type :: mccut_result
    real(dp),allocatable :: value(:,:)
  contains
    procedure :: nstatistics => mccut_nstatistics
    procedure :: nsu => mccut_nsu
  end type mccut_result
contains
  function evalmcmod(model,nsv,nsu,seed) result(res)
    type(mcmodel),intent(in)::model
    integer,intent(in)::nsv,nsu
    integer,intent(in),optional::seed
    type(mc)::res
    if(.not.associated(model%evaluate))error stop 'evalmcmod: model callback is not associated'
    if(present(seed))call seed_random(seed)
    call model%evaluate(nsv,nsu,res)
  end function evalmcmod

  function evalmccut(model,nsv,nsu,seed) result(res)
    type(mcmodelcut),intent(in)::model
    integer,intent(in)::nsv,nsu
    integer,intent(in),optional::seed
    type(mc)::res
    if(.not.associated(model%evaluate))then
      error stop 'evalmccut: full-model callback is not associated; use evalmccut_reduce for column mode'
    end if
    if(present(seed))call seed_random(seed)
    if(associated(model%setup))call model%setup(nsv,nsu)
    call model%evaluate(nsv,nsu,res)
  end function evalmccut

  function evalmccut_reduce(model,nsv,nsu,reducer,seed) result(res)
    type(mcmodelcut),intent(in)::model
    integer,intent(in)::nsv,nsu
    procedure(cut_reduce_callback)::reducer
    integer,intent(in),optional::seed
    type(mccut_result)::res
    type(mc)::column
    real(dp),allocatable::stat(:)
    integer::j,nstat

    if(nsv<=0.or.nsu<=0)error stop 'evalmccut_reduce: nsv and nsu must be positive'
    if(.not.associated(model%evaluate_column)) &
      error stop 'evalmccut_reduce: column callback is not associated'
    if(present(seed))call seed_random(seed)
    if(associated(model%setup))call model%setup(nsv,nsu)

    do j=1,nsu
      call model%evaluate_column(nsv,j,column)
      call reducer(column,stat)
      if(.not.allocated(stat))error stop 'evalmccut_reduce: reducer did not allocate its statistic vector'
      if(j==1)then
        nstat=size(stat)
        if(nstat<=0)error stop 'evalmccut_reduce: reducer returned no statistics'
        allocate(res%value(nstat,nsu))
      else if(size(stat)/=nstat)then
        error stop 'evalmccut_reduce: reducer output size changed between columns'
      end if
      res%value(:,j)=stat
      deallocate(stat)
    end do
  end function evalmccut_reduce

  integer function mccut_nstatistics(self) result(n)
    class(mccut_result),intent(in)::self
    if(allocated(self%value))then
      n=size(self%value,1)
    else
      n=0
    end if
  end function mccut_nstatistics

  integer function mccut_nsu(self) result(n)
    class(mccut_result),intent(in)::self
    if(allocated(self%value))then
      n=size(self%value,2)
    else
      n=0
    end if
  end function mccut_nsu
end module mc2d_model
