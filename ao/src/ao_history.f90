! SPDX-License-Identifier: GPL-3.0-only
module ao_history_mod
  use ao_kinds, only : dp
  use ao_types, only : ao_history
  implicit none
  private
  public :: history_initialize, history_append, history_best_index
  public :: history_first_index, history_total_seconds
contains
  subroutine history_initialize(h, parameter, value)
    type(ao_history), intent(out) :: h
    real(dp), intent(in) :: parameter(:), value
    integer, parameter :: initial_capacity = 32
    integer :: npar
    npar=size(parameter); h%capacity=initial_capacity; h%n=1
    allocate(h%iteration(h%capacity),h%value(h%capacity),h%seconds(h%capacity))
    allocate(h%parameter(npar,h%capacity),h%active(npar,h%capacity))
    h%iteration=0; h%value=0.0_dp; h%seconds=0.0_dp; h%parameter=0.0_dp; h%active=.false.
    h%iteration(1)=0; h%value(1)=value; h%parameter(:,1)=parameter
  end subroutine history_initialize

  subroutine grow_history(h)
    type(ao_history), intent(inout) :: h
    integer :: newcap, npar
    integer, allocatable :: iteration(:)
    real(dp), allocatable :: value(:), parameter(:,:), seconds(:)
    logical, allocatable :: active(:,:)
    newcap=max(2*h%capacity,h%capacity+16); npar=size(h%parameter,1)
    allocate(iteration(newcap),value(newcap),seconds(newcap))
    allocate(parameter(npar,newcap),active(npar,newcap))
    iteration=0; value=0.0_dp; seconds=0.0_dp; parameter=0.0_dp; active=.false.
    iteration(1:h%n)=h%iteration(1:h%n); value(1:h%n)=h%value(1:h%n)
    seconds(1:h%n)=h%seconds(1:h%n); parameter(:,1:h%n)=h%parameter(:,1:h%n)
    active(:,1:h%n)=h%active(:,1:h%n)
    call move_alloc(iteration,h%iteration); call move_alloc(value,h%value)
    call move_alloc(seconds,h%seconds); call move_alloc(parameter,h%parameter)
    call move_alloc(active,h%active); h%capacity=newcap
  end subroutine grow_history

  subroutine history_append(h, iteration, value, parameter, block, seconds)
    type(ao_history), intent(inout) :: h
    integer, intent(in) :: iteration, block(:)
    real(dp), intent(in) :: value, parameter(:), seconds
    integer :: k
    if(h%n>=h%capacity) call grow_history(h)
    h%n=h%n+1; k=h%n
    h%iteration(k)=iteration; h%value(k)=value; h%parameter(:,k)=parameter
    h%active(:,k)=.false.; h%active(block,k)=.true.; h%seconds(k)=seconds
  end subroutine history_append

  integer function history_best_index(h, minimize) result(idx)
    type(ao_history), intent(in) :: h
    logical, intent(in) :: minimize
    if(minimize) then
      idx=minloc(h%value(1:h%n),dim=1)
    else
      idx=maxloc(h%value(1:h%n),dim=1)
    end if
  end function history_best_index

  integer function history_first_index(h, iteration) result(idx)
    type(ao_history), intent(in) :: h
    integer, intent(in) :: iteration
    integer :: k
    idx=0
    do k=1,h%n
      if(h%iteration(k)==iteration) then; idx=k; return; end if
    end do
  end function history_first_index

  real(dp) function history_total_seconds(h) result(value)
    type(ao_history), intent(in) :: h
    value=sum(h%seconds(1:h%n))
  end function history_total_seconds
end module ao_history_mod
