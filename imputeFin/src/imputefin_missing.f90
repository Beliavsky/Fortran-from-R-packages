! SPDX-License-Identifier: GPL-3.0-only
module imputefin_missing
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use imputefin_kinds, only : dp
  implicit none
  private
  public :: missing_blocks, is_inner_na, any_inner_na, collect_indices, trim_observed_range

  type, public :: missing_blocks
    integer :: n = 0
    integer :: n_block = 0
    integer, allocatable :: obs(:)
    integer, allocatable :: miss(:)
    integer, allocatable :: first(:)
    integer, allocatable :: last(:)
    integer, allocatable :: length(:)
  contains
    procedure :: build => build_missing_blocks
  end type missing_blocks
contains
  elemental logical function is_missing(x)
    real(dp), intent(in) :: x
    is_missing = ieee_is_nan(x)
  end function is_missing

  function is_inner_na(y) result(mask)
    real(dp), intent(in) :: y(:)
    logical, allocatable :: mask(:)
    integer :: i, first_obs, last_obs
    allocate(mask(size(y))); mask=.false.
    first_obs=0; last_obs=0
    do i=1,size(y)
      if (.not. is_missing(y(i))) then; first_obs=i; exit; end if
    end do
    do i=size(y),1,-1
      if (.not. is_missing(y(i))) then; last_obs=i; exit; end if
    end do
    if(first_obs>0 .and. last_obs>first_obs) then
      do i=first_obs,last_obs; mask(i)=is_missing(y(i)); end do
    end if
  end function is_inner_na

  logical function any_inner_na(y)
    real(dp), intent(in) :: y(:)
    logical, allocatable :: m(:)
    m=is_inner_na(y); any_inner_na=any(m)
  end function any_inner_na

  subroutine collect_indices(mask, idx)
    logical, intent(in) :: mask(:)
    integer, allocatable, intent(out) :: idx(:)
    integer :: i,k
    allocate(idx(count(mask))); k=0
    do i=1,size(mask)
      if(mask(i)) then; k=k+1; idx(k)=i; end if
    end do
  end subroutine collect_indices

  subroutine trim_observed_range(y, first_obs, last_obs, status)
    real(dp), intent(in) :: y(:)
    integer, intent(out) :: first_obs,last_obs,status
    integer :: i
    first_obs=0; last_obs=0; status=0
    do i=1,size(y)
      if(.not.is_missing(y(i))) then; first_obs=i; exit; end if
    end do
    do i=size(y),1,-1
      if(.not.is_missing(y(i))) then; last_obs=i; exit; end if
    end do
    if(first_obs==0 .or. last_obs==0) status=1
  end subroutine trim_observed_range

  subroutine build_missing_blocks(self,y)
    class(missing_blocks), intent(inout) :: self
    real(dp), intent(in) :: y(:)
    logical, allocatable :: missmask(:), obsmask(:)
    integer :: i,b,start
    self%n=size(y)
    missmask=is_inner_na(y)
    obsmask=.not.[(is_missing(y(i)),i=1,size(y))]
    call collect_indices(obsmask,self%obs)
    call collect_indices(missmask,self%miss)
    self%n_block=0; i=1
    do while(i<=size(y))
      if(missmask(i)) then
        self%n_block=self%n_block+1
        do while(i<=size(y) .and. missmask(i)); i=i+1; end do
      else
        i=i+1
      end if
    end do
    allocate(self%first(self%n_block),self%last(self%n_block),self%length(self%n_block))
    b=0; i=1
    do while(i<=size(y))
      if(missmask(i)) then
        b=b+1; start=i
        do while(i<=size(y) .and. missmask(i)); i=i+1; end do
        self%first(b)=start; self%last(b)=i-1; self%length(b)=i-start
      else
        i=i+1
      end if
    end do
  end subroutine build_missing_blocks
end module imputefin_missing
