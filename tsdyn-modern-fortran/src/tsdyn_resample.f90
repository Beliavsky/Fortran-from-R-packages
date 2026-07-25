! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_resample
  use tsdyn_kinds, only: dp
  use tsdyn_utils, only: random_normal
  implicit none
  private
  public :: resample_vector, block_resample_matrix
contains
  subroutine resample_vector(x, method, y, info, block_length)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: method
    real(dp), allocatable, intent(out) :: y(:)
    integer, intent(out) :: info
    integer, intent(in), optional :: block_length
    integer :: n, i, j, b, start, pos
    real(dp) :: u, w
    character(len=16) :: m
    n=size(x); allocate(y(n)); y=0.0_dp
    if(n<1)then;info=-1;return;end if
    m=adjustl(method)
    select case(trim(m))
    case('resample','iid')
      do i=1,n
        call random_number(u);j=min(n,1+int(u*real(n,dp)));y(i)=x(j)
      end do
    case('block')
      b=max(1,nint(sqrt(real(n,dp))));if(present(block_length))b=max(1,min(n,block_length))
      pos=1
      do while(pos<=n)
        call random_number(u);start=min(n,1+int(u*real(n,dp)))
        do j=0,b-1
          if(pos>n)exit
          y(pos)=x(1+mod(start-1+j,n));pos=pos+1
        end do
      end do
    case('wild1','rademacher')
      do i=1,n
        call random_number(u);w=merge(1.0_dp,-1.0_dp,u>=0.5_dp);y(i)=x(i)*w
      end do
    case('wild2','mammen')
      do i=1,n
        call random_number(u)
        if(u<(sqrt(5.0_dp)+1.0_dp)/(2.0_dp*sqrt(5.0_dp)))then
          w=(1.0_dp-sqrt(5.0_dp))/2.0_dp
        else
          w=(1.0_dp+sqrt(5.0_dp))/2.0_dp
        end if
        y(i)=x(i)*w
      end do
    case('normal')
      do i=1,n;y(i)=x(i)*random_normal();end do
    case('check','identity')
      y=x
    case default
      info=-2;return
    end select
    info=0
  end subroutine resample_vector

  subroutine block_resample_matrix(x, block_length, y, info)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: block_length
    real(dp), allocatable, intent(out) :: y(:,:)
    integer, intent(out) :: info
    integer :: n, k, pos, start, j
    real(dp) :: u
    n=size(x,1);k=size(x,2);allocate(y(n,k));y=0.0_dp
    if(n<1.or.k<1.or.block_length<1)then;info=-1;return;end if
    pos=1
    do while(pos<=n)
      call random_number(u);start=min(n,1+int(u*real(n,dp)))
      do j=0,block_length-1
        if(pos>n)exit
        y(pos,:)=x(1+mod(start-1+j,n),:);pos=pos+1
      end do
    end do
    info=0
  end subroutine block_resample_matrix
end module tsdyn_resample
