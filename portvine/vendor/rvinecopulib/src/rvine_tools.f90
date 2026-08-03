! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_tools
  use rvine_kinds, only : dp
  implicit none
  private
  public :: pseudo_observations, empirical_cdf, is_valid_permutation
contains
  subroutine pseudo_observations(x,u)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: u(:,:)
    integer :: i,j,k,n,d,nless,nequal
    n=size(x,2); d=size(x,1)
    do i=1,d
      do j=1,n
        nless=0; nequal=0
        do k=1,n
          if (x(i,k)<x(i,j)) nless=nless+1
          if (.not.(x(i,k)<x(i,j)) .and. .not.(x(i,k)>x(i,j))) nequal=nequal+1
        end do
        u(i,j)=(real(nless,dp)+0.5_dp*real(nequal+1,dp))/real(n+1,dp)
      end do
    end do
  end subroutine pseudo_observations

  subroutine empirical_cdf(x,p)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: p(:)
    integer :: i,j,n,count
    n=size(x)
    do i=1,n
      count=0
      do j=1,n
        if (x(j)<=x(i)) count=count+1
      end do
      p(i)=real(count,dp)/real(n,dp)
    end do
  end subroutine empirical_cdf

  pure logical function is_valid_permutation(order) result(ok)
    integer, intent(in) :: order(:)
    logical :: seen(size(order))
    integer :: i
    seen=.false.; ok=.true.
    do i=1,size(order)
      if (order(i)<1 .or. order(i)>size(order)) then
        ok=.false.; return
      end if
      if (seen(order(i))) then
        ok=.false.; return
      end if
      seen(order(i))=.true.
    end do
  end function is_valid_permutation
end module rvine_tools
