! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_returns
  use qrmtools_kinds, only : dp
  implicit none
  private
  integer, parameter, public :: return_logarithmic=1, return_simple=2, return_difference=3
  public :: compute_returns, invert_returns
contains
  function compute_returns(prices,method) result(values)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in), optional :: method
    real(dp), allocatable :: values(:,:)
    integer :: m
    m=return_logarithmic; if(present(method))m=method
    allocate(values(max(size(prices,1)-1,0),size(prices,2)))
    select case(m)
    case(return_logarithmic); values=log(prices(2:,:)/prices(:size(prices,1)-1,:))
    case(return_simple); values=prices(2:,:)/prices(:size(prices,1)-1,:)-1.0_dp
    case(return_difference); values=prices(2:,:)-prices(:size(prices,1)-1,:)
    end select
  end function compute_returns

  function invert_returns(values,start,method) result(prices)
    real(dp), intent(in) :: values(:,:),start(:)
    integer, intent(in), optional :: method
    real(dp), allocatable :: prices(:,:)
    integer :: m,i
    m=return_logarithmic; if(present(method))m=method
    allocate(prices(size(values,1)+1,size(values,2))); prices(1,:)=start
    do i=1,size(values,1)
      select case(m)
      case(return_logarithmic); prices(i+1,:)=prices(i,:)*exp(values(i,:))
      case(return_simple); prices(i+1,:)=prices(i,:)*(1.0_dp+values(i,:))
      case(return_difference); prices(i+1,:)=prices(i,:)+values(i,:)
      end select
    end do
  end function invert_returns
end module qrmtools_returns
