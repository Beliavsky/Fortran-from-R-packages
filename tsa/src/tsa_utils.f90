! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_utils
  use tsa_kinds, only : dp
  use tseries_linalg, only : least_squares, solve_linear
  implicit none
  private
  public :: mean_value, variance_n, sd_n, difference_series, seasonal_difference
  public :: lag_vector, convolution_filter, recursive_filter, sort_with_index
  public :: normal_quantile, polynomial_convolution, arma_to_ma, linear_fit
  public :: build_lag_matrix, sample_with_replacement, quantile_sorted

contains

  pure real(dp) function mean_value(x) result(m)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      m = 0.0_dp
    else
      m = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure real(dp) function variance_n(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      v = 0.0_dp
      return
    end if
    m = mean_value(x)
    v = sum((x-m)**2) / real(size(x), dp)
  end function variance_n

  pure real(dp) function sd_n(x) result(s)
    real(dp), intent(in) :: x(:)
    s = sqrt(max(0.0_dp, variance_n(x)))
  end function sd_n

  function difference_series(x, d) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: d
    real(dp), allocatable :: y(:), tmp(:), nxt(:)
    integer :: k
    allocate(tmp(size(x)))
    tmp = x
    do k = 1, max(0,d)
      if (size(tmp) <= 1) then
        allocate(y(0))
        return
      end if
      allocate(nxt(size(tmp)-1))
      nxt = tmp(2:) - tmp(:size(tmp)-1)
      call move_alloc(nxt,tmp)
    end do
    call move_alloc(tmp,y)
  end function difference_series

  function seasonal_difference(x, d, period) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: d, period
    real(dp), allocatable :: y(:), tmp(:), nxt(:)
    integer :: k, s
    s=max(1,period)
    allocate(tmp(size(x)))
    tmp=x
    do k=1,max(0,d)
      if(size(tmp)<=s) then
        allocate(y(0))
        return
      end if
      allocate(nxt(size(tmp)-s))
      nxt=tmp(s+1:)-tmp(:size(tmp)-s)
      call move_alloc(nxt,tmp)
    end do
    call move_alloc(tmp,y)
  end function seasonal_difference

  subroutine lag_vector(x, lag, y, missing_value)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    real(dp), intent(out) :: y(:)
    real(dp), intent(in), optional :: missing_value
    real(dp) :: mv
    integer :: n
    n=size(x)
    mv=0.0_dp
    if(present(missing_value)) mv=missing_value
    y=mv
    if(lag==0) then
      y=x
    else if(lag>0 .and. lag<n) then
      y(lag+1:)=x(:n-lag)
    end if
  end subroutine lag_vector

  subroutine convolution_filter(x, coeff, y)
    real(dp), intent(in) :: x(:), coeff(:)
    real(dp), intent(out) :: y(:)
    integer :: i,j
    y=0.0_dp
    do i=1,size(x)
      do j=1,min(size(coeff),i)
        y(i)=y(i)+coeff(j)*x(i-j+1)
      end do
    end do
  end subroutine convolution_filter

  subroutine recursive_filter(x, coeff, y, init)
    ! R-style recursive filter: y_t = x_t + sum_j coeff(j) y_{t-j}
    real(dp), intent(in) :: x(:), coeff(:)
    real(dp), intent(out) :: y(:)
    real(dp), intent(in), optional :: init(:)
    integer :: i,j,k
    real(dp) :: prev
    y=0.0_dp
    do i=1,size(x)
      y(i)=x(i)
      do j=1,size(coeff)
        k=i-j
        if(k>=1) then
          prev=y(k)
        else if (present(init)) then
          if (size(init) >= j-i+1) then
            prev = init(size(init)-(j-i))
          else
            prev = 0.0_dp
          end if
        else
          prev = 0.0_dp
        end if
        y(i)=y(i)+coeff(j)*prev
      end do
    end do
  end subroutine recursive_filter

  subroutine sort_with_index(x, xs, idx)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: xs(:)
    integer, intent(out) :: idx(:)
    integer :: i,j,ti
    real(dp) :: key
    xs=x
    idx=[(i,i=1,size(x))]
    do i=2,size(x)
      key=xs(i)
      ti=idx(i)
      j=i-1
      do while (j >= 1)
        if (xs(j) <= key) exit
        xs(j+1) = xs(j)
        idx(j+1) = idx(j)
        j = j-1
      end do
      xs(j+1)=key
      idx(j+1)=ti
    end do
  end subroutine sort_with_index

  pure real(dp) function normal_quantile(p) result(x)
    ! Peter J. Acklam rational approximation; adequate for statistical cutoffs.
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6)=[-3.969683028665376e1_dp,2.209460984245205e2_dp, &
      -2.759285104469687e2_dp,1.383577518672690e2_dp,-3.066479806614716e1_dp,2.506628277459239_dp]
    real(dp), parameter :: b(5)=[-5.447609879822406e1_dp,1.615858368580409e2_dp, &
      -1.556989798598866e2_dp,6.680131188771972e1_dp,-1.328068155288572e1_dp]
    real(dp), parameter :: c(6)=[-7.784894002430293e-3_dp,-3.223964580411365e-1_dp, &
      -2.400758277161838_dp,-2.549732539343734_dp,4.374664141464968_dp,2.938163982698783_dp]
    real(dp), parameter :: d(4)=[7.784695709041462e-3_dp,3.224671290700398e-1_dp, &
      2.445134137142996_dp,3.754408661907416_dp]
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q,r
    if(p<=0.0_dp) then
    x=-huge(1.0_dp)
    return
    end if
    if(p>=1.0_dp) then
    x= huge(1.0_dp)
    return
    end if
    if(p<plow) then
      q=sqrt(-2.0_dp*log(p))
      x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if(p<=phigh) then
      q=p-0.5_dp
      r=q*q
      x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q=sqrt(-2.0_dp*log(1.0_dp-p))
      x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  function polynomial_convolution(a,b) result(c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp), allocatable :: c(:)
    integer :: i,j
    allocate(c(size(a)+size(b)-1))
    c=0.0_dp
    do i=1,size(a)
    do j=1,size(b)
    c(i+j-1)=c(i+j-1)+a(i)*b(j)
    end do
    end do
  end function polynomial_convolution

  subroutine arma_to_ma(ar, ma, lag_max, psi)
    real(dp), intent(in) :: ar(:),ma(:)
    integer, intent(in) :: lag_max
    real(dp), intent(out) :: psi(:)
    integer :: k,j
    psi=0.0_dp
    do k=1,min(lag_max,size(psi))
      if(k<=size(ma)) psi(k)=ma(k)
      do j=1,min(k,size(ar))
        if(k-j==0) then
          psi(k)=psi(k)+ar(j)
        else
          psi(k)=psi(k)+ar(j)*psi(k-j)
        end if
      end do
    end do
  end subroutine arma_to_ma

  subroutine linear_fit(x,y,beta,resid,status,intercept)
    real(dp), intent(in) :: x(:,:),y(:)
    real(dp), allocatable, intent(out) :: beta(:),resid(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: intercept
    logical :: inc
    real(dp), allocatable :: design(:,:)
    integer :: n,p
    inc=.true.
    if(present(intercept)) inc=intercept
    n=size(x,1)
    p=size(x,2)+merge(1,0,inc)
    allocate(design(n,p),beta(p),resid(n))
    if(inc) then
      design(:,1)=1.0_dp
      if(size(x,2)>0) design(:,2:)=x
    else if(size(x,2)>0) then
      design=x
    end if
    call least_squares(design,y,beta,residuals=resid,status=status)
  end subroutine linear_fit

  subroutine build_lag_matrix(x,p,start,y,design,intercept)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p,start
    real(dp), allocatable, intent(out) :: y(:),design(:,:)
    logical, intent(in), optional :: intercept
    logical :: inc
    integer :: n,i,j,cols
    inc=.true.
    if(present(intercept)) inc=intercept
    n=size(x)-start+1
    cols=p+merge(1,0,inc)
    allocate(y(n),design(n,cols))
    do i=1,n
      y(i)=x(start+i-1)
      if(inc) design(i,1)=1.0_dp
      do j=1,p
        design(i,j+merge(1,0,inc))=x(start+i-1-j)
      end do
    end do
  end subroutine build_lag_matrix

  subroutine sample_with_replacement(x,y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp) :: u
    integer :: i,j
    do i=1,size(y)
      call random_number(u)
      j=min(size(x),1+int(u*real(size(x),dp)))
      y(i)=x(j)
    end do
  end subroutine sample_with_replacement

  real(dp) function quantile_sorted(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: s(:)
    integer, allocatable :: idx(:)
    real(dp) :: h,w
    integer :: lo,hi
    allocate(s(size(x)),idx(size(x)))
    call sort_with_index(x,s,idx)
    if(size(x)==1) then
    q=s(1)
    return
    end if
    h=1.0_dp+(real(size(x)-1,dp))*max(0.0_dp,min(1.0_dp,p))
    lo=floor(h)
    hi=ceiling(h)
    w=h-real(lo,dp)
    q=(1.0_dp-w)*s(lo)+w*s(hi)
  end function quantile_sorted
end module tsa_utils
