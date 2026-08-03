! SPDX-License-Identifier: Artistic-2.0
module ecd_options
  use ecd_kinds, only : dp, ecd_ok, ecd_invalid, ecd_no_convergence
  use ecd_math, only : normal_cdf, brent_root, nan_dp
  implicit none
  private

  public :: bs_option_price, bs_call_price, bs_put_price, bs_implied_volatility
  public :: polyfit_option, option_intrinsic_value

contains

  pure elemental function option_intrinsic_value(strike, spot, option_type) result(v)
    real(dp), intent(in) :: strike, spot
    character(len=*), intent(in) :: option_type
    real(dp) :: v
    if (option_type(1:1)=='p' .or. option_type(1:1)=='P') then
      v=max(strike-spot,0.0_dp)
    else
      v=max(spot-strike,0.0_dp)
    end if
  end function option_intrinsic_value

  pure elemental function bs_option_price(volatility, strike, spot, ttm, rate, dividend, option_type) result(v)
    real(dp), intent(in) :: volatility, strike, spot, ttm
    real(dp), intent(in), optional :: rate, dividend
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v,r,q,sqrt_t,d1,d2,sd,kd,sgn
    character(len=1) :: ot
    r=0.0_dp; q=0.0_dp; ot='c'
    if(present(rate))r=rate
    if(present(dividend))q=dividend
    if(present(option_type))ot=option_type(1:1)
    if(strike<=0.0_dp .or. spot<=0.0_dp .or. ttm<0.0_dp .or. volatility<0.0_dp) then
      v=nan_dp(); return
    end if
    sd=spot*exp(-q*ttm); kd=strike*exp(-r*ttm)
    sgn=merge(-1.0_dp,1.0_dp,ot=='p'.or.ot=='P')
    if(ttm==0.0_dp .or. volatility==0.0_dp) then
      v=max(sgn*(sd-kd),0.0_dp); return
    end if
    sqrt_t=sqrt(ttm)
    d1=(log(spot/strike)+(r-q+0.5_dp*volatility**2)*ttm)/(volatility*sqrt_t)
    d2=d1-volatility*sqrt_t
    v=sgn*(sd*normal_cdf(sgn*d1)-kd*normal_cdf(sgn*d2))
  end function bs_option_price

  pure elemental function bs_call_price(volatility,strike,spot,ttm,rate,dividend) result(v)
    real(dp), intent(in) :: volatility,strike,spot,ttm
    real(dp), intent(in), optional :: rate,dividend
    real(dp) :: v,r,q
    r=0.0_dp; q=0.0_dp
    if(present(rate))r=rate
    if(present(dividend))q=dividend
    v=bs_option_price(volatility,strike,spot,ttm,r,q,'c')
  end function bs_call_price

  pure elemental function bs_put_price(volatility,strike,spot,ttm,rate,dividend) result(v)
    real(dp), intent(in) :: volatility,strike,spot,ttm
    real(dp), intent(in), optional :: rate,dividend
    real(dp) :: v,r,q
    r=0.0_dp; q=0.0_dp
    if(present(rate))r=rate
    if(present(dividend))q=dividend
    v=bs_option_price(volatility,strike,spot,ttm,r,q,'p')
  end function bs_put_price

  function bs_implied_volatility(price,strike,spot,ttm,rate,dividend,option_type,status) result(volatility)
    real(dp), intent(in) :: price,strike,spot,ttm
    real(dp), intent(in), optional :: rate,dividend
    character(len=*), intent(in), optional :: option_type
    integer, intent(out), optional :: status
    real(dp) :: volatility,r,q,lo,hi,flo,fhi
    integer :: st
    character(len=1) :: ot
    r=0.0_dp; q=0.0_dp; ot='c'
    if(present(rate))r=rate
    if(present(dividend))q=dividend
    if(present(option_type))ot=option_type(1:1)
    if(present(status))status=ecd_ok
    if(price<0.0_dp .or. strike<=0.0_dp .or. spot<=0.0_dp .or. ttm<=0.0_dp) then
      volatility=nan_dp(); if(present(status))status=ecd_invalid; return
    end if
    lo=1.0e-8_dp; hi=4.0_dp
    flo=objective(lo); fhi=objective(hi)
    do while(flo*fhi>0.0_dp .and. hi<64.0_dp)
      hi=2.0_dp*hi; fhi=objective(hi)
    end do
    if(flo*fhi>0.0_dp) then
      volatility=nan_dp(); if(present(status))status=ecd_no_convergence; return
    end if
    volatility=brent_root(objective,lo,hi,1.0e-10_dp,300,st)
    if(present(status))status=st
  contains
    function objective(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y=bs_option_price(x,strike,spot,ttm,r,q,ot)-price
    end function objective
  end function bs_implied_volatility

  subroutine polyfit_option(k,price,k_cusp,k_new,price_new,degree_left,degree_right,status)
    real(dp), intent(in) :: k(:),price(:),k_cusp,k_new(:)
    real(dp), intent(out) :: price_new(:)
    integer, intent(in), optional :: degree_left,degree_right
    integer, intent(out), optional :: status
    integer :: dl,dr,nl,nr,i,st1,st2
    real(dp), allocatable :: xl(:),yl(:),xr(:),yr(:),cl(:),cr(:)
    dl=6; dr=6
    if(present(degree_left))dl=degree_left
    if(present(degree_right))dr=degree_right
    if(present(status))status=ecd_ok
    if(size(k)/=size(price) .or. size(price_new)/=size(k_new) .or. any(price<=0.0_dp)) then
      price_new=nan_dp(); if(present(status))status=ecd_invalid; return
    end if
    nl=count(k<k_cusp); nr=size(k)-nl
    if(nl<dl+1 .or. nr<dr+1) then
      price_new=nan_dp(); if(present(status))status=ecd_invalid; return
    end if
    allocate(xl(nl),yl(nl),xr(nr),yr(nr),cl(dl+1),cr(dr+1))
    nl=0; nr=0
    do i=1,size(k)
      if(k(i)<k_cusp) then
        nl=nl+1; xl(nl)=k(i); yl(nl)=log(price(i))
      else
        nr=nr+1; xr(nr)=k(i); yr(nr)=log(price(i))
      end if
    end do
    call polynomial_least_squares(xl,yl,dl,cl,st1)
    call polynomial_least_squares(xr,yr,dr,cr,st2)
    if(st1/=ecd_ok .or. st2/=ecd_ok) then
      price_new=nan_dp(); if(present(status))status=ecd_no_convergence; return
    end if
    do i=1,size(k_new)
      if(k_new(i)<k_cusp) then
        price_new(i)=exp(polyval(cl,k_new(i)))
      else
        price_new(i)=exp(polyval(cr,k_new(i)))
      end if
    end do
  end subroutine polyfit_option

  subroutine polynomial_least_squares(x,y,degree,c,status)
    real(dp), intent(in) :: x(:),y(:)
    integer, intent(in) :: degree
    real(dp), intent(out) :: c(:)
    integer, intent(out) :: status
    real(dp), allocatable :: a(:,:),b(:),powers(:)
    integer :: n,p,i,j,l
    n=size(x); p=degree+1
    allocate(a(p,p),b(p),powers(0:2*degree))
    a=0.0_dp; b=0.0_dp
    do i=1,n
      powers(0)=1.0_dp
      do j=1,2*degree; powers(j)=powers(j-1)*x(i); end do
      do j=1,p
        b(j)=b(j)+y(i)*powers(j-1)
        do l=1,p
          a(j,l)=a(j,l)+powers(j+l-2)
        end do
      end do
    end do
    call solve_dense(a,b,c,status)
  end subroutine polynomial_least_squares

  subroutine solve_dense(a,b,x,status)
    real(dp), intent(inout) :: a(:,:)
    real(dp), intent(inout) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    integer :: n,i,j,k,p
    real(dp) :: pivot,tmp,factor
    n=size(b); status=ecd_ok
    do k=1,n
      p=k
      do i=k+1,n
        if(abs(a(i,k))>abs(a(p,k)))p=i
      end do
      if(abs(a(p,k))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
        status=ecd_no_convergence; x=nan_dp(); return
      end if
      if(p/=k) then
        do j=k,n
          tmp=a(k,j); a(k,j)=a(p,j); a(p,j)=tmp
        end do
        tmp=b(k); b(k)=b(p); b(p)=tmp
      end if
      pivot=a(k,k)
      do i=k+1,n
        factor=a(i,k)/pivot
        a(i,k:n)=a(i,k:n)-factor*a(k,k:n)
        b(i)=b(i)-factor*b(k)
      end do
    end do
    do i=n,1,-1
      x(i)=(b(i)-dot_product(a(i,i+1:n),x(i+1:n)))/a(i,i)
    end do
  end subroutine solve_dense

  pure function polyval(c,x) result(y)
    real(dp), intent(in) :: c(:),x
    real(dp) :: y
    integer :: i
    y=c(size(c))
    do i=size(c)-1,1,-1; y=y*x+c(i); end do
  end function polyval

end module ecd_options
