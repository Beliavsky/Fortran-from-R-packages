! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_bonds
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use jrvfinance_kinds, only: dp
  use jrvfinance_types, only: date_t, bond_cashflows, root_result, JRV_OK, &
    JRV_INVALID_ARGUMENT, JRV_SIZE_MISMATCH
  use jrvfinance_dates, only: edate, year_fraction, daycount_actual, &
    operator(<=), operator(>), operator(>=)
  use jrvfinance_cashflows, only: equiv_rate, duration
  use jrvfinance_roots, only: irr_solve
  implicit none
  private
  public :: coupons_n, coupons_next, coupons_prev, coupons_dates
  public :: bond_tcf, bond_price, bond_yield, bond_duration
  public :: bond_prices, bond_yields, bond_durations

  type :: bond_yield_context
    real(dp), allocatable :: time(:), cashflow(:)
    real(dp) :: price, accrued, freq
  end type bond_yield_context
contains
  integer function coupons_n(settle,mature,freq,status) result(ncoupon)
    type(date_t),intent(in)::settle,mature
    integer,intent(in),optional::freq
    integer,intent(out),optional::status
    integer::f,n,m
    f=2;if(present(freq))f=freq
    if(f<=0.or.mod(12,f)/=0.or.settle>=mature)then
      ncoupon=0;if(present(status))status=JRV_INVALID_ARGUMENT;return
    end if
    n=int(real(f,dp)*real(date_days(settle,mature),dp)/365.25_dp)
    m=-12/f
    do while(edate(mature,n*m)<=settle);n=n-1;end do
    do while(edate(mature,(n+1)*m)>settle);n=n+1;end do
    ncoupon=n+1
    if(present(status))status=JRV_OK
  end function coupons_n

  function coupons_next(settle,mature,freq,status) result(d)
    type(date_t),intent(in)::settle,mature
    integer,intent(in),optional::freq
    integer,intent(out),optional::status
    type(date_t)::d
    integer::f,n,st
    f=2;if(present(freq))f=freq;n=coupons_n(settle,mature,f,st)
    if(st/=JRV_OK)then;d=date_t();if(present(status))status=st;return;end if
    d=edate(mature,(-12/f)*(n-1));if(present(status))status=JRV_OK
  end function coupons_next

  function coupons_prev(settle,mature,freq,status) result(d)
    type(date_t),intent(in)::settle,mature
    integer,intent(in),optional::freq
    integer,intent(out),optional::status
    type(date_t)::d
    integer::f,n,st
    f=2;if(present(freq))f=freq;n=coupons_n(settle,mature,f,st)
    if(st/=JRV_OK)then;d=date_t();if(present(status))status=st;return;end if
    d=edate(mature,(-12/f)*n);if(present(status))status=JRV_OK
  end function coupons_prev

  function coupons_dates(settle,mature,freq,status) result(dates)
    type(date_t),intent(in)::settle,mature
    integer,intent(in),optional::freq
    integer,intent(out),optional::status
    type(date_t),allocatable::dates(:)
    integer::f,n,i,st,m
    f=2;if(present(freq))f=freq;n=coupons_n(settle,mature,f,st)
    if(st/=JRV_OK)then;allocate(dates(0));if(present(status))status=st;return;end if
    allocate(dates(n));m=-12/f
    do i=1,n;dates(i)=edate(mature,m*(n-i));end do
    if(present(status))status=JRV_OK
  end function coupons_dates

  function bond_tcf(settle,mature,coupon,freq,convention,redemption_value) result(tcf)
    type(date_t),intent(in)::settle,mature
    real(dp),intent(in)::coupon
    integer,intent(in),optional::freq
    character(len=*),intent(in),optional::convention
    real(dp),intent(in),optional::redemption_value
    type(bond_cashflows)::tcf
    type(date_t)::nextc,prevc
    integer::f,n,i,st
    real(dp)::redemption,start
    character(len=16)::conv
    f=2;redemption=100.0_dp;conv='30/360'
    if(present(freq))f=freq;if(present(redemption_value))redemption=redemption_value
    if(present(convention))conv=convention
    n=coupons_n(settle,mature,f,st)
    if(st/=JRV_OK)then;allocate(tcf%time(0),tcf%cashflow(0));tcf%status=st;return;end if
    nextc=coupons_next(settle,mature,f);prevc=coupons_prev(settle,mature,f)
    tcf%accrued=100.0_dp*coupon*year_fraction(prevc,settle,prevc,nextc,f,conv)
    start=year_fraction(settle,nextc,prevc,nextc,f,conv)
    allocate(tcf%time(n),tcf%cashflow(n))
    do i=1,n;tcf%time(i)=start+real(i-1,dp)/real(f,dp);end do
    tcf%cashflow=coupon*100.0_dp/real(f,dp)
    tcf%cashflow(n)=redemption+coupon*100.0_dp/real(f,dp)
    tcf%status=JRV_OK
  end function bond_tcf

  function bond_price(settle,mature,coupon,yield_rate,freq,convention,comp_freq,redemption_value,status) result(price)
    type(date_t),intent(in)::settle,mature
    real(dp),intent(in)::coupon,yield_rate
    integer,intent(in),optional::freq
    character(len=*),intent(in),optional::convention
    real(dp),intent(in),optional::comp_freq,redemption_value
    integer,intent(out),optional::status
    real(dp)::price,y,compf
    integer::f
    type(bond_cashflows)::tcf
    f=2;compf=real(f,dp);if(present(freq))f=freq;if(.not.present(comp_freq))compf=real(f,dp)
    if(present(comp_freq))compf=comp_freq
    tcf=bond_tcf(settle,mature,coupon,f,convention,redemption_value)
    if(tcf%status/=JRV_OK)then;price=nan();if(present(status))status=tcf%status;return;end if
    if(size(tcf%time)==1)then
      price=tcf%cashflow(1)/(1.0_dp+yield_rate*tcf%time(1))-tcf%accrued
    else
      y=yield_rate;if(abs(compf-real(f,dp))>epsilon(1.0_dp))y=equiv_rate(y,compf,real(f,dp))
      price=sum(tcf%cashflow*(1.0_dp+y/real(f,dp))**(-real(f,dp)*tcf%time))-tcf%accrued
    end if
    if(present(status))status=JRV_OK
  end function bond_price

  function bond_yield(settle,mature,coupon,price,freq,convention,comp_freq,redemption_value,status) result(yield_rate)
    type(date_t),intent(in)::settle,mature
    real(dp),intent(in)::coupon,price
    integer,intent(in),optional::freq
    character(len=*),intent(in),optional::convention
    real(dp),intent(in),optional::comp_freq,redemption_value
    integer,intent(out),optional::status
    real(dp)::yield_rate,compf
    integer::f
    type(bond_cashflows)::tcf
    type(bond_yield_context)::ctx
    type(root_result)::rr
    f=2;if(present(freq))f=freq;compf=real(f,dp);if(present(comp_freq))compf=comp_freq
    tcf=bond_tcf(settle,mature,coupon,f,convention,redemption_value)
    if(tcf%status/=JRV_OK)then;yield_rate=nan();if(present(status))status=tcf%status;return;end if
    if(size(tcf%time)==1)then
      yield_rate=(tcf%cashflow(1)/(price+tcf%accrued)-1.0_dp)/tcf%time(1)
      if(present(status))status=JRV_OK;return
    end if
    ctx%time=tcf%time;ctx%cashflow=tcf%cashflow;ctx%price=price;ctx%accrued=tcf%accrued;ctx%freq=real(f,dp)
    rr=irr_solve(bond_yield_callback,ctx)
    if(rr%status==JRV_OK)then;yield_rate=equiv_rate(rr%root,real(f,dp),compf);else;yield_rate=nan();end if
    if(present(status))status=rr%status
  end function bond_yield

  subroutine bond_yield_callback(r,context,value,gradient)
    real(dp),intent(in)::r
    class(*),intent(in)::context
    real(dp),intent(out)::value,gradient
    real(dp),allocatable::df(:)
    select type(ctx=>context)
    type is(bond_yield_context)
      if(1.0_dp+r/ctx%freq<=0.0_dp)then;value=huge(1.0_dp);gradient=-huge(1.0_dp);return;end if
      allocate(df(size(ctx%time)));df=(1.0_dp+r/ctx%freq)**(-ctx%freq*ctx%time)
      value=sum(ctx%cashflow*df)-ctx%accrued-ctx%price
      gradient=-sum(ctx%cashflow*df*ctx%time)/(1.0_dp+r/ctx%freq)
    class default
      value=nan();gradient=nan()
    end select
  end subroutine bond_yield_callback

  function bond_duration(settle, mature, coupon, yield_rate, freq, convention, &
      modified, comp_freq, redemption_value, status) result(value)
    type(date_t),intent(in)::settle,mature
    real(dp),intent(in)::coupon,yield_rate
    integer,intent(in),optional::freq
    character(len=*),intent(in),optional::convention
    logical,intent(in),optional::modified
    real(dp),intent(in),optional::comp_freq,redemption_value
    integer,intent(out),optional::status
    real(dp)::value,compf
    integer::f
    logical::modif
    type(bond_cashflows)::tcf
    f=2;if(present(freq))f=freq;compf=real(f,dp);if(present(comp_freq))compf=comp_freq
    modif=.false.;if(present(modified))modif=modified
    tcf=bond_tcf(settle,mature,coupon,f,convention,redemption_value)
    if(tcf%status/=JRV_OK)then;value=nan();if(present(status))status=tcf%status;return;end if
    if(size(tcf%time)==1)then
      value=real(daycount_actual(settle,mature),dp)/365.0_dp
    else
      value=duration(tcf%cashflow,yield_rate,real(f,dp),compf,tcf%time,.false.,modif)
    end if
    if(present(status))status=JRV_OK
  end function bond_duration

  subroutine bond_prices(settle,mature,coupon,yield_rate,prices,freq,convention,comp_freq,redemption_value,status)
    type(date_t),intent(in)::settle(:),mature(:)
    real(dp),intent(in)::coupon(:),yield_rate(:)
    real(dp),allocatable,intent(out)::prices(:)
    integer,intent(in),optional::freq(:)
    character(len=*),intent(in),optional::convention
    real(dp),intent(in),optional::comp_freq(:),redemption_value(:)
    integer,intent(out),optional::status
    integer::n,i,f,st
    real(dp)::cf,rv
    n=size(settle);allocate(prices(n))
    if (size(mature) /= n .or. size(coupon) /= n .or. &
        size(yield_rate) /= n) then
      prices = nan()
      if (present(status)) status = JRV_SIZE_MISMATCH
      return
    end if
    do i=1,n
      f=2;if(present(freq))f=freq(min(i,size(freq)))
      cf=real(f,dp);if(present(comp_freq))cf=comp_freq(min(i,size(comp_freq)))
      rv=100.0_dp;if(present(redemption_value))rv=redemption_value(min(i,size(redemption_value)))
      prices(i)=bond_price(settle(i),mature(i),coupon(i),yield_rate(i),f,convention,cf,rv,st)
      if(st/=JRV_OK)then;if(present(status))status=st;return;end if
    end do
    if(present(status))status=JRV_OK
  end subroutine bond_prices

  subroutine bond_yields(settle,mature,coupon,price,yields,freq,convention,comp_freq,redemption_value,status)
    type(date_t),intent(in)::settle(:),mature(:)
    real(dp),intent(in)::coupon(:),price(:)
    real(dp),allocatable,intent(out)::yields(:)
    integer,intent(in),optional::freq(:)
    character(len=*),intent(in),optional::convention
    real(dp),intent(in),optional::comp_freq(:),redemption_value(:)
    integer,intent(out),optional::status
    integer::n,i,f,st
    real(dp)::cf,rv
    n=size(settle);allocate(yields(n))
    if (size(mature) /= n .or. size(coupon) /= n .or. size(price) /= n) then
      yields = nan()
      if (present(status)) status = JRV_SIZE_MISMATCH
      return
    end if
    do i=1,n
      f=2;if(present(freq))f=freq(min(i,size(freq)))
      cf=real(f,dp);if(present(comp_freq))cf=comp_freq(min(i,size(comp_freq)))
      rv=100.0_dp;if(present(redemption_value))rv=redemption_value(min(i,size(redemption_value)))
      yields(i)=bond_yield(settle(i),mature(i),coupon(i),price(i),f,convention,cf,rv,st)
      if(st/=JRV_OK)then;if(present(status))status=st;return;end if
    end do
    if(present(status))status=JRV_OK
  end subroutine bond_yields

  subroutine bond_durations(settle, mature, coupon, yield_rate, durations, &
      freq, convention, modified, comp_freq, redemption_value, status)
    type(date_t),intent(in)::settle(:),mature(:)
    real(dp),intent(in)::coupon(:),yield_rate(:)
    real(dp),allocatable,intent(out)::durations(:)
    integer,intent(in),optional::freq(:)
    character(len=*),intent(in),optional::convention
    logical,intent(in),optional::modified
    real(dp),intent(in),optional::comp_freq(:),redemption_value(:)
    integer,intent(out),optional::status
    integer::n,i,f,st
    real(dp)::cf,rv
    logical::modif
    n=size(settle);allocate(durations(n));modif=.false.;if(present(modified))modif=modified
    if (size(mature) /= n .or. size(coupon) /= n .or. &
        size(yield_rate) /= n) then
      durations = nan()
      if (present(status)) status = JRV_SIZE_MISMATCH
      return
    end if
    do i=1,n
      f=2;if(present(freq))f=freq(min(i,size(freq)))
      cf=real(f,dp);if(present(comp_freq))cf=comp_freq(min(i,size(comp_freq)))
      rv=100.0_dp;if(present(redemption_value))rv=redemption_value(min(i,size(redemption_value)))
      durations(i)=bond_duration(settle(i),mature(i),coupon(i),yield_rate(i),f,convention,modif,cf,rv,st)
      if(st/=JRV_OK)then;if(present(status))status=st;return;end if
    end do
    if(present(status))status=JRV_OK
  end subroutine bond_durations

  pure integer function date_days(a,b) result(days)
    use jrvfinance_dates, only: daycount_actual
    type(date_t),intent(in)::a,b
    days=daycount_actual(a,b)
  end function date_days
  pure real(dp) function nan() result(x)
    x=ieee_value(0.0_dp,ieee_quiet_nan)
  end function nan
end module jrvfinance_bonds
