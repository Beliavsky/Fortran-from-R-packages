! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_bonds
  use rq_kinds, only: dp
  use rq_curves, only: discount_curve_t
  implicit none
  private
  public :: bond_result, zero_price_by_yield, zero_yield_by_price
  public :: fixed_rate_bond_from_yield, fixed_rate_bond_from_curve
  public :: fixed_rate_bond_yield, floating_rate_bond_from_curves
  public :: cashflow_present_value

  type :: bond_result
    real(dp) :: npv=0.0_dp
    real(dp) :: clean_price=0.0_dp
    real(dp) :: dirty_price=0.0_dp
    real(dp) :: accrued=0.0_dp
    real(dp) :: yield=0.0_dp
    real(dp) :: duration=0.0_dp
    real(dp) :: modified_duration=0.0_dp
    real(dp) :: convexity=0.0_dp
    real(dp) :: bps=0.0_dp
    real(dp),allocatable :: cashflow_times(:)
    real(dp),allocatable :: cashflows(:)
  end type bond_result
contains
  pure real(dp) function discount_from_yield(y,t,frequency,compounding) result(df)
    real(dp),intent(in)::y,t
    integer,intent(in)::frequency
    character(len=*),intent(in)::compounding
    if(index(lower(compounding),'continuous')==1) then
      df=exp(-y*t)
    else if(index(lower(compounding),'simple')==1) then
      df=1.0_dp/(1.0_dp+y*t)
    else
      df=(1.0_dp+y/real(max(1,frequency),dp))**(-real(max(1,frequency),dp)*t)
    end if
  end function discount_from_yield

  pure real(dp) function zero_price_by_yield(yield,face,maturity,frequency,compounding) result(price)
    real(dp),intent(in)::yield,face,maturity
    integer,intent(in),optional::frequency
    character(len=*),intent(in),optional::compounding
    integer::f
    character(len=16)::c
    f=1; c='compounded'; if(present(frequency)) f=frequency; if(present(compounding)) c=compounding
    price=face*discount_from_yield(yield,maturity,f,c)
  end function zero_price_by_yield

  subroutine zero_yield_by_price(price,face,maturity,yield,status,frequency,compounding)
    real(dp),intent(in)::price,face,maturity
    real(dp),intent(out)::yield
    integer,intent(out)::status
    integer,intent(in),optional::frequency
    character(len=*),intent(in),optional::compounding
    integer::f
    character(len=16)::c
    f=1; c='compounded'; if(present(frequency)) f=frequency; if(present(compounding)) c=compounding
    if(price<=0.0_dp.or.face<=0.0_dp.or.maturity<=0.0_dp) then; status=1; yield=0.0_dp; return; end if
    if(index(lower(c),'continuous')==1) then
      yield=-log(price/face)/maturity
    else if(index(lower(c),'simple')==1) then
      yield=(face/price-1.0_dp)/maturity
    else
      yield=real(f,dp)*((face/price)**(1.0_dp/(real(f,dp)*maturity))-1.0_dp)
    end if
    status=0
  end subroutine zero_yield_by_price

  pure real(dp) function cashflow_present_value(times,cashflows,curve) result(pv)
    real(dp),intent(in)::times(:),cashflows(:)
    type(discount_curve_t),intent(in)::curve
    integer::i
    pv=0.0_dp
    do i=1,min(size(times),size(cashflows)); pv=pv+cashflows(i)*curve%discount(times(i)); end do
  end function cashflow_present_value

  subroutine fixed_rate_bond_from_yield(face,coupon_rate,maturity,frequency,yield,result,redemption, &
                                        settlement_time,compounding)
    real(dp),intent(in)::face,coupon_rate,maturity,yield
    integer,intent(in)::frequency
    type(bond_result),intent(out)::result
    real(dp),intent(in),optional::redemption,settlement_time
    character(len=*),intent(in),optional::compounding
    integer::n,i,f
    real(dp)::red,settle,coupon,pv,dur_num,conv_num,dy,pup,pdn
    character(len=16)::comp
    f=max(1,frequency); red=face; settle=0.0_dp; comp='compounded'
    if(present(redemption)) red=redemption
    if(present(settlement_time)) settle=settlement_time
    if(present(compounding)) comp=compounding
    n=max(1,nint((maturity-settle)*real(f,dp)))
    allocate(result%cashflow_times(n),result%cashflows(n)); coupon=face*coupon_rate/real(f,dp)
    do i=1,n
      result%cashflow_times(i)=settle+real(i,dp)/real(f,dp)
      result%cashflows(i)=coupon
    end do
    result%cashflows(n)=result%cashflows(n)+red
    pv=0.0_dp; dur_num=0.0_dp; conv_num=0.0_dp
    do i=1,n
      pv=pv+result%cashflows(i)*discount_from_yield(yield,result%cashflow_times(i)-settle,f,comp)
      dur_num=dur_num+(result%cashflow_times(i)-settle)*result%cashflows(i)* &
                      discount_from_yield(yield,result%cashflow_times(i)-settle,f,comp)
      conv_num=conv_num+(result%cashflow_times(i)-settle)**2*result%cashflows(i)* &
                        discount_from_yield(yield,result%cashflow_times(i)-settle,f,comp)
    end do
    result%npv=pv; result%dirty_price=100.0_dp*pv/face; result%clean_price=result%dirty_price
    result%yield=yield; result%duration=dur_num/pv
    if(index(lower(comp),'compounded')==1) then
      result%modified_duration=result%duration/(1.0_dp+yield/real(f,dp))
    else
      result%modified_duration=result%duration
    end if
    result%convexity=conv_num/pv
    dy=1.0e-4_dp
    pup=fixed_price_only(face,coupon_rate,maturity,f,yield+dy,red,settle,comp)
    pdn=fixed_price_only(face,coupon_rate,maturity,f,yield-dy,red,settle,comp)
    result%bps=(pdn-pup)*0.5_dp
  end subroutine fixed_rate_bond_from_yield

  subroutine fixed_rate_bond_from_curve(face,coupon_rate,maturity,frequency,curve,result,redemption,settlement_time)
    real(dp),intent(in)::face,coupon_rate,maturity
    integer,intent(in)::frequency
    type(discount_curve_t),intent(in)::curve
    type(bond_result),intent(out)::result
    real(dp),intent(in),optional::redemption,settlement_time
    integer::n,i,f,status
    real(dp)::red,settle,coupon,pv,dur_num
    f=max(1,frequency); red=face; settle=0.0_dp
    if(present(redemption)) red=redemption
    if(present(settlement_time)) settle=settlement_time
    n=max(1,nint((maturity-settle)*real(f,dp)))
    allocate(result%cashflow_times(n),result%cashflows(n)); coupon=face*coupon_rate/real(f,dp)
    do i=1,n
      result%cashflow_times(i)=settle+real(i,dp)/real(f,dp); result%cashflows(i)=coupon
    end do
    result%cashflows(n)=result%cashflows(n)+red
    pv=0.0_dp; dur_num=0.0_dp
    do i=1,n
      pv=pv+result%cashflows(i)*curve%discount(result%cashflow_times(i))
      dur_num=dur_num+(result%cashflow_times(i)-settle)*result%cashflows(i)*curve%discount(result%cashflow_times(i))
    end do
    result%npv=pv; result%dirty_price=100.0_dp*pv/face; result%clean_price=result%dirty_price
    result%duration=dur_num/pv
    call fixed_rate_bond_yield(face,coupon_rate,maturity,f,pv,result%yield,status,red,settle)
    result%modified_duration=result%duration/(1.0_dp+result%yield/real(f,dp))
  end subroutine fixed_rate_bond_from_curve

  subroutine fixed_rate_bond_yield(face,coupon_rate,maturity,frequency,price,yield,status,redemption,settlement_time)
    real(dp),intent(in)::face,coupon_rate,maturity,price
    integer,intent(in)::frequency
    real(dp),intent(out)::yield
    integer,intent(out)::status
    real(dp),intent(in),optional::redemption,settlement_time
    real(dp)::red,settle,lo,hi,mid,flo,fmid
    integer::i
    red=face; settle=0.0_dp; if(present(redemption)) red=redemption; if(present(settlement_time)) settle=settlement_time
    lo=-0.95_dp*real(max(1,frequency),dp); hi=5.0_dp
    flo=fixed_price_only(face,coupon_rate,maturity,frequency,lo,red,settle,'compounded')-price
    if(flo*(fixed_price_only(face,coupon_rate,maturity,frequency,hi,red,settle,'compounded')-price)>0.0_dp) then
      status=1; yield=0.0_dp; return
    end if
    do i=1,200
      mid=0.5_dp*(lo+hi); fmid=fixed_price_only(face,coupon_rate,maturity,frequency,mid,red,settle,'compounded')-price
      if(abs(fmid)<1.0e-10_dp) exit
      if(fmid*flo<=0.0_dp) then; hi=mid; else; lo=mid; flo=fmid; end if
    end do
    yield=mid; status=0
  end subroutine fixed_rate_bond_yield

  subroutine floating_rate_bond_from_curves(face,maturity,frequency,forward_curve,discount_curve,result, &
                                            gearing,spread,cap,floor,redemption)
    real(dp),intent(in)::face,maturity
    integer,intent(in)::frequency
    type(discount_curve_t),intent(in)::forward_curve,discount_curve
    type(bond_result),intent(out)::result
    real(dp),intent(in),optional::gearing,spread,cap,floor,redemption
    integer::n,i,f
    real(dp)::g,s,c,fl,red,t0,t1,rate,cf,pv
    f=max(1,frequency); g=1.0_dp; s=0.0_dp; c=huge(1.0_dp); fl=-huge(1.0_dp); red=face
    if(present(gearing)) g=gearing; if(present(spread)) s=spread; if(present(cap)) c=cap
    if(present(floor)) fl=floor; if(present(redemption)) red=redemption
    n=max(1,nint(maturity*real(f,dp))); allocate(result%cashflow_times(n),result%cashflows(n)); pv=0.0_dp
    t0=0.0_dp
    do i=1,n
      t1=real(i,dp)/real(f,dp); rate=g*forward_curve%forward_rate(t0,t1)+s; rate=min(max(rate,fl),c)
      cf=face*rate*(t1-t0); if(i==n) cf=cf+red
      result%cashflow_times(i)=t1; result%cashflows(i)=cf; pv=pv+cf*discount_curve%discount(t1); t0=t1
    end do
    result%npv=pv; result%dirty_price=100.0_dp*pv/face; result%clean_price=result%dirty_price
  end subroutine floating_rate_bond_from_curves

  pure real(dp) function fixed_price_only(face,coupon_rate,maturity,frequency,yield,redemption,settlement,compounding) result(price)
    real(dp),intent(in)::face,coupon_rate,maturity,yield,redemption,settlement
    integer,intent(in)::frequency
    character(len=*),intent(in)::compounding
    integer::n,i
    real(dp)::coupon,t
    n=max(1,nint((maturity-settlement)*real(frequency,dp))); coupon=face*coupon_rate/real(frequency,dp); price=0.0_dp
    do i=1,n
      t=real(i,dp)/real(frequency,dp); price=price+coupon*discount_from_yield(yield,t,frequency,compounding)
    end do
    price=price+redemption*discount_from_yield(yield,real(n,dp)/real(frequency,dp),frequency,compounding)
  end function fixed_price_only

  pure function lower(text) result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    integer::i,c
    do i=1,len(text); c=iachar(text(i:i)); if(c>=65.and.c<=90) then; out(i:i)=achar(c+32); else; out(i:i)=text(i:i); end if; end do
  end function lower
end module rq_bonds
