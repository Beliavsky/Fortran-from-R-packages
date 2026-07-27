! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_hull_white
  use rq_kinds, only: dp
  use rq_math, only: normal_cdf
  use rq_curves, only: discount_curve_t
  implicit none
  private
  public :: hull_white_b, hull_white_discount_bond, hull_white_bond_option
  public :: hull_white_swaption, hull_white_caplet
  public :: hull_white_calibration_result, calibrate_hull_white_caplets

  type :: hull_white_calibration_result
    real(dp) :: a=0.1_dp
    real(dp) :: sigma=0.01_dp
    real(dp) :: rmse=0.0_dp
    integer :: iterations=0
    integer :: status=1
    real(dp),allocatable :: fitted_prices(:)
  end type hull_white_calibration_result
contains
  pure elemental real(dp) function hull_white_b(a,t_now,maturity) result(b)
    real(dp),intent(in)::a,t_now,maturity
    if(abs(a)<1.0e-12_dp) then
      b=maturity-t_now
    else
      b=(1.0_dp-exp(-a*(maturity-t_now)))/a
    end if
  end function hull_white_b

  pure real(dp) function hull_white_discount_bond(curve,a,sigma,t_now,maturity,r_t) &
      result(price)
    type(discount_curve_t),intent(in)::curve
    real(dp),intent(in)::a,sigma,t_now,maturity,r_t
    real(dp)::b,forward_t,variance,a_factor,dt
    if(maturity<=t_now) then
      price=1.0_dp
      return
    end if
    b=hull_white_b(a,t_now,maturity)
    dt=max(1.0e-5_dp,min(1.0e-3_dp,0.1_dp*max(t_now,1.0e-3_dp)))
    if(t_now<=dt) then
      forward_t=curve%zero_rate(max(dt,1.0e-5_dp))
    else
      forward_t=log(curve%discount(t_now-dt)/curve%discount(t_now+dt))/(2.0_dp*dt)
    end if
    if(abs(a)<1.0e-12_dp) then
      variance=0.5_dp*sigma*sigma*t_now*b*b
    else
      variance=sigma*sigma*(1.0_dp-exp(-2.0_dp*a*t_now))*b*b/(4.0_dp*a)
    end if
    a_factor=curve%discount(maturity)/curve%discount(t_now)*exp(b*forward_t-variance)
    price=a_factor*exp(-b*r_t)
  end function hull_white_discount_bond

  pure real(dp) function hull_white_bond_option(option_type,curve,a,sigma, &
      expiry,bond_maturity,strike) result(value)
    character(len=*),intent(in)::option_type
    type(discount_curve_t),intent(in)::curve
    real(dp),intent(in)::a,sigma,expiry,bond_maturity,strike
    real(dp)::p0t,p0s,sigmap,h
    p0t=curve%discount(bond_maturity)
    p0s=curve%discount(expiry)
    if(abs(a)<1.0e-12_dp) then
      sigmap=sigma*sqrt(expiry)*(bond_maturity-expiry)
    else
      sigmap=sigma*sqrt((1.0_dp-exp(-2.0_dp*a*expiry))/(2.0_dp*a))* &
        hull_white_b(a,expiry,bond_maturity)
    end if
    if(sigmap<=1.0e-14_dp) then
      if(is_call(option_type)) then
        value=max(p0t-strike*p0s,0.0_dp)
      else
        value=max(strike*p0s-p0t,0.0_dp)
      end if
      return
    end if
    h=log(p0t/(strike*p0s))/sigmap+0.5_dp*sigmap
    if(is_call(option_type)) then
      value=p0t*normal_cdf(h)-strike*p0s*normal_cdf(h-sigmap)
    else
      value=strike*p0s*normal_cdf(-h+sigmap)-p0t*normal_cdf(-h)
    end if
  end function hull_white_bond_option

  subroutine hull_white_swaption(option_type,curve,a,sigma,expiry,payment_times, &
                                 accruals,fixed_rate,notional,value,status)
    character(len=*),intent(in)::option_type
    type(discount_curve_t),intent(in)::curve
    real(dp),intent(in)::a,sigma,expiry,payment_times(:),accruals(:)
    real(dp),intent(in)::fixed_rate,notional
    real(dp),intent(out)::value
    integer,intent(out)::status
    real(dp),allocatable::cashflows(:),strikes(:)
    real(dp)::lo,hi,mid,flo,fhi,fmid,rstar
    integer::i,n
    n=min(size(payment_times),size(accruals))
    if(n<1) then
      value=0.0_dp
      status=1
      return
    end if
    allocate(cashflows(n),strikes(n))
    cashflows=fixed_rate*accruals(1:n)
    cashflows(n)=cashflows(n)+1.0_dp
    lo=-2.0_dp
    hi=2.0_dp
    flo=bond_sum(lo)-1.0_dp
    fhi=bond_sum(hi)-1.0_dp
    if(flo*fhi>0.0_dp) then
      value=0.0_dp
      status=2
      return
    end if
    mid=0.0_dp
    do i=1,200
      mid=0.5_dp*(lo+hi)
      fmid=bond_sum(mid)-1.0_dp
      if(abs(fmid)<1.0e-12_dp) exit
      if(fmid*flo<=0.0_dp) then
        hi=mid
      else
        lo=mid
        flo=fmid
      end if
    end do
    rstar=mid
    do i=1,n
      strikes(i)=hull_white_discount_bond(curve,a,sigma,expiry, &
                                          payment_times(i),rstar)
    end do
    value=0.0_dp
    if(index(lower(option_type),'receiver')==1) then
      do i=1,n
        value=value+cashflows(i)*hull_white_bond_option('call',curve,a, &
          sigma,expiry,payment_times(i),strikes(i))
      end do
    else
      do i=1,n
        value=value+cashflows(i)*hull_white_bond_option('put',curve,a, &
          sigma,expiry,payment_times(i),strikes(i))
      end do
    end if
    value=notional*value
    status=0
  contains
    real(dp) function bond_sum(r) result(s)
      real(dp),intent(in)::r
      integer::j
      s=0.0_dp
      do j=1,n
        s=s+cashflows(j)*hull_white_discount_bond(curve,a,sigma,expiry, &
                                                  payment_times(j),r)
      end do
    end function bond_sum
  end subroutine hull_white_swaption

  pure real(dp) function hull_white_caplet(curve,a,sigma,start_time,end_time, &
      strike,notional) result(value)
    type(discount_curve_t),intent(in)::curve
    real(dp),intent(in)::a,sigma,start_time,end_time,strike,notional
    real(dp)::delta,bond_strike
    delta=end_time-start_time
    bond_strike=1.0_dp/(1.0_dp+strike*delta)
    value=notional*(1.0_dp+strike*delta)* &
      hull_white_bond_option('put',curve,a,sigma,start_time,end_time,bond_strike)
  end function hull_white_caplet

  subroutine calibrate_hull_white_caplets(curve,start_times,end_times,strikes, &
      market_prices,notionals,result,max_iter)
    type(discount_curve_t),intent(in)::curve
    real(dp),intent(in)::start_times(:),end_times(:),strikes(:)
    real(dp),intent(in)::market_prices(:),notionals(:)
    type(hull_white_calibration_result),intent(out)::result
    integer,intent(in),optional::max_iter
    real(dp)::p(2),trial(2),steps(2),best,score
    integer::iter,j,nmax,n
    n=min(size(start_times),size(end_times),size(strikes), &
          size(market_prices),size(notionals))
    p=[0.1_dp,0.01_dp]
    steps=[0.05_dp,0.005_dp]
    nmax=400
    if(present(max_iter)) nmax=max_iter
    best=sse(p)
    do iter=1,nmax
      do j=1,2
        trial=p
        trial(j)=max(1.0e-5_dp,p(j)+steps(j))
        score=sse(trial)
        if(score<best) then
          p=trial
          best=score
          cycle
        end if
        trial=p
        trial(j)=max(1.0e-5_dp,p(j)-steps(j))
        score=sse(trial)
        if(score<best) then
          p=trial
          best=score
        else
          steps(j)=0.75_dp*steps(j)
        end if
      end do
      if(maxval(steps)<1.0e-8_dp) exit
    end do
    result%a=p(1)
    result%sigma=p(2)
    result%iterations=iter
    allocate(result%fitted_prices(n))
    do j=1,n
      result%fitted_prices(j)=hull_white_caplet(curve,p(1),p(2), &
        start_times(j),end_times(j),strikes(j),notionals(j))
    end do
    result%rmse=sqrt(sum((result%fitted_prices-market_prices(1:n))**2)/ &
      real(n,dp))
    result%status=0
  contains
    real(dp) function sse(x) result(v)
      real(dp),intent(in)::x(:)
      integer::k
      v=0.0_dp
      do k=1,n
        v=v+(hull_white_caplet(curve,x(1),x(2),start_times(k), &
          end_times(k),strikes(k),notionals(k))-market_prices(k))**2
      end do
    end function sse
  end subroutine calibrate_hull_white_caplets

  pure logical function is_call(option_type) result(ok)
    character(len=*),intent(in)::option_type
    ok=index(lower(option_type),'call')==1
  end function is_call

  pure function lower(text) result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    integer::i,c
    do i=1,len(text)
      c=iachar(text(i:i))
      if(c>=65.and.c<=90) then
        out(i:i)=achar(c+32)
      else
        out(i:i)=text(i:i)
      end if
    end do
  end function lower
end module rq_hull_white
