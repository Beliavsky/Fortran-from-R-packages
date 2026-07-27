! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_options
  use rq_kinds, only: dp
  use rq_math, only: normal_pdf, normal_cdf, random_normal, seed_random
  implicit none
  private
  public :: option_result, european_option, american_option, binary_option
  public :: barrier_option, geometric_asian_option, arithmetic_asian_mc
  public :: european_implied_volatility, american_implied_volatility
  public :: european_option_array

  type :: option_result
    real(dp) :: value=0.0_dp
    real(dp) :: delta=0.0_dp
    real(dp) :: gamma=0.0_dp
    real(dp) :: vega=0.0_dp
    real(dp) :: theta=0.0_dp
    real(dp) :: rho=0.0_dp
    real(dp) :: dividend_rho=0.0_dp
  end type option_result
contains
  pure function adjusted_spot(spot,rate,div_times,div_amounts) result(sadj)
    real(dp),intent(in)::spot,rate
    real(dp),intent(in),optional::div_times(:),div_amounts(:)
    real(dp)::sadj
    integer::i
    sadj=spot
    if(present(div_times).and.present(div_amounts)) then
      do i=1,min(size(div_times),size(div_amounts))
        if(div_times(i)>0.0_dp) then
          sadj=sadj-div_amounts(i)*exp(-rate*div_times(i))
        end if
      end do
    end if
  end function adjusted_spot

  pure function european_option(option_type,spot,strike,dividend_yield,rate, &
                                 maturity,volatility,discrete_dividend_times, &
                                 discrete_dividends) result(res)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strike,dividend_yield,rate,maturity,volatility
    real(dp),intent(in),optional::discrete_dividend_times(:)
    real(dp),intent(in),optional::discrete_dividends(:)
    type(option_result)::res
    real(dp)::s,dfq,dfr,sigroot,d1,d2
    s=adjusted_spot(spot,rate,discrete_dividend_times,discrete_dividends)
    if(maturity<=0.0_dp .or. volatility<=0.0_dp .or. &
       s<=0.0_dp .or. strike<=0.0_dp) then
      if(is_call(option_type)) then
        res%value=max(s-strike,0.0_dp)
        if(s>strike) res%delta=1.0_dp
      else
        res%value=max(strike-s,0.0_dp)
        if(s<strike) res%delta=-1.0_dp
      end if
      return
    end if
    dfq=exp(-dividend_yield*maturity)
    dfr=exp(-rate*maturity)
    sigroot=volatility*sqrt(maturity)
    d1=(log(s/strike)+(rate-dividend_yield+ &
       0.5_dp*volatility**2)*maturity)/sigroot
    d2=d1-sigroot
    if(is_call(option_type)) then
      res%value=s*dfq*normal_cdf(d1)-strike*dfr*normal_cdf(d2)
      res%delta=dfq*normal_cdf(d1)
      res%theta=-s*dfq*normal_pdf(d1)*volatility/(2.0_dp*sqrt(maturity)) &
                -rate*strike*dfr*normal_cdf(d2) &
                +dividend_yield*s*dfq*normal_cdf(d1)
      res%rho=strike*maturity*dfr*normal_cdf(d2)
      res%dividend_rho=-s*maturity*dfq*normal_cdf(d1)
    else
      res%value=strike*dfr*normal_cdf(-d2)-s*dfq*normal_cdf(-d1)
      res%delta=-dfq*normal_cdf(-d1)
      res%theta=-s*dfq*normal_pdf(d1)*volatility/(2.0_dp*sqrt(maturity)) &
                +rate*strike*dfr*normal_cdf(-d2) &
                -dividend_yield*s*dfq*normal_cdf(-d1)
      res%rho=-strike*maturity*dfr*normal_cdf(-d2)
      res%dividend_rho=s*maturity*dfq*normal_cdf(-d1)
    end if
    res%gamma=dfq*normal_pdf(d1)/(s*sigroot)
    res%vega=s*dfq*normal_pdf(d1)*sqrt(maturity)
  end function european_option

  function american_option(option_type,spot,strike,dividend_yield,rate, &
                           maturity,volatility,steps, &
                           discrete_dividend_times,discrete_dividends) result(res)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strike,dividend_yield,rate,maturity,volatility
    integer,intent(in),optional::steps
    real(dp),intent(in),optional::discrete_dividend_times(:)
    real(dp),intent(in),optional::discrete_dividends(:)
    type(option_result)::res
    integer::n
    real(dp)::s,h,vup,vdn
    n=500
    if(present(steps)) n=max(2,steps)
    s=adjusted_spot(spot,rate,discrete_dividend_times,discrete_dividends)
    res%value=american_value_only(option_type,s,strike,dividend_yield, &
                                  rate,maturity,volatility,n)
    h=max(1.0e-4_dp,1.0e-4_dp*spot)
    vup=american_value_only(option_type,s+h,strike,dividend_yield, &
                            rate,maturity,volatility,n)
    vdn=american_value_only(option_type,max(s-h,1.0e-12_dp),strike, &
                            dividend_yield,rate,maturity,volatility,n)
    res%delta=(vup-vdn)/(2.0_dp*h)
    res%gamma=(vup-2.0_dp*res%value+vdn)/(h*h)
  end function american_option

  function american_value_only(option_type,spot,strike,q,r,t,sigma,n) &
      result(value)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma
    integer,intent(in)::n
    real(dp)::value,dt,u,d,p,disc,snode
    real(dp),allocatable::v(:)
    integer::i,j
    dt=t/real(n,dp)
    u=exp(sigma*sqrt(dt))
    d=1.0_dp/u
    p=(exp((r-q)*dt)-d)/(u-d)
    p=min(max(p,0.0_dp),1.0_dp)
    disc=exp(-r*dt)
    allocate(v(0:n))
    do j=0,n
      snode=spot*u**j*d**(n-j)
      v(j)=payoff(option_type,snode,strike)
    end do
    do i=n-1,0,-1
      do j=0,i
        snode=spot*u**j*d**(i-j)
        v(j)=max(payoff(option_type,snode,strike), &
                 disc*(p*v(j+1)+(1.0_dp-p)*v(j)))
      end do
    end do
    value=v(0)
  end function american_value_only

  function binary_option(binary_type,option_type,exercise_type,spot,strike, &
                         q,r,t,sigma,cash_payoff,steps) result(res)
    character(len=*),intent(in)::binary_type,option_type,exercise_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma
    real(dp),intent(in),optional::cash_payoff
    integer,intent(in),optional::steps
    type(option_result)::res
    real(dp)::cash,sigroot,d1,d2
    integer::n
    cash=1.0_dp
    if(present(cash_payoff)) cash=cash_payoff
    if(index(lower(exercise_type),'amer')==1) then
      n=500
      if(present(steps)) n=max(2,steps)
      res%value=binary_binomial(binary_type,option_type,spot,strike, &
                                q,r,t,sigma,cash,n)
      return
    end if
    sigroot=sigma*sqrt(t)
    d1=(log(spot/strike)+(r-q+0.5_dp*sigma*sigma)*t)/sigroot
    d2=d1-sigroot
    if(index(lower(binary_type),'asset')==1) then
      if(is_call(option_type)) then
        res%value=spot*exp(-q*t)*normal_cdf(d1)
      else
        res%value=spot*exp(-q*t)*normal_cdf(-d1)
      end if
    else
      if(is_call(option_type)) then
        res%value=cash*exp(-r*t)*normal_cdf(d2)
      else
        res%value=cash*exp(-r*t)*normal_cdf(-d2)
      end if
    end if
  end function binary_option

  function binary_binomial(binary_type,option_type,spot,strike,q,r,t, &
                           sigma,cash,n) result(value)
    character(len=*),intent(in)::binary_type,option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma,cash
    integer,intent(in)::n
    real(dp)::value,dt,u,d,p,disc,snode,imm
    real(dp),allocatable::v(:)
    integer::i,j
    dt=t/real(n,dp)
    u=exp(sigma*sqrt(dt))
    d=1.0_dp/u
    p=min(max((exp((r-q)*dt)-d)/(u-d),0.0_dp),1.0_dp)
    disc=exp(-r*dt)
    allocate(v(0:n))
    do j=0,n
      snode=spot*u**j*d**(n-j)
      v(j)=binary_payoff(binary_type,option_type,snode,strike,cash)
    end do
    do i=n-1,0,-1
      do j=0,i
        snode=spot*u**j*d**(i-j)
        imm=binary_payoff(binary_type,option_type,snode,strike,cash)
        v(j)=max(imm,disc*(p*v(j+1)+(1.0_dp-p)*v(j)))
      end do
    end do
    value=v(0)
  end function binary_binomial

  function barrier_option(barrier_type,option_type,spot,strike,q,r,t,sigma, &
                          barrier,rebate,steps) result(res)
    character(len=*),intent(in)::barrier_type,option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma,barrier,rebate
    integer,intent(in),optional::steps
    type(option_result)::res,tmp
    integer::n
    real(dp)::out_zero
    logical::knock_in
    n=1000
    if(present(steps)) n=max(20,steps)
    knock_in=index(lower(barrier_type),'in')>0
    if(knock_in) then
      tmp=european_option(option_type,spot,strike,q,r,t,sigma)
      out_zero=barrier_out_binomial(replace_in_with_out(barrier_type), &
               option_type,spot,strike,q,r,t,sigma,barrier,0.0_dp,n)
      res%value=max(tmp%value-out_zero,0.0_dp)
      if(barrier_already_hit(barrier_type,spot,barrier)) then
        res%value=tmp%value
      end if
    else
      res%value=barrier_out_binomial(barrier_type,option_type,spot, &
                strike,q,r,t,sigma,barrier,rebate,n)
    end if
  end function barrier_option

  function barrier_out_binomial(barrier_type,option_type,spot,strike,q,r,t, &
                                sigma,barrier,rebate,n) result(value)
    character(len=*),intent(in)::barrier_type,option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma,barrier,rebate
    integer,intent(in)::n
    real(dp)::value,dt,u,d,p,disc,snode
    real(dp),allocatable::v(:)
    integer::i,j
    logical::down,hit
    down=index(lower(barrier_type),'down')==1
    hit=(down.and.spot<=barrier).or.((.not.down).and.spot>=barrier)
    if(hit) then
      value=rebate
      return
    end if
    dt=t/real(n,dp)
    u=exp(sigma*sqrt(dt))
    d=1.0_dp/u
    p=min(max((exp((r-q)*dt)-d)/(u-d),0.0_dp),1.0_dp)
    disc=exp(-r*dt)
    allocate(v(0:n))
    do j=0,n
      snode=spot*u**j*d**(n-j)
      hit=(down.and.snode<=barrier).or.((.not.down).and.snode>=barrier)
      if(hit) then
        v(j)=rebate
      else
        v(j)=payoff(option_type,snode,strike)
      end if
    end do
    do i=n-1,0,-1
      do j=0,i
        snode=spot*u**j*d**(i-j)
        hit=(down.and.snode<=barrier).or.((.not.down).and.snode>=barrier)
        if(hit) then
          v(j)=rebate
        else
          v(j)=disc*(p*v(j+1)+(1.0_dp-p)*v(j))
        end if
      end do
    end do
    value=v(0)
  end function barrier_out_binomial

  pure function geometric_asian_option(option_type,spot,strike,q,r,t,sigma) &
      result(res)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma
    type(option_result)::res
    real(dp)::varg,sg,fg,d1,d2,disc
    varg=sigma*sigma*t/3.0_dp
    sg=sqrt(varg)
    fg=spot*exp(0.5_dp*(r-q)*t-sigma*sigma*t/12.0_dp)
    disc=exp(-r*t)
    d1=(log(fg/strike)+0.5_dp*varg)/sg
    d2=d1-sg
    if(is_call(option_type)) then
      res%value=disc*(fg*normal_cdf(d1)-strike*normal_cdf(d2))
    else
      res%value=disc*(strike*normal_cdf(-d2)-fg*normal_cdf(-d1))
    end if
  end function geometric_asian_option

  subroutine arithmetic_asian_mc(option_type,spot,strike,q,r,t,sigma, &
                                 n_fixings,n_paths,value,std_error,seed)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strike,q,r,t,sigma
    integer,intent(in)::n_fixings,n_paths
    real(dp),intent(out)::value,std_error
    integer,intent(in),optional::seed
    integer::i,j,np
    real(dp)::dt,s1,s2,a1,a2,z,p1,p2,sumv,sum2,disc
    if(present(seed)) call seed_random(seed)
    np=max(2,n_paths/2)
    dt=t/real(n_fixings,dp)
    sumv=0.0_dp
    sum2=0.0_dp
    disc=exp(-r*t)
    do i=1,np
      s1=spot
      s2=spot
      a1=0.0_dp
      a2=0.0_dp
      do j=1,n_fixings
        z=random_normal()
        s1=s1*exp((r-q-0.5_dp*sigma*sigma)*dt+sigma*sqrt(dt)*z)
        s2=s2*exp((r-q-0.5_dp*sigma*sigma)*dt-sigma*sqrt(dt)*z)
        a1=a1+s1
        a2=a2+s2
      end do
      p1=disc*payoff(option_type,a1/real(n_fixings,dp),strike)
      p2=disc*payoff(option_type,a2/real(n_fixings,dp),strike)
      sumv=sumv+p1+p2
      sum2=sum2+p1*p1+p2*p2
    end do
    value=sumv/real(2*np,dp)
    std_error=sqrt(max(sum2/real(2*np,dp)-value*value,0.0_dp)/ &
                   real(2*np,dp))
  end subroutine arithmetic_asian_mc

  subroutine european_implied_volatility(option_type,target,spot,strike,q,r, &
                                         t,volatility,status)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::target,spot,strike,q,r,t
    real(dp),intent(out)::volatility
    integer,intent(out)::status
    real(dp)::lo,hi,mid,flo,fhi,fmid
    integer::i
    type(option_result)::tmp
    lo=1.0e-8_dp
    hi=5.0_dp
    tmp=european_option(option_type,spot,strike,q,r,t,lo)
    flo=tmp%value-target
    tmp=european_option(option_type,spot,strike,q,r,t,hi)
    fhi=tmp%value-target
    if(fhi*flo>0.0_dp) then
      status=1
      volatility=0.0_dp
      return
    end if
    mid=0.5_dp*(lo+hi)
    do i=1,150
      mid=0.5_dp*(lo+hi)
      tmp=european_option(option_type,spot,strike,q,r,t,mid)
      fmid=tmp%value-target
      if(abs(fmid)<1.0e-10_dp) exit
      if(fmid*flo<=0.0_dp) then
        hi=mid
      else
        lo=mid
        flo=fmid
      end if
    end do
    volatility=mid
    status=0
  end subroutine european_implied_volatility

  subroutine american_implied_volatility(option_type,target,spot,strike,q,r, &
                                         t,volatility,status,steps)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::target,spot,strike,q,r,t
    real(dp),intent(out)::volatility
    integer,intent(out)::status
    integer,intent(in),optional::steps
    real(dp)::lo,hi,mid,flo,fhi,fmid
    integer::i,n
    n=500
    if(present(steps)) n=max(2,steps)
    lo=1.0e-5_dp
    hi=5.0_dp
    flo=american_value_only(option_type,spot,strike,q,r,t,lo,n)-target
    fhi=american_value_only(option_type,spot,strike,q,r,t,hi,n)-target
    if(fhi*flo>0.0_dp) then
      status=1
      volatility=0.0_dp
      return
    end if
    mid=0.5_dp*(lo+hi)
    do i=1,100
      mid=0.5_dp*(lo+hi)
      fmid=american_value_only(option_type,spot,strike,q,r,t,mid,n)-target
      if(abs(fmid)<1.0e-8_dp) exit
      if(fmid*flo<=0.0_dp) then
        hi=mid
      else
        lo=mid
        flo=fmid
      end if
    end do
    volatility=mid
    status=0
  end subroutine american_implied_volatility

  subroutine european_option_array(option_type,spot,strikes,q,r,maturities, &
                                   volatility,values)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::spot,strikes(:),q,r,maturities(:),volatility
    real(dp),intent(out)::values(size(strikes),size(maturities))
    integer::i,j
    type(option_result)::tmp
    do j=1,size(maturities)
      do i=1,size(strikes)
        tmp=european_option(option_type,spot,strikes(i),q,r, &
                            maturities(j),volatility)
        values(i,j)=tmp%value
      end do
    end do
  end subroutine european_option_array

  pure logical function is_call(option_type) result(ok)
    character(len=*),intent(in)::option_type
    ok=index(lower(option_type),'call')==1
  end function is_call

  pure real(dp) function payoff(option_type,s,k) result(v)
    character(len=*),intent(in)::option_type
    real(dp),intent(in)::s,k
    if(is_call(option_type)) then
      v=max(s-k,0.0_dp)
    else
      v=max(k-s,0.0_dp)
    end if
  end function payoff

  pure real(dp) function binary_payoff(binary_type,option_type,s,k,cash) &
      result(v)
    character(len=*),intent(in)::binary_type,option_type
    real(dp),intent(in)::s,k,cash
    logical::hit
    if(is_call(option_type)) then
      hit=s>k
    else
      hit=s<k
    end if
    if(.not.hit) then
      v=0.0_dp
    else if(index(lower(binary_type),'asset')==1) then
      v=s
    else
      v=cash
    end if
  end function binary_payoff

  pure logical function barrier_already_hit(barrier_type,s,b) result(hit)
    character(len=*),intent(in)::barrier_type
    real(dp),intent(in)::s,b
    if(index(lower(barrier_type),'down')==1) then
      hit=s<=b
    else
      hit=s>=b
    end if
  end function barrier_already_hit

  pure function replace_in_with_out(text) result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    out=lower(text)
    if(index(out,'downin')==1) out='downout'
    if(index(out,'upin')==1) out='upout'
  end function replace_in_with_out

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
end module rq_options
