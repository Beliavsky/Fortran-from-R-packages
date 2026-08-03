! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_options
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use jrvfinance_kinds, only: dp
  use jrvfinance_types, only: black_scholes_result, root_result, JRV_OK, &
    JRV_INVALID_ARGUMENT
  use jrvfinance_roots, only: newton_raphson_root
  implicit none
  private
  public :: gen_bs, gen_bs_implied, gen_bs_implied_guess

  real(dp), parameter :: pi = acos(-1.0_dp)
  type :: implied_context
    real(dp) :: spot, strike, time, target
  end type implied_context
contains
  function gen_bs(s,strike,rate,sigma,time,div_yield) result(res)
    real(dp),intent(in)::s,strike,rate,sigma,time
    real(dp),intent(in),optional::div_yield
    type(black_scholes_result)::res
    real(dp)::q,g,sqrt_t,a,b,c,pdf
    q=0.0_dp;if(present(div_yield))q=div_yield
    if(s<0.0_dp.or.strike<0.0_dp.or.sigma<0.0_dp.or.time<0.0_dp)then
      res%status=JRV_INVALID_ARGUMENT;return
    end if
    if(time<=epsilon(1.0_dp))then
      res%call=max(s-strike,0.0_dp);res%put=max(strike-s,0.0_dp)
      if(s>strike)then;res%call_delta=1.0_dp;else if(s<strike)then;res%put_delta=-1.0_dp;end if
      res%status=JRV_OK;return
    end if
    if(s<=tiny(1.0_dp).or.strike<=tiny(1.0_dp).or.sigma<=tiny(1.0_dp))then
      call deterministic_option(s,strike,rate,time,q,res);return
    end if
    g=rate-q;sqrt_t=sqrt(time)
    res%d1=(log(s/strike)+(g+0.5_dp*sigma*sigma)*time)/(sigma*sqrt_t)
    res%d2=res%d1-sigma*sqrt_t
    res%nd1=normal_cdf(res%d1);res%nd2=normal_cdf(res%d2)
    res%nminusd1=normal_cdf(-res%d1);res%nminusd2=normal_cdf(-res%d2)
    res%call=exp(-rate*time)*(s*exp(g*time)*res%nd1-strike*res%nd2)
    res%put=exp(-rate*time)*(-s*exp(g*time)*res%nminusd1+strike*res%nminusd2)
    res%call_delta=exp(-q*time)*res%nd1
    res%put_delta=-exp(-q*time)*res%nminusd1
    pdf=normal_pdf(res%d1)
    a=-s*pdf*sigma*exp(-q*time)/(2.0_dp*sqrt_t)
    b=rate*strike*exp(-rate*time)*res%nd2
    c=q*s*res%nd1*exp(-q*time)
    res%call_theta=a-b+c
    b=rate*strike*exp(-rate*time)*res%nminusd2
    c=q*s*res%nminusd1*exp(-q*time)
    res%put_theta=a+b-c
    res%gamma=pdf*exp(-q*time)/(s*sigma*sqrt_t)
    res%vega=s*sqrt_t*pdf*exp(-q*time)
    res%call_rho=strike*time*exp(-rate*time)*res%nd2
    res%put_rho=-strike*time*exp(-rate*time)*res%nminusd2
    res%call_probability=res%nd2;res%put_probability=res%nminusd2
    res%status=JRV_OK
  end function gen_bs

  subroutine deterministic_option(s,strike,rate,time,q,res)
    real(dp),intent(in)::s,strike,rate,time,q
    type(black_scholes_result),intent(out)::res
    real(dp)::ds,dk
    ds=s*exp(-q*time);dk=strike*exp(-rate*time)
    res%call=max(ds-dk,0.0_dp);res%put=max(dk-ds,0.0_dp)
    if(ds>dk)then
      res%call_delta=exp(-q*time);res%call_probability=1.0_dp;res%d1=huge(1.0_dp);res%d2=huge(1.0_dp)
      res%nd1=1.0_dp;res%nd2=1.0_dp
    else if(ds<dk)then
      res%put_delta=-exp(-q*time);res%put_probability=1.0_dp;res%d1=-huge(1.0_dp);res%d2=-huge(1.0_dp)
      res%nminusd1=1.0_dp;res%nminusd2=1.0_dp
    else
      res%d1=0.0_dp;res%d2=0.0_dp;res%nd1=0.5_dp;res%nd2=0.5_dp
      res%nminusd1=0.5_dp;res%nminusd2=0.5_dp
    end if
    res%status=JRV_OK
  end subroutine deterministic_option

  function gen_bs_implied(s,strike,rate,price,time,div_yield,put_option,toler,max_iter,convergence,status) result(sigma)
    real(dp),intent(in)::s,strike,rate,price,time
    real(dp),intent(in),optional::div_yield,toler,convergence
    logical,intent(in),optional::put_option
    integer,intent(in),optional::max_iter
    integer,intent(out),optional::status
    real(dp)::sigma,q,target,spot0,strike0,sminus,splus,guess,upper,ftol,xtol
    logical::is_put
    integer::nmax,i
    type(implied_context)::ctx
    type(root_result)::rr
    q=0.0_dp;is_put=.false.;ftol=1e-6_dp;xtol=1e-8_dp;nmax=100
    if(present(div_yield))q=div_yield;if(present(put_option))is_put=put_option
    if(present(toler))ftol=toler;if(present(convergence))xtol=convergence;if(present(max_iter))nmax=max_iter
    if (s < 0.0_dp .or. strike < 0.0_dp .or. time <= 0.0_dp .or. &
        price < 0.0_dp) then
      sigma = nan()
      if (present(status)) status = JRV_INVALID_ARGUMENT
      return
    end if
    strike0=strike*exp(-rate*time);spot0=s*exp(-q*time);target=price
    if(is_put)target=target+(spot0-strike0)
    sminus=spot0-strike0;splus=spot0+strike0
    if (target < max(sminus, 0.0_dp) .or. target > spot0) then
      sigma = nan()
      if (present(status)) status = JRV_INVALID_ARGUMENT
      return
    end if
    if (abs(target-max(sminus, 0.0_dp)) <= ftol .or. target <= ftol) then
      sigma = 0.0_dp
      if (present(status)) status = JRV_OK
      return
    end if
    if(strike0<=tiny(1.0_dp))then;sigma=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
    guess=gen_bs_implied_guess(target,sminus,splus,time);upper=1.0_dp
    ctx%spot=spot0;ctx%strike=strike0;ctx%time=time;ctx%target=target
    do i=1,60
      if(call_price_zero_rate(spot0,strike0,upper,time)>=target)exit
      upper=upper*2.0_dp
    end do
    rr=newton_raphson_root(implied_callback,ctx,max(0.0_dp,min(guess,upper)),0.0_dp,upper,nmax,ftol,xtol)
    if(rr%status==JRV_OK)then;sigma=rr%root;else;sigma=bisect_implied(ctx,0.0_dp,upper,ftol,xtol,nmax,rr%status);end if
    if(present(status))status=rr%status
  end function gen_bs_implied

  pure real(dp) function gen_bs_implied_guess(price,sminusx,splusx,time) result(sigma)
    real(dp),intent(in)::price,sminusx,splusx,time
    real(dp)::root2pi,h,temp,radical,sigma_root_t
    root2pi=sqrt(2.0_dp*pi);h=0.5_dp*splusx;temp=price-0.5_dp*sminusx
    radical=temp*temp-sminusx*sminusx/pi
    if (radical < 0.0_dp) then
      sigma_root_t = (root2pi/h)*temp
    else
      sigma_root_t = (root2pi/h)*(temp/2.0_dp+sqrt(radical))
    end if
    sigma=sigma_root_t/sqrt(time)
  end function gen_bs_implied_guess

  subroutine implied_callback(x,context,value,gradient)
    real(dp),intent(in)::x
    class(*),intent(in)::context
    real(dp),intent(out)::value,gradient
    real(dp)::d1
    select type(ctx=>context)
    type is(implied_context)
      value=call_price_zero_rate(ctx%spot,ctx%strike,x,ctx%time)-ctx%target
      if(x<=0.0_dp)then
        gradient=0.0_dp
      else
        d1=(log(ctx%spot/ctx%strike)+0.5_dp*x*x*ctx%time)/(x*sqrt(ctx%time))
        gradient=ctx%spot*sqrt(ctx%time)*normal_pdf(d1)
      end if
    class default
      value=nan();gradient=nan()
    end select
  end subroutine implied_callback

  function bisect_implied(ctx,lo0,hi0,ftol,xtol,nmax,status) result(root)
    type(implied_context),intent(in)::ctx
    real(dp),intent(in)::lo0,hi0,ftol,xtol
    integer,intent(in)::nmax
    integer,intent(out)::status
    real(dp)::root,lo,hi,fm
    integer::i
    lo=lo0;hi=hi0
    do i=1,max(200,nmax)
      root=0.5_dp*(lo+hi);fm=call_price_zero_rate(ctx%spot,ctx%strike,root,ctx%time)-ctx%target
      if(abs(fm)<=ftol.or.abs(hi-lo)<=xtol*max(1.0_dp,root))then;status=JRV_OK;return;end if
      if(fm<0.0_dp)then;lo=root;else;hi=root;end if
    end do
    status=2
  end function bisect_implied

  pure real(dp) function call_price_zero_rate(s,strike,sigma,time) result(value)
    real(dp),intent(in)::s,strike,sigma,time
    real(dp)::d1,d2
    if(sigma<=0.0_dp)then;value=max(s-strike,0.0_dp);return;end if
    d1=(log(s/strike)+0.5_dp*sigma*sigma*time)/(sigma*sqrt(time));d2=d1-sigma*sqrt(time)
    value=s*normal_cdf(d1)-strike*normal_cdf(d2)
  end function call_price_zero_rate
  pure real(dp) function normal_cdf(x) result(value)
    real(dp),intent(in)::x
    value=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf
  pure real(dp) function normal_pdf(x) result(value)
    real(dp),intent(in)::x
    value=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf
  pure real(dp) function nan() result(x)
    x=ieee_value(0.0_dp,ieee_quiet_nan)
  end function nan
end module jrvfinance_options
