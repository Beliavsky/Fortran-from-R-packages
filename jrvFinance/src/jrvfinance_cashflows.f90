! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_cashflows
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_quiet_nan
  use jrvfinance_kinds, only: dp
  use jrvfinance_types, only: annuity_breakup_result, root_result, JRV_OK, &
    JRV_INVALID_ARGUMENT, JRV_NO_ROOT
  use jrvfinance_roots, only: irr_solve
  implicit none
  private
  public :: equiv_rate, npv, irr, duration
  public :: annuity_pv, annuity_fv, annuity_instalment
  public :: annuity_periods, annuity_rate, annuity_instalment_breakup

  type :: irr_context
    real(dp), allocatable :: cf(:), time(:)
  end type irr_context
  type :: annuity_rate_context
    real(dp) :: n_periods, instalment, pv, fv, terminal_payment
    logical :: immediate_start
  end type annuity_rate_context
contains
  pure real(dp) function equiv_rate(rate, from_freq, to_freq) result(converted)
    real(dp), intent(in) :: rate
    real(dp), intent(in), optional :: from_freq, to_freq
    real(dp) :: from_f, to_f, cc_rate
    from_f=1.0_dp; to_f=1.0_dp
    if(present(from_freq)) from_f=from_freq
    if(present(to_freq)) to_f=to_freq
    if(from_f<=0.0_dp) then
      cc_rate=rate
    else
      cc_rate=log(1.0_dp+rate/from_f)*from_f
    end if
    if(to_f<=0.0_dp) then
      converted=cc_rate
    else
      converted=(exp(cc_rate/to_f)-1.0_dp)*to_f
    end if
  end function equiv_rate

  function npv(cf, rate, cf_freq, comp_freq, cf_t, immediate_start, status) result(value)
    real(dp), intent(in) :: cf(:), rate
    real(dp), intent(in), optional :: cf_freq, comp_freq, cf_t(:)
    logical, intent(in), optional :: immediate_start
    integer, intent(out), optional :: status
    real(dp) :: value, cff, compf, start, cc_rate
    real(dp), allocatable :: time(:)
    integer :: i
    logical :: immediate
    cff=1.0_dp; compf=1.0_dp; immediate=.false.
    if(present(cf_freq)) cff=cf_freq
    if(present(comp_freq)) compf=comp_freq
    if(present(immediate_start)) immediate=immediate_start
    if(cff<=0.0_dp .or. size(cf)==0) then
      value=nan(); if(present(status)) status=JRV_INVALID_ARGUMENT; return
    end if
    allocate(time(size(cf)))
    if(present(cf_t)) then
      if(size(cf_t)/=size(cf)) then
        value=nan(); if(present(status)) status=JRV_INVALID_ARGUMENT; return
      end if
      time=cf_t
    else
      start=merge(0.0_dp,1.0_dp/cff,immediate)
      do i=1,size(cf); time(i)=start+real(i-1,dp)/cff; end do
    end if
    cc_rate=equiv_rate(rate,compf,0.0_dp)
    value=sum(cf*exp(-cc_rate*time))
    if(present(status)) status=JRV_OK
  end function npv

  function irr(cf, interval, cf_freq, comp_freq, cf_t, r_guess, toler, convergence, &
      max_iter, method, status) result(rate)
    real(dp), intent(in) :: cf(:)
    real(dp), intent(in), optional :: interval(2), cf_freq, comp_freq, cf_t(:), r_guess, toler, convergence
    integer, intent(in), optional :: max_iter
    character(len=*), intent(in), optional :: method
    integer, intent(out), optional :: status
    real(dp) :: rate, cff, compf
    type(irr_context) :: ctx
    type(root_result) :: rr
    integer :: i
    cff=1.0_dp; compf=1.0_dp
    if(present(cf_freq)) cff=cf_freq
    if(present(comp_freq)) compf=comp_freq
    if(size(cf)==0 .or. cff<=0.0_dp .or. maxval(cf)*minval(cf)>=0.0_dp) then
      rate=nan(); if(present(status)) status=JRV_NO_ROOT; return
    end if
    allocate(ctx%cf(size(cf)),ctx%time(size(cf))); ctx%cf=cf
    if(present(cf_t)) then
      if(size(cf_t)/=size(cf)) then; rate=nan(); if(present(status)) status=JRV_INVALID_ARGUMENT; return; end if
      ctx%time=cf_t
    else
      do i=1,size(cf); ctx%time(i)=real(i-1,dp)/cff; end do
    end if
    if(present(interval)) then
      if(present(r_guess)) then
        rr=call_irr_solve(ctx, interval=interval, r_guess=r_guess, &
          toler=toler, convergence=convergence, max_iter=max_iter, method=method)
      else
        rr=call_irr_solve(ctx,interval=interval,toler=toler,convergence=convergence,max_iter=max_iter,method=method)
      end if
    else if(present(r_guess)) then
      rr=call_irr_solve(ctx,r_guess=r_guess,toler=toler,convergence=convergence,max_iter=max_iter,method=method)
    else
      rr=call_irr_solve(ctx,toler=toler,convergence=convergence,max_iter=max_iter,method=method)
    end if
    if(rr%status==JRV_OK) then
      rate=equiv_rate(rr%root,0.0_dp,compf)
    else
      rate=nan()
    end if
    if(present(status)) status=rr%status
  end function irr

  function call_irr_solve(ctx,interval,r_guess,toler,convergence,max_iter,method) result(rr)
    type(irr_context), intent(in) :: ctx
    real(dp),intent(in),optional::interval(2),r_guess,toler,convergence
    integer,intent(in),optional::max_iter
    character(len=*),intent(in),optional::method
    type(root_result)::rr
    real(dp)::iv(2),rg,ft,xt
    integer::mi
    character(len=16)::md
    iv=[-1.0_dp,sqrt(real(huge(1),dp))]; if(present(interval)) iv=interval
    rg=merge(0.5_dp*(iv(1)+iv(2)),0.0_dp,present(interval)); if(present(r_guess)) rg=r_guess
    ft=1e-6_dp; if(present(toler)) ft=toler
    xt=1e-8_dp; if(present(convergence)) xt=convergence
    mi=100; if(present(max_iter)) mi=max_iter
    md='default'; if(present(method)) md=method
    rr=irr_solve(irr_callback,ctx,iv,rg,ft,xt,mi,md)
  end function call_irr_solve

  subroutine irr_callback(x,context,value,gradient)
    real(dp),intent(in)::x
    class(*),intent(in)::context
    real(dp),intent(out)::value,gradient
    real(dp),allocatable::df(:)
    select type(ctx=>context)
    type is(irr_context)
      allocate(df(size(ctx%cf))); df=exp(-x*ctx%time)
      value=sum(ctx%cf*df); gradient=-sum(ctx%cf*df*ctx%time)
    class default
      value=nan(); gradient=nan()
    end select
  end subroutine irr_callback

  function duration(cf,rate,cf_freq,comp_freq,cf_t,immediate_start,modified,status) result(value)
    real(dp),intent(in)::cf(:),rate
    real(dp),intent(in),optional::cf_freq,comp_freq,cf_t(:)
    logical,intent(in),optional::immediate_start,modified
    integer,intent(out),optional::status
    real(dp)::value,cff,compf,start,cc_rate,den
    real(dp),allocatable::time(:),df(:)
    integer::i
    logical::imm,modif
    cff=1.0_dp;compf=1.0_dp;imm=.false.;modif=.false.
    if(present(cf_freq)) cff=cf_freq
    if(present(comp_freq)) compf=comp_freq
    if(present(immediate_start)) imm=immediate_start
    if(present(modified)) modif=modified
    if(size(cf)==0.or.cff<=0.0_dp) then;value=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
    allocate(time(size(cf)),df(size(cf)))
    if(present(cf_t)) then
      if(size(cf_t)/=size(cf)) then;value=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
      time=cf_t
    else
      start=merge(0.0_dp,1.0_dp/cff,imm);do i=1,size(cf);time(i)=start+real(i-1,dp)/cff;end do
    end if
    cc_rate=equiv_rate(rate,compf,0.0_dp);df=exp(-cc_rate*time);den=sum(cf*df)
    if(abs(den)<=tiny(1.0_dp)) then;value=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
    value=sum(cf*df*time)/den
    if(modif) then
      if(compf<=0.0_dp) then
        value=value
      else
        value=value/(1.0_dp+rate/compf)
      end if
    end if
    if(present(status))status=JRV_OK
  end function duration

  function annuity_pv(rate,n_periods,instalment,terminal_payment,immediate_start,cf_freq,comp_freq,status) result(value)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::n_periods,instalment,terminal_payment,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(out),optional::status
    real(dp)::value,n,inst,term,cff,compf,r,df,adjust
    logical::imm
    n=ieee_value(0.0_dp,ieee_positive_inf);inst=1.0_dp;term=0.0_dp;cff=1.0_dp;compf=1.0_dp;imm=.false.
    if(present(n_periods))n=n_periods;if(present(instalment))inst=instalment
    if(present(terminal_payment))term=terminal_payment;if(present(cf_freq))cff=cf_freq
    if(present(comp_freq))compf=comp_freq;if(present(immediate_start))imm=immediate_start
    if(cff<=0.0_dp.or.n<0.0_dp) then;value=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
    if(abs(rate)<=epsilon(1.0_dp)) then;value=n*inst+term;if(present(status))status=JRV_OK;return;end if
    r=equiv_rate(rate,compf,cff)/cff;df=(1.0_dp+r)**(-n);adjust=merge(1.0_dp+r,1.0_dp,imm)
    value=adjust*(inst*(1.0_dp-df)/r+term*df);if(present(status))status=JRV_OK
  end function annuity_pv

  function annuity_fv(rate,n_periods,instalment,terminal_payment,immediate_start,cf_freq,comp_freq,status) result(value)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::n_periods,instalment,terminal_payment,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(out),optional::status
    real(dp)::value,n,inst,term,cff,compf,r,growth,adjust
    logical::imm
    n=ieee_value(0.0_dp,ieee_positive_inf);inst=1.0_dp;term=0.0_dp;cff=1.0_dp;compf=1.0_dp;imm=.false.
    if(present(n_periods))n=n_periods;if(present(instalment))inst=instalment
    if(present(terminal_payment))term=terminal_payment;if(present(cf_freq))cff=cf_freq
    if(present(comp_freq))compf=comp_freq;if(present(immediate_start))imm=immediate_start
    if(abs(rate)<=epsilon(1.0_dp)) then;value=n*inst;if(present(status))status=JRV_OK;return;end if
    r=equiv_rate(rate,compf,cff)/cff;growth=(1.0_dp+r)**n;adjust=merge(1.0_dp+r,1.0_dp,imm)
    value=adjust*inst*(growth-1.0_dp)/r+term;if(present(status))status=JRV_OK
  end function annuity_fv

  function annuity_instalment(rate, n_periods, pv, fv, terminal_payment, &
      immediate_start, cf_freq, comp_freq, status) result(value)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::n_periods,pv,fv,terminal_payment,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(out),optional::status
    real(dp)::value,n,pv0,fv0,term,cff,compf,r,df,lhs,adjust
    logical::imm
    n=ieee_value(0.0_dp,ieee_positive_inf);pv0=1.0_dp;fv0=0.0_dp;term=0.0_dp;cff=1.0_dp;compf=1.0_dp;imm=.false.
    if(present(n_periods))n=n_periods;if(present(pv))pv0=pv;if(present(fv))fv0=fv
    if(present(terminal_payment))term=terminal_payment;if(present(cf_freq))cff=cf_freq
    if(present(comp_freq))compf=comp_freq;if(present(immediate_start))imm=immediate_start
    r=equiv_rate(rate,compf,cff)/cff;df=(1.0_dp+r)**(-n);lhs=pv0+(fv0-term)*df
    if (abs(rate) <= epsilon(1.0_dp)) then
      value = lhs/n
    else
      adjust = merge(1.0_dp+r, 1.0_dp, imm)
      value = r*lhs/(adjust*(1.0_dp-df))
    end if
    if(present(status))status=JRV_OK
  end function annuity_instalment

  function annuity_periods(rate, instalment, pv, fv, terminal_payment, &
      immediate_start, cf_freq, comp_freq, round_digits, status) result(value)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::instalment,pv,fv,terminal_payment,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(in),optional::round_digits
    integer,intent(out),optional::status
    real(dp)::value,inst,pv0,fv0,term,cff,compf,r,df,n,tol
    logical::imm
    integer::digits,nearest
    inst=1.0_dp;pv0=1.0_dp;fv0=0.0_dp;term=0.0_dp;cff=1.0_dp;compf=1.0_dp;imm=.false.;digits=3
    if(present(instalment))inst=instalment;if(present(pv))pv0=pv;if(present(fv))fv0=fv
    if(present(terminal_payment))term=terminal_payment;if(present(cf_freq))cff=cf_freq
    if(present(comp_freq))compf=comp_freq;if(present(immediate_start))imm=immediate_start
    if(present(round_digits))digits=round_digits
    if(abs(rate)<=epsilon(1.0_dp)) then;value=(pv0+fv0-term)/inst;if(present(status))status=JRV_OK;return;end if
    r=equiv_rate(rate,compf,cff)/cff;if(imm)pv0=pv0-inst
    df=(inst/r-pv0)/(fv0-term+inst/r)
    if(df<=0.0_dp.or.1.0_dp+r<=0.0_dp) then;value=nan();if(present(status))status=JRV_INVALID_ARGUMENT;return;end if
    n=-log(df)/log(1.0_dp+r)+merge(1.0_dp,0.0_dp,imm)
    nearest=nint(n);tol=0.5_dp*10.0_dp**(-digits)
    value=merge(real(nearest,dp),n,abs(n-real(nearest,dp))<tol)
    if(present(status))status=JRV_OK
  end function annuity_periods

  function annuity_rate(n_periods, instalment, pv, fv, terminal_payment, &
      immediate_start, cf_freq, comp_freq, status) result(value)
    real(dp),intent(in),optional::n_periods,instalment,pv,fv,terminal_payment,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(out),optional::status
    real(dp)::value,n,cff,compf
    type(annuity_rate_context)::ctx
    type(root_result)::rr
    n=ieee_value(0.0_dp,ieee_positive_inf);ctx%instalment=1.0_dp;ctx%pv=1.0_dp;ctx%fv=0.0_dp
    ctx%terminal_payment=0.0_dp;ctx%immediate_start=.false.;cff=1.0_dp;compf=1.0_dp
    if(present(n_periods))n=n_periods;if(present(instalment))ctx%instalment=instalment
    if(present(pv))ctx%pv=pv;if(present(fv))ctx%fv=fv;if(present(terminal_payment))ctx%terminal_payment=terminal_payment
    if (present(immediate_start)) ctx%immediate_start = immediate_start
    if (present(cf_freq)) cff = cf_freq
    if (present(comp_freq)) compf = comp_freq
    if(.not.ieee_is_finite(n)) then;value=ctx%pv/ctx%instalment;if(present(status))status=JRV_OK;return;end if
    ctx%n_periods=n
    rr=irr_solve(annuity_rate_callback,ctx)
    if(rr%status==JRV_OK)then;value=equiv_rate(rr%root*cff,cff,compf);else;value=nan();end if
    if(present(status))status=rr%status
  end function annuity_rate

  subroutine annuity_rate_callback(r,context,value,gradient)
    real(dp),intent(in)::r
    class(*),intent(in)::context
    real(dp),intent(out)::value,gradient
    real(dp)::df,dfg,af,afg,adjust
    select type(ctx=>context)
    type is(annuity_rate_context)
      df=(1.0_dp+r)**(-ctx%n_periods);dfg=-ctx%n_periods*df/(1.0_dp+r)
      adjust=merge(1.0_dp+r,1.0_dp,ctx%immediate_start)
      if(abs(r)<=sqrt(epsilon(1.0_dp))) then
        af=ctx%n_periods;afg=-ctx%n_periods*(ctx%n_periods+1.0_dp)/2.0_dp
      else
        af=(1.0_dp-df)/r;afg=(-(1.0_dp-df)-r*dfg)/r**2
      end if
      value=adjust*(ctx%instalment*af+(ctx%terminal_payment-ctx%fv)*df)-ctx%pv
      gradient=adjust*(ctx%instalment*afg+(ctx%terminal_payment-ctx%fv)*dfg)
    class default
      value=nan();gradient=nan()
    end select
  end subroutine annuity_rate_callback

  function annuity_instalment_breakup(rate,n_periods,pv,immediate_start,cf_freq,comp_freq,period_no) result(res)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::n_periods,pv,cf_freq,comp_freq
    logical,intent(in),optional::immediate_start
    integer,intent(in),optional::period_no
    type(annuity_breakup_result)::res
    real(dp)::n,pv0,cff,compf,inst,r,growth
    logical::imm
    integer::period
    n=ieee_value(0.0_dp,ieee_positive_inf);pv0=1.0_dp;cff=1.0_dp;compf=1.0_dp;imm=.false.;period=1
    if(present(n_periods))n=n_periods;if(present(pv))pv0=pv;if(present(cf_freq))cff=cf_freq
    if (present(comp_freq)) compf = comp_freq
    if (present(immediate_start)) imm = immediate_start
    if (present(period_no)) period = period_no
    if(period<1)then;res%status=JRV_INVALID_ARGUMENT;return;end if
    inst=annuity_instalment(rate,n,pv0,0.0_dp,0.0_dp,imm,cff,compf)
    r=equiv_rate(rate,compf,cff)/cff;growth=(1.0_dp+r)**real(period-1,dp)
    res%opening_principal=pv0*growth-annuity_fv(rate,real(period-1,dp),inst,0.0_dp,imm,cff,compf)
    res%interest_part=res%opening_principal*r
    res%principal_part=inst-res%interest_part
    res%closing_principal=res%opening_principal+res%interest_part-inst
    res%status=JRV_OK
  end function annuity_instalment_breakup

  pure real(dp) function nan() result(x)
    x=ieee_value(0.0_dp,ieee_quiet_nan)
  end function nan
end module jrvfinance_cashflows
