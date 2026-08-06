! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_backtest_module
  use ghyp_kinds, only : dp
  use tsd_distributions, only : qdist
  use tsd_math, only : regularized_gamma_p
  use tsgarch_types
  use tsgarch_model, only : filter_garch
  use tsgarch_fit_module, only : estimate_garch
  implicit none
  private
  public :: backtest_var, coverage_tests
contains
  function backtest_var(y,spec,start_index,probability,refit_every,window,vreg,options) result(out)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    integer,intent(in)::start_index
    real(dp),intent(in),optional::probability
    integer,intent(in),optional::refit_every,window
    real(dp),intent(in),optional::vreg(:,:)
    type(fit_options),intent(in),optional::options
    type(backtest_result)::out
    type(garch_fit)::fit
    type(garch_filter_result)::filtered
    type(fit_options)::opt
    real(dp)::pr,z
    integer::step,win,nout,k,t,lo,last_fit
    if(start_index<max(21,max(spec%p,spec%q)+3).or.start_index>=size(y))then
    out%message='invalid backtest start'
    return
    end if
    pr=0.01_dp
    if(present(probability))pr=probability
    step=1
    if(present(refit_every))step=max(1,refit_every)
    win=0
    if(present(window))win=max(0,window)
    if(pr<=0.0_dp.or.pr>=1.0_dp)then
    out%message='probability must lie in (0,1)'
    return
    end if
    if(present(vreg))then
    if(size(vreg,1)/=size(y))then
    out%message='variance regressor row count mismatch'
    return
    end if
    end if
    opt=fit_options()
    opt%compute_inference=.false.
    if(present(options))opt=options
    opt%compute_inference=.false.
    nout=size(y)-start_index
    allocate(out%actual(nout),out%mean(nout),out%sigma(nout),out%value_at_risk(nout),out%exceedance(nout))
    out%actual=y(start_index+1:size(y))
    out%mean=0.0_dp
    out%sigma=0.0_dp
    out%value_at_risk=0.0_dp
    out%exceedance=.false.
    last_fit=-huge(1)
    do k=1,nout
      t=start_index+k
      if(k==1.or.k-last_fit>=step)then
        lo=1
        if(win>0)lo=max(1,t-win)
        if(present(vreg))then
        fit=estimate_garch(y(lo:t-1),spec,vreg=vreg(lo:t-1,:),options=opt)
        else
        fit=estimate_garch(y(lo:t-1),spec,options=opt)
        end if
        if(fit%filtered%status/=tsg_success)then
        out%message='backtest estimation failed: '//trim(fit%message)
        return
        end if
        last_fit=k
        out%refits=out%refits+1
      end if
      lo=1
      if(win>0)lo=max(1,t-win)
      if(present(vreg))then
        filtered=filter_garch([y(lo:t-1),fit%parameters%mu],fit%spec,fit%parameters,vreg(lo:t,:))
      else
        filtered=filter_garch([y(lo:t-1),fit%parameters%mu],fit%spec,fit%parameters)
      end if
      if(filtered%status/=tsg_success)then
      out%message='backtest filtering failed'
      return
      end if
      z=qdist(fit%spec%distribution,pr,fit%parameters%dist)
      out%mean(k)=fit%parameters%mu
      out%sigma(k)=filtered%sigma(size(filtered%sigma))
      out%value_at_risk(k)=out%mean(k)+out%sigma(k)*z
      out%exceedance(k)=out%actual(k)<out%value_at_risk(k)
    end do
    out%expected_coverage=pr
    out%coverage=real(count(out%exceedance),dp)/real(nout,dp)
    call coverage_tests(out%exceedance,pr,out%kupiec_statistic,out%kupiec_pvalue,out%independence_statistic,&
      out%independence_pvalue,out%conditional_coverage_statistic,out%conditional_coverage_pvalue)
    out%status=tsg_success
    out%message='ok'
  end function backtest_var

  subroutine coverage_tests(exceedance,probability,kupiec,kupiec_p,independence,independence_p,conditional,conditional_p)
    logical,intent(in)::exceedance(:)
    real(dp),intent(in)::probability
    real(dp),intent(out)::kupiec,kupiec_p,independence,independence_p,conditional,conditional_p
    integer::n,x,i,n00,n01,n10,n11
    real(dp)::phat,p01,p11,pall,ll0,ll1,tiny_p
    n=size(exceedance)
    x=count(exceedance)
    tiny_p=1.0e-14_dp
    phat=min(max(real(x,dp)/real(max(1,n),dp),tiny_p),1.0_dp-tiny_p)
    kupiec=-2.0_dp*((real(n-x,dp)*log(max(1.0_dp-probability,tiny_p))+real(x,dp)*log(max(probability,tiny_p)))-&
      (real(n-x,dp)*log(1.0_dp-phat)+real(x,dp)*log(phat)))
    kupiec=max(kupiec,0.0_dp)
    kupiec_p=1.0_dp-regularized_gamma_p(0.5_dp,0.5_dp*kupiec)
    n00=0
    n01=0
    n10=0
    n11=0
    do i=2,n
      if(.not.exceedance(i-1).and..not.exceedance(i))n00=n00+1
      if(.not.exceedance(i-1).and.exceedance(i))n01=n01+1
      if(exceedance(i-1).and..not.exceedance(i))n10=n10+1
      if(exceedance(i-1).and.exceedance(i))n11=n11+1
    end do
    p01=real(n01,dp)/real(max(1,n00+n01),dp)
    p11=real(n11,dp)/real(max(1,n10+n11),dp)
    pall=real(n01+n11,dp)/real(max(1,n00+n01+n10+n11),dp)
    p01=min(max(p01,tiny_p),1.0_dp-tiny_p)
    p11=min(max(p11,tiny_p),1.0_dp-tiny_p)
    pall=min(max(pall,tiny_p),1.0_dp-tiny_p)
    ll0=real(n00+n10,dp)*log(1.0_dp-pall)+real(n01+n11,dp)*log(pall)
    ll1=real(n00,dp)*log(1.0_dp-p01)+real(n01,dp)*log(p01)+real(n10,dp)*log(1.0_dp-p11)+real(n11,dp)*log(p11)
    independence=max(0.0_dp,-2.0_dp*(ll0-ll1))
    independence_p=1.0_dp-regularized_gamma_p(0.5_dp,0.5_dp*independence)
    conditional=kupiec+independence
    conditional_p=exp(-0.5_dp*conditional)
  end subroutine coverage_tests
end module tsgarch_backtest_module
