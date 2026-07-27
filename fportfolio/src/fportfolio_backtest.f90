! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_backtest
  use fportfolio_kinds, only: dp
  use fportfolio_types, only: linear_constraints, optimizer_result, backtest_result
  use fportfolio_statistics, only: sample_estimator
  use fportfolio_optimization, only: minvariance_portfolio, tangency_portfolio, risk_parity_portfolio, &
    maximum_diversification_portfolio, minimum_cvar_portfolio
  use fportfolio_risk, only: maximum_drawdown, historical_var, historical_es, drawdown_at_risk, conditional_drawdown_at_risk
  implicit none
  private
  public :: run_backtest, rolling_sigma, rolling_var, rolling_cvar, rolling_dar, rolling_cdar, &
            exponential_weight_smoother, net_performance
contains
  subroutine run_backtest(data,window,rebalance_every,strategy,con,result,risk_free,alpha, &
      smoothing_lambda,transaction_cost,annualization)
    real(dp),intent(in)::data(:,:)
    integer,intent(in)::window,rebalance_every
    character(len=*),intent(in)::strategy
    type(linear_constraints),intent(in)::con
    type(backtest_result),intent(out)::result
    real(dp),intent(in),optional::risk_free,alpha,smoothing_lambda,transaction_cost,annualization
    real(dp),allocatable::mu(:),sigma(:,:),w(:),wnew(:)
    type(optimizer_result)::opt
    real(dp)::rf,a,lam,tc,ann,cost,meanr,sdr
    integer::n,p,t,k,status
    n=size(data,1);p=size(data,2)
    if(window<2 .or. window>=n)error stop "run_backtest: invalid window"
    rf=0.0_dp;if(present(risk_free))rf=risk_free
    a=0.05_dp;if(present(alpha))a=alpha
    lam=0.0_dp;if(present(smoothing_lambda))lam=max(0.0_dp,min(1.0_dp,smoothing_lambda))
    tc=0.0_dp;if(present(transaction_cost))tc=max(0.0_dp,transaction_cost)
    ann=252.0_dp;if(present(annualization))ann=annualization
    allocate(result%portfolio_returns(n-window),result%wealth(n-window),result%turnover(n-window))
    allocate(result%transaction_costs(n-window),result%weights(p,n-window),w(p),wnew(p))
    w=1.0_dp/real(p,dp);result%wealth=1.0_dp
    do t=window+1,n
      k=t-window
      if(k==1 .or. mod(k-1,max(1,rebalance_every))==0)then
        call sample_estimator(data(t-window:t-1,:),mu,sigma)
        select case(trim(adjustl(strategy)))
        case("minvariance","minimum_variance")
          call minvariance_portfolio(mu,sigma,con,opt)
        case("tangency","maxratio")
          call tangency_portfolio(mu,sigma,rf,con,opt)
        case("riskparity","risk_parity")
          call risk_parity_portfolio(mu,sigma,con,opt)
        case("maxdiversification","maximum_diversification")
          call maximum_diversification_portfolio(mu,sigma,con,opt)
        case("cvar","minimum_cvar")
          call minimum_cvar_portfolio(data(t-window:t-1,:),mu,a,con,opt,max_iter=30000)
        case default
          opt%converged=.true.;allocate(opt%weights(p));opt%weights=1.0_dp/real(p,dp)
        end select
        if(opt%converged)then
          wnew=lam*w+(1.0_dp-lam)*opt%weights
        else
          wnew=w
        end if
        result%turnover(k)=sum(abs(wnew-w))
        cost=tc*result%turnover(k)
        result%transaction_costs(k)=cost
        w=wnew
      else
        result%turnover(k)=0.0_dp;result%transaction_costs(k)=0.0_dp;cost=0.0_dp
      end if
      result%weights(:,k)=w
      result%portfolio_returns(k)=dot_product(data(t,:),w)-cost
      if(k==1)then
        result%wealth(k)=1.0_dp+result%portfolio_returns(k)
      else
        result%wealth(k)=result%wealth(k-1)*(1.0_dp+result%portfolio_returns(k))
      end if
    end do
    result%total_return=result%wealth(n-window)-1.0_dp
    meanr=sum(result%portfolio_returns)/real(n-window,dp)
    sdr=sqrt(sum((result%portfolio_returns-meanr)**2)/real(max(1,n-window-1),dp))
    result%annualized_return=(1.0_dp+result%total_return)**(ann/real(n-window,dp))-1.0_dp
    result%annualized_volatility=sdr*sqrt(ann)
    if(sdr>0.0_dp)result%sharpe=meanr/sdr*sqrt(ann)
    result%max_drawdown=maximum_drawdown(result%portfolio_returns)
    status=0
  end subroutine run_backtest

  function exponential_weight_smoother(weights,lambda,double_smoothing,initial_weights) result(smoothed)
    real(dp),intent(in)::weights(:,:),lambda
    logical,intent(in),optional::double_smoothing
    real(dp),intent(in),optional::initial_weights(:)
    real(dp)::smoothed(size(weights,1),size(weights,2)),level(size(weights,1)),trend(size(weights,1))
    logical::dbl
    integer::t
    dbl=.false.;if(present(double_smoothing))dbl=double_smoothing
    if(present(initial_weights))then;level=initial_weights;else;level=weights(:,1);end if
    trend=0.0_dp
    do t=1,size(weights,2)
      if(dbl)then
        trend=lambda*trend+(1.0_dp-lambda)*(weights(:,t)-level)
        level=lambda*level+(1.0_dp-lambda)*weights(:,t)
        smoothed(:,t)=level+trend
      else
        level=lambda*level+(1.0_dp-lambda)*weights(:,t)
        smoothed(:,t)=level
      end if
      if(abs(sum(smoothed(:,t)))>tiny(1.0_dp)) smoothed(:,t)=smoothed(:,t)/sum(smoothed(:,t))
    end do
  end function exponential_weight_smoother

  subroutine net_performance(gross_returns,turnover,proportional_cost,fixed_cost,net_returns,wealth)
    real(dp),intent(in)::gross_returns(:),turnover(:),proportional_cost,fixed_cost
    real(dp),allocatable,intent(out)::net_returns(:),wealth(:)
    integer::i,n
    n=size(gross_returns);allocate(net_returns(n),wealth(n))
    net_returns=gross_returns-proportional_cost*turnover
    where(turnover>0.0_dp)net_returns=net_returns-fixed_cost
    do i=1,n
      if(i==1)then;wealth(i)=1.0_dp+net_returns(i);else;wealth(i)=wealth(i-1)*(1.0_dp+net_returns(i));end if
    end do
  end subroutine net_performance

  function rolling_sigma(returns,width) result(out)
    real(dp),intent(in)::returns(:)
    integer,intent(in)::width
    real(dp)::out(size(returns)),m
    integer::i
    out=0.0_dp
    do i=width,size(returns)
      m=sum(returns(i-width+1:i))/real(width,dp)
      out(i)=sqrt(sum((returns(i-width+1:i)-m)**2)/real(max(1,width-1),dp))
    end do
  end function rolling_sigma

  function rolling_var(returns,width,alpha) result(out)
    real(dp),intent(in)::returns(:),alpha
    integer,intent(in)::width
    real(dp)::out(size(returns)),x(width,1),w(1)
    integer::i
    out=0.0_dp;w=1.0_dp
    do i=width,size(returns);x(:,1)=returns(i-width+1:i);out(i)=historical_var(x,w,alpha);end do
  end function rolling_var

  function rolling_cvar(returns,width,alpha) result(out)
    real(dp),intent(in)::returns(:),alpha
    integer,intent(in)::width
    real(dp)::out(size(returns)),x(width,1),w(1)
    integer::i
    out=0.0_dp;w=1.0_dp
    do i=width,size(returns);x(:,1)=returns(i-width+1:i);out(i)=historical_es(x,w,alpha);end do
  end function rolling_cvar

  function rolling_dar(returns,width,alpha) result(out)
    real(dp),intent(in)::returns(:),alpha
    integer,intent(in)::width
    real(dp)::out(size(returns))
    integer::i
    out=0.0_dp
    do i=width,size(returns);out(i)=drawdown_at_risk(returns(i-width+1:i),alpha);end do
  end function rolling_dar

  function rolling_cdar(returns,width,alpha) result(out)
    real(dp),intent(in)::returns(:),alpha
    integer,intent(in)::width
    real(dp)::out(size(returns))
    integer::i
    out=0.0_dp
    do i=width,size(returns);out(i)=conditional_drawdown_at_risk(returns(i-width+1:i),alpha);end do
  end function rolling_cdar
end module fportfolio_backtest
