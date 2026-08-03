! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_performance
  use ptr_kinds, only : dp
  use ptr_types, only : backtest_result, performance_result
  use ptr_utils, only : is_finite, finite_mean, finite_sd, percentile, nan_dp
  implicit none
  private
  public :: calculate_drawdown_series, calculate_annualized_return
  public :: analyze_performance, perf_metrics, benchmark_statistics
  public :: calculate_recovery_time, create_regime_buckets

contains

  subroutine calculate_drawdown_series(values, drawdowns)
    real(dp), intent(in) :: values(:)
    real(dp), allocatable, intent(out) :: drawdowns(:)
    real(dp) :: peak
    integer :: i
    allocate(drawdowns(size(values))); drawdowns = nan_dp()
    peak = -huge(1.0_dp)
    do i = 1, size(values)
      if (.not. is_finite(values(i))) cycle
      peak = max(peak, values(i))
      if (peak > 0.0_dp) drawdowns(i) = values(i) / peak - 1.0_dp
    end do
  end subroutine calculate_drawdown_series

  real(dp) function calculate_annualized_return(values, frequency)
    real(dp), intent(in) :: values(:), frequency
    integer :: first, last, nperiod
    first = 0; last = 0
    do nperiod = 1, size(values)
      if (is_finite(values(nperiod)) .and. values(nperiod) > 0.0_dp) then
        first = nperiod; exit
      end if
    end do
    do nperiod = size(values), 1, -1
      if (is_finite(values(nperiod)) .and. values(nperiod) > 0.0_dp) then
        last = nperiod; exit
      end if
    end do
    if (first == 0 .or. last <= first) then
      calculate_annualized_return = nan_dp()
    else
      nperiod = last - first
      calculate_annualized_return = (values(last)/values(first))**(frequency/real(nperiod,dp)) - 1.0_dp
    end if
  end function calculate_annualized_return

  subroutine perf_metrics(returns, frequency, result, values, turnover, confidence)
    real(dp), intent(in) :: returns(:), frequency
    type(performance_result), intent(out) :: result
    real(dp), intent(in), optional :: values(:), turnover(:), confidence
    real(dp), allocatable :: neg(:), dd(:)
    real(dp) :: mu, sdv, dsd, conf, threshold, gains, losses
    integer :: i, nneg, nfinite
    conf = 0.95_dp; if (present(confidence)) conf = confidence
    mu = finite_mean(returns); sdv = finite_sd(returns)
    result%annualized_return = mu * frequency
    result%annualized_volatility = sdv * sqrt(frequency)
    if (is_finite(sdv) .and. sdv > 0.0_dp) result%sharpe = mu/sdv*sqrt(frequency)
    nneg = count(is_finite(returns) .and. returns < 0.0_dp)
    allocate(neg(max(1,nneg)))
    nneg = 0
    do i=1,size(returns)
      if(is_finite(returns(i)) .and. returns(i)<0.0_dp)then
        nneg=nneg+1;neg(nneg)=returns(i)
      end if
    end do
    if(nneg>1)then
      dsd=sqrt(sum(neg(:nneg)**2)/real(nneg,dp))
      if(dsd>0.0_dp)result%sortino=mu/dsd*sqrt(frequency)
    end if
    result%var = -percentile(returns, 1.0_dp-conf)
    threshold = percentile(returns, 1.0_dp-conf)
    if(is_finite(threshold))then
      nneg=0;losses=0.0_dp
      do i=1,size(returns)
        if(is_finite(returns(i)).and.returns(i)<=threshold)then
          losses=losses+returns(i);nneg=nneg+1
        end if
      end do
      if(nneg>0)result%cvar=-losses/real(nneg,dp)
    end if
    gains=0.0_dp;losses=0.0_dp;nfinite=0
    do i=1,size(returns)
      if(.not.is_finite(returns(i)))cycle
      nfinite=nfinite+1
      if(returns(i)>0.0_dp)then;gains=gains+returns(i);else;losses=losses-returns(i);end if
    end do
    if(losses>0.0_dp)result%omega=gains/losses
    if(nfinite>0)result%win_rate=real(count(is_finite(returns).and.returns>0.0_dp),dp)/real(nfinite,dp)
    if(present(values))then
      call calculate_drawdown_series(values,dd)
      if(any(is_finite(dd)))result%max_drawdown=-minval(dd,mask=is_finite(dd))
      result%total_return=values(size(values))/values(1)-1.0_dp
      result%annualized_return=calculate_annualized_return(values,frequency)
      if(result%max_drawdown>0.0_dp)result%calmar=result%annualized_return/result%max_drawdown
    end if
    if(present(turnover))result%turnover=finite_mean(turnover)
  end subroutine perf_metrics

  subroutine analyze_performance(backtest, frequency, result, confidence)
    type(backtest_result), intent(in) :: backtest
    real(dp), intent(in) :: frequency
    type(performance_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence
    call perf_metrics(backtest%returns, frequency, result, backtest%portfolio_value, &
      backtest%turnover, confidence)
  end subroutine analyze_performance

  subroutine benchmark_statistics(portfolio_returns, benchmark_returns, frequency, &
      beta, alpha, correlation, tracking_error, information_ratio)
    real(dp),intent(in)::portfolio_returns(:),benchmark_returns(:),frequency
    real(dp),intent(out)::beta,alpha,correlation,tracking_error,information_ratio
    real(dp)::mp,mb,cpb,vb,sp,sb,active_mu
    integer::i,n
    mp=0.0_dp;mb=0.0_dp;n=0
    do i=1,min(size(portfolio_returns),size(benchmark_returns))
      if(is_finite(portfolio_returns(i)).and.is_finite(benchmark_returns(i)))then
        mp=mp+portfolio_returns(i);mb=mb+benchmark_returns(i);n=n+1
      end if
    end do
    beta=nan_dp();alpha=nan_dp();correlation=nan_dp();tracking_error=nan_dp();information_ratio=nan_dp()
    if(n<2)return
    mp=mp/real(n,dp);mb=mb/real(n,dp);cpb=0.0_dp;vb=0.0_dp;sp=0.0_dp;sb=0.0_dp
    do i=1,min(size(portfolio_returns),size(benchmark_returns))
      if(is_finite(portfolio_returns(i)).and.is_finite(benchmark_returns(i)))then
        cpb=cpb+(portfolio_returns(i)-mp)*(benchmark_returns(i)-mb)
        vb=vb+(benchmark_returns(i)-mb)**2
        sp=sp+(portfolio_returns(i)-mp)**2;sb=sb+(benchmark_returns(i)-mb)**2
      end if
    end do
    if(vb>0.0_dp)beta=cpb/vb
    if(is_finite(beta))alpha=(mp-beta*mb)*frequency
    if(sp>0.0_dp.and.sb>0.0_dp)correlation=cpb/sqrt(sp*sb)
    tracking_error=finite_sd(portfolio_returns(:min(size(portfolio_returns),size(benchmark_returns)))- &
      benchmark_returns(:min(size(portfolio_returns),size(benchmark_returns))))*sqrt(frequency)
    active_mu=mp-mb
    if(is_finite(tracking_error).and.tracking_error>0.0_dp)information_ratio=active_mu*frequency/tracking_error
  end subroutine benchmark_statistics

  integer function calculate_recovery_time(drawdowns)
    real(dp),intent(in)::drawdowns(:)
    integer::i,trough
    real(dp)::mindd
    calculate_recovery_time=-1;trough=0;mindd=0.0_dp
    do i=1,size(drawdowns)
      if(is_finite(drawdowns(i)).and.drawdowns(i)<mindd)then;mindd=drawdowns(i);trough=i;end if
    end do
    if(trough==0)then;calculate_recovery_time=0;return;end if
    do i=trough+1,size(drawdowns)
      if(is_finite(drawdowns(i)).and.drawdowns(i)>=-1.0e-12_dp)then
        calculate_recovery_time=i-trough;return
      end if
    end do
  end function calculate_recovery_time

  subroutine create_regime_buckets(indicator, breakpoints, buckets, use_percentiles)
    real(dp),intent(in)::indicator(:),breakpoints(:)
    integer,allocatable,intent(out)::buckets(:)
    logical,intent(in),optional::use_percentiles
    logical::pct
    real(dp),allocatable::cuts(:)
    integer::i,j
    pct=.false.;if(present(use_percentiles))pct=use_percentiles
    allocate(buckets(size(indicator)),cuts(size(breakpoints)));buckets=0
    if(pct)then
      do j=1,size(breakpoints);cuts(j)=percentile(indicator,breakpoints(j));end do
    else;cuts=breakpoints;end if
    do i=1,size(indicator)
      if(.not.is_finite(indicator(i)))cycle
      buckets(i)=1
      do j=1,size(cuts)
        if(indicator(i)>cuts(j))buckets(i)=j+1
      end do
    end do
  end subroutine create_regime_buckets

end module ptr_performance
