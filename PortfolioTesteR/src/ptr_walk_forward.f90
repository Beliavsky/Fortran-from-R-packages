! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_walk_forward
  use ptr_kinds, only : dp
  use ptr_types, only : backtest_result, grid_result, walk_forward_result
  use ptr_utils, only : nan_dp
  use ptr_backtest, only : run_backtest
  use ptr_optimization, only : strategy_builder, metric_callback, run_param_grid, metric_sharpe
  implicit none
  private
  public :: make_walk_forward_splits, run_walk_forward, stitch_returns

contains

  subroutine make_walk_forward_splits(n, is_periods, oos_periods, step, &
      is_start, is_end, oos_start, oos_end)
    integer,intent(in)::n,is_periods,oos_periods,step
    integer,allocatable,intent(out)::is_start(:),is_end(:),oos_start(:),oos_end(:)
    integer::nw,s,k
    nw=0;s=1
    do while(s+is_periods+oos_periods-1<=n);nw=nw+1;s=s+step;end do
    allocate(is_start(nw),is_end(nw),oos_start(nw),oos_end(nw));s=1
    do k=1,nw
      is_start(k)=s;is_end(k)=s+is_periods-1
      oos_start(k)=is_end(k)+1;oos_end(k)=oos_start(k)+oos_periods-1
      s=s+step
    end do
  end subroutine make_walk_forward_splits

  subroutine run_walk_forward(prices,params,builder,is_periods,oos_periods,step,result, &
      metric,warmup_periods,initial_capital,cost_bps,frequency)
    real(dp),intent(in)::prices(:,:),params(:,:)
    procedure(strategy_builder)::builder
    integer,intent(in)::is_periods,oos_periods,step
    type(walk_forward_result),intent(out)::result
    procedure(metric_callback),optional::metric
    integer,intent(in),optional::warmup_periods
    real(dp),intent(in),optional::initial_capital,cost_bps,frequency
    type(grid_result)::grid
    type(backtest_result)::bt
    real(dp),allocatable::is_prices(:,:),oos_prices(:,:),weights(:,:),agg(:)
    integer,allocatable::counts(:)
    real(dp)::capital,cost,freq,score
    integer::warm,nw,k,status,ext_start,local_t,global_t
    capital=100000.0_dp;if(present(initial_capital))capital=initial_capital
    cost=0.0_dp;if(present(cost_bps))cost=cost_bps
    freq=52.0_dp;if(present(frequency))freq=frequency
    warm=0;if(present(warmup_periods))warm=max(0,warmup_periods)
    call make_walk_forward_splits(size(prices,1),is_periods,oos_periods,step, &
      result%is_start,result%is_end,result%oos_start,result%oos_end)
    nw=size(result%is_start)
    allocate(result%best_index(nw),result%chosen_params(size(params,1),nw))
    allocate(result%is_score(nw),result%oos_score(nw),result%oos_return(nw))
    result%best_index=0;result%chosen_params=nan_dp();result%is_score=nan_dp()
    result%oos_score=nan_dp();result%oos_return=nan_dp()
    allocate(agg(size(prices,1)),counts(size(prices,1)));agg=0.0_dp;counts=0
    do k=1,nw
      is_prices=prices(result%is_start(k):result%is_end(k),:)
      if(present(metric))then
        call run_param_grid(is_prices,params,builder,grid,metric,capital,cost,freq)
      else
        call run_param_grid(is_prices,params,builder,grid,initial_capital=capital, &
          cost_bps=cost,frequency=freq)
      end if
      if(grid%best_index==0)cycle
      result%best_index(k)=grid%best_index
      result%chosen_params(:,k)=params(:,grid%best_index)
      result%is_score(k)=grid%best_score
      ext_start=max(1,result%oos_start(k)-warm)
      oos_prices=prices(ext_start:result%oos_end(k),:)
      call builder(oos_prices,params(:,grid%best_index),weights,status)
      if(status/=0.or.any(shape(weights)/=shape(oos_prices)))cycle
      if(warm>0)then
        weights=weights(result%oos_start(k)-ext_start+1:,:)
        oos_prices=prices(result%oos_start(k):result%oos_end(k),:)
      end if
      call run_backtest(oos_prices,weights,capital,bt,cost_bps=cost,frequency=freq)
      if(present(metric))then;score=metric(bt);else;score=metric_sharpe(bt);end if
      result%oos_score(k)=score;result%oos_return(k)=bt%total_return
      do local_t=2,size(bt%returns)
        global_t=result%oos_start(k)+local_t-1
        if(counts(global_t)==0)then
          agg(global_t)=bt%returns(local_t)
        else
          agg(global_t)=(1.0_dp+agg(global_t))*(1.0_dp+bt%returns(local_t))-1.0_dp
        end if
        counts(global_t)=counts(global_t)+1
      end do
    end do
    call stitch_returns(agg,counts,capital,result%stitched_index,result%stitched_value)
  end subroutine run_walk_forward

  subroutine stitch_returns(returns,counts,initial_value,indices,values)
    real(dp),intent(in)::returns(:),initial_value
    integer,intent(in)::counts(:)
    integer,allocatable,intent(out)::indices(:)
    real(dp),allocatable,intent(out)::values(:)
    integer::n,i,k
    n=count(counts>0)
    if(n==0)then;allocate(indices(0),values(0));return;end if
    allocate(indices(n+1),values(n+1));k=1
    do i=1,size(counts);if(counts(i)>0)then;indices(1)=max(1,i-1);exit;end if;end do
    values(1)=initial_value
    do i=1,size(counts)
      if(counts(i)>0)then;k=k+1;indices(k)=i;values(k)=values(k-1)*(1.0_dp+returns(i));end if
    end do
  end subroutine stitch_returns

end module ptr_walk_forward
