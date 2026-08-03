! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_backtest
  use ptr_kinds, only : dp
  use ptr_types, only : backtest_result
  use ptr_utils, only : is_finite, normalize_nonnegative, finite_mean, finite_sd
  use ptr_performance, only : calculate_drawdown_series, calculate_annualized_return
  implicit none
  private
  public :: run_backtest, portfolio_returns_from_weights

contains

  subroutine run_backtest(prices, target_weights, initial_capital, result, &
      cost_bps, integer_shares, stop_loss, frequency)
    real(dp), intent(in) :: prices(:,:), target_weights(:,:), initial_capital
    type(backtest_result), intent(out) :: result
    real(dp), intent(in), optional :: cost_bps, stop_loss, frequency
    logical, intent(in), optional :: integer_shares
    real(dp), allocatable :: shares(:), desired(:), roww(:), entry(:), &
      position_values(:), previous_weights(:), drawdowns(:)
    real(dp) :: cash, value, fee_rate, stop, freq, trade, proceeds, cost
    real(dp) :: max_buy, sdv, mu
    integer :: nt, na, t, j
    logical :: whole

    nt = size(prices,1); na = size(prices,2)
    allocate(result%portfolio_value(nt), result%returns(nt), result%cash(nt))
    allocate(result%shares(nt,na), result%executed_weights(nt,na), result%turnover(nt))
    result%portfolio_value = initial_capital
    result%returns = 0.0_dp; result%cash = initial_capital
    result%shares = 0.0_dp; result%executed_weights = 0.0_dp; result%turnover = 0.0_dp
    result%initial_capital = initial_capital
    result%n_transactions = 0; result%first_trade = 0; result%bankrupt = .false.
    if (any(shape(prices) /= shape(target_weights)) .or. nt == 0 .or. na == 0) then
      result%bankrupt = .true.; return
    end if

    fee_rate = 0.0_dp; if (present(cost_bps)) fee_rate = max(0.0_dp,cost_bps)*1.0e-4_dp
    whole = .true.; if (present(integer_shares)) whole = integer_shares
    stop = 0.0_dp; if (present(stop_loss)) stop = max(0.0_dp,min(0.999_dp,stop_loss))
    freq = 52.0_dp; if (present(frequency)) freq = frequency
    allocate(shares(na), desired(na), roww(na), entry(na), position_values(na), previous_weights(na))
    shares = 0.0_dp; entry = 0.0_dp; previous_weights = 0.0_dp; cash = initial_capital

    do t = 1, nt
      do j = 1, na
        if (shares(j) > 0.0_dp .and. stop > 0.0_dp .and. entry(j) > 0.0_dp .and. &
            is_finite(prices(t,j)) .and. prices(t,j) <= entry(j)*(1.0_dp-stop)) then
          proceeds = shares(j)*prices(t,j)
          cash = cash + proceeds*(1.0_dp-fee_rate)
          shares(j) = 0.0_dp; entry(j) = 0.0_dp
          result%n_transactions = result%n_transactions + 1
        end if
      end do

      position_values = 0.0_dp
      do j=1,na
        if(is_finite(prices(t,j)).and.prices(t,j)>0.0_dp) position_values(j)=shares(j)*prices(t,j)
      end do
      value = cash + sum(position_values)
      if (.not. is_finite(value) .or. value <= 0.0_dp) then
        result%bankrupt = .true.
        if(t>1)then
          result%portfolio_value(t:)=result%portfolio_value(t-1)
          result%cash(t:)=cash
          result%shares(t:,:)=spread(shares,1,nt-t+1)
        end if
        exit
      end if

      roww = target_weights(t,:)
      where(.not.is_finite(roww).or.roww<0.0_dp)roww=0.0_dp
      if(sum(roww)>1.0_dp)call normalize_nonnegative(roww)
      desired = shares
      do j=1,na
        if(is_finite(prices(t,j)).and.prices(t,j)>0.0_dp)then
          desired(j)=value*roww(j)/prices(t,j)
          if(whole)desired(j)=floor(desired(j))
        else
          desired(j)=shares(j)
        end if
      end do

      do j=1,na
        trade=desired(j)-shares(j)
        if(trade<0.0_dp.and.is_finite(prices(t,j)).and.prices(t,j)>0.0_dp)then
          proceeds=(-trade)*prices(t,j)
          cash=cash+proceeds*(1.0_dp-fee_rate)
          shares(j)=desired(j)
          result%n_transactions=result%n_transactions+1
          if(shares(j)<=0.0_dp)entry(j)=0.0_dp
        end if
      end do

      do j=1,na
        trade=desired(j)-shares(j)
        if(trade>0.0_dp.and.is_finite(prices(t,j)).and.prices(t,j)>0.0_dp)then
          cost=trade*prices(t,j)*(1.0_dp+fee_rate)
          if(cost>cash)then
            max_buy=cash/(prices(t,j)*(1.0_dp+fee_rate))
            if(whole)max_buy=floor(max_buy)
            trade=max(0.0_dp,min(trade,max_buy))
            cost=trade*prices(t,j)*(1.0_dp+fee_rate)
          end if
          if(trade>0.0_dp)then
            cash=cash-cost;shares(j)=shares(j)+trade;entry(j)=prices(t,j)
            result%n_transactions=result%n_transactions+1
            if(result%first_trade==0)result%first_trade=t
          end if
        end if
      end do

      position_values=0.0_dp
      do j=1,na
        if(is_finite(prices(t,j)).and.prices(t,j)>0.0_dp) position_values(j)=shares(j)*prices(t,j)
      end do
      value=cash+sum(position_values)
      result%portfolio_value(t)=value;result%cash(t)=cash;result%shares(t,:)=shares
      if(value>0.0_dp)result%executed_weights(t,:)=position_values/value
      result%turnover(t)=0.5_dp*sum(abs(result%executed_weights(t,:)-previous_weights))
      previous_weights=result%executed_weights(t,:)
      if (t > 1) then
        if (result%portfolio_value(t-1) > 0.0_dp) then
          result%returns(t)=result%portfolio_value(t)/result%portfolio_value(t-1)-1.0_dp
        end if
      end if
    end do

    result%total_return=result%portfolio_value(nt)/initial_capital-1.0_dp
    result%annualized_return=calculate_annualized_return(result%portfolio_value,freq)
    mu=finite_mean(result%returns(2:));sdv=finite_sd(result%returns(2:))
    if(is_finite(sdv).and.sdv>0.0_dp)result%sharpe=mu/sdv*sqrt(freq)
    if(is_finite(sdv))result%annualized_volatility=sdv*sqrt(freq)
    call calculate_drawdown_series(result%portfolio_value,drawdowns)
    if(any(is_finite(drawdowns)))result%max_drawdown=-minval(drawdowns,mask=is_finite(drawdowns))
  end subroutine run_backtest

  subroutine portfolio_returns_from_weights(prices,weights,returns,cost_bps)
    real(dp),intent(in)::prices(:,:),weights(:,:)
    real(dp),allocatable,intent(out)::returns(:)
    real(dp),intent(in),optional::cost_bps
    real(dp)::fee,turn
    integer::t,j
    fee=0.0_dp;if(present(cost_bps))fee=max(0.0_dp,cost_bps)*1.0e-4_dp
    allocate(returns(size(prices,1)));returns=0.0_dp
    do t=2,size(prices,1)
      do j=1,size(prices,2)
        if(is_finite(prices(t,j)).and.is_finite(prices(t-1,j)).and.prices(t-1,j)>0.0_dp)then
          returns(t)=returns(t)+weights(t-1,j)*(prices(t,j)/prices(t-1,j)-1.0_dp)
        end if
      end do
      turn=0.5_dp*sum(abs(weights(t,:)-weights(t-1,:)))
      returns(t)=returns(t)-fee*turn
    end do
  end subroutine portfolio_returns_from_weights

end module ptr_backtest
