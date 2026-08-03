! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_strategies
  use ptr_kinds, only : dp
  use ptr_indicators, only : calc_momentum, calc_rsi, calc_rolling_volatility
  use ptr_filters, only : filter_top_n
  use ptr_weighting, only : weight_equally, weight_by_rank, weight_by_signal
  use ptr_utils, only : nan_dp, is_finite
  implicit none
  private
  public :: momentum_top_n_strategy, rsi_reversion_strategy
  public :: volatility_adjusted_momentum_strategy

contains

  subroutine momentum_top_n_strategy(prices,params,weights,status)
    real(dp),intent(in)::prices(:,:),params(:)
    real(dp),allocatable,intent(out)::weights(:,:)
    integer,intent(out)::status
    real(dp),allocatable::momentum(:,:),selection(:,:)
    integer::lookback,n_top,method
    status=1
    if(size(params)<2)return
    lookback=max(1,nint(params(1)));n_top=max(1,nint(params(2)));method=0
    if(size(params)>=3)method=nint(params(3))
    call calc_momentum(prices,lookback,momentum)
    call filter_top_n(momentum,n_top,selection)
    select case(method)
    case(1);call weight_by_rank(selection,momentum,weights)
    case(2);call weight_by_signal(selection,momentum,weights)
    case default;call weight_equally(selection,weights)
    end select
    status=0
  end subroutine momentum_top_n_strategy

  subroutine rsi_reversion_strategy(prices,params,weights,status)
    real(dp),intent(in)::prices(:,:),params(:)
    real(dp),allocatable,intent(out)::weights(:,:)
    integer,intent(out)::status
    real(dp),allocatable::rsi(:,:),selection(:,:)
    integer::period,n_top
    status=1
    if(size(params)<2)return
    period=max(2,nint(params(1)));n_top=max(1,nint(params(2)))
    call calc_rsi(prices,period,rsi)
    call filter_top_n(rsi,n_top,selection,ascending=.true.)
    call weight_equally(selection,weights)
    status=0
  end subroutine rsi_reversion_strategy

  subroutine volatility_adjusted_momentum_strategy(prices,params,weights,status)
    real(dp),intent(in)::prices(:,:),params(:)
    real(dp),allocatable,intent(out)::weights(:,:)
    integer,intent(out)::status
    real(dp),allocatable::momentum(:,:),volatility(:,:),score(:,:),selection(:,:)
    integer::mom_lookback,vol_lookback,n_top,t,j
    status=1
    if(size(params)<3)return
    mom_lookback=max(1,nint(params(1)));vol_lookback=max(2,nint(params(2)))
    n_top=max(1,nint(params(3)))
    call calc_momentum(prices,mom_lookback,momentum)
    call calc_rolling_volatility(prices,vol_lookback,volatility)
    allocate(score(size(prices,1),size(prices,2)));score=nan_dp()
    do t=1,size(score,1)
      do j=1,size(score,2)
        if(is_finite(momentum(t,j)).and.is_finite(volatility(t,j)).and.volatility(t,j)>0.0_dp)then
          score(t,j)=momentum(t,j)/volatility(t,j)
        end if
      end do
    end do
    call filter_top_n(score,n_top,selection)
    call weight_by_rank(selection,score,weights)
    status=0
  end subroutine volatility_adjusted_momentum_strategy

end module ptr_strategies
