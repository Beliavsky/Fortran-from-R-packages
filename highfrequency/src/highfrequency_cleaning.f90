! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_cleaning
  use highfrequency_kinds, only: dp
  use highfrequency_data, only: previous_tick
  implicit none
  private
  public :: no_zero_prices_mask, no_zero_quotes_mask
  public :: nonnegative_spread_mask, maximum_spread_mask
  public :: price_outlier_mask, match_trades_quotes, spread_prices
  public :: business_time_groups

contains

  pure function no_zero_prices_mask(price) result(keep)
    real(dp),intent(in)::price(:)
    logical::keep(size(price))
    keep=price>0.0_dp
  end function no_zero_prices_mask

  pure function no_zero_quotes_mask(bid,offer,bid_size,offer_size) result(keep)
    real(dp),intent(in)::bid(:),offer(:),bid_size(:),offer_size(:)
    logical::keep(size(bid))
    keep=bid>0.0_dp .and. offer>0.0_dp .and. bid_size>0.0_dp .and. offer_size>0.0_dp
  end function no_zero_quotes_mask

  pure function nonnegative_spread_mask(bid,offer) result(keep)
    real(dp),intent(in)::bid(:),offer(:)
    logical::keep(size(bid))
    keep=offer>=bid
  end function nonnegative_spread_mask

  pure function maximum_spread_mask(bid,offer,max_relative_spread) result(keep)
    real(dp),intent(in)::bid(:),offer(:),max_relative_spread
    logical::keep(size(bid))
    real(dp)::mid(size(bid))
    mid=0.5_dp*(bid+offer)
    keep=mid>0.0_dp
    where(keep)keep=(offer-bid)/mid<=max_relative_spread
  end function maximum_spread_mask

  function price_outlier_mask(price,window,scale) result(keep)
    real(dp),intent(in)::price(:),scale
    integer,intent(in)::window
    logical::keep(size(price))
    real(dp),allocatable::log_price(:),change(:),local(:)
    real(dp)::med,mad
    integer::i,left,right,n,m
    n=size(price)
    keep=.true.
    if(n<3 .or. window<1 .or. any(price<=0.0_dp))return
    allocate(log_price(n),change(n),local(2*window+1))
    log_price=log(price)
    change=0.0_dp
    change(2:)=abs(log_price(2:)-log_price(:n-1))
    do i=2,n
      left=max(2,i-window)
      right=min(n,i+window)
      m=right-left+1
      local(:m)=change(left:right)
      call sort_values(local(:m))
      med=median_sorted(local(:m))
      local(:m)=abs(local(:m)-med)
      call sort_values(local(:m))
      mad=median_sorted(local(:m))
      if(mad>0.0_dp)keep(i)=change(i)<=med+scale*1.4826_dp*mad
    end do
  contains
    subroutine sort_values(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::ii,jj
      do ii=2,size(x)
        key=x(ii);jj=ii-1
        do while(jj>=1)
          if(x(jj)<=key)exit
          x(jj+1)=x(jj);jj=jj-1
        end do
        x(jj+1)=key
      end do
    end subroutine sort_values
    pure real(dp) function median_sorted(x) result(value)
      real(dp),intent(in)::x(:)
      integer::nn
      nn=size(x)
      if(mod(nn,2)==1)then
        value=x((nn+1)/2)
      else
        value=0.5_dp*(x(nn/2)+x(nn/2+1))
      end if
    end function median_sorted
  end function price_outlier_mask

  subroutine match_trades_quotes(trade_times,quote_times,bid,offer,matched_bid,matched_offer,available)
    integer,intent(in)::trade_times(:),quote_times(:)
    real(dp),intent(in)::bid(:),offer(:)
    real(dp),intent(out)::matched_bid(size(trade_times)),matched_offer(size(trade_times))
    logical,intent(out),optional::available(size(trade_times))
    logical::a1(size(trade_times)),a2(size(trade_times))
    call previous_tick(quote_times,bid,trade_times,matched_bid,a1)
    call previous_tick(quote_times,offer,trade_times,matched_offer,a2)
    if(present(available))available=a1.and.a2
  end subroutine match_trades_quotes

  pure function spread_prices(bid,offer) result(midpoint)
    real(dp),intent(in)::bid(:),offer(:)
    real(dp)::midpoint(size(bid))
    midpoint=0.5_dp*(bid+offer)
  end function spread_prices

  subroutine business_time_groups(activity,target,group,n_groups)
    real(dp),intent(in)::activity(:),target
    integer,intent(out)::group(size(activity)),n_groups
    real(dp)::cumulative
    integer::i
    group=0
    n_groups=0
    cumulative=0.0_dp
    if(target<=0.0_dp)return
    do i=1,size(activity)
      if(n_groups==0) n_groups=1
      group(i)=n_groups
      cumulative=cumulative+max(0.0_dp,activity(i))
      if(cumulative>=target .and. i<size(activity))then
        n_groups=n_groups+1
        cumulative=0.0_dp
      end if
    end do
  end subroutine business_time_groups

end module highfrequency_cleaning
