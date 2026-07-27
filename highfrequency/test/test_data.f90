! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program test_data
  use highfrequency, only: dp, make_returns, aggregate_last, previous_tick
  use highfrequency, only: refresh_time_pair, trade_direction, liquidity_measures
  use highfrequency, only: liquidity_result, make_ohlcv
  use highfrequency, only: no_zero_prices_mask, nonnegative_spread_mask
  implicit none
  real(dp) :: prices(4), expected(3), values(5,1), target(4)
  real(dp) :: bid(5), offer(5), trade_price(5), volume(5), bid_size(5), offer_size(5)
  real(dp), allocatable :: returns(:), out_values(:,:), x1(:), x2(:), ohlcv(:,:)
  integer :: times(5), source_times(3), target_times(4), t1(3), t2(3)
  integer, allocatable :: out_times(:), refresh_times(:), directions(:)
  integer :: nout
  logical :: available(4), keep5(5)
  type(liquidity_result) :: liq

  prices=[100.0_dp,101.0_dp,100.5_dp,102.0_dp]
  expected=log(prices(2:)/prices(:3))
  returns=make_returns(prices)
  call assert_array(returns,expected,1.0e-14_dp)

  times=[1,2,7,8,11]
  values(:,1)=[10.0_dp,11.0_dp,12.0_dp,13.0_dp,14.0_dp]
  call aggregate_last(times,values,5,out_times,out_values,nout)
  if(nout/=3)error stop 1
  call assert_array(out_values(:,1),[11.0_dp,13.0_dp,14.0_dp],0.0_dp)

  source_times=[1,4,9]
  target_times=[0,2,5,10]
  call previous_tick(source_times,[2.0_dp,3.0_dp,4.0_dp],target_times,target,available)
  call assert_array(target,[0.0_dp,2.0_dp,3.0_dp,4.0_dp],0.0_dp)
  if(any(available .neqv. [.false.,.true.,.true.,.true.]))error stop 1

  t1=[1,4,8];t2=[2,5,7]
  call refresh_time_pair(t1,[100.0_dp,101.0_dp,102.0_dp],t2,[50.0_dp,51.0_dp,52.0_dp], &
    refresh_times,x1,x2,nout)
  if(nout<2)error stop 1

  trade_price=[100.0_dp,100.2_dp,100.2_dp,99.8_dp,100.1_dp]
  bid=[99.9_dp,100.0_dp,100.1_dp,99.7_dp,100.0_dp]
  offer=[100.1_dp,100.2_dp,100.3_dp,99.9_dp,100.2_dp]
  allocate(directions(5))
  call trade_direction(trade_price,bid,offer,directions)
  if(any(abs(directions)/=1))error stop 1
  volume=10.0_dp
  bid_size=[100.0_dp,110.0_dp,120.0_dp,130.0_dp,140.0_dp]
  offer_size=[90.0_dp,100.0_dp,110.0_dp,120.0_dp,130.0_dp]
  call liquidity_measures(trade_price,volume,bid,offer,bid_size,offer_size,1,liq)
  if(size(liq%direction)/=5)error stop 1
  call assert_close(liq%quoted_spread(1),0.2_dp,1.0e-12_dp)

  keep5=no_zero_prices_mask(trade_price)
  if(.not.all(keep5))error stop 1
  keep5=nonnegative_spread_mask(bid,offer)
  if(.not.all(keep5))error stop 1

  call make_ohlcv(times,trade_price,volume,5,out_times,ohlcv,nout)
  if(nout/=3 .or. size(ohlcv,2)/=5)error stop 1
  print '(a)', 'test_data: PASS'
contains
  subroutine assert_close(actual,expected,tolerance)
    real(dp),intent(in)::actual,expected,tolerance
    if(abs(actual-expected)>tolerance)error stop 1
  end subroutine assert_close
  subroutine assert_array(actual,expected,tolerance)
    real(dp),intent(in)::actual(:),expected(:),tolerance
    if(size(actual)/=size(expected))error stop 1
    if(maxval(abs(actual-expected))>tolerance)error stop 1
  end subroutine assert_array
end program test_data
