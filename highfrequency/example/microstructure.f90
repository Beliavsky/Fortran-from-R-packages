! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program microstructure
  use highfrequency, only: dp, lead_lag, lead_lag_result
  use highfrequency, only: liquidity_measures, liquidity_result
  implicit none
  integer :: times(6), lags(5), i
  real(dp) :: p1(6),p2(6),bid(6),offer(6),volume(6),bid_size(6),offer_size(6)
  type(lead_lag_result) :: ll
  type(liquidity_result) :: liq

  times=[0,1,2,3,4,5]
  do i=1,6
    p1(i)=100.0_dp*exp(0.001_dp*real(i-1,dp))
  end do
  p2(1)=p1(1)
  p2(2:)=p1(:5)
  lags=[-2,-1,0,1,2]
  ll=lead_lag(times,p1,times,p2,lags)
  print '(a,i0)', 'estimated lag: ',ll%optimal_lag

  bid=p1-0.01_dp
  offer=p1+0.01_dp
  volume=100.0_dp
  bid_size=500.0_dp
  offer_size=450.0_dp
  call liquidity_measures(p1,volume,bid,offer,bid_size,offer_size,1,liq)
  print '(a,f10.6)', 'mean effective spread: ',sum(liq%effective_spread)/real(size(p1),dp)
end program microstructure
