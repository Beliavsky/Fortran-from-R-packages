! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program test_microstructure
  use highfrequency, only: dp, lead_lag, lead_lag_result
  use highfrequency, only: spot_volatility, spot_drift, drift_burst_statistic
  use highfrequency, only: remedi, rhy_cov, rtsvar
  implicit none
  integer :: t(8), lags(5), i
  real(dp) :: p1(8), p2(8), eval(3), times(7), returns(7)
  real(dp), allocatable :: vol(:), drift(:), burst(:), rem(:)
  type(lead_lag_result) :: ll

  t=[0,1,2,3,4,5,6,7]
  do i=1,8
    p1(i)=100.0_dp*exp(0.002_dp*real(i-1,dp)+0.001_dp*sin(real(i,dp)))
  end do
  p2(1)=p1(1)
  p2(2:8)=p1(1:7)
  lags=[-2,-1,0,1,2]
  ll=lead_lag(t,p1,t,p2,lags)
  if(size(ll%contrast)/=5)error stop 1
  if(abs(ll%optimal_lag)>2)error stop 1
  if(rhy_cov(t,p1,t,p1)<=0.0_dp)error stop 1
  if(rtsvar(p1,2,1)<0.0_dp)error stop 1

  times=real(t(2:8),dp)
  returns=log(p1(2:8)/p1(1:7))
  eval=[2.0_dp,4.0_dp,6.0_dp]
  vol=spot_volatility(times,returns,eval,2.0_dp,'gaussian')
  drift=spot_drift(times,returns,eval,2.0_dp,'mean','gaussian')
  burst=drift_burst_statistic(times,returns,eval,2.0_dp,'gaussian')
  if(any(vol<0.0_dp))error stop 1
  if(any(abs(burst)>100.0_dp))error stop 1
  rem=remedi(p1,1,[0,1])
  if(size(rem)/=2)error stop 1
  print '(a)', 'test_microstructure: PASS'
end program test_microstructure
