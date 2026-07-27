! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program highfrequency_demo
  use highfrequency, only: dp, make_returns, rrvar, rbpvar, rskew, rkurt
  use highfrequency, only: fit_har, har_model, har_forecast
  implicit none
  real(dp) :: prices(9), daily_rv(30), forecast
  real(dp), allocatable :: returns(:)
  type(har_model) :: model
  integer :: i

  prices=[100.0_dp,100.2_dp,100.1_dp,100.5_dp,100.4_dp,100.8_dp,100.7_dp,101.0_dp,101.2_dp]
  returns=make_returns(prices)
  print '(a,f12.8)', 'realized variance: ',rrvar(returns)
  print '(a,f12.8)', 'bipower variation: ',rbpvar(returns)
  print '(a,f12.8)', 'realized skewness: ',rskew(returns)
  print '(a,f12.8)', 'realized kurtosis: ',rkurt(returns)

  daily_rv(1)=0.0001_dp
  do i=2,size(daily_rv)
    daily_rv(i)=0.00002_dp+0.75_dp*daily_rv(i-1)+0.00001_dp*(1.0_dp+sin(real(i,dp)))
  end do
  model=fit_har(daily_rv,[1,5,10])
  forecast=har_forecast(model,daily_rv,[1,5,10])
  print '(a,es14.6)', 'next HAR variance forecast: ',forecast
end program highfrequency_demo
