program rolling_backtest
   use quarks
   implicit none
   real(dp) :: returns(180), losses(30)
   type(rollcast_result) :: forecasts
   type(coverage_result) :: coverage
   type(traffic_result) :: traffic
   integer :: i

   do i = 1, size(returns)
      returns(i) = 0.009_dp * sin(0.27_dp * real(i, dp)) + &
         0.004_dp * cos(0.091_dp * real(i, dp))
   end do
   forecasts = rollcast(returns, p=0.95_dp, method=method_age, lambda=0.98_dp, &
      nout=30, nwin=100)
   losses = -forecasts%xout
   coverage = cvgtest(losses, forecasts%var, forecasts%p)
   traffic = trftest(losses, forecasts%var, forecasts%p)

   print '(a,i0)', 'violations: ', coverage%violations
   print '(a,f10.6)', 'Kupiec p-value: ', coverage%p_uc
   print '(a,f10.6)', 'independence p-value: ', coverage%p_ind
   print '(a,f10.6)', 'traffic-light cumulative probability: ', &
      traffic%cumulative_probability
end program rolling_backtest
