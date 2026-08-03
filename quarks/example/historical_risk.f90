program historical_risk
   use quarks
   implicit none
   real(dp) :: returns(12)
   type(risk_result) :: plain, age

   returns = [0.010_dp, -0.008_dp, 0.004_dp, -0.020_dp, 0.013_dp, -0.005_dp, &
      0.006_dp, -0.012_dp, 0.009_dp, -0.018_dp, 0.003_dp, -0.007_dp]
   plain = hs(returns, p=0.95_dp, method=method_plain)
   age = hs(returns, p=0.95_dp, method=method_age, lambda=0.98_dp)

   print '(a,f12.6)', 'plain VaR: ', plain%var
   print '(a,f12.6)', 'plain ES : ', plain%es
   print '(a,f12.6)', 'age VaR  : ', age%var
   print '(a,f12.6)', 'age ES   : ', age%es
end program historical_risk
