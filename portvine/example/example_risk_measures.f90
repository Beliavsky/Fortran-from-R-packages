program example_risk_measures
   use portvine, only : dp, est_var, est_es, risk_es_mean
   implicit none
   real(dp) :: sample(10),alpha(2),value(2)
   sample=[-0.08_dp,-0.04_dp,-0.03_dp,-0.01_dp,0.0_dp,0.01_dp,0.02_dp,0.03_dp,0.05_dp,0.08_dp]
   alpha=[0.05_dp,0.10_dp]
   call est_var(sample,alpha,value)
   print '(a,2f10.5)', 'VaR:     ',value
   call est_es(sample,alpha,value,risk_es_mean)
   print '(a,2f10.5)', 'mean ES: ',value
end program example_risk_measures
