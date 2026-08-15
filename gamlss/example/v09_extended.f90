program v09_extended
   use gamlss
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_GA
   implicit none
   real(dp) :: eta(3),sigma_b
   type(marginal_prediction_result_t) :: exact,quantile,conditional

   eta=[log(1.0_dp),log(2.0_dp),log(4.0_dp)]
   sigma_b=0.6_dp
   call marginal_predict_eta(eta,sigma_b,GAMLSS_GA,1,conditional,MARGINAL_NONE)
   call marginal_predict_eta(eta,sigma_b,GAMLSS_GA,1,exact,MARGINAL_INTEGRATE)
   call marginal_predict_eta(eta,sigma_b,GAMLSS_GA,1,quantile,MARGINAL_QFUNCTION)

   write(*,'(a,3f10.4)') 'Conditional means:          ',conditional%fitted
   write(*,'(a,3f10.4)') 'Marginal means (integrate): ',exact%fitted
   write(*,'(a,3f10.4)') 'Marginal means (qfunction): ',quantile%fitted
end program v09_extended
