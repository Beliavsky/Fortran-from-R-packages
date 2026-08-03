program volatility_weighting
   use iso_fortran_env, only : int64
   use quarks
   implicit none
   real(dp) :: returns(120)
   type(risk_result) :: weighted, filtered
   type(rng_state) :: rng
   integer :: i

   do i = 1, size(returns)
      returns(i) = 0.008_dp * sin(0.19_dp * real(i, dp)) + &
         0.003_dp * cos(0.047_dp * real(i, dp))
   end do
   weighted = vwhs(returns, p=0.975_dp, model=volatility_ewma, lambda=0.94_dp)
   call seed_rng(rng, 20260802_int64)
   filtered = fhs(returns, p=0.975_dp, model=volatility_ewma, lambda=0.94_dp, &
      nboot=20000, rng=rng)

   print '(a,f12.6)', 'volatility-weighted VaR: ', weighted%var
   print '(a,f12.6)', 'volatility-weighted ES : ', weighted%es
   print '(a,f12.6)', 'filtered VaR           : ', filtered%var
   print '(a,f12.6)', 'filtered ES            : ', filtered%es
end program volatility_weighting
