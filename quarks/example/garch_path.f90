program garch_path
   use quarks
   implicit none
   real(dp) :: returns(160)
   type(risk_result) :: result
   integer :: i

   do i = 1, size(returns)
      returns(i) = (0.004_dp + 0.00002_dp * real(i, dp)) * &
         sin(0.31_dp * real(i, dp))
   end do
   result = vwhs(returns, p=0.95_dp, model=volatility_garch, &
      garch_max_iterations=400)
   print '(a,f12.6)', 'GARCH-weighted VaR: ', result%var
   print '(a,f12.6)', 'GARCH-weighted ES : ', result%es
   print '(a,l1)', 'EWMA fallback used: ', result%used_fallback
end program garch_path
