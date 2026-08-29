program lmtest_example
   use lmtest, only : dp, lm_result, test_result, lm_fit, &
      breusch_pagan_test, durbin_watson_test
   implicit none
   integer, parameter :: n = 60
   real(dp) :: x(n,3), y(n), t
   integer :: i
   type(lm_result) :: fit
   type(test_result) :: bp, dw

   do i = 1, n
      t = real(i,dp)
      x(i,:) = [1.0_dp, t/real(n,dp), sin(0.2_dp*t)]
      y(i) = 1.0_dp + 0.8_dp*x(i,2) - 0.4_dp*x(i,3) + &
         0.15_dp*sin(0.8_dp*t)
   end do

   fit = lm_fit(x,y)
   bp = breusch_pagan_test(x,y,x)
   dw = durbin_watson_test(x,y,exact=.false.)

   print '(a,3f12.6)', 'coefficients: ', fit%beta
   print '(a,f12.6,a,es12.4)', 'Breusch-Pagan: ', bp%statistic, &
      '  p=', bp%p_value
   print '(a,f12.6,a,es12.4)', 'Durbin-Watson: ', dw%statistic, &
      '  p=', dw%p_value
end program lmtest_example
