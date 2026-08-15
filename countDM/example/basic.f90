program basic
   use countdm
   implicit none
   integer, allocatable :: x(:)
   type(bell_closed_result_t) :: closed
   type(mle_result_t) :: bt

   call get_data_sbirth(x)
   closed = bell_mle_closed(x)
   bt = mle_bt(x, 0.3_dp, 1.0_dp)

   print '(a,f10.6)', 'Closed-form Bell theta: ', closed%theta
   print '(a,2f10.6)', 'Bell-Touchard lambda, theta: ', bt%estimate
   print '(a,f12.4)', 'Bell-Touchard AIC: ', bt%aic
   print '(a,es14.6)', 'P_BT(X=2; lambda=2, theta=2): ', dbellt(2, 2.0_dp, 2.0_dp)
end program basic
