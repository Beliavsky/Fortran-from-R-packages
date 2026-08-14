program basic_isotone
   use isotone
   implicit none
   real(dp) :: y(8)
   type(gpava_result) :: fit
   y=[3.0_dp,1.0_dp,2.0_dp,5.0_dp,4.0_dp,7.0_dp,6.0_dp,8.0_dp]
   call gpava_fit(y,fit)
   print '(a,8f8.3)', 'observed: ', y
   print '(a,8f8.3)', 'isotonic: ', fit%x
end program
