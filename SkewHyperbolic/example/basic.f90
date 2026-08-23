program basic
  use skewhyperbolic
  implicit none
  real(dp) :: x(5)
  call rskewhyp(x,0.0_dp,1.0_dp,0.5_dp,10.0_dp)
  print '(a,es14.6)','density at zero = ',dskewhyp(0.0_dp,0.0_dp,1.0_dp,0.5_dp,10.0_dp)
  print '(a,5(1x,f9.4))','draws:',x
end program
