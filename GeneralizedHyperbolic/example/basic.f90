program basic
  use generalized_hyperbolic
  implicit none
  real(dp) :: x(5)

  print '(a,f12.8)', 'NIG density at zero: ', dnig(0.0_dp,0.0_dp,1.0_dp,2.0_dp,0.0_dp)
  call rhyperb(x,0.0_dp,1.0_dp,2.0_dp,0.2_dp)
  print '(a,5f10.4)', 'Hyperbolic draws: ', x
end program basic
