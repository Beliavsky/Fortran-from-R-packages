program basic
  use betafunctions
  implicit none
  type(moment_set) :: m
  real(dp) :: x

  x = beta4_pdf(0.5_dp, 0.0_dp, 1.0_dp, 5.0_dp, 3.0_dp)
  print '(a,f12.8)', 'beta4 pdf at 0.5 = ', x
  call beta_moments(5.0_dp, 3.0_dp, 0.0_dp, 1.0_dp, 4, m)
  print '(a,f12.8)', 'mean             = ', m%raw(1)
  print '(a,f12.8)', 'variance         = ', m%central(2)
end program basic
