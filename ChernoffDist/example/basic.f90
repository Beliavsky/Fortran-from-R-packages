program basic
  use chernoffdist, only: dp, dchern, pchern, qchern
  implicit none

  print '(a,f14.10)', 'density at 0 = ', dchern(0.0_dp)
  print '(a,f14.10)', 'cdf at 1     = ', pchern(1.0_dp)
  print '(a,f14.10)', 'q(.90)       = ', qchern(0.90_dp)
end program basic
