program basic_sgt
  use sgt
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  implicit none
  real(dp) :: inf
  inf = ieee_value(0.0_dp, ieee_positive_inf)
  print '(a,f12.8)', 'density = ', dsgt(0.8_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)
  print '(a,f12.8)', 'cdf     = ', psgt(0.8_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)
  print '(a,f12.8)', 'q(.37)  = ', qsgt(0.37_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)
  print '(a,f12.8)', 'normal  = ', dsgt(1.2_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, inf)
end program basic_sgt
