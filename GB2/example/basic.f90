program basic
  use gb2, only : dp, dgb2, pgb2, qgb2, moment_gb2, main_gb2
  implicit none
  real(dp) :: ind(6)
  real(dp), parameter :: a=2.3_dp, b=4.2_dp, p=1.7_dp, q=3.4_dp

  print '(a,f12.8)', 'density at 3.1 = ', dgb2(3.1_dp,a,b,p,q)
  print '(a,f12.8)', 'cdf at 3.1     = ', pgb2(3.1_dp,a,b,p,q)
  print '(a,f12.8)', 'median         = ', qgb2(0.5_dp,a,b,p,q)
  print '(a,f12.8)', 'mean           = ', moment_gb2(1.0_dp,a,b,p,q)

  call main_gb2(0.6_dp,a,b,p,q,ind)
  print '(a,6(1x,f12.6))', 'indicators =', ind
end program basic
