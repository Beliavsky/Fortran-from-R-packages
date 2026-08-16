program basic
  use mbbefd, only : dp, dmbbefd, pmbbefd, qmbbefd, ecmbbefd, dr_fit_result, fit_dr, rmbbefd
  implicit none
  real(dp) :: x(500)
  type(dr_fit_result) :: fit
  print '(a,f12.8)', 'density at 0.4 = ', dmbbefd(0.4_dp,0.5_dp,0.3_dp)
  print '(a,f12.8)', 'CDF at 0.4     = ', pmbbefd(0.4_dp,0.5_dp,0.3_dp)
  print '(a,f12.8)', 'Q(0.25)        = ', qmbbefd(0.25_dp,0.5_dp,0.3_dp)
  print '(a,f12.8)', 'exposure G(.5) = ', ecmbbefd(0.5_dp,0.5_dp,0.3_dp)
  call rmbbefd(x,0.5_dp,0.3_dp)
  call fit_dr(x,'mbbefd',fit,method='tlmme')
  print '(a,2f12.6)', 'TLMMe estimate = ', fit%estimate
end program basic
