program demo
  use mev, only: dp, dgev, pgev, qgev, rgp, gpd_fit, mev_fit_result
  implicit none
  real(dp) :: x(1000), q
  type(mev_fit_result) :: fit

  q = qgev(0.99_dp, loc=0.0_dp, scale=1.0_dp, shape=0.1_dp)
  print '(a,f12.6)', 'GEV 0.99 quantile: ', q
  print '(a,es14.6)', 'density at quantile: ', dgev(q,0.0_dp,1.0_dp,0.1_dp)
  print '(a,f12.8)', 'CDF check: ', pgev(q,0.0_dp,1.0_dp,0.1_dp)

  call rgp(size(x), x, scale=1.5_dp, shape=0.15_dp)
  call gpd_fit(x, fit)
  print '(a,2f12.6)', 'GPD fit (scale, shape): ', fit%estimate
end program demo
