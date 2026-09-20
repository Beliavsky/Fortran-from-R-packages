program basic_stats
  use hmisc
  implicit none
  real(dp) :: x(5), w(5)
  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 8.0_dp]
  w = [1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp]
  print '(a,f8.4)', 'Gini mean difference: ', gini_md(x)
  print '(a,f8.4)', 'Weighted mean:        ', weighted_mean(x,w)
  print '(a,f8.4)', 'Weighted median:      ', weighted_quantile(x,w,0.5_dp,.false.)
  print '(a,f8.4)', 'Pseudomedian:         ', pseudomedian(x)
end program basic_stats
