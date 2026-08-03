program example_weighted
  use corpcor, only : dp, moments_result, wt_moments
  implicit none
  real(dp) :: x(5, 3), w(5)
  type(moments_result) :: moments

  x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
    2.0_dp, 1.0_dp, 4.0_dp, 6.0_dp, 7.0_dp, &
    5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], [5, 3])
  w = [0.10_dp, 0.15_dp, 0.20_dp, 0.25_dp, 0.30_dp]
  moments = wt_moments(x, w)
  print '(a,*(f10.5,1x))', 'weighted means:     ', moments%mean
  print '(a,*(f10.5,1x))', 'weighted variances: ', moments%variance
end program example_weighted
