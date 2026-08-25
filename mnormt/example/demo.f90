program demo
  use mnormt
  implicit none
  real(dp) :: mean(2), s(2,2), x(2), p
  real(dp) :: lower(2), upper(2)
  type(probability_result) :: pr

  mean = [0.0_dp, 0.0_dp]
  s = reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2,2])
  x = [0.5_dp, -0.25_dp]
  p = pmnorm(x, mean, s)
  print '(a,f12.8)', 'P(X <= x) = ', p

  lower = [-1.0_dp, -2.0_dp]
  upper = [0.5_dp, 1.5_dp]
  pr = sadmvt_prob(5.0_dp, lower, upper, mean, s)
  print '(a,f12.8)', 't rectangle probability = ', pr%value
end program demo
