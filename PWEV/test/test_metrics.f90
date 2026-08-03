program test_metrics
  use pwev, only : dp, pwev_metric_vector
  implicit none
  real(dp) :: actual(3), predicted(3), metric(9), expected(9)
  actual = [1.0_dp, 2.0_dp, 4.0_dp]
  predicted = [1.0_dp, 3.0_dp, 2.0_dp]
  expected = [1.2909944487358056_dp, 0.3333333333333333_dp, 1.0_dp, &
    1.0350983390135313_dp, 1.0_dp, 0.3384788485297446_dp, 0.9_dp, &
    0.3555555555555556_dp, 0.1071428571428571_dp]
  call pwev_metric_vector(actual, predicted, metric)
  if (maxval(abs(metric - expected)) > 1.0e-12_dp) error stop 'metric reference values failed'
  print *, 'test_metrics: PASS'
end program test_metrics
