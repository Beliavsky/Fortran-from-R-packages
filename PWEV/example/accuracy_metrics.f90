program accuracy_metrics
  use pwev
  implicit none
  real(dp) :: actual(5), predicted(5), metric(9)
  actual = [10.0_dp, 12.0_dp, 11.0_dp, 15.0_dp, 14.0_dp]
  predicted = [9.5_dp, 12.5_dp, 10.8_dp, 14.2_dp, 14.4_dp]
  call pwev_metric_vector(actual, predicted, metric)
  print '(a,f10.5)', 'RMSE: ', metric(1)
  print '(a,f10.5)', 'MAPE: ', metric(2)
  print '(a,f10.5)', 'R squared: ', metric(9)
end program accuracy_metrics
