! SPDX-License-Identifier: MIT
program jump_tests_example
  use jumptest, only : dp, statp_result, adjp_result, jumptestday, jumptestperiod
  implicit none

  real(dp) :: returns(12), return_matrix(12, 3)
  type(statp_result) :: daily
  type(adjp_result) :: period

  returns = [0.001_dp, -0.002_dp, 0.0015_dp, 0.0005_dp, -0.001_dp, &
    0.012_dp, -0.0015_dp, 0.0008_dp, -0.0004_dp, 0.0011_dp, &
    -0.0007_dp, 0.0006_dp]
  return_matrix(:, 1) = returns
  return_matrix(:, 2) = 0.7_dp*returns
  return_matrix(:, 3) = returns
  return_matrix(6, 3) = 0.025_dp

  call jumptestday(returns, daily, 'BNS')
  print '(a,f10.5,a,f10.6)', 'BNS statistic = ', daily%stat, ', p = ', daily%pvalue

  call jumptestperiod(return_matrix, period, 'Amed')
  print '(a)', 'Amed tests by period:'
  print '(3(f10.5,1x))', period%stat
  print '(a)', 'BH-adjusted p-values:'
  print '(3(f10.6,1x))', period%adjp
end program jump_tests_example
