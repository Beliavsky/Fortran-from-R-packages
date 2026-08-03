! SPDX-License-Identifier: AGPL-3.0-or-later
program demo_ren
  use ren, only : dp, analysis_options, analysis_result, perform_analysis, ren_success
  implicit none
  integer, parameter :: months = 8, days = 10, assets = 5, n = months * days
  real(dp) :: returns(n, assets)
  integer :: month(n), date(n)
  type(analysis_options) :: options
  type(analysis_result) :: result
  integer :: i, j, k
  k = 0
  do i = 1, months
    do j = 1, days
      k = k + 1
      month(k) = i
      date(k) = 20200100 + 100 * (i - 1) + j
      returns(k, :) = [(0.05_dp + 0.25_dp * sin(0.13_dp * real(k * j, dp)) + &
        0.02_dp * real(j, dp), j=1,assets)]
    end do
  end do
  options%cluster_repetitions = 10
  options%stochastic_samples = 25
  call perform_analysis(returns, month, date, result, options)
  if (result%status /= ren_success) error stop 'REN demo failed'
  print '(a)', 'method                 turnover%       sharpe    volatility%       maxdd%'
  do j = 1, size(result%method)
    print '(a20,4(1x,f13.6))', trim(result%method(j)), result%turnover_mean(j), &
      result%sharpe_ratio(j), result%volatility(j), result%max_drawdown(j)
  end do
end program demo_ren
