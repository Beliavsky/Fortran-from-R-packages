! SPDX-License-Identifier: AGPL-3.0-or-later
program test_analysis
  use ren, only : dp, analysis_options, analysis_result, perform_analysis, ren_success, ren_method_count
  implicit none
  integer, parameter :: nmonth = 8, per_month = 8, n = nmonth * per_month, p = 5
  real(dp) :: x(n, p)
  integer :: month(n), date(n)
  type(analysis_options) :: options
  type(analysis_result) :: result
  integer :: i, j, row
  row = 0
  do i = 1, nmonth
    do j = 1, per_month
      row = row + 1
      month(row) = i
      date(row) = 20200100 + min(28, j) + 100 * (i - 1)
      x(row, :) = [(0.2_dp * sin(0.11_dp * real(row * j, dp)) + &
        0.05_dp * real(j, dp) + 0.01_dp * real(i, dp), j=1,p)]
    end do
  end do
  options%cluster_repetitions = 3
  options%stochastic_samples = 5
  options%random_seed = 91
  call perform_analysis(x, month, date, result, options)
  if (result%status /= ren_success) error stop 'analysis status'
  if (size(result%turnover, 2) /= ren_method_count) error stop 'analysis methods'
  if (size(result%gross_returns, 1) /= 2 * per_month) error stop 'analysis return rows'
  do j = 1, ren_method_count
    if (abs(sum(result%weights(1, :, j)) - 1.0_dp) > 1.0e-6_dp) error stop 'analysis weight sum'
  end do
  print '(a)', 'test_analysis: PASS'
end program test_analysis
