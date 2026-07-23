! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program example_autocorrelation
  use timsac, only: dp, autocorrelation, autocorrelation_result
  implicit none

  real(dp) :: x(10)
  integer :: lag
  type(autocorrelation_result) :: result

  x = [1.0_dp, 2.0_dp, 1.5_dp, 3.0_dp, 2.5_dp, &
       4.0_dp, 3.5_dp, 5.0_dp, 4.5_dp, 6.0_dp]

  result = autocorrelation(x, max_lag=4)

  print '(a,f10.5)', 'mean = ', result%mean
  print '(a)', 'lag   covariance   correlation'
  do lag = lbound(result%covariance, 1), ubound(result%covariance, 1)
    print '(i3,2f14.6)', lag, result%covariance(lag), result%correlation(lag)
  end do
end program example_autocorrelation
