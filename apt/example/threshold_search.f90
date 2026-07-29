! SPDX-License-Identifier: GPL-2.0-or-later
program threshold_search
  use apt, only : dp, apt_tar, ci_tar_threshold_result, ci_tar_threshold
  implicit none
  integer, parameter :: n = 60
  real(dp) :: x(n), y(n), e
  integer :: i
  type(ci_tar_threshold_result) :: result
  x(1) = 5.0_dp
  e = 0.5_dp
  y(1) = 1.0_dp + 1.8_dp*x(1) + e
  do i = 2, n
    x(i) = x(i-1) + 0.1_dp + 0.15_dp*sin(0.4_dp*real(i,dp))
    if (e >= 0.08_dp) then
      e = 0.72_dp*e + 0.04_dp*cos(0.7_dp*real(i,dp))
    else
      e = 0.45_dp*e + 0.04_dp*cos(0.7_dp*real(i,dp))
    end if
    y(i) = 1.0_dp + 1.8_dp*x(i) + e
  end do
  call ci_tar_threshold(y, x, result, apt_tar, lag=1, trim_fraction=0.15_dp)
  print '(a,f10.5)', 'Selected TAR threshold: ', result%threshold
  print '(a,f12.7)', 'Minimum SSE: ', result%minimum_sse
  print '(a,i0)', 'Thresholds examined: ', size(result%path_threshold)
end program threshold_search
