! SPDX-License-Identifier: GPL-2.0-or-later
program demo_segmented
  use segmented
  implicit none
  integer, parameter :: n = 100
  real(dp) :: y(n), x(n, 2), z(n, 1), psi(1), value
  real(dp), allocatable :: slopes(:), slope_se(:)
  type(segmented_result) :: fit
  type(test_result) :: test
  integer :: i, status
  do i = 1, n
    value = 10.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    z(i, 1) = value
    y(i) = 0.5_dp + 0.35_dp * value + 1.2_dp * max(0.0_dp, value - 4.5_dp) + &
        0.05_dp * sin(real(i, dp))
  end do
  psi = 6.0_dp
  call segmented_lm(y, x, z, psi, fit)
  call segment_slopes(fit, 2, slopes, slope_se, status)
  call davies_test(y, x, z(:, 1), 30, test)
  write(*, '(a,f8.4)') 'breakpoint = ', fit%breakpoints(1)
  write(*, '(a,2f9.4)') 'segment slopes = ', slopes
  write(*, '(a,es12.4)') 'Davies-grid p-value = ', test%p_value
end program demo_segmented
