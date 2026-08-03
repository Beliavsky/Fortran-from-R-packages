! SPDX-License-Identifier: GPL-2.0-or-later
program segmented_linear_example
  use segmented
  implicit none
  integer, parameter :: n = 120
  real(dp) :: y(n), x(n, 2), z(n, 1), psi(1), value
  type(segmented_result) :: fit
  integer :: i
  do i = 1, n
    value = 10.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    z(i, 1) = value
    y(i) = 1.0_dp + 0.5_dp * value + 1.4_dp * max(0.0_dp, value - 4.0_dp) + &
        0.04_dp * sin(real(i, dp))
  end do
  psi = 6.0_dp
  call segmented_lm(y, x, z, psi, fit)
  write(*, '(a,f9.4)') 'estimated breakpoint: ', fit%breakpoints(1)
  write(*, '(a,3f10.5)') 'coefficients: ', fit%coefficients
end program segmented_linear_example
