! SPDX-License-Identifier: GPL-2.0-or-later
program stepmented_example
  use segmented
  implicit none
  integer, parameter :: n = 100
  real(dp) :: y(n), x(n, 1), z(n, 1), psi(1), value
  type(segmented_result) :: fit
  integer :: i
  do i = 1, n
    value = real(i - 1, dp) / real(n - 1, dp)
    x(i, 1) = 1.0_dp
    z(i, 1) = value
    y(i) = 2.0_dp + merge(1.5_dp, 0.0_dp, value > 0.62_dp) + &
        0.02_dp * cos(real(i, dp))
  end do
  psi = 0.4_dp
  call stepmented_lm(y, x, z, psi, fit)
  write(*, '(a,f9.5)') 'estimated change point: ', fit%breakpoints(1)
  write(*, '(a,2f10.5)') 'level and jump: ', fit%coefficients
end program stepmented_example
