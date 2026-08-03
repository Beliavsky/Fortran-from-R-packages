! SPDX-License-Identifier: GPL-2.0-or-later
program segmented_glm_example
  use segmented
  implicit none
  integer, parameter :: n = 150
  real(dp) :: y(n), x(n, 2), z(n, 1), psi(1), value, eta
  type(segmented_result) :: fit
  integer :: i
  do i = 1, n
    value = 8.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    z(i, 1) = value
    eta = -1.5_dp + 0.25_dp * value + 0.8_dp * max(0.0_dp, value - 3.2_dp)
    y(i) = 1.0_dp / (1.0_dp + exp(-eta))
  end do
  psi = 5.0_dp
  call segmented_glm(y, x, z, psi, FAMILY_BINOMIAL, fit)
  write(*, '(a,f9.4)') 'estimated logistic breakpoint: ', fit%breakpoints(1)
  write(*, '(a,3f10.5)') 'linear-predictor coefficients: ', fit%coefficients
end program segmented_glm_example
