! SPDX-License-Identifier: GPL-2.0-or-later
program segmented_lme_example
  use segmented
  implicit none
  integer, parameter :: ng = 8, ni = 9, n = ng * ni
  real(dp) :: y(n), x(n, 2), z(n, 1), zr(n, 1), psi(1), value, b0
  integer :: group(n), g, i, k
  type(segmented_lme_result) :: fit
  type(segmented_lme_options) :: options
  type(segmented_control) :: control
  k = 0
  do g = 1, ng
    b0 = 0.15_dp * sin(real(g, dp))
    do i = 1, ni
      k = k + 1
      value = 7.0_dp * real(i - 1, dp) / real(ni - 1, dp)
      x(k, :) = [1.0_dp, value]
      z(k, 1) = value
      zr(k, 1) = 1.0_dp
      group(k) = g
      y(k) = 1.0_dp + b0 + 0.3_dp * value + &
          0.9_dp * max(0.0_dp, value - 3.0_dp)
    end do
  end do
  psi = 4.5_dp
  options = segmented_lme_options()
  options%control%max_iter = 80
  control = segmented_control()
  control%max_iter = 12
  call segmented_lme(y, x, z, psi, zr, group, fit, options, control)
  write(*, '(a,f9.4)') 'mixed-effects breakpoint: ', fit%breakpoints(1)
  write(*, '(a,3f10.5)') 'fixed effects: ', fit%fit%beta
end program segmented_lme_example
