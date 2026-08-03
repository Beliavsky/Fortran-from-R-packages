! SPDX-License-Identifier: GPL-2.0-or-later
program test_segmented_lme
  use segmented
  implicit none
  integer, parameter :: groups = 10, observations = 10, n = groups * observations
  real(dp) :: y(n), x(n, 2), z(n, 1), random_design(n, 1), psi(1)
  real(dp) :: value, random_intercept
  integer :: group(n), g, i, k
  type(segmented_lme_result) :: fit
  type(segmented_lme_options) :: options
  type(segmented_control) :: control

  k = 0
  do g = 1, groups
    random_intercept = 0.18_dp * sin(real(g, dp))
    do i = 1, observations
      k = k + 1
      value = 8.0_dp * real(i - 1, dp) / real(observations - 1, dp)
      x(k, :) = [1.0_dp, value]
      z(k, 1) = value
      random_design(k, 1) = 1.0_dp
      group(k) = g
      y(k) = 1.0_dp + random_intercept + 0.4_dp * value + &
          max(0.0_dp, value - 3.7_dp) + 0.03_dp * sin(real(k, dp))
    end do
  end do
  psi = 5.5_dp
  control = segmented_control()
  control%max_iter = 12
  options = segmented_lme_options()
  options%control%max_iter = 80
  options%control%max_outer = 10
  call fit_segmented_lme(y, x, z, psi, random_design, group, fit, options, control)
  call check(fit%status == SEG_SUCCESS, 'segmented LME status')
  call check(fit%converged, 'segmented LME convergence')
  call check(abs(fit%breakpoints(1) - 3.7_dp) < 0.04_dp, 'mixed breakpoint')
  call check(maxval(abs(fit%fit%beta - [1.0_dp, 0.4_dp, 1.0_dp])) < 0.04_dp, &
      'mixed fixed effects')
  call check(size(fit%fit%random_effects, 1) == groups, 'BLUP group count')
  call check(fit%breakpoint_se(1) > 0.0_dp, 'mixed breakpoint standard error')
  print *, 'test_segmented_lme: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine check
end program test_segmented_lme
