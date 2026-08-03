! SPDX-License-Identifier: GPL-2.0-or-later
program test_stepmented_glm
  use segmented
  implicit none
  integer, parameter :: n = 160
  real(dp) :: y(n), yb(n), yp(n), x(n, 2), xi(n, 1), z(n, 1), psi(1)
  real(dp) :: value, probability, lambda
  type(segmented_result) :: fit_step, fit_binom, fit_poisson
  integer :: i

  do i = 1, n
    value = 8.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    xi(i, 1) = 1.0_dp
    z(i, 1) = value
    y(i) = 2.0_dp + merge(3.0_dp, 0.0_dp, value > 5.0_dp) + &
        0.03_dp * sin(real(i, dp))
    probability = 1.0_dp / (1.0_dp + exp(-(-1.8_dp + 0.35_dp * value + &
        0.9_dp * max(0.0_dp, value - 3.5_dp))))
    yb(i) = min(0.999_dp, max(0.001_dp, probability + 0.015_dp * sin(real(i, dp))))
    lambda = exp(0.2_dp + 0.12_dp * value + 0.3_dp * max(0.0_dp, value - 4.0_dp))
    yp(i) = lambda
  end do
  psi = 3.0_dp
  call fit_stepmented_lm(y, xi, z, psi, fit_step)
  call check(fit_step%status == SEG_SUCCESS, 'stepmented LM status')
  call check(abs(fit_step%breakpoints(1) - 5.0_dp) < 0.06_dp, 'step breakpoint')
  call check(maxval(abs(fit_step%coefficients - [2.0_dp, 3.0_dp])) < 0.01_dp, &
      'step coefficients')

  psi = 5.5_dp
  call fit_segmented_glm(yb, x, z, psi, FAMILY_BINOMIAL, fit_binom)
  call check(fit_binom%status == SEG_SUCCESS, 'binomial segmented status')
  call check(abs(fit_binom%breakpoints(1) - 3.5_dp) < 0.08_dp, &
      'binomial breakpoint')
  call check(maxval(abs(fit_binom%coefficients - [-1.8_dp, 0.35_dp, 0.9_dp])) &
      < 0.03_dp, 'binomial coefficients')

  psi = 5.0_dp
  call fit_segmented_glm(yp, x, z, psi, FAMILY_POISSON, fit_poisson)
  call check(fit_poisson%status == SEG_SUCCESS, 'Poisson segmented status')
  call check(abs(fit_poisson%breakpoints(1) - 4.0_dp) < 0.08_dp, &
      'Poisson breakpoint')
  call check(maxval(abs(fit_poisson%coefficients - [0.2_dp, 0.12_dp, 0.3_dp])) &
      < 0.03_dp, 'Poisson coefficients')
  print *, 'test_stepmented_glm: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine check
end program test_stepmented_glm
