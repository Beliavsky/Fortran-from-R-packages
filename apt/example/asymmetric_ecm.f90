! SPDX-License-Identifier: GPL-2.0-or-later
program asymmetric_ecm
  use apt, only : dp, apt_mtar, ecm_fit_result, ecm_asymmetry_test_result, &
    ecm_asymmetric_fit, ecm_asymmetry_tests
  implicit none
  integer, parameter :: n = 80
  real(dp) :: x(n), y(n), z(n), dx, previous_dz, phi
  integer :: i
  type(ecm_fit_result) :: fit
  type(ecm_asymmetry_test_result) :: tests

  x(1) = 10.0_dp
  z(1) = 0.8_dp
  previous_dz = 0.0_dp
  do i = 2, n
    dx = 0.12_dp + 0.22_dp*sin(0.37_dp*real(i,dp)) + &
      0.08_dp*cos(0.11_dp*real(i,dp))
    x(i) = x(i-1) + dx
    if (z(i-1) >= 0.0_dp) then
      phi = -0.18_dp
    else
      phi = -0.42_dp
    end if
    z(i) = z(i-1) + phi*z(i-1) + 0.22_dp*previous_dz + &
      0.07_dp*sin(0.73_dp*real(i,dp))
    previous_dz = z(i)-z(i-1)
  end do
  do i = 1, n
    y(i) = 2.5_dp + 1.35_dp*x(i) + z(i) + &
      0.03_dp*cos(0.21_dp*real(i,dp))
  end do

  call ecm_asymmetric_fit(y, x, fit, lag=2, split=.true., &
    model=apt_mtar, threshold=0.0_dp)
  if (fit%status /= 0) error stop 'asymmetric ECM fit failed'
  call ecm_asymmetry_tests(fit, tests)
  if (tests%status /= 0) error stop 'asymmetry tests failed'
  print '(a,*(f11.6,1x))', 'Equation-y coefficients: ', fit%equation_y%coefficients
  print '(a,f10.4,a,f10.6)', 'ECT symmetry F = ', &
    tests%tests(1)%equation_y%f_statistic, ', p = ', &
    tests%tests(1)%equation_y%p_value
end program asymmetric_ecm
