! SPDX-License-Identifier: GPL-2.0-or-later
program apt_demo
  use apt, only : dp, apt_mtar, ci_tar_fit_result, ci_tar_threshold_result, &
    ecm_fit_result, ecm_diagnostics_result, ci_tar_fit, ci_tar_threshold, &
    ecm_asymmetric_fit, ecm_diagnostics
  implicit none
  integer, parameter :: n = 80
  real(dp) :: x(n), y(n), z(n), dx, dzlag, phi, previous_dz
  integer :: t
  type(ci_tar_threshold_result) :: threshold_search
  type(ci_tar_fit_result) :: cointegration
  type(ecm_fit_result) :: ecm
  type(ecm_diagnostics_result) :: diagnostics

  x(1) = 10.0_dp
  z(1) = 0.8_dp
  previous_dz = 0.0_dp
  do t = 2, n
    dx = 0.12_dp + 0.22_dp*sin(0.37_dp*real(t,dp)) + &
      0.08_dp*cos(0.11_dp*real(t,dp))
    x(t) = x(t-1) + dx
    dzlag = previous_dz
    if (z(t-1) >= 0.0_dp) then
      phi = -0.18_dp
    else
      phi = -0.42_dp
    end if
    z(t) = z(t-1) + phi*z(t-1) + 0.22_dp*dzlag + &
      0.07_dp*sin(0.73_dp*real(t,dp))
    previous_dz = z(t) - z(t-1)
  end do
  do t = 1, n
    y(t) = 2.5_dp + 1.35_dp*x(t) + z(t) + &
      0.03_dp*cos(0.21_dp*real(t,dp))
  end do

  call ci_tar_threshold(y, x, threshold_search, apt_mtar, lag=1)
  call ci_tar_fit(y, x, cointegration, apt_mtar, lag=1, &
    threshold=threshold_search%threshold)
  call ecm_asymmetric_fit(y, x, ecm, lag=2, split=.true., &
    model=apt_mtar, threshold=threshold_search%threshold)
  call ecm_diagnostics(ecm, diagnostics)

  print '(a,f12.6)', 'Selected MTAR threshold: ', threshold_search%threshold
  print '(a,*(f12.6,1x))', 'Threshold-regression coefficients: ', &
    cointegration%threshold_regression%coefficients
  print '(a,f12.6)', 'No-cointegration F statistic: ', &
    cointegration%no_cointegration_test%f_statistic
  print '(a,f12.6)', 'Long-run symmetry F statistic: ', &
    cointegration%symmetry_test%f_statistic
  print '(a,*(f12.6,1x))', 'ECM equation-x coefficients: ', &
    ecm%equation_x%coefficients
  print '(a,f12.6)', 'Equation-x Durbin-Watson: ', &
    diagnostics%equation_x%durbin_watson
end program apt_demo
