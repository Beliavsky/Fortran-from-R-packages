! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
! This program is free software under GNU GPL version 2 only.
program attribution_example
  use pa
  implicit none
  integer, parameter :: n = 9
  integer :: categorical(n,1), levels(1), status
  real(dp) :: numeric(n,2), wb(n), wp(n), beta(5), ret(n)
  real(dp), allocatable :: x(:, :), values(:)
  integer, allocatable :: group_start(:), group_end(:)
  type(regression_period_result) :: model

  categorical(:,1) = [1,2,3,1,2,3,1,2,3]
  levels = [3]
  numeric(:,1) = [-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,-0.8_dp,0.2_dp,0.7_dp,-0.3_dp]
  numeric(:,2) = [0.4_dp,-0.2_dp,0.1_dp,0.6_dp,-0.5_dp,0.3_dp,-0.7_dp,0.2_dp,0.8_dp]
  wb = [0.10_dp,0.12_dp,0.08_dp,0.13_dp,0.09_dp,0.11_dp,0.12_dp,0.15_dp,0.10_dp]
  wp = [0.12_dp,0.10_dp,0.08_dp,0.10_dp,0.12_dp,0.10_dp,0.14_dp,0.12_dp,0.12_dp]
  beta = [0.01_dp,0.02_dp,-0.005_dp,0.003_dp,-0.002_dp]

  call build_design_matrix(categorical, levels, numeric, x, group_start, group_end, status)
  ret = matmul(x, beta)
  call fit_regression_period(ret, x, wb, wp, model)
  call summarize_regression_period(model, group_start, group_end, values, status)

  print '(a,5f13.7)', 'factor returns: ', model%coefficients
  print '(a,5f13.7)', 'active exposure:', model%active_exposure
  print '(a,7f13.7)', 'attribution:    ', values
end program attribution_example
