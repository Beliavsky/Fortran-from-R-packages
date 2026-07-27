! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

program simulate_and_fit
  use betategarch, only : dp, tegarch_parameters, tegarch_filter_result, tegarch_fit_result, &
    tegarch_simulate, tegarch_fit, tegarch_forecast, tegarch_bic_per_observation
  implicit none

  type(tegarch_parameters) :: truth
  type(tegarch_filter_result) :: simulated
  type(tegarch_fit_result) :: fitted
  real(dp), allocatable :: sigma_forecast(:), stdev_forecast(:)
  integer :: i

  truth%components = 1
  truth%asym = .true.
  truth%skewed = .true.
  truth%omega = 0.01_dp
  truth%phi1 = 0.90_dp
  truth%kappa1 = 0.08_dp
  truth%kappastar = 0.03_dp
  truth%df = 9.0_dp
  truth%skew = 0.85_dp

  call tegarch_simulate(500, truth, simulated, seed=123)
  call tegarch_fit(simulated%y, 1, .true., .true., fitted, compute_hessian=.false., &
    max_iterations=1200, tolerance=1.0e-6_dp)

  print '(a)', "Beta-Skew-t-EGARCH one-component demonstration"
  print '(a,1x,f12.6)', "log-likelihood:", fitted%log_likelihood
  print '(a,1x,f12.6)', "BIC per observation:", &
    tegarch_bic_per_observation(fitted%log_likelihood, size(simulated%y), size(fitted%free_parameters))
  print '(a,1x,i0)', "optimizer convergence code:", fitted%convergence
  print '(a,1x,i0)', "iterations:", fitted%iterations
  print '(a,1x,i0)', "evaluations:", fitted%evaluations
  print '(a)', "free parameter estimates:"
  do i = 1, size(fitted%free_parameters)
    print '(i3,1x,es18.10)', i, fitted%free_parameters(i)
  end do

  allocate(sigma_forecast(5), stdev_forecast(5))
  call tegarch_forecast(simulated%y, fitted%parameters, 5, sigma_forecast, stdev_forecast, &
    n_sim=20000, seed=456)
  print '(a)', "five-step conditional standard-deviation forecast:"
  do i = 1, 5
    print '(i3,1x,es18.10)', i, stdev_forecast(i)
  end do
end program simulate_and_fit
