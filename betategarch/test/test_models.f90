! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

program test_models
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use betategarch, only : dp, tegarch_parameters, tegarch_filter_result, tegarch_fit_result, &
    tegarch_default_parameters, tegarch_params_to_free, tegarch_simulate, tegarch_filter, &
    tegarch_loglik, tegarch_fit, tegarch_forecast, skew_t_variance, &
    tegarch_bic_per_observation, tegarch_standard_errors
  implicit none

  type(tegarch_parameters) :: p_sym, p_full, p_two, default_p
  type(tegarch_filter_result) :: sim_sym, sim_full, sim_two, filtered
  type(tegarch_fit_result) :: fit_sym, fit_full, fit_two
  real(dp), allocatable :: initial(:), sigma(:), stdev(:), standard_errors(:)
  real(dp) :: initial_loglik, expected_one_step, bic
  integer :: se_info

  p_sym%components = 1
  p_sym%asym = .false.
  p_sym%skewed = .false.
  p_sym%omega = 0.01_dp
  p_sym%phi1 = 0.86_dp
  p_sym%kappa1 = 0.07_dp
  p_sym%kappastar = 0.0_dp
  p_sym%df = 8.0_dp
  p_sym%skew = 1.0_dp
  call tegarch_simulate(260, p_sym, sim_sym, seed=7101)
  default_p = tegarch_default_parameters(1, .false., .false.)
  allocate(initial(4))
  call tegarch_params_to_free(default_p, initial)
  initial_loglik = tegarch_loglik(sim_sym%y, default_p)
  call tegarch_fit(sim_sym%y, 1, .false., .false., fit_sym, initial=initial, &
    compute_hessian=.true., max_iterations=700, tolerance=2.0e-6_dp)
  call assert_fit_improved(fit_sym, initial_loglik, "symmetric one-component fit")
  if (.not. fit_sym%hessian_available) error stop "Hessian was not computed"
  if (maxval(abs(fit_sym%hessian - transpose(fit_sym%hessian))) > 1.0e-8_dp) then
    error stop "Hessian is not symmetric"
  end if
  if (.not. fit_sym%covariance_available) error stop "Covariance matrix inversion failed"
  if (.not. all(ieee_is_finite(fit_sym%covariance))) error stop "Covariance contains non-finite values"
  allocate(standard_errors(size(fit_sym%free_parameters)))
  call tegarch_standard_errors(fit_sym, standard_errors, se_info)
  if (se_info /= 0 .or. any(standard_errors < 0.0_dp)) error stop "standard error extraction failed"
  bic = tegarch_bic_per_observation(fit_sym%log_likelihood, size(sim_sym%y), size(fit_sym%free_parameters))
  if (.not. ieee_is_finite(bic)) error stop "BIC calculation failed"
  deallocate(standard_errors)
  deallocate(initial)

  p_full%components = 1
  p_full%asym = .true.
  p_full%skewed = .true.
  p_full%omega = 0.015_dp
  p_full%phi1 = 0.90_dp
  p_full%kappa1 = 0.08_dp
  p_full%kappastar = 0.035_dp
  p_full%df = 9.0_dp
  p_full%skew = 0.85_dp
  call tegarch_simulate(220, p_full, sim_full, seed=7102)
  default_p = tegarch_default_parameters(1, .true., .true.)
  allocate(initial(6))
  call tegarch_params_to_free(default_p, initial)
  initial_loglik = tegarch_loglik(sim_full%y, default_p)
  call tegarch_fit(sim_full%y, 1, .true., .true., fit_full, initial=initial, &
    compute_hessian=.false., max_iterations=900, tolerance=2.0e-6_dp)
  call assert_fit_improved(fit_full, initial_loglik, "full one-component fit")
  deallocate(initial)

  p_two%components = 2
  p_two%asym = .true.
  p_two%skewed = .true.
  p_two%omega = 0.01_dp
  p_two%phi1 = 0.92_dp
  p_two%phi2 = 0.65_dp
  p_two%kappa1 = 0.02_dp
  p_two%kappa2 = 0.05_dp
  p_two%kappastar = 0.02_dp
  p_two%df = 10.0_dp
  p_two%skew = 0.90_dp
  call tegarch_simulate(220, p_two, sim_two, seed=7103)
  default_p = tegarch_default_parameters(2, .true., .true.)
  allocate(initial(8))
  call tegarch_params_to_free(default_p, initial)
  initial_loglik = tegarch_loglik(sim_two%y, default_p)
  call tegarch_fit(sim_two%y, 2, .true., .true., fit_two, initial=initial, &
    compute_hessian=.false., max_iterations=1200, tolerance=3.0e-6_dp)
  call assert_fit_improved(fit_two, initial_loglik, "two-component fit")
  deallocate(initial)

  allocate(sigma(5), stdev(5))
  call tegarch_forecast(sim_full%y, p_full, 5, sigma, stdev, n_sim=15000, seed=7201)
  if (.not. all(ieee_is_finite(sigma)) .or. any(sigma <= 0.0_dp)) error stop "invalid one-component forecast"
  if (.not. all(ieee_is_finite(stdev)) .or. any(stdev <= 0.0_dp)) error stop "invalid one-component stdev forecast"
  call tegarch_filter(sim_full%y, p_full, filtered)
  expected_one_step = exp(p_full%omega + p_full%phi1*filtered%lambda1_dagger(size(sim_full%y)) + &
    p_full%kappa1*last_score(sim_full%y, filtered%lambda(size(sim_full%y)), p_full) + &
    p_full%kappastar*sign_real(-sim_full%y(size(sim_full%y)))*(last_score(sim_full%y, &
    filtered%lambda(size(sim_full%y)), p_full) + 1.0_dp))
  call assert_close(sigma(1), expected_one_step, 2.0e-13_dp, "one-step forecast")
  call assert_close(stdev(1), sigma(1)*sqrt(skew_t_variance(p_full%df, p_full%skew)), &
    2.0e-13_dp, "one-step standard deviation")

  call tegarch_forecast(sim_two%y, p_two, 5, sigma, stdev, n_sim=15000, seed=7202)
  if (.not. all(ieee_is_finite(sigma)) .or. any(sigma <= 0.0_dp)) error stop "invalid two-component forecast"
  if (.not. all(ieee_is_finite(stdev)) .or. any(stdev <= 0.0_dp)) error stop "invalid two-component stdev forecast"

  print '(a)', "Fitting, Hessian, covariance, and forecasting tests passed."

contains

  function last_score(y, lambda_last, p) result(value)
    real(dp), intent(in) :: y(:), lambda_last
    type(tegarch_parameters), intent(in) :: p
    real(dp) :: value
    real(dp) :: mu, ylast, arg

    mu = skew_t_mean_local(p%df, p%skew)
    ylast = y(size(y))
    arg = ylast + mu*exp(lambda_last)
    value = (p%df + 1.0_dp)*(ylast*ylast + ylast*mu*exp(lambda_last))/ &
      (p%df*exp(2.0_dp*lambda_last)*p%skew**(2*sign_integer(arg)) + arg*arg) - 1.0_dp
  end function last_score

  pure function skew_t_mean_local(df, skew) result(value)
    real(dp), intent(in) :: df, skew
    real(dp) :: value
    real(dp), parameter :: pi_local = acos(-1.0_dp)

    value = (skew - 1.0_dp/skew)*sqrt(df)*exp(log_gamma(0.5_dp*(df - 1.0_dp)) + &
      log_gamma(0.5_dp) - log_gamma(0.5_dp*df))/pi_local
  end function skew_t_mean_local

  pure function sign_integer(x) result(value)
    real(dp), intent(in) :: x
    integer :: value

    if (x > 0.0_dp) then
      value = 1
    else if (x < 0.0_dp) then
      value = -1
    else
      value = 0
    end if
  end function sign_integer

  pure function sign_real(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = real(sign_integer(x), dp)
  end function sign_real

  subroutine assert_fit_improved(fit, initial_value, label)
    type(tegarch_fit_result), intent(in) :: fit
    real(dp), intent(in) :: initial_value
    character(len=*), intent(in) :: label

    if (.not. ieee_is_finite(fit%log_likelihood)) error stop trim(label)//": non-finite log-likelihood"
    if (fit%log_likelihood + 1.0e-8_dp < initial_value) then
      print '(a,1x,es16.8)', trim(label)//" initial:", initial_value
      print '(a,1x,es16.8)', trim(label)//" fitted:", fit%log_likelihood
      error stop trim(label)//": optimizer did not improve"
    end if
    if (fit%parameters%df <= 2.0_dp .or. fit%parameters%skew <= 0.0_dp) then
      error stop trim(label)//": invalid fitted distribution parameters"
    end if
    if (abs(fit%parameters%phi1) >= 1.0_dp) error stop trim(label)//": invalid phi1"
    if (fit%parameters%components == 2 .and. abs(fit%parameters%phi2) >= 1.0_dp) then
      error stop trim(label)//": invalid phi2"
    end if
  end subroutine assert_fit_improved

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      print '(a,1x,es24.16)', trim(label)//" actual:", actual
      print '(a,1x,es24.16)', trim(label)//" expected:", expected
      error stop "assert_close failed"
    end if
  end subroutine assert_close

end program test_models
