! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module tegarch_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use betategarch_kinds, only : dp
  use betategarch_math, only : pi, signum, invert_matrix, numerical_hessian
  use betategarch_rng, only : set_random_seed
  use skew_t_mod, only : skew_t_random, skew_t_mean, skew_t_variance
  use bounded_nelder_mead_mod, only : bounded_nelder_mead
  implicit none
  private

  type, public :: tegarch_parameters
    integer :: components = 1
    logical :: asym = .true.
    logical :: skewed = .true.
    real(dp) :: omega = 0.02_dp
    real(dp) :: phi1 = 0.95_dp
    real(dp) :: phi2 = 0.90_dp
    real(dp) :: kappa1 = 0.05_dp
    real(dp) :: kappa2 = 0.01_dp
    real(dp) :: kappastar = 0.01_dp
    real(dp) :: df = 10.0_dp
    real(dp) :: skew = 0.98_dp
  end type tegarch_parameters

  type, public :: tegarch_filter_result
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: stdev(:)
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: lambda1_dagger(:)
    real(dp), allocatable :: lambda2_dagger(:)
    real(dp), allocatable :: score(:)
    real(dp), allocatable :: epsilon(:)
    real(dp), allocatable :: residual_std(:)
  end type tegarch_filter_result

  type, public :: tegarch_fit_result
    type(tegarch_parameters) :: parameters
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: convergence = 1
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: hessian_available = .false.
    logical :: covariance_available = .false.
    real(dp), allocatable :: free_parameters(:)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: covariance(:, :)
  end type tegarch_fit_result

  type :: tegarch_fit_context
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: lambda_initial(:)
    integer :: components = 1
    logical :: asym = .true.
    logical :: skewed = .true.
    logical :: has_lambda_initial = .false.
  end type tegarch_fit_context

  public :: tegarch_filter, tegarch_simulate, tegarch_loglik
  public :: tegarch_fit, tegarch_forecast
  public :: tegarch_bic_per_observation, tegarch_standard_errors
  public :: tegarch_default_parameters, tegarch_default_bounds
  public :: tegarch_free_parameter_count, tegarch_params_to_free, tegarch_free_to_params

contains

  pure function tegarch_free_parameter_count(components, asym, skewed) result(n)
    integer, intent(in) :: components
    logical, intent(in) :: asym, skewed
    integer :: n

    if (components == 1) then
      n = 4
      if (asym) n = n + 1
      if (skewed) n = n + 1
    else if (components == 2) then
      n = 7
      if (skewed) n = n + 1
    else
      n = 0
    end if
  end function tegarch_free_parameter_count

  pure function tegarch_default_parameters(components, asym, skewed) result(params)
    integer, intent(in) :: components
    logical, intent(in) :: asym, skewed
    type(tegarch_parameters) :: params

    params%components = components
    params%asym = asym
    params%skewed = skewed
    if (components == 1) then
      params%omega = 0.02_dp
      params%phi1 = 0.95_dp
      params%kappa1 = 0.05_dp
      params%kappastar = merge(0.01_dp, 0.0_dp, asym)
      params%df = 10.0_dp
      params%skew = merge(0.98_dp, 1.0_dp, skewed)
      params%phi2 = 0.0_dp
      params%kappa2 = 0.0_dp
    else
      params%omega = 0.02_dp
      params%phi1 = 0.95_dp
      params%phi2 = 0.90_dp
      params%kappa1 = 0.001_dp
      params%kappa2 = 0.01_dp
      params%kappastar = 0.005_dp
      params%df = 10.0_dp
      params%skew = merge(0.98_dp, 1.0_dp, skewed)
    end if
  end function tegarch_default_parameters

  subroutine tegarch_default_bounds(components, asym, skewed, lower, upper)
    integer, intent(in) :: components
    logical, intent(in) :: asym, skewed
    real(dp), allocatable, intent(out) :: lower(:), upper(:)

    integer :: n, j
    real(dp) :: inf, eps

    n = tegarch_free_parameter_count(components, asym, skewed)
    allocate(lower(n), upper(n))
    inf = huge(1.0_dp)
    eps = epsilon(1.0_dp)
    lower = -inf
    upper = inf

    if (components == 1) then
      lower(2) = -1.0_dp + eps
      upper(2) = 1.0_dp - eps
      j = 4
      if (asym) j = 5
      lower(j) = 2.0_dp + eps
      if (skewed) lower(j+1) = eps
    else
      lower(2:3) = -1.0_dp + eps
      upper(2:3) = 1.0_dp - eps
      lower(7) = 2.0_dp + eps
      if (skewed) lower(8) = eps
    end if
  end subroutine tegarch_default_bounds

  pure subroutine tegarch_params_to_free(params, free)
    type(tegarch_parameters), intent(in) :: params
    real(dp), intent(out) :: free(:)

    integer :: j

    j = 0
    if (params%components == 1) then
      j = j + 1; free(j) = params%omega
      j = j + 1; free(j) = params%phi1
      j = j + 1; free(j) = params%kappa1
      if (params%asym) then
        j = j + 1; free(j) = params%kappastar
      end if
      j = j + 1; free(j) = params%df
      if (params%skewed) then
        j = j + 1; free(j) = params%skew
      end if
    else
      j = j + 1; free(j) = params%omega
      j = j + 1; free(j) = params%phi1
      j = j + 1; free(j) = params%phi2
      j = j + 1; free(j) = params%kappa1
      j = j + 1; free(j) = params%kappa2
      j = j + 1; free(j) = params%kappastar
      j = j + 1; free(j) = params%df
      if (params%skewed) then
        j = j + 1; free(j) = params%skew
      end if
    end if
  end subroutine tegarch_params_to_free

  pure subroutine tegarch_free_to_params(free, components, asym, skewed, params)
    real(dp), intent(in) :: free(:)
    integer, intent(in) :: components
    logical, intent(in) :: asym, skewed
    type(tegarch_parameters), intent(out) :: params

    integer :: j

    params = tegarch_default_parameters(components, asym, skewed)
    j = 0
    if (components == 1) then
      j = j + 1; params%omega = free(j)
      j = j + 1; params%phi1 = free(j)
      j = j + 1; params%kappa1 = free(j)
      if (asym) then
        j = j + 1; params%kappastar = free(j)
      else
        params%kappastar = 0.0_dp
      end if
      j = j + 1; params%df = free(j)
      if (skewed) then
        j = j + 1; params%skew = free(j)
      else
        params%skew = 1.0_dp
      end if
      params%phi2 = 0.0_dp
      params%kappa2 = 0.0_dp
    else
      j = j + 1; params%omega = free(j)
      j = j + 1; params%phi1 = free(j)
      j = j + 1; params%phi2 = free(j)
      j = j + 1; params%kappa1 = free(j)
      j = j + 1; params%kappa2 = free(j)
      j = j + 1; params%kappastar = free(j)
      j = j + 1; params%df = free(j)
      if (skewed) then
        j = j + 1; params%skew = free(j)
      else
        params%skew = 1.0_dp
      end if
    end if
  end subroutine tegarch_free_to_params

  pure function valid_parameters(params) result(ok)
    type(tegarch_parameters), intent(in) :: params
    logical :: ok

    ok = params%components == 1 .or. params%components == 2
    ok = ok .and. abs(params%phi1) < 1.0_dp
    if (params%components == 2) ok = ok .and. abs(params%phi2) < 1.0_dp .and. params%asym
    ok = ok .and. params%df > 2.0_dp .and. params%skew > 0.0_dp
  end function valid_parameters

  subroutine allocate_filter(result, n)
    type(tegarch_filter_result), intent(out) :: result
    integer, intent(in) :: n

    allocate(result%y(n), result%sigma(n), result%stdev(n), result%lambda(n), &
      result%lambda1_dagger(n), result%lambda2_dagger(n), result%score(n), &
      result%epsilon(n), result%residual_std(n))
    result%y = 0.0_dp
    result%sigma = 0.0_dp
    result%stdev = 0.0_dp
    result%lambda = 0.0_dp
    result%lambda1_dagger = 0.0_dp
    result%lambda2_dagger = 0.0_dp
    result%score = 0.0_dp
    result%epsilon = 0.0_dp
    result%residual_std = 0.0_dp
  end subroutine allocate_filter

  subroutine tegarch_filter(y, params, result, lambda_initial)
    real(dp), intent(in) :: y(:)
    type(tegarch_parameters), intent(in) :: params
    type(tegarch_filter_result), intent(out) :: result
    real(dp), intent(in), optional :: lambda_initial(:)

    real(dp) :: mu_eps, sd_eps, sign_arg, skew_term, df_plus_one
    integer :: n, i

    if (.not. valid_parameters(params)) error stop "tegarch_filter: invalid parameters"
    n = size(y)
    if (n < 1) error stop "tegarch_filter: empty series"
    call allocate_filter(result, n)
    result%y = y

    mu_eps = skew_t_mean(params%df, params%skew)
    sd_eps = sqrt(skew_t_variance(params%df, params%skew))
    df_plus_one = params%df + 1.0_dp

    if (params%components == 1) then
      result%lambda(1) = params%omega
      result%lambda1_dagger(1) = 0.0_dp
      if (present(lambda_initial)) then
        if (size(lambda_initial) >= 1) result%lambda(1) = lambda_initial(1)
        if (size(lambda_initial) >= 2) result%lambda1_dagger(1) = lambda_initial(2)
      end if

      do i = 1, n - 1
        sign_arg = y(i) + mu_eps*exp(result%lambda(i))
        skew_term = params%skew**(2*signum(sign_arg))
        result%score(i) = df_plus_one*sign_arg*y(i)/(params%df*exp(2.0_dp*result%lambda(i))* &
          skew_term + sign_arg*sign_arg) - 1.0_dp
        result%lambda1_dagger(i+1) = params%phi1*result%lambda1_dagger(i) + &
          params%kappa1*result%score(i) + params%kappastar*real(signum(-y(i)), dp)*(result%score(i) + 1.0_dp)
        result%lambda(i+1) = params%omega + result%lambda1_dagger(i+1)
      end do
    else
      result%lambda(1) = params%omega
      result%lambda1_dagger(1) = 0.0_dp
      result%lambda2_dagger(1) = 0.0_dp
      if (present(lambda_initial)) then
        if (size(lambda_initial) >= 1) result%lambda(1) = lambda_initial(1)
        if (size(lambda_initial) >= 2) result%lambda1_dagger(1) = lambda_initial(2)
        if (size(lambda_initial) >= 3) result%lambda2_dagger(1) = lambda_initial(3)
      end if

      do i = 1, n - 1
        sign_arg = y(i) + mu_eps*exp(result%lambda(i))
        skew_term = params%skew**(2*signum(sign_arg))
        result%score(i) = df_plus_one*sign_arg*y(i)/(params%df*exp(2.0_dp*result%lambda(i))* &
          skew_term + sign_arg*sign_arg) - 1.0_dp
        result%lambda1_dagger(i+1) = params%phi1*result%lambda1_dagger(i) + params%kappa1*result%score(i)
        result%lambda2_dagger(i+1) = params%phi2*result%lambda2_dagger(i) + &
          params%kappa2*result%score(i) + params%kappastar*real(signum(-y(i)), dp)*(result%score(i) + 1.0_dp)
        result%lambda(i+1) = params%omega + result%lambda1_dagger(i+1) + result%lambda2_dagger(i+1)
      end do
    end if

    result%sigma = exp(result%lambda)
    result%stdev = result%sigma*sd_eps
    result%epsilon = y/result%sigma
    result%residual_std = result%epsilon/sd_eps
  end subroutine tegarch_filter

  function tegarch_loglik(y, params, lambda_initial, penalty) result(log_likelihood)
    real(dp), intent(in) :: y(:)
    type(tegarch_parameters), intent(in) :: params
    real(dp), intent(in), optional :: lambda_initial(:)
    real(dp), intent(in), optional :: penalty
    real(dp) :: log_likelihood

    type(tegarch_filter_result) :: filtered
    real(dp) :: term1, term2, term3, yterm, denom, penalty_value, mu_eps
    integer :: i, n

    penalty_value = -1.0e100_dp
    if (present(penalty)) penalty_value = penalty
    if (.not. valid_parameters(params) .or. size(y) < 1) then
      log_likelihood = penalty_value
      return
    end if

    call tegarch_filter(y, params, filtered, lambda_initial)
    n = size(y)
    term1 = real(n, dp)*(log(2.0_dp) - log(params%skew + 1.0_dp/params%skew) + &
      log_gamma(0.5_dp*(params%df + 1.0_dp)) - log_gamma(0.5_dp*params%df) - &
      0.5_dp*log(pi*params%df))
    term2 = sum(filtered%lambda)
    mu_eps = skew_t_mean(params%df, params%skew)
    term3 = 0.0_dp
    do i = 1, n
      yterm = y(i) + mu_eps*exp(filtered%lambda(i))
      denom = params%skew**(2*signum(yterm))*params%df*exp(2.0_dp*filtered%lambda(i))
      term3 = term3 + 0.5_dp*(params%df + 1.0_dp)*log(1.0_dp + yterm*yterm/denom)
    end do
    log_likelihood = term1 - term2 - term3
    if (.not. ieee_is_finite(log_likelihood)) log_likelihood = penalty_value
  end function tegarch_loglik

  subroutine tegarch_simulate(n, params, result, seed, lambda_initial)
    integer, intent(in) :: n
    type(tegarch_parameters), intent(in) :: params
    type(tegarch_filter_result), intent(out) :: result
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: lambda_initial(:)

    real(dp), allocatable :: innovation(:)
    real(dp) :: mu_eps, sd_eps, eps2, score_value, leverage_value
    integer :: i

    if (.not. valid_parameters(params)) error stop "tegarch_simulate: invalid parameters"
    if (n < 1) error stop "tegarch_simulate: n must be positive"
    if (present(seed)) call set_random_seed(seed)
    call allocate_filter(result, n)
    allocate(innovation(n))
    call skew_t_random(innovation, params%df, params%skew)
    mu_eps = skew_t_mean(params%df, params%skew)
    sd_eps = sqrt(skew_t_variance(params%df, params%skew))

    if (params%components == 1) then
      result%lambda1_dagger(1) = 0.0_dp
      result%lambda(1) = params%omega
      if (present(lambda_initial)) then
        if (size(lambda_initial) >= 1) result%lambda(1) = lambda_initial(1)
        if (size(lambda_initial) >= 2) result%lambda1_dagger(1) = lambda_initial(2)
      end if
      do i = 1, n - 1
        eps2 = innovation(i)*innovation(i)
        score_value = (params%df + 1.0_dp)*(eps2 - mu_eps*innovation(i))/ &
          (params%df*params%skew**(2*signum(innovation(i))) + eps2) - 1.0_dp
        leverage_value = real(signum(mu_eps - innovation(i)), dp)*(score_value + 1.0_dp)
        result%score(i) = score_value
        result%lambda1_dagger(i+1) = params%phi1*result%lambda1_dagger(i) + &
          params%kappa1*score_value + params%kappastar*leverage_value
        result%lambda(i+1) = params%omega + result%lambda1_dagger(i+1)
      end do
    else
      result%lambda1_dagger(1) = 0.0_dp
      result%lambda2_dagger(1) = 0.0_dp
      result%lambda(1) = params%omega
      if (present(lambda_initial)) then
        if (size(lambda_initial) >= 1) result%lambda(1) = lambda_initial(1)
        if (size(lambda_initial) >= 2) result%lambda1_dagger(1) = lambda_initial(2)
        if (size(lambda_initial) >= 3) result%lambda2_dagger(1) = lambda_initial(3)
      end if
      do i = 1, n - 1
        eps2 = innovation(i)*innovation(i)
        score_value = (params%df + 1.0_dp)*(eps2 - mu_eps*innovation(i))/ &
          (params%df*params%skew**(2*signum(innovation(i))) + eps2) - 1.0_dp
        leverage_value = real(signum(mu_eps - innovation(i)), dp)*(score_value + 1.0_dp)
        result%score(i) = score_value
        result%lambda1_dagger(i+1) = params%phi1*result%lambda1_dagger(i) + params%kappa1*score_value
        result%lambda2_dagger(i+1) = params%phi2*result%lambda2_dagger(i) + &
          params%kappa2*score_value + params%kappastar*leverage_value
        result%lambda(i+1) = params%omega + result%lambda1_dagger(i+1) + result%lambda2_dagger(i+1)
      end do
    end if

    result%sigma = exp(result%lambda)
    result%stdev = result%sigma*sd_eps
    result%epsilon = innovation - mu_eps
    result%y = result%sigma*result%epsilon
    result%residual_std = result%epsilon/sd_eps
  end subroutine tegarch_simulate

  subroutine tegarch_fit(y, components, asym, skewed, fit, initial, lower, upper, &
      lambda_initial, compute_hessian, max_iterations, tolerance)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: components
    logical, intent(in) :: asym, skewed
    type(tegarch_fit_result), intent(out) :: fit
    real(dp), intent(in), optional :: initial(:), lower(:), upper(:), lambda_initial(:)
    logical, intent(in), optional :: compute_hessian
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance

    type(tegarch_parameters) :: start_params
    type(tegarch_fit_context) :: context
    real(dp), allocatable :: start(:), lo(:), hi(:), solution(:), neg_hessian(:, :)
    real(dp) :: objective_value, tol
    integer :: npar, max_iter, hinfo, invinfo
    logical :: want_hessian

    if (components /= 1 .and. components /= 2) error stop "tegarch_fit: components must be 1 or 2"
    if (components == 2 .and. .not. asym) error stop "tegarch_fit: asymmetry is required for two components"
    if (size(y) < 2) error stop "tegarch_fit: at least two observations are required"

    npar = tegarch_free_parameter_count(components, asym, skewed)
    allocate(start(npar), solution(npar))
    start_params = tegarch_default_parameters(components, asym, skewed)
    call tegarch_params_to_free(start_params, start)
    if (present(initial)) then
      if (size(initial) /= npar) error stop "tegarch_fit: initial has wrong size"
      start = initial
    end if

    call tegarch_default_bounds(components, asym, skewed, lo, hi)
    if (present(lower)) then
      if (size(lower) /= npar) error stop "tegarch_fit: lower has wrong size"
      lo = lower
    end if
    if (present(upper)) then
      if (size(upper) /= npar) error stop "tegarch_fit: upper has wrong size"
      hi = upper
    end if

    max_iter = 2000
    if (present(max_iterations)) max_iter = max_iterations
    tol = 1.0e-7_dp
    if (present(tolerance)) tol = tolerance
    context%components = components
    context%asym = asym
    context%skewed = skewed
    context%y = y
    if (present(lambda_initial)) then
      context%lambda_initial = lambda_initial
      context%has_lambda_initial = .true.
    end if
    call bounded_nelder_mead(negative_loglik_context, context, start, lo, hi, solution, objective_value, &
      fit%convergence, fit%iterations, fit%evaluations, max_iterations=max_iter, tolerance=tol)

    allocate(fit%free_parameters(npar))
    fit%free_parameters = solution
    call tegarch_free_to_params(solution, components, asym, skewed, fit%parameters)
    fit%log_likelihood = -objective_value

    want_hessian = .true.
    if (present(compute_hessian)) want_hessian = compute_hessian
    if (want_hessian) then
      allocate(fit%hessian(npar, npar), fit%covariance(npar, npar), neg_hessian(npar, npar))
      call numerical_hessian(loglik_context, solution, context, fit%hessian, hinfo)
      fit%hessian_available = hinfo == 0
      if (hinfo == 0) then
        neg_hessian = -fit%hessian
        call invert_matrix(neg_hessian, fit%covariance, invinfo)
        fit%covariance_available = invinfo == 0
        if (invinfo /= 0) fit%covariance = 0.0_dp
      end if
    end if

  end subroutine tegarch_fit

  function negative_loglik_context(free, context) result(value)
    real(dp), intent(in) :: free(:)
    class(*), intent(in) :: context
    real(dp) :: value

    type(tegarch_parameters) :: params

    select type (ctx => context)
    type is (tegarch_fit_context)
      call tegarch_free_to_params(free, ctx%components, ctx%asym, ctx%skewed, params)
      if (ctx%has_lambda_initial) then
        value = -tegarch_loglik(ctx%y, params, lambda_initial=ctx%lambda_initial)
      else
        value = -tegarch_loglik(ctx%y, params)
      end if
    class default
      value = huge(1.0_dp)
    end select
  end function negative_loglik_context

  function loglik_context(free, context) result(value)
    real(dp), intent(in) :: free(:)
    class(*), intent(in) :: context
    real(dp) :: value

    value = -negative_loglik_context(free, context)
  end function loglik_context

  subroutine tegarch_forecast(y, params, n_ahead, sigma, stdev, n_sim, seed, lambda_initial)
    real(dp), intent(in) :: y(:)
    type(tegarch_parameters), intent(in) :: params
    integer, intent(in) :: n_ahead
    real(dp), intent(out) :: sigma(:), stdev(:)
    integer, intent(in), optional :: n_sim, seed
    real(dp), intent(in), optional :: lambda_initial(:)

    type(tegarch_filter_result) :: filtered
    real(dp), allocatable :: innovation(:), g1(:), g2(:)
    real(dp) :: mu_eps, sd_eps, y_last, lambda_last, dagger1_last, dagger2_last
    real(dp) :: u_last, g1_last, g2_last, factor, eps2, exponent_power
    integer :: simulations, h, j

    if (n_ahead < 1 .or. size(sigma) /= n_ahead .or. size(stdev) /= n_ahead) then
      error stop "tegarch_forecast: inconsistent horizon"
    end if
    simulations = 10000
    if (present(n_sim)) simulations = n_sim
    if (simulations < 1) error stop "tegarch_forecast: n_sim must be positive"
    if (present(seed)) call set_random_seed(seed)

    if (present(lambda_initial)) then
      call tegarch_filter(y, params, filtered, lambda_initial=lambda_initial)
    else
      call tegarch_filter(y, params, filtered)
    end if
    y_last = filtered%y(size(y))
    lambda_last = filtered%lambda(size(y))
    dagger1_last = filtered%lambda1_dagger(size(y))
    dagger2_last = filtered%lambda2_dagger(size(y))
    mu_eps = skew_t_mean(params%df, params%skew)
    sd_eps = sqrt(skew_t_variance(params%df, params%skew))
    u_last = (params%df + 1.0_dp)*(y_last*y_last + y_last*mu_eps*exp(lambda_last))/ &
      (params%df*exp(2.0_dp*lambda_last)*params%skew**(2*signum(y_last + mu_eps*exp(lambda_last))) + &
      (y_last + mu_eps*exp(lambda_last))**2) - 1.0_dp

    if (params%components == 1) then
      g1_last = params%kappa1*u_last + params%kappastar*real(signum(-y_last), dp)*(u_last + 1.0_dp)
      sigma(1) = exp(params%omega + params%phi1*dagger1_last + g1_last)
      if (n_ahead > 1) then
        allocate(innovation(simulations), g1(simulations))
        call skew_t_random(innovation, params%df, params%skew)
        do j = 1, simulations
          eps2 = innovation(j)*innovation(j)
          g1(j) = params%kappa1*((params%df + 1.0_dp)*(eps2 - mu_eps*innovation(j))/ &
            (params%df*params%skew**(2*signum(innovation(j))) + eps2) - 1.0_dp)
          g1(j) = g1(j) + params%kappastar*real(signum(mu_eps - innovation(j)), dp)* &
            (((params%df + 1.0_dp)*(eps2 - mu_eps*innovation(j))/ &
            (params%df*params%skew**(2*signum(innovation(j))) + eps2) - 1.0_dp) + 1.0_dp)
        end do
        do h = 2, n_ahead
          factor = exp(params%phi1**(h-1)*g1_last)
          do j = 1, h - 1
            exponent_power = params%phi1**(h-1-j)
            factor = factor*sum(exp(exponent_power*g1))/real(simulations, dp)
          end do
          sigma(h) = exp(params%omega + params%phi1**h*dagger1_last)*factor
        end do
      end if
    else
      g1_last = params%kappa1*u_last
      g2_last = params%kappa2*u_last + params%kappastar*real(signum(-y_last), dp)*(u_last + 1.0_dp)
      sigma(1) = exp(params%omega + params%phi1*dagger1_last + params%phi2*dagger2_last + g1_last + g2_last)
      if (n_ahead > 1) then
        allocate(innovation(simulations), g1(simulations), g2(simulations))
        call skew_t_random(innovation, params%df, params%skew)
        do j = 1, simulations
          eps2 = innovation(j)*innovation(j)
          u_last = (params%df + 1.0_dp)*(eps2 - mu_eps*innovation(j))/ &
            (params%df*params%skew**(2*signum(innovation(j))) + eps2) - 1.0_dp
          g1(j) = params%kappa1*u_last
          g2(j) = params%kappa2*u_last + params%kappastar*real(signum(mu_eps - innovation(j)), dp)*(u_last + 1.0_dp)
        end do
        do h = 2, n_ahead
          factor = exp(params%phi1**(h-1)*g1_last + params%phi2**(h-1)*g2_last)
          do j = 1, h - 1
            factor = factor*sum(exp(params%phi1**(h-1-j)*g1 + params%phi2**(h-1-j)*g2))/ &
              real(simulations, dp)
          end do
          sigma(h) = exp(params%omega + params%phi1**h*dagger1_last + params%phi2**h*dagger2_last)*factor
        end do
      end if
    end if
    stdev = sigma*sd_eps
  end subroutine tegarch_forecast

  pure function tegarch_bic_per_observation(log_likelihood, nobs, npar) result(value)
    real(dp), intent(in) :: log_likelihood
    integer, intent(in) :: nobs, npar
    real(dp) :: value

    if (nobs < 1 .or. npar < 0) then
      value = huge(1.0_dp)
    else
      value = (-2.0_dp*log_likelihood + real(npar, dp)*log(real(nobs, dp)))/real(nobs, dp)
    end if
  end function tegarch_bic_per_observation

  subroutine tegarch_standard_errors(fit, standard_errors, info)
    type(tegarch_fit_result), intent(in) :: fit
    real(dp), intent(out) :: standard_errors(:)
    integer, intent(out) :: info

    integer :: i, n

    info = 0
    n = size(standard_errors)
    if (.not. fit%covariance_available .or. .not. allocated(fit%covariance)) then
      info = 1
      standard_errors = 0.0_dp
      return
    end if
    if (size(fit%covariance, 1) /= n .or. size(fit%covariance, 2) /= n) then
      info = -1
      standard_errors = 0.0_dp
      return
    end if
    do i = 1, n
      if (fit%covariance(i, i) >= 0.0_dp) then
        standard_errors(i) = sqrt(fit%covariance(i, i))
      else
        standard_errors(i) = 0.0_dp
        info = 2
      end if
    end do
  end subroutine tegarch_standard_errors

end module tegarch_mod
