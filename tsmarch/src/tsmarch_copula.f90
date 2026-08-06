! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_copula
  use ghyp_kinds, only : dp
  use tsd_types, only : distribution_parameters
  use tsd_distributions, only : pdist, qdist, dstd
  use tsd_optimize, only : nelder_mead
  use tsgarch, only : garch_spec, garch_fit, fit_options
  use tsmarch_types
  use tsmarch_dcc, only : dcc_filter, dcc_constant_filter, fit_marginal_garch, forecast_dcc
  implicit none
  private
  public :: probability_transform, copula_transform
  public :: copula_filter, estimate_copula, forecast_copula

  type :: copula_objective_context
    real(dp), allocatable :: uniforms(:, :)
    type(copula_spec) :: spec
  end type copula_objective_context

contains

  function probability_transform(fits) result(u)
    type(garch_fit), intent(in) :: fits(:)
    real(dp), allocatable :: u(:, :)
    integer :: n, m, i, j
    type(distribution_parameters) :: p
    m = size(fits)
    if (m == 0) then
      allocate(u(0, 0))
      return
    end if
    n = fits(1)%filtered%nobs
    allocate(u(n, m))
    do j = 1, m
      p = fits(j)%parameters%dist
      p%mu = 0.0_dp
      p%sigma = 1.0_dp
      do i = 1, n
        u(i, j) = pdist(fits(j)%spec%distribution, fits(j)%filtered%standardized_residuals(i), p)
        u(i, j) = min(max(u(i, j), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
      end do
    end do
  end function probability_transform

  function copula_transform(u, distribution, shape_parameter) result(z)
    real(dp), intent(in) :: u(:, :)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in), optional :: shape_parameter
    real(dp), allocatable :: z(:, :)
    type(distribution_parameters) :: p
    real(dp) :: shapev
    integer :: i, j
    shapev = 8.0_dp
    if (present(shape_parameter)) shapev = shape_parameter
    allocate(z(size(u, 1), size(u, 2)))
    p%mu = 0.0_dp
    p%sigma = 1.0_dp
    p%shape = shapev
    do j = 1, size(u, 2)
      do i = 1, size(u, 1)
        if (trim(lower_string(distribution)) == 'student' .or. &
            trim(lower_string(distribution)) == 'mvt') then
          z(i, j) = qdist('std', min(max(u(i, j), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp), p)
        else
          z(i, j) = qdist('norm', min(max(u(i, j), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp), p)
        end if
      end do
    end do
  end function copula_transform

  function copula_filter(u, specification, parameters) result(out)
    real(dp), intent(in) :: u(:, :)
    type(copula_spec), intent(in) :: specification
    type(dcc_parameters), intent(in) :: parameters
    type(dcc_filter_result) :: out
    type(dcc_spec) :: ds
    real(dp), allocatable :: z(:, :), sigma(:, :)
    real(dp) :: correction
    integer :: i, j, n, m

    n = size(u, 1)
    m = size(u, 2)
    if (n < 2 .or. m < 2 .or. any(u <= 0.0_dp) .or. any(u >= 1.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'copula uniforms must lie strictly between zero and one'
      return
    end if
    z = copula_transform(u, specification%distribution, parameters%shape)
    allocate(sigma(n, m))
    sigma = 1.0_dp
    ds%distribution = merge('mvt         ', 'mvn         ', &
      trim(lower_string(specification%distribution)) == 'student' .or. &
      trim(lower_string(specification%distribution)) == 'mvt')
    ds%alpha_order = specification%alpha_order
    ds%gamma_order = specification%gamma_order
    ds%beta_order = specification%beta_order
    ds%constant_correlation = specification%constant_correlation
    if (specification%constant_correlation) then
      out = dcc_constant_filter(z, sigma, ds%distribution, 1, parameters%shape)
    else
      out = dcc_filter(z, sigma, ds, parameters)
    end if
    if (out%status /= tsm_success) return
    if (trim(ds%distribution) == 'mvt') then
      do i = 1, n
        correction = 0.0_dp
        do j = 1, m
          correction = correction + log(max(dstd(z(i, j), 0.0_dp, 1.0_dp, parameters%shape, .false.), tiny(1.0_dp)))
        end do
        out%dcc_loglik_vector(i) = out%dcc_loglik_vector(i) + correction
        out%loglik_vector(i) = out%dcc_loglik_vector(i)
      end do
      out%dcc_log_likelihood = -sum(out%dcc_loglik_vector)
      out%log_likelihood = out%dcc_log_likelihood
    else
      out%loglik_vector = out%dcc_loglik_vector
      out%log_likelihood = out%dcc_log_likelihood
    end if
  end function copula_filter

  function estimate_copula(data, marginal_spec, specification, options) result(fit)
    real(dp), intent(in) :: data(:, :)
    type(garch_spec), intent(in) :: marginal_spec
    type(copula_spec), intent(in) :: specification
    type(fit_options), intent(in), optional :: options
    type(copula_fit) :: fit
    type(fit_options) :: opt
    type(copula_objective_context) :: context
    real(dp), allocatable :: sigma(:, :), z(:, :), x(:)
    type(garch_fit), allocatable :: marginal_fits(:)
    real(dp) :: fval
    integer :: status, iterations, j, ncoef

    if (size(data, 1) < 20 .or. size(data, 2) < 2) then
      fit%status = tsm_invalid_argument
      fit%message = 'copula estimation requires at least 20 observations and two series'
      return
    end if
    opt%compute_inference = .false.
    if (present(options)) opt = options
    allocate(marginal_fits(size(data, 2)), z(size(data, 1), size(data, 2)), sigma(size(data, 1), size(data, 2)))
    call fit_marginal_garch(data, marginal_spec, marginal_fits, z, sigma, opt)
    call move_alloc(marginal_fits, fit%marginals)
    do j = 1, size(fit%marginals)
      if (fit%marginals(j)%status /= 0) then
        fit%status = tsm_numerical_failure
        fit%message = 'one or more marginal GARCH fits failed'
        return
      end if
    end do
    fit%uniforms = probability_transform(fit%marginals)
    fit%spec = specification
    fit%parameters = default_parameters(specification)
    if (.not. specification%constant_correlation) then
      context%uniforms = fit%uniforms
      context%spec = specification
      x = to_unconstrained(fit%parameters, specification)
      call nelder_mead(copula_objective, context, x, fval, status, iterations, &
        max_iterations=opt%max_iterations, tolerance=opt%tolerance, scale=opt%simplex_scale)
      fit%parameters = from_unconstrained(x, specification)
    end if
    fit%filtered = copula_filter(fit%uniforms, specification, fit%parameters)
    fit%log_likelihood = fit%filtered%log_likelihood
    fit%status = fit%filtered%status
    fit%message = fit%filtered%message
    ncoef = specification%alpha_order + specification%gamma_order + specification%beta_order
    if (ncoef < 0) fit%message = fit%message
  end function estimate_copula


  function forecast_copula(data, fit, horizon, paths, seed) result(out)
    use ghyp_kinds, only : i8
    real(dp), intent(in) :: data(:, :)
    type(copula_fit), intent(in) :: fit
    integer, intent(in) :: horizon
    integer, intent(in), optional :: paths
    integer(i8), intent(in), optional :: seed
    type(dcc_forecast) :: out
    type(dcc_fit) :: temporary
    temporary%parameters = fit%parameters
    temporary%marginals = fit%marginals
    temporary%filtered = fit%filtered
    temporary%spec%distribution = merge('mvt         ', 'mvn         ', &
      trim(lower_string(fit%spec%distribution)) == 'student' .or. &
      trim(lower_string(fit%spec%distribution)) == 'mvt')
    temporary%spec%alpha_order = fit%spec%alpha_order
    temporary%spec%gamma_order = fit%spec%gamma_order
    temporary%spec%beta_order = fit%spec%beta_order
    temporary%spec%constant_correlation = fit%spec%constant_correlation
    temporary%status = fit%status
    temporary%message = fit%message
    if (present(paths) .and. present(seed)) then
      out = forecast_dcc(data, temporary, horizon, paths=paths, seed=seed)
    else if (present(paths)) then
      out = forecast_dcc(data, temporary, horizon, paths=paths)
    else if (present(seed)) then
      out = forecast_dcc(data, temporary, horizon, seed=seed)
    else
      out = forecast_dcc(data, temporary, horizon)
    end if
  end function forecast_copula

  function copula_objective(x, generic_context) result(value)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout) :: generic_context
    real(dp) :: value
    type(dcc_parameters) :: p
    type(dcc_filter_result) :: f
    select type (context => generic_context)
    type is (copula_objective_context)
      p = from_unconstrained(x, context%spec)
      f = copula_filter(context%uniforms, context%spec, p)
      if (f%status == tsm_success) then
        value = -f%log_likelihood
      else
        value = huge(1.0_dp) / 100.0_dp
      end if
    class default
      value = huge(1.0_dp) / 100.0_dp
    end select
  end function copula_objective

  function default_parameters(spec) result(p)
    type(copula_spec), intent(in) :: spec
    type(dcc_parameters) :: p
    allocate(p%alpha(spec%alpha_order), p%gamma(spec%gamma_order), p%beta(spec%beta_order))
    if (spec%alpha_order > 0) p%alpha = 0.04_dp / real(spec%alpha_order, dp)
    if (spec%gamma_order > 0) p%gamma = 0.02_dp / real(spec%gamma_order, dp)
    if (spec%beta_order > 0) p%beta = 0.90_dp / real(spec%beta_order, dp)
    p%shape = 8.0_dp
  end function default_parameters

  function to_unconstrained(p, spec) result(x)
    type(dcc_parameters), intent(in) :: p
    type(copula_spec), intent(in) :: spec
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: slack
    integer :: n, pos
    n = spec%alpha_order + spec%gamma_order + spec%beta_order
    allocate(x(n + 1), w(n))
    pos = 0
    if (spec%alpha_order > 0) then
      w(pos + 1:pos + spec%alpha_order) = p%alpha
      pos = pos + spec%alpha_order
    end if
    if (spec%gamma_order > 0) then
      w(pos + 1:pos + spec%gamma_order) = p%gamma
      pos = pos + spec%gamma_order
    end if
    if (spec%beta_order > 0) w(pos + 1:pos + spec%beta_order) = p%beta
    slack = max(1.0e-8_dp, 0.999_dp - sum(w))
    x(1:n) = log(max(w, 1.0e-8_dp) / slack)
    x(n + 1) = log(max(p%shape - 2.01_dp, 1.0e-6_dp))
  end function to_unconstrained

  function from_unconstrained(x, spec) result(p)
    real(dp), intent(in) :: x(:)
    type(copula_spec), intent(in) :: spec
    type(dcc_parameters) :: p
    real(dp), allocatable :: w(:)
    integer :: n, pos
    n = spec%alpha_order + spec%gamma_order + spec%beta_order
    allocate(w(n))
    w = exp(max(min(x(1:n), 50.0_dp), -50.0_dp))
    w = 0.995_dp * w / (1.0_dp + sum(w))
    allocate(p%alpha(spec%alpha_order), p%gamma(spec%gamma_order), p%beta(spec%beta_order))
    pos = 0
    if (spec%alpha_order > 0) then
      p%alpha = w(pos + 1:pos + spec%alpha_order)
      pos = pos + spec%alpha_order
    end if
    if (spec%gamma_order > 0) then
      p%gamma = w(pos + 1:pos + spec%gamma_order)
      pos = pos + spec%gamma_order
    end if
    if (spec%beta_order > 0) p%beta = w(pos + 1:pos + spec%beta_order)
    if (trim(lower_string(spec%distribution)) == 'student' .or. trim(lower_string(spec%distribution)) == 'mvt') then
      p%shape = 2.01_dp + exp(min(x(n + 1), 8.0_dp))
    else
      p%shape = 8.0_dp
    end if
  end function from_unconstrained

  function lower_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code
    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lower_string

end module tsmarch_copula
