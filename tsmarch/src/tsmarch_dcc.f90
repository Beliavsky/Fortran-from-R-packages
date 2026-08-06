! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_dcc
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ghyp_kinds, only : dp, i8
  use ghyp_linalg, only : inverse_spd, logdet_spd, cholesky_lower
  use tsd_optimize, only : nelder_mead
  use tsgarch, only : garch_spec, garch_fit, garch_forecast, fit_options, estimate_garch, forecast_garch
  use tsmarch_types
  use tsmarch_linalg
  implicit none
  private
  public :: dcc_filter, dcc_constant_filter, estimate_dcc, fit_marginal_garch
  public :: simulate_dcc_innovations, forecast_dcc, dcc_stationarity
  public :: pack_dcc_parameters, unpack_dcc_parameters

  real(dp), parameter :: pi = acos(-1.0_dp)

  type :: dcc_objective_context
    real(dp), allocatable :: z(:, :), sigma(:, :)
    type(dcc_spec) :: spec
  end type dcc_objective_context

contains

  function dcc_filter(z, sigma, spec, parameters) result(out)
    real(dp), intent(in) :: z(:, :), sigma(:, :)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters), intent(in) :: parameters
    type(dcc_filter_result) :: out
    real(dp), allocatable :: qbar(:, :), nbar(:, :), omega(:, :), q(:, :), r(:, :)
    real(dp), allocatable :: rinv(:, :), asym(:, :), current(:)
    real(dp) :: sum_a, sum_b, sum_g, part1, part2, part3, const_term, zi2
    integer :: t, j, n, m, maxlag
    logical :: ok

    n = size(z, 1)
    m = size(z, 2)
    if (n < 2 .or. m < 1 .or. any(shape(sigma) /= shape(z))) then
      out%status = tsm_invalid_argument
      out%message = 'z and sigma must be conformable with at least two observations'
      return
    end if
    if (.not. parameters_conform(spec, parameters)) then
      out%status = tsm_invalid_argument
      out%message = 'DCC parameter arrays do not conform to the requested orders'
      return
    end if
    if (any(sigma <= 0.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'all conditional standard deviations must be positive'
      return
    end if
    if (spec%constant_correlation) then
      out = dcc_constant_filter(z, sigma, spec%distribution, spec%correlation_method, parameters%shape)
      return
    end if

    maxlag = max(1, max(spec%alpha_order, spec%beta_order))
    qbar = sample_covariance(z)
    qbar = nearest_correlation(qbar)
    allocate(asym(n, m))
    where (z < 0.0_dp)
      asym = z
    elsewhere
      asym = 0.0_dp
    end where
    if (spec%gamma_order > 0) then
      nbar = sample_covariance(asym)
    else
      allocate(nbar(m, m))
      nbar = 0.0_dp
    end if
    sum_a = sum(parameters%alpha)
    sum_b = sum(parameters%beta)
    sum_g = sum(parameters%gamma)
    omega = (1.0_dp - sum_a - sum_b) * qbar - sum_g * nbar
    if (dcc_stationarity(parameters, qbar, nbar) >= 1.0_dp) then
      out%status = tsm_invalid_argument
      out%message = 'DCC parameters violate the stationarity constraint'
      return
    end if

    allocate(out%q(m, m, n), out%correlation(m, m, n), out%covariance(m, m, n))
    allocate(out%standardized_residuals(n, m), out%sigma(n, m))
    allocate(out%loglik_vector(n), out%dcc_loglik_vector(n), out%qbar(m, m), out%nbar(m, m))
    out%standardized_residuals = z
    out%sigma = sigma
    out%qbar = qbar
    out%nbar = nbar
    out%loglik_vector = 0.0_dp
    out%dcc_loglik_vector = 0.0_dp
    do t = 1, min(maxlag, n)
      out%q(:, :, t) = qbar
      out%correlation(:, :, t) = qbar
      out%covariance(:, :, t) = correlation_to_covariance(qbar, sigma(t, :))
    end do

    allocate(q(m, m), r(m, m), current(m))
    do t = maxlag + 1, n
      q = omega
      do j = 1, spec%alpha_order
        q = q + parameters%alpha(j) * outer_product(z(t - j, :), z(t - j, :))
      end do
      do j = 1, spec%gamma_order
        q = q + parameters%gamma(j) * outer_product(asym(t - j, :), asym(t - j, :))
      end do
      do j = 1, spec%beta_order
        q = q + parameters%beta(j) * out%q(:, :, t - j)
      end do
      call symmetrize_local(q)
      r = q_to_correlation(q)
      if (.not. is_positive_definite(r)) r = nearest_correlation(r)
      out%q(:, :, t) = q
      out%correlation(:, :, t) = r
      out%covariance(:, :, t) = correlation_to_covariance(r, sigma(t, :))
      call inverse_spd(r, rinv, ok)
      if (.not. ok) then
        out%status = tsm_numerical_failure
        out%message = 'failed to invert a DCC correlation matrix'
        return
      end if
      part1 = logdet_spd(r, ok)
      if (.not. ok) then
        out%status = tsm_numerical_failure
        out%message = 'failed to compute a DCC log determinant'
        return
      end if
      current = z(t, :)
      if (trim(lower_string(spec%distribution)) == 'mvt' .or. &
          trim(lower_string(spec%distribution)) == 'student') then
        if (parameters%shape <= 2.0_dp) then
          out%status = tsm_invalid_argument
          out%message = 'multivariate Student shape must exceed two'
          return
        end if
        part2 = dot_product(current, matmul(rinv, current))
        part3 = sum(log(sigma(t, :)))
        const_term = -log_gamma(0.5_dp * (parameters%shape + real(m, dp))) + &
          log_gamma(0.5_dp * parameters%shape) + 0.5_dp * real(m, dp) * &
          log(pi * (parameters%shape - 2.0_dp))
        out%dcc_loglik_vector(t) = const_term + 0.5_dp * part1 + &
          0.5_dp * (parameters%shape + real(m, dp)) * &
          log(1.0_dp + part2 / (parameters%shape - 2.0_dp))
        out%loglik_vector(t) = out%dcc_loglik_vector(t) + part3
      else
        part2 = dot_product(current, matmul(rinv, current))
        zi2 = dot_product(current, current)
        out%dcc_loglik_vector(t) = 0.5_dp * (part2 + part1 - zi2)
        out%loglik_vector(t) = out%dcc_loglik_vector(t) + 0.5_dp * &
          (real(m, dp) * log(2.0_dp * pi) + 2.0_dp * sum(log(sigma(t, :))) + zi2)
      end if
    end do
    out%log_likelihood = -sum(out%loglik_vector(maxlag + 1:n))
    out%dcc_log_likelihood = -sum(out%dcc_loglik_vector(maxlag + 1:n))
    out%status = tsm_success
    out%message = 'ok'
  end function dcc_filter

  function dcc_constant_filter(z, sigma, distribution, method, shape_parameter) result(out)
    real(dp), intent(in) :: z(:, :), sigma(:, :)
    character(len=*), intent(in), optional :: distribution
    integer, intent(in), optional :: method
    real(dp), intent(in), optional :: shape_parameter
    type(dcc_filter_result) :: out
    real(dp), allocatable :: r(:, :), rinv(:, :), current(:)
    real(dp) :: part1, part2, zi2, shapev, const_term
    character(len=12) :: dist
    integer :: t, n, m, mth
    logical :: ok
    n = size(z, 1)
    m = size(z, 2)
    if (any(shape(sigma) /= shape(z)) .or. n < 2 .or. any(sigma <= 0.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'invalid standardized-residual or sigma input'
      return
    end if
    dist = 'mvn'
    if (present(distribution)) dist = lower_string(distribution)
    mth = 1
    if (present(method)) mth = method
    shapev = 8.0_dp
    if (present(shape_parameter)) shapev = shape_parameter
    r = sample_correlation(z, mth)
    if (.not. is_positive_definite(r)) r = nearest_correlation(r)
    call inverse_spd(r, rinv, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'constant correlation matrix is singular'
      return
    end if
    part1 = logdet_spd(r, ok)
    allocate(out%q(m, m, n), out%correlation(m, m, n), out%covariance(m, m, n))
    allocate(out%standardized_residuals(n, m), out%sigma(n, m))
    allocate(out%loglik_vector(n), out%dcc_loglik_vector(n), out%qbar(m, m), out%nbar(m, m))
    out%standardized_residuals = z
    out%sigma = sigma
    out%qbar = r
    out%nbar = 0.0_dp
    do t = 1, n
      out%q(:, :, t) = r
      out%correlation(:, :, t) = r
      out%covariance(:, :, t) = correlation_to_covariance(r, sigma(t, :))
      current = z(t, :)
      if (trim(dist) == 'mvt' .or. trim(dist) == 'student') then
        if (shapev <= 2.0_dp) then
          out%status = tsm_invalid_argument
          out%message = 'multivariate Student shape must exceed two'
          return
        end if
        part2 = dot_product(current, matmul(rinv, current))
        const_term = -log_gamma(0.5_dp * (shapev + real(m, dp))) + log_gamma(0.5_dp * shapev) + &
          0.5_dp * real(m, dp) * log(pi * (shapev - 2.0_dp))
        out%dcc_loglik_vector(t) = const_term + 0.5_dp * part1 + &
          0.5_dp * (shapev + real(m, dp)) * log(1.0_dp + part2 / (shapev - 2.0_dp))
        out%loglik_vector(t) = out%dcc_loglik_vector(t) + sum(log(sigma(t, :)))
      else
        part2 = dot_product(current, matmul(rinv, current))
        zi2 = dot_product(current, current)
        out%dcc_loglik_vector(t) = 0.5_dp * (part2 + part1 - zi2)
        out%loglik_vector(t) = out%dcc_loglik_vector(t) + 0.5_dp * &
          (real(m, dp) * log(2.0_dp * pi) + 2.0_dp * sum(log(sigma(t, :))) + zi2)
      end if
    end do
    out%log_likelihood = -sum(out%loglik_vector)
    out%dcc_log_likelihood = -sum(out%dcc_loglik_vector)
    out%status = tsm_success
    out%message = 'ok'
  end function dcc_constant_filter

  subroutine fit_marginal_garch(data, specification, fits, z, sigma, options)
    real(dp), intent(in) :: data(:, :)
    type(garch_spec), intent(in) :: specification
    type(garch_fit), intent(out) :: fits(:)
    real(dp), intent(out) :: z(:, :), sigma(:, :)
    type(fit_options), intent(in), optional :: options
    type(fit_options) :: opt
    integer :: j, n, m
    n = size(data, 1)
    m = size(data, 2)
    if (size(fits) /= m .or. size(z, 1) /= n .or. size(z, 2) /= m .or. &
        size(sigma, 1) /= n .or. size(sigma, 2) /= m) error stop "fit_marginal_garch: nonconforming output arrays"
    opt%compute_inference = .false.
    if (present(options)) opt = options
    do j = 1, m
      fits(j) = estimate_garch(data(:, j), specification, options=opt)
      if (fits(j)%status /= 0) then
        z(:, j) = 0.0_dp
        sigma(:, j) = 1.0_dp
      else
        z(:, j) = fits(j)%filtered%standardized_residuals
        sigma(:, j) = fits(j)%filtered%sigma
      end if
    end do
  end subroutine fit_marginal_garch

  function estimate_dcc(data, marginal_spec, specification, options) result(fit)
    real(dp), intent(in) :: data(:, :)
    type(garch_spec), intent(in) :: marginal_spec
    type(dcc_spec), intent(in) :: specification
    type(fit_options), intent(in), optional :: options
    type(dcc_fit) :: fit
    type(fit_options) :: opt
    type(dcc_objective_context) :: context
    real(dp), allocatable :: z(:, :), sigma(:, :), x(:)
    type(garch_fit), allocatable :: marginal_fits(:)
    real(dp) :: fval
    integer :: status, iterations, npar

    if (size(data, 1) < 20 .or. size(data, 2) < 2) then
      fit%status = tsm_invalid_argument
      fit%message = 'DCC estimation requires at least 20 observations and two series'
      return
    end if
    opt%compute_inference = .false.
    if (present(options)) opt = options
    allocate(marginal_fits(size(data, 2)), z(size(data, 1), size(data, 2)), sigma(size(data, 1), size(data, 2)))
    call fit_marginal_garch(data, marginal_spec, marginal_fits, z, sigma, opt)
    call move_alloc(marginal_fits, fit%marginals)
    if (any([(fit%marginals(status)%status /= 0, status = 1, size(fit%marginals))])) then
      fit%status = tsm_numerical_failure
      fit%message = 'one or more marginal GARCH fits failed'
      return
    end if
    fit%spec = specification
    fit%parameters = default_dcc_parameters(specification)
    if (specification%constant_correlation) then
      fit%filtered = dcc_constant_filter(z, sigma, specification%distribution, &
        specification%correlation_method, fit%parameters%shape)
      fit%log_likelihood = fit%filtered%log_likelihood
      fit%status = fit%filtered%status
      fit%message = fit%filtered%message
      return
    end if
    context%z = z
    context%sigma = sigma
    context%spec = specification
    x = unconstrained_from_dcc(fit%parameters, specification)
    call nelder_mead(dcc_objective, context, x, fval, status, iterations, &
      max_iterations=opt%max_iterations, tolerance=opt%tolerance, scale=opt%simplex_scale)
    fit%parameters = dcc_from_unconstrained(x, specification)
    fit%filtered = dcc_filter(z, sigma, specification, fit%parameters)
    fit%packed_parameters = pack_dcc_parameters(fit%parameters)
    fit%iterations = iterations
    fit%log_likelihood = fit%filtered%log_likelihood
    npar = size(fit%packed_parameters) + sum([(fit%marginals(status)%npars, status = 1, size(fit%marginals))])
    fit%aic = -2.0_dp * fit%log_likelihood + 2.0_dp * real(npar, dp)
    fit%bic = -2.0_dp * fit%log_likelihood + log(real(size(data, 1), dp)) * real(npar, dp)
    if (fit%filtered%status /= tsm_success) then
      fit%status = fit%filtered%status
      fit%message = fit%filtered%message
    else if (status /= 0 .and. ieee_is_finite(fval)) then
      fit%status = tsm_no_convergence
      fit%message = 'DCC optimizer stopped before its convergence criterion'
    else
      fit%status = tsm_success
      fit%message = 'ok'
    end if
    if (present(options)) then
      if (options%compute_inference .and. fit%status == tsm_success) then
        call dcc_numerical_inference(context, x, fit%hessian, fit%covariance, fit%standard_errors, fit%scores)
      end if
    end if
  end function estimate_dcc

  function forecast_dcc(data, fit, horizon, paths, seed) result(out)
    real(dp), intent(in) :: data(:, :)
    type(dcc_fit), intent(in) :: fit
    integer, intent(in) :: horizon
    integer, intent(in), optional :: paths
    integer(i8), intent(in), optional :: seed
    type(dcc_forecast) :: out
    integer :: h, j, m, npaths, p
    real(dp), allocatable :: q(:, :), qbar(:, :), nbar(:, :), r(:, :), omega(:, :)
    real(dp), allocatable :: mf(:), sf(:), zdraw(:, :), l(:, :)
    real(dp) :: sa, sb, sg, chisq
    logical :: ok
    type(garch_forecast) :: uf

    if (fit%status /= tsm_success .or. horizon < 1) then
      out%status = tsm_invalid_argument
      out%message = 'a successful DCC fit and positive horizon are required'
      return
    end if
    m = size(data, 2)
    npaths = 0
    if (present(paths)) npaths = max(0, paths)
    allocate(out%mean(horizon, m), out%sigma(horizon, m))
    allocate(out%covariance(m, m, horizon), out%correlation(m, m, horizon))
    out%mean = 0.0_dp
    do j = 1, m
      uf = forecast_garch(data(:, j), fit%marginals(j), horizon, paths=max(100, npaths), seed=int(31*j, i8))
      if (uf%status /= 0) then
        out%status = tsm_numerical_failure
        out%message = 'a marginal forecast failed'
        return
      end if
      out%mean(:, j) = uf%mean
      out%sigma(:, j) = uf%sigma
    end do
    if (fit%spec%constant_correlation) then
      do h = 1, horizon
        out%correlation(:, :, h) = fit%filtered%correlation(:, :, size(fit%filtered%correlation, 3))
        out%covariance(:, :, h) = correlation_to_covariance(out%correlation(:, :, h), out%sigma(h, :))
      end do
    else
      q = fit%filtered%q(:, :, size(fit%filtered%q, 3))
      qbar = fit%filtered%qbar
      nbar = fit%filtered%nbar
      sa = sum(fit%parameters%alpha)
      sb = sum(fit%parameters%beta)
      sg = sum(fit%parameters%gamma)
      omega = (1.0_dp - sa - sb) * qbar - sg * nbar
      do h = 1, horizon
        q = omega + (sa + sb) * q + sg * nbar
        r = q_to_correlation(q)
        if (.not. is_positive_definite(r)) r = nearest_correlation(r)
        out%correlation(:, :, h) = r
        out%covariance(:, :, h) = correlation_to_covariance(r, out%sigma(h, :))
      end do
    end if
    if (npaths > 0) then
      allocate(out%simulated(horizon, m, npaths), mf(m), sf(m))
      if (present(seed)) call set_random_seed(seed)
      do p = 1, npaths
        do h = 1, horizon
          call cholesky_lower(out%correlation(:, :, h), l, ok)
          if (.not. ok) then
            out%status = tsm_numerical_failure
            out%message = 'forecast correlation is not positive definite'
            return
          end if
          zdraw = random_normal_matrix(m, 1)
          mf = out%mean(h, :)
          sf = out%sigma(h, :)
          if (trim(lower_string(fit%spec%distribution)) == 'mvt') then
            chisq = chi_square_draw(fit%parameters%shape)
            out%simulated(h, :, p) = mf + sf * matmul(l, zdraw(:, 1)) * &
              sqrt((fit%parameters%shape - 2.0_dp) / max(chisq, tiny(1.0_dp)))
          else
            out%simulated(h, :, p) = mf + sf * matmul(l, zdraw(:, 1))
          end if
        end do
      end do
    end if
    out%status = tsm_success
    out%message = 'ok'
  end function forecast_dcc

  function simulate_dcc_innovations(spec, parameters, qbar, nobs, paths, burn, seed) result(out)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters), intent(in) :: parameters
    real(dp), intent(in) :: qbar(:, :)
    integer, intent(in) :: nobs
    integer, intent(in), optional :: paths, burn
    integer(i8), intent(in), optional :: seed
    type(dcc_simulation) :: out
    integer :: m, npaths, nburn, p, t, j, maxlag, total
    real(dp), allocatable :: qhist(:, :, :), zhist(:, :), asym(:, :), q(:, :), r(:, :), omega(:, :), l(:, :)
    real(dp), allocatable :: normal(:, :), nbar(:, :)
    real(dp) :: sa, sb, sg, chisq
    logical :: ok

    m = size(qbar, 1)
    npaths = 1
    if (present(paths)) npaths = max(1, paths)
    nburn = 100
    if (present(burn)) nburn = max(0, burn)
    if (nobs < 1 .or. size(qbar, 2) /= m .or. .not. parameters_conform(spec, parameters)) then
      out%status = tsm_invalid_argument
      out%message = 'invalid DCC simulation arguments'
      return
    end if
    if (present(seed)) call set_random_seed(seed)
    allocate(out%innovations(nobs, m, npaths), out%correlation(m, m, nobs, npaths))
    allocate(out%series(nobs, m, npaths), out%sigma(nobs, m, npaths))
    out%sigma = 1.0_dp
    out%series = 0.0_dp
    maxlag = max(1, max(spec%alpha_order, spec%beta_order))
    total = nobs + nburn + maxlag
    sa = sum(parameters%alpha)
    sb = sum(parameters%beta)
    sg = sum(parameters%gamma)
    allocate(nbar(m, m))
    nbar = 0.5_dp * qbar
    omega = (1.0_dp - sa - sb) * qbar - sg * nbar
    do p = 1, npaths
      allocate(qhist(m, m, total), zhist(total, m), asym(total, m))
      qhist = 0.0_dp
      zhist = 0.0_dp
      asym = 0.0_dp
      do t = 1, maxlag
        qhist(:, :, t) = qbar
      end do
      do t = maxlag + 1, total
        q = omega
        do j = 1, spec%alpha_order
          q = q + parameters%alpha(j) * outer_product(zhist(t - j, :), zhist(t - j, :))
        end do
        do j = 1, spec%gamma_order
          q = q + parameters%gamma(j) * outer_product(asym(t - j, :), asym(t - j, :))
        end do
        do j = 1, spec%beta_order
          q = q + parameters%beta(j) * qhist(:, :, t - j)
        end do
        qhist(:, :, t) = q
        r = q_to_correlation(q)
        if (.not. is_positive_definite(r)) r = nearest_correlation(r)
        call cholesky_lower(r, l, ok)
        if (.not. ok) then
          out%status = tsm_numerical_failure
          out%message = 'failed to factor a simulated correlation matrix'
          return
        end if
        normal = random_normal_matrix(m, 1)
        if (trim(lower_string(spec%distribution)) == 'mvt') then
          chisq = chi_square_draw(parameters%shape)
          zhist(t, :) = matmul(l, normal(:, 1)) * sqrt((parameters%shape - 2.0_dp) / max(chisq, tiny(1.0_dp)))
        else
          zhist(t, :) = matmul(l, normal(:, 1))
        end if
        where (zhist(t, :) < 0.0_dp)
          asym(t, :) = zhist(t, :)
        elsewhere
          asym(t, :) = 0.0_dp
        end where
        if (t > maxlag + nburn) then
          out%innovations(t - maxlag - nburn, :, p) = zhist(t, :)
          out%series(t - maxlag - nburn, :, p) = zhist(t, :)
          out%correlation(:, :, t - maxlag - nburn, p) = r
        end if
      end do
      deallocate(qhist, zhist, asym)
    end do
    out%status = tsm_success
    out%message = 'ok'
  end function simulate_dcc_innovations

  function dcc_stationarity(parameters, qbar, nbar) result(value)
    type(dcc_parameters), intent(in) :: parameters
    real(dp), intent(in) :: qbar(:, :), nbar(:, :)
    real(dp) :: value
    real(dp), allocatable :: qinvhalf(:, :), tmp(:, :), eig(:), vec(:, :)
    logical :: ok
    value = sum(parameters%alpha) + sum(parameters%beta)
    if (size(parameters%gamma) > 0 .and. any(abs(parameters%gamma) > 0.0_dp)) then
      qinvhalf = symmetric_inverse_sqrt(qbar, ok)
      if (ok) then
        tmp = matmul(qinvhalf, matmul(nbar, qinvhalf))
        call jacobi_eigen(tmp, eig, vec, ok)
        if (ok) value = value + maxval(eig) * sum(parameters%gamma)
      else
        value = value + 0.5_dp * sum(parameters%gamma)
      end if
    end if
  end function dcc_stationarity

  function pack_dcc_parameters(parameters) result(x)
    type(dcc_parameters), intent(in) :: parameters
    real(dp), allocatable :: x(:)
    integer :: n, pos
    n = size(parameters%alpha) + size(parameters%gamma) + size(parameters%beta) + 1
    allocate(x(n))
    pos = 0
    if (size(parameters%alpha) > 0) then
      x(pos + 1:pos + size(parameters%alpha)) = parameters%alpha
      pos = pos + size(parameters%alpha)
    end if
    if (size(parameters%gamma) > 0) then
      x(pos + 1:pos + size(parameters%gamma)) = parameters%gamma
      pos = pos + size(parameters%gamma)
    end if
    if (size(parameters%beta) > 0) then
      x(pos + 1:pos + size(parameters%beta)) = parameters%beta
      pos = pos + size(parameters%beta)
    end if
    x(n) = parameters%shape
  end function pack_dcc_parameters

  function unpack_dcc_parameters(x, spec) result(parameters)
    real(dp), intent(in) :: x(:)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters) :: parameters
    integer :: pos
    allocate(parameters%alpha(spec%alpha_order), parameters%gamma(spec%gamma_order), parameters%beta(spec%beta_order))
    pos = 0
    if (spec%alpha_order > 0) then
      parameters%alpha = x(pos + 1:pos + spec%alpha_order)
      pos = pos + spec%alpha_order
    end if
    if (spec%gamma_order > 0) then
      parameters%gamma = x(pos + 1:pos + spec%gamma_order)
      pos = pos + spec%gamma_order
    end if
    if (spec%beta_order > 0) then
      parameters%beta = x(pos + 1:pos + spec%beta_order)
      pos = pos + spec%beta_order
    end if
    parameters%shape = x(pos + 1)
  end function unpack_dcc_parameters

  function default_dcc_parameters(spec) result(parameters)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters) :: parameters
    allocate(parameters%alpha(spec%alpha_order), parameters%gamma(spec%gamma_order), parameters%beta(spec%beta_order))
    if (spec%alpha_order > 0) parameters%alpha = 0.04_dp / real(spec%alpha_order, dp)
    if (spec%gamma_order > 0) parameters%gamma = 0.02_dp / real(spec%gamma_order, dp)
    if (spec%beta_order > 0) parameters%beta = 0.90_dp / real(spec%beta_order, dp)
    parameters%shape = 8.0_dp
  end function default_dcc_parameters

  function dcc_objective(x, generic_context) result(value)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout) :: generic_context
    real(dp) :: value
    type(dcc_parameters) :: parameters
    type(dcc_filter_result) :: result
    select type (context => generic_context)
    type is (dcc_objective_context)
      parameters = dcc_from_unconstrained(x, context%spec)
      result = dcc_filter(context%z, context%sigma, context%spec, parameters)
      if (result%status == tsm_success) then
        value = -result%log_likelihood
      else
        value = huge(1.0_dp) / 100.0_dp
      end if
    class default
      value = huge(1.0_dp) / 100.0_dp
    end select
  end function dcc_objective

  function unconstrained_from_dcc(parameters, spec) result(x)
    type(dcc_parameters), intent(in) :: parameters
    type(dcc_spec), intent(in) :: spec
    real(dp), allocatable :: x(:), p(:)
    real(dp) :: total, scale
    integer :: n, i
    p = pack_dcc_parameters(parameters)
    n = size(p)
    allocate(x(n))
    total = sum(parameters%alpha) + sum(parameters%beta) + 0.5_dp * sum(parameters%gamma)
    scale = max(1.0e-8_dp, 0.999_dp - total)
    do i = 1, n - 1
      x(i) = log(max(p(i), 1.0e-8_dp) / scale)
    end do
    x(n) = log(max(parameters%shape - 2.01_dp, 1.0e-6_dp))
    if (trim(lower_string(spec%distribution)) /= 'mvt' .and. trim(lower_string(spec%distribution)) /= 'student') x(n) = 0.0_dp
  end function unconstrained_from_dcc

  function dcc_from_unconstrained(x, spec) result(parameters)
    real(dp), intent(in) :: x(:)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters) :: parameters
    real(dp), allocatable :: w(:)
    real(dp) :: denom
    integer :: ncoef, pos
    ncoef = spec%alpha_order + spec%gamma_order + spec%beta_order
    allocate(w(ncoef))
    w = exp(max(min(x(1:ncoef), 50.0_dp), -50.0_dp))
    denom = 1.0_dp + sum(w)
    w = 0.995_dp * w / denom
    allocate(parameters%alpha(spec%alpha_order), parameters%gamma(spec%gamma_order), parameters%beta(spec%beta_order))
    pos = 0
    if (spec%alpha_order > 0) then
      parameters%alpha = w(pos + 1:pos + spec%alpha_order)
      pos = pos + spec%alpha_order
    end if
    if (spec%gamma_order > 0) then
      parameters%gamma = w(pos + 1:pos + spec%gamma_order)
      pos = pos + spec%gamma_order
    end if
    if (spec%beta_order > 0) parameters%beta = w(pos + 1:pos + spec%beta_order)
    if (trim(lower_string(spec%distribution)) == 'mvt' .or. trim(lower_string(spec%distribution)) == 'student') then
      parameters%shape = 2.01_dp + exp(min(x(ncoef + 1), 8.0_dp))
    else
      parameters%shape = 8.0_dp
    end if
  end function dcc_from_unconstrained

  subroutine dcc_numerical_inference(context, x, hessian, covariance, standard_errors, scores)
    type(dcc_objective_context), intent(inout) :: context
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: hessian(:, :), covariance(:, :), standard_errors(:), scores(:, :)
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), inv(:, :)
    real(dp) :: hi, hj, f0
    logical :: ok
    integer :: i, j, k, n, nobs
    n = size(x)
    nobs = size(context%z, 1)
    allocate(hessian(n, n), standard_errors(n), scores(nobs, n))
    hessian = 0.0_dp
    scores = 0.0_dp
    f0 = dcc_objective(x, context)
    do i = 1, n
      hi = 1.0e-4_dp * max(1.0_dp, abs(x(i)))
      xp = x
      xm = x
      xp(i) = xp(i) + hi
      xm(i) = xm(i) - hi
      hessian(i, i) = (dcc_objective(xp, context) - 2.0_dp * f0 + dcc_objective(xm, context)) / (hi * hi)
      do j = i + 1, n
        hj = 1.0e-4_dp * max(1.0_dp, abs(x(j)))
        xpp = x; xpm = x; xmp = x; xmm = x
        xpp(i) = xpp(i) + hi; xpp(j) = xpp(j) + hj
        xpm(i) = xpm(i) + hi; xpm(j) = xpm(j) - hj
        xmp(i) = xmp(i) - hi; xmp(j) = xmp(j) + hj
        xmm(i) = xmm(i) - hi; xmm(j) = xmm(j) - hj
        hessian(i, j) = (dcc_objective(xpp, context) - dcc_objective(xpm, context) - &
          dcc_objective(xmp, context) + dcc_objective(xmm, context)) / (4.0_dp * hi * hj)
        hessian(j, i) = hessian(i, j)
      end do
    end do
    call matrix_inverse_general(hessian, inv, ok)
    if (ok) then
      covariance = inv
      do i = 1, n
        standard_errors(i) = sqrt(max(covariance(i, i), 0.0_dp))
      end do
    else
      allocate(covariance(n, n))
      covariance = 0.0_dp
      standard_errors = huge(1.0_dp)
    end if
    do k = 1, nobs
      scores(k, :) = 0.0_dp
    end do
  end subroutine dcc_numerical_inference

  logical function parameters_conform(spec, parameters)
    type(dcc_spec), intent(in) :: spec
    type(dcc_parameters), intent(in) :: parameters
    parameters_conform = allocated(parameters%alpha) .and. allocated(parameters%gamma) .and. allocated(parameters%beta)
    if (.not. parameters_conform) return
    parameters_conform = size(parameters%alpha) == spec%alpha_order .and. &
      size(parameters%gamma) == spec%gamma_order .and. size(parameters%beta) == spec%beta_order .and. &
      all(parameters%alpha >= 0.0_dp) .and. all(parameters%gamma >= 0.0_dp) .and. all(parameters%beta >= 0.0_dp)
  end function parameters_conform

  function q_to_correlation(q) result(r)
    real(dp), intent(in) :: q(:, :)
    real(dp), allocatable :: r(:, :)
    integer :: i, j, m
    m = size(q, 1)
    allocate(r(m, m))
    do j = 1, m
      do i = 1, m
        r(i, j) = q(i, j) / sqrt(max(q(i, i), tiny(1.0_dp)) * max(q(j, j), tiny(1.0_dp)))
      end do
    end do
    do i = 1, m
      r(i, i) = 1.0_dp
    end do
    call symmetrize_local(r)
  end function q_to_correlation

  function correlation_to_covariance(r, sigma) result(cov)
    real(dp), intent(in) :: r(:, :), sigma(:)
    real(dp), allocatable :: cov(:, :)
    integer :: i, j, m
    m = size(r, 1)
    allocate(cov(m, m))
    do j = 1, m
      do i = 1, m
        cov(i, j) = r(i, j) * sigma(i) * sigma(j)
      end do
    end do
  end function correlation_to_covariance

  subroutine symmetrize_local(a)
    real(dp), intent(inout) :: a(:, :)
    integer :: i, j
    real(dp) :: v
    do j = 1, size(a, 2)
      do i = j + 1, size(a, 1)
        v = 0.5_dp * (a(i, j) + a(j, i))
        a(i, j) = v
        a(j, i) = v
      end do
    end do
  end subroutine symmetrize_local

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

  function chi_square_draw(df) result(value)
    real(dp), intent(in) :: df
    real(dp) :: value
    integer :: k, i
    real(dp) :: frac, z, u, sumsq
    real(dp), allocatable :: normal(:, :)
    k = int(floor(df))
    frac = df - real(k, dp)
    normal = random_normal_matrix(max(1, k), 1)
    sumsq = sum(normal(:, 1) ** 2)
    if (frac > 1.0e-12_dp) then
      call random_number(u)
      u = max(u, tiny(1.0_dp))
      z = -2.0_dp * log(u)
      sumsq = sumsq + frac * z
    end if
    value = max(sumsq, tiny(1.0_dp))
    if (k == 0) then
      i = 0
    end if
  end function chi_square_draw

end module tsmarch_dcc
