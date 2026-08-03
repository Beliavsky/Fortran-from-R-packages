! SPDX-License-Identifier: MIT
module mfgarch_components
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use mfgarch_kinds, only : dp
  use mfgarch_math, only : finite_value, pi_dp, sample_variance, variance_ignore_nan
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument, &
    mfgarch_dimension_error, mfgarch_numerical_error
  use mfgarch_types, only : mfgarch_model
  implicit none
  private

  public :: beta_weights, low_frequency_log_tau, build_tau
  public :: calculate_g, likelihood_contributions, log_likelihood
  public :: forecast_garch, forecast_tau, variance_ratio
  public :: model_parameter_count, model_to_raw, raw_to_model
  public :: model_parameters, parameter_names, valid_model

contains

  subroutine beta_weights(k, w1, w2, weights, status)
    integer, intent(in) :: k
    real(dp), intent(in) :: w1, w2
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp) :: x, total
    integer :: j

    status = mfgarch_success
    if (k < 0 .or. w1 <= 0.0_dp .or. w2 <= 0.0_dp) then
      status = mfgarch_invalid_argument
      allocate(weights(0))
      return
    end if
    allocate(weights(k))
    if (k == 0) return
    do j = 1, k
      x = real(j, dp) / real(k + 1, dp)
      weights(j) = x**(w1 - 1.0_dp) * (1.0_dp - x)**(w2 - 1.0_dp)
    end do
    total = sum(weights)
    if (.not. finite_value(total) .or. total <= tiny(1.0_dp)) then
      status = mfgarch_numerical_error
      weights = 0.0_dp
      return
    end if
    weights = weights / total
  end subroutine beta_weights

  subroutine low_frequency_log_tau(covariate, k, m, theta, w1, w2, log_tau, status)
    real(dp), intent(in) :: covariate(:)
    integer, intent(in) :: k
    real(dp), intent(in) :: m, theta, w1, w2
    real(dp), allocatable, intent(out) :: log_tau(:)
    integer, intent(out) :: status
    real(dp), allocatable :: weights(:)
    real(dp) :: nan_value
    integer :: t, j

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    call beta_weights(k, w1, w2, weights, status)
    if (status /= mfgarch_success) then
      allocate(log_tau(0))
      return
    end if
    if (size(covariate) <= k) then
      status = mfgarch_invalid_argument
      allocate(log_tau(0))
      return
    end if
    allocate(log_tau(size(covariate)))
    log_tau = nan_value
    do t = k + 1, size(covariate)
      log_tau(t) = m
      do j = 1, k
        log_tau(t) = log_tau(t) + theta * weights(j) * covariate(t-j)
      end do
    end do
  end subroutine low_frequency_log_tau

  subroutine build_tau(model, period, covariate, tau, status, period_two, covariate_two)
    type(mfgarch_model), intent(in) :: model
    integer, intent(in) :: period(:)
    real(dp), intent(in), optional :: covariate(:)
    real(dp), allocatable, intent(out) :: tau(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: period_two(:)
    real(dp), intent(in), optional :: covariate_two(:)
    real(dp), allocatable :: log_tau_one(:), log_tau_two(:)
    real(dp) :: log_value, nan_value
    integer :: i

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    status = mfgarch_success
    allocate(tau(size(period)))
    tau = nan_value

    if (model%k == 0) then
      tau = exp(model%m)
      return
    end if
    if (.not. present(covariate)) then
      status = mfgarch_invalid_argument
      return
    end if
    if (any(period < 1) .or. any(period > size(covariate))) then
      status = mfgarch_invalid_argument
      return
    end if
    call low_frequency_log_tau(covariate, model%k, model%m, model%theta, &
      model%w1, model%w2, log_tau_one, status)
    if (status /= mfgarch_success) return

    if (model%has_second) then
      if (.not. present(period_two) .or. .not. present(covariate_two)) then
        status = mfgarch_invalid_argument
        return
      end if
      if (size(period_two) /= size(period)) then
        status = mfgarch_dimension_error
        return
      end if
      if (any(period_two < 1) .or. any(period_two > size(covariate_two))) then
        status = mfgarch_invalid_argument
        return
      end if
      call low_frequency_log_tau(covariate_two, model%k_two, 0.0_dp, &
        model%theta_two, model%w1_two, model%w2_two, log_tau_two, status)
      if (status /= mfgarch_success) return
    end if

    do i = 1, size(period)
      if (.not. finite_value(log_tau_one(period(i)))) cycle
      log_value = log_tau_one(period(i))
      if (model%has_second) then
        if (.not. finite_value(log_tau_two(period_two(i)))) cycle
        log_value = log_value + log_tau_two(period_two(i))
      end if
      if (log_value > log(huge(1.0_dp)) - 2.0_dp .or. &
          log_value < log(tiny(1.0_dp)) + 2.0_dp) cycle
      tau(i) = exp(log_value)
    end do
  end subroutine build_tau

  subroutine calculate_g(standardized_returns, alpha, beta, gamma, g0, g, status)
    real(dp), intent(in) :: standardized_returns(:)
    real(dp), intent(in) :: alpha, beta, gamma, g0
    real(dp), allocatable, intent(out) :: g(:)
    integer, intent(out) :: status
    real(dp) :: omega, shock, nan_value
    integer :: i

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    status = mfgarch_success
    allocate(g(size(standardized_returns)))
    g = nan_value
    omega = 1.0_dp - alpha - beta - 0.5_dp * gamma
    if (alpha < 0.0_dp .or. beta < 0.0_dp .or. gamma < 0.0_dp .or. &
        omega <= 0.0_dp .or. g0 <= 0.0_dp) then
      status = mfgarch_invalid_argument
      return
    end if

    if (size(standardized_returns) > 0) then
      if (finite_value(standardized_returns(1))) g(1) = g0
    end if
    do i = 2, size(standardized_returns)
      if (.not. finite_value(standardized_returns(i))) cycle
      if (.not. finite_value(standardized_returns(i-1)) .or. &
          .not. finite_value(g(i-1))) then
        g(i) = g0
      else
        shock = alpha * standardized_returns(i-1)**2
        if (standardized_returns(i-1) < 0.0_dp) then
          shock = shock + gamma * standardized_returns(i-1)**2
        end if
        g(i) = omega + shock + beta * g(i-1)
      end if
      if (.not. finite_value(g(i)) .or. g(i) <= 0.0_dp) then
        status = mfgarch_numerical_error
        return
      end if
    end do
  end subroutine calculate_g

  subroutine likelihood_contributions(model, returns, period, contributions, tau, g, &
      residuals, status, covariate, period_two, covariate_two, log_variance_only)
    type(mfgarch_model), intent(in) :: model
    real(dp), intent(in) :: returns(:)
    integer, intent(in) :: period(:)
    real(dp), allocatable, intent(out) :: contributions(:)
    real(dp), allocatable, intent(out) :: tau(:), g(:), residuals(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: covariate(:), covariate_two(:)
    integer, intent(in), optional :: period_two(:)
    logical, intent(in), optional :: log_variance_only
    real(dp), allocatable :: standardized(:)
    real(dp) :: g0, variance_value, nan_value
    logical :: variance_only
    integer :: i

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    status = mfgarch_success
    if (size(returns) /= size(period) .or. size(returns) < 2) then
      status = mfgarch_dimension_error
      allocate(contributions(0), tau(0), g(0), residuals(0))
      return
    end if
    variance_only = .false.
    if (present(log_variance_only)) variance_only = log_variance_only

    if (model%k == 0) then
      call build_tau(model, period, tau=tau, status=status)
    else if (model%has_second) then
      if (.not. present(covariate) .or. .not. present(period_two) .or. &
          .not. present(covariate_two)) then
        status = mfgarch_invalid_argument
        allocate(contributions(0), g(0), residuals(0))
        return
      end if
      call build_tau(model, period, covariate, tau, status, period_two, covariate_two)
    else
      if (.not. present(covariate)) then
        status = mfgarch_invalid_argument
        allocate(contributions(0), g(0), residuals(0))
        return
      end if
      call build_tau(model, period, covariate, tau, status)
    end if
    if (status /= mfgarch_success) then
      allocate(contributions(0), g(0), residuals(0))
      return
    end if

    allocate(standardized(size(returns)), contributions(size(returns)), residuals(size(returns)))
    standardized = nan_value
    contributions = nan_value
    residuals = nan_value
    do i = 1, size(returns)
      if (finite_value(tau(i)) .and. tau(i) > 0.0_dp) then
        standardized(i) = (returns(i) - model%mu) / sqrt(tau(i))
      end if
    end do
    g0 = sample_variance(returns)
    if (g0 <= tiny(1.0_dp)) g0 = 1.0_dp
    call calculate_g(standardized, model%alpha, model%beta, merge(model%gamma, 0.0_dp, &
      model%asymmetric), g0, g, status)
    if (status /= mfgarch_success) return

    do i = 1, size(returns)
      if (.not. finite_value(tau(i)) .or. .not. finite_value(g(i))) cycle
      variance_value = tau(i) * g(i)
      if (variance_value <= 0.0_dp .or. .not. finite_value(variance_value)) then
        status = mfgarch_numerical_error
        return
      end if
      residuals(i) = (returns(i) - model%mu) / sqrt(variance_value)
      if (variance_only) then
        contributions(i) = log(variance_value)
      else
        contributions(i) = -0.5_dp * (log(2.0_dp * pi_dp) + log(variance_value) + residuals(i)**2)
      end if
    end do
  end subroutine likelihood_contributions

  function log_likelihood(model, returns, period, status, covariate, period_two, &
      covariate_two) result(value)
    type(mfgarch_model), intent(in) :: model
    real(dp), intent(in) :: returns(:)
    integer, intent(in) :: period(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: covariate(:), covariate_two(:)
    integer, intent(in), optional :: period_two(:)
    real(dp) :: value
    real(dp), allocatable :: c(:), tau(:), g(:), residuals(:)
    integer :: i, n

    if (model%k == 0) then
      call likelihood_contributions(model, returns, period, c, tau, g, residuals, status)
    else if (model%has_second) then
      if (.not. present(covariate) .or. .not. present(period_two) .or. &
          .not. present(covariate_two)) then
        status = mfgarch_invalid_argument
        value = -huge(1.0_dp)
        return
      end if
      call likelihood_contributions(model, returns, period, c, tau, g, residuals, &
        status, covariate, period_two, covariate_two)
    else
      if (.not. present(covariate)) then
        status = mfgarch_invalid_argument
        value = -huge(1.0_dp)
        return
      end if
      call likelihood_contributions(model, returns, period, c, tau, g, residuals, &
        status, covariate)
    end if
    if (status /= mfgarch_success) then
      value = -huge(1.0_dp)
      return
    end if
    value = 0.0_dp
    n = 0
    do i = 1, size(c)
      if (finite_value(c(i))) then
        value = value + c(i)
        n = n + 1
      end if
    end do
    if (n == 0 .or. .not. finite_value(value)) then
      status = mfgarch_numerical_error
      value = -huge(1.0_dp)
    end if
  end function log_likelihood

  pure function forecast_garch(omega, alpha, beta, gamma, g, return_value, steps_ahead) result(value)
    real(dp), intent(in) :: omega, alpha, beta, gamma, g, return_value
    integer, intent(in) :: steps_ahead
    real(dp) :: value, persistence, unconditional, first_variance

    persistence = alpha + beta + 0.5_dp * gamma
    unconditional = omega / max(1.0_dp - persistence, tiny(1.0_dp))
    first_variance = omega + (alpha + merge(gamma, 0.0_dp, return_value < 0.0_dp)) * &
      return_value**2 + beta * g
    if (steps_ahead <= 1) then
      value = first_variance
    else
      value = unconditional + persistence**(steps_ahead - 1) * (first_variance - unconditional)
    end if
  end function forecast_garch

  function forecast_tau(model, covariate, status, covariate_two) result(value)
    type(mfgarch_model), intent(in) :: model
    real(dp), intent(in), optional :: covariate(:), covariate_two(:)
    integer, intent(out) :: status
    real(dp) :: value, log_value
    real(dp), allocatable :: weights(:)
    integer :: j

    status = mfgarch_success
    if (model%k == 0) then
      value = exp(model%m)
      return
    end if
    if (.not. present(covariate) .or. size(covariate) < model%k) then
      status = mfgarch_invalid_argument
      value = 0.0_dp
      return
    end if
    call beta_weights(model%k, model%w1, model%w2, weights, status)
    if (status /= mfgarch_success) then
      value = 0.0_dp
      return
    end if
    log_value = model%m
    do j = 1, model%k
      log_value = log_value + model%theta * weights(j) * covariate(size(covariate) + 1 - j)
    end do
    if (model%has_second) then
      if (.not. present(covariate_two) .or. size(covariate_two) < model%k_two) then
        status = mfgarch_invalid_argument
        value = 0.0_dp
        return
      end if
      call beta_weights(model%k_two, model%w1_two, model%w2_two, weights, status)
      if (status /= mfgarch_success) then
        value = 0.0_dp
        return
      end if
      do j = 1, model%k_two
        log_value = log_value + model%theta_two * weights(j) * &
          covariate_two(size(covariate_two) + 1 - j)
      end do
    end if
    if (log_value > log(huge(1.0_dp)) - 2.0_dp) then
      status = mfgarch_numerical_error
      value = huge(1.0_dp)
    else
      value = exp(log_value)
    end if
  end function forecast_tau

  function variance_ratio(tau, g, groups, status) result(value)
    real(dp), intent(in) :: tau(:), g(:)
    integer, intent(in) :: groups(:)
    integer, intent(out) :: status
    real(dp) :: value
    real(dp), allocatable :: mean_tau(:), mean_total(:), log_tau(:), log_total(:)
    integer, allocatable :: counts(:)
    integer :: i, ng, nvalid

    status = mfgarch_success
    value = 0.0_dp
    if (size(tau) /= size(g) .or. size(tau) /= size(groups) .or. size(tau) == 0) then
      status = mfgarch_dimension_error
      return
    end if
    if (any(groups < 1)) then
      status = mfgarch_invalid_argument
      return
    end if
    ng = maxval(groups)
    allocate(mean_tau(ng), mean_total(ng), counts(ng), log_tau(ng), log_total(ng))
    mean_tau = 0.0_dp
    mean_total = 0.0_dp
    counts = 0
    do i = 1, size(tau)
      if (finite_value(tau(i)) .and. finite_value(g(i)) .and. tau(i) > 0.0_dp .and. g(i) > 0.0_dp) then
        mean_tau(groups(i)) = mean_tau(groups(i)) + tau(i)
        mean_total(groups(i)) = mean_total(groups(i)) + tau(i) * g(i)
        counts(groups(i)) = counts(groups(i)) + 1
      end if
    end do
    log_tau = ieee_value(0.0_dp, ieee_quiet_nan)
    log_total = ieee_value(0.0_dp, ieee_quiet_nan)
    nvalid = 0
    do i = 1, ng
      if (counts(i) > 0) then
        log_tau(i) = log(mean_tau(i) / real(counts(i), dp))
        log_total(i) = log(mean_total(i) / real(counts(i), dp))
        nvalid = nvalid + 1
      end if
    end do
    if (nvalid <= 1 .or. variance_ignore_nan(log_total) <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = 100.0_dp * variance_ignore_nan(log_tau) / variance_ignore_nan(log_total)
    end if
  end function variance_ratio

  pure logical function valid_model(model)
    type(mfgarch_model), intent(in) :: model
    real(dp) :: persistence

    persistence = model%alpha + model%beta
    if (model%asymmetric) persistence = persistence + 0.5_dp * model%gamma
    valid_model = model%alpha >= 0.0_dp .and. model%beta >= 0.0_dp .and. &
      model%gamma >= 0.0_dp .and. persistence < 1.0_dp .and. &
      model%w1 >= 1.0_dp .and. model%w2 >= 1.0_dp .and. model%k >= 0
    if (model%has_second) valid_model = valid_model .and. model%k_two >= 1 .and. &
      model%w1_two >= 1.0_dp .and. model%w2_two >= 1.0_dp
  end function valid_model

  pure integer function model_parameter_count(model) result(n)
    type(mfgarch_model), intent(in) :: model

    n = 4
    if (model%asymmetric) n = n + 1
    if (model%k > 0) then
      n = n + 1
      if (model%k > 1) then
        n = n + 1
        if (model%unrestricted_weights) n = n + 1
      end if
    end if
    if (model%has_second) then
      n = n + 1
      if (model%k_two > 1) n = n + 1
    end if
  end function model_parameter_count

  subroutine model_to_raw(model, raw, status)
    type(mfgarch_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: raw(:)
    integer, intent(out) :: status
    real(dp) :: intercept, floor_value
    integer :: i

    status = mfgarch_success
    if (.not. valid_model(model)) then
      status = mfgarch_invalid_argument
      allocate(raw(0))
      return
    end if
    floor_value = 1.0e-10_dp
    allocate(raw(model_parameter_count(model)))
    i = 1
    raw(i) = model%mu
    i = i + 1
    intercept = max(floor_value, 1.0_dp - model%alpha - model%beta - &
      merge(0.5_dp * model%gamma, 0.0_dp, model%asymmetric))
    raw(i) = log(max(model%alpha, floor_value) / intercept)
    i = i + 1
    raw(i) = log(max(model%beta, floor_value) / intercept)
    i = i + 1
    if (model%asymmetric) then
      raw(i) = log(max(0.5_dp * model%gamma, floor_value) / intercept)
      i = i + 1
    end if
    raw(i) = model%m
    i = i + 1
    if (model%k > 0) then
      raw(i) = model%theta
      i = i + 1
      if (model%k > 1) then
        if (model%unrestricted_weights) then
          raw(i) = log(max(model%w1 - 1.0_dp, 1.0e-8_dp))
          i = i + 1
        end if
        raw(i) = log(max(model%w2 - 1.0_dp, 1.0e-8_dp))
        i = i + 1
      end if
    end if
    if (model%has_second) then
      raw(i) = model%theta_two
      i = i + 1
      if (model%k_two > 1) raw(i) = log(max(model%w2_two - 1.0_dp, 1.0e-8_dp))
    end if
  end subroutine model_to_raw

  subroutine raw_to_model(raw, template, model, status)
    real(dp), intent(in) :: raw(:)
    type(mfgarch_model), intent(in) :: template
    type(mfgarch_model), intent(out) :: model
    integer, intent(out) :: status
    real(dp) :: ea, eb, eg, denominator, scale, max_raw
    integer :: i

    status = mfgarch_success
    if (size(raw) /= model_parameter_count(template)) then
      status = mfgarch_dimension_error
      model = template
      return
    end if
    model = template
    i = 1
    model%mu = raw(i)
    i = i + 1
    if (template%asymmetric) then
      max_raw = max(0.0_dp, max(raw(i), max(raw(i+1), raw(i+2))))
      ea = exp(max(-700.0_dp, min(700.0_dp, raw(i) - max_raw)))
      eb = exp(max(-700.0_dp, min(700.0_dp, raw(i+1) - max_raw)))
      eg = exp(max(-700.0_dp, min(700.0_dp, raw(i+2) - max_raw)))
      denominator = exp(-max_raw) + ea + eb + eg
      scale = 1.0_dp
      model%alpha = scale * ea / denominator
      model%beta = scale * eb / denominator
      model%gamma = 2.0_dp * scale * eg / denominator
      i = i + 3
    else
      max_raw = max(0.0_dp, max(raw(i), raw(i+1)))
      ea = exp(max(-700.0_dp, min(700.0_dp, raw(i) - max_raw)))
      eb = exp(max(-700.0_dp, min(700.0_dp, raw(i+1) - max_raw)))
      denominator = exp(-max_raw) + ea + eb
      scale = 1.0_dp
      model%alpha = scale * ea / denominator
      model%beta = scale * eb / denominator
      model%gamma = 0.0_dp
      i = i + 2
    end if
    model%m = raw(i)
    i = i + 1
    if (model%k > 0) then
      model%theta = raw(i)
      i = i + 1
      if (model%k > 1) then
        if (model%unrestricted_weights) then
          model%w1 = 1.0_dp + exp(max(-30.0_dp, min(30.0_dp, raw(i))))
          i = i + 1
        else
          model%w1 = 1.0_dp
        end if
        model%w2 = 1.0_dp + exp(max(-30.0_dp, min(30.0_dp, raw(i))))
        i = i + 1
      else
        model%w1 = 1.0_dp
        model%w2 = 1.0_dp
      end if
    end if
    if (model%has_second) then
      model%theta_two = raw(i)
      i = i + 1
      model%w1_two = 1.0_dp
      if (model%k_two > 1) then
        model%w2_two = 1.0_dp + exp(max(-30.0_dp, min(30.0_dp, raw(i))))
      else
        model%w2_two = 1.0_dp
      end if
    end if
  end subroutine raw_to_model

  subroutine model_parameters(model, parameters)
    type(mfgarch_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: parameters(:)
    integer :: i

    allocate(parameters(model_parameter_count(model)))
    i = 1
    parameters(i) = model%mu
    i = i + 1
    parameters(i) = model%alpha
    i = i + 1
    parameters(i) = model%beta
    i = i + 1
    if (model%asymmetric) then
      parameters(i) = model%gamma
      i = i + 1
    end if
    parameters(i) = model%m
    i = i + 1
    if (model%k > 0) then
      parameters(i) = model%theta
      i = i + 1
      if (model%k > 1) then
        if (model%unrestricted_weights) then
          parameters(i) = model%w1
          i = i + 1
        end if
        parameters(i) = model%w2
        i = i + 1
      end if
    end if
    if (model%has_second) then
      parameters(i) = model%theta_two
      i = i + 1
      if (model%k_two > 1) parameters(i) = model%w2_two
    end if
  end subroutine model_parameters

  subroutine parameter_names(model, names)
    type(mfgarch_model), intent(in) :: model
    character(len=16), allocatable, intent(out) :: names(:)
    integer :: i

    allocate(names(model_parameter_count(model)))
    i = 1
    names(i) = 'mu'
    i = i + 1
    names(i) = 'alpha'
    i = i + 1
    names(i) = 'beta'
    i = i + 1
    if (model%asymmetric) then
      names(i) = 'gamma'
      i = i + 1
    end if
    names(i) = 'm'
    i = i + 1
    if (model%k > 0) then
      names(i) = 'theta'
      i = i + 1
      if (model%k > 1) then
        if (model%unrestricted_weights) then
          names(i) = 'w1'
          i = i + 1
        end if
        names(i) = 'w2'
        i = i + 1
      end if
    end if
    if (model%has_second) then
      names(i) = 'theta_two'
      i = i + 1
      if (model%k_two > 1) names(i) = 'w2_two'
    end if
  end subroutine parameter_names

end module mfgarch_components
