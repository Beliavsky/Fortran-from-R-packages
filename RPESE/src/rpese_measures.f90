! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_measures
  use rpese_kinds, only : dp
  use rpese_types, only : rpese_options, rpese_success, rpese_invalid_argument, &
    rpese_numerical_failure, rpese_unknown_estimator
  use rpeif_stats, only : mean_value, sample_sd, quantile_type7, lower_partial_moment, &
    upper_partial_moment, lower_string
  use rpeif_robust, only : robust_location_scale
  implicit none
  private
  public :: point_estimate, supported_measure, canonical_measure
  public :: mean_measure, sd_measure, semisd_measure, var_measure, es_measure
  public :: sharpe_measure, sortino_measure, downside_sharpe_measure
  public :: es_ratio_measure, var_ratio_measure, rachev_ratio_measure
  public :: robust_mean_measure, lpm_measure, omega_ratio_measure, upm_measure
contains
  function canonical_measure(name_in) result(name)
    character(len=*), intent(in) :: name_in
    character(len=:), allocatable :: name
    character(len=:), allocatable :: raw

    raw = trim(lower_string(adjustl(name_in)))
    select case (raw)
    case ('mean')
      name = 'mean'
    case ('sd', 'standarddeviation', 'standard_deviation')
      name = 'sd'
    case ('semisd', 'semi_sd', 'semistandarddeviation')
      name = 'semisd'
    case ('var', 'valueatrisk', 'value_at_risk')
      name = 'var'
    case ('es', 'expectedshortfall', 'expected_shortfall')
      name = 'es'
    case ('sr', 'sharpe', 'sharperatio', 'sharpe_ratio')
      name = 'sr'
    case ('sor', 'sortino', 'sortinoratio', 'sortino_ratio')
      name = 'sor'
    case ('dsr', 'downside_sharpe', 'downside_sharpe_ratio')
      name = 'dsr'
    case ('esratio', 'es_ratio', 'expected_shortfall_ratio')
      name = 'esratio'
    case ('varratio', 'var_ratio', 'value_at_risk_ratio')
      name = 'varratio'
    case ('rachevratio', 'rachev_ratio')
      name = 'rachevratio'
    case ('robmean', 'robustmean', 'robust_mean')
      name = 'robmean'
    case ('lpm', 'lowerpartialmoment', 'lower_partial_moment')
      name = 'lpm'
    case ('omegaratio', 'omega', 'omega_ratio')
      name = 'omegaratio'
    case default
      name = raw
    end select
  end function canonical_measure

  logical function supported_measure(name) result(ok)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: key
    key = canonical_measure(name)
    select case (key)
    case ('mean', 'sd', 'semisd', 'var', 'es', 'sr', 'sor', 'dsr', &
          'esratio', 'varratio', 'rachevratio', 'robmean', 'lpm', 'omegaratio')
      ok = .true.
    case default
      ok = .false.
    end select
  end function supported_measure

  subroutine point_estimate(returns, estimator, value, options, status, message)
    real(dp), intent(in) :: returns(:)
    character(len=*), intent(in) :: estimator
    real(dp), intent(out) :: value
    type(rpese_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    type(rpese_options) :: opts
    character(len=:), allocatable :: key
    integer :: stat
    character(len=160) :: msg

    opts = rpese_options()
    if (present(options)) opts = options
    key = canonical_measure(estimator)
    value = 0.0_dp
    stat = rpese_success
    msg = ''

    if (size(returns) < 2) then
      stat = rpese_invalid_argument
      msg = 'At least two observations are required.'
    else if (.not. supported_measure(key)) then
      stat = rpese_unknown_estimator
      msg = 'Unknown risk or performance measure.'
    else if (opts%alpha <= 0.0_dp .or. opts%alpha >= 1.0_dp .or. &
             opts%beta <= 0.0_dp .or. opts%beta >= 1.0_dp) then
      stat = rpese_invalid_argument
      msg = 'Tail probabilities must lie strictly between zero and one.'
    else if (opts%lpm_order < 1) then
      stat = rpese_invalid_argument
      msg = 'The lower partial moment order must be positive.'
    else
      select case (key)
      case ('mean')
        value = mean_measure(returns)
      case ('sd')
        value = sd_measure(returns)
      case ('semisd')
        value = semisd_measure(returns)
      case ('var')
        value = var_measure(returns, opts%alpha)
      case ('es')
        value = es_measure(returns, opts%alpha, stat)
      case ('sr')
        value = sharpe_measure(returns, opts%risk_free, stat)
      case ('sor')
        value = sortino_measure(returns, opts%risk_free, opts%threshold_constant, &
          opts%sortino_threshold, stat)
      case ('dsr')
        value = downside_sharpe_measure(returns, opts%risk_free, opts%source_compatibility, stat)
      case ('esratio')
        value = es_ratio_measure(returns, opts%alpha, opts%risk_free, stat)
      case ('varratio')
        value = var_ratio_measure(returns, opts%alpha, opts%risk_free, stat)
      case ('rachevratio')
        value = rachev_ratio_measure(returns, opts%alpha, opts%beta, stat)
      case ('robmean')
        value = robust_mean_measure(returns, opts%robust_family, opts%robust_efficiency, stat)
      case ('lpm')
        value = lpm_measure(returns, opts%threshold_constant, opts%lpm_order)
      case ('omegaratio')
        value = omega_ratio_measure(returns, opts%threshold_constant, stat)
      end select
      if (stat /= rpese_success) msg = 'The requested measure is undefined for these data.'
    end if

    if (present(status)) status = stat
    if (present(message)) message = trim(msg)
  end subroutine point_estimate

  real(dp) function mean_measure(x) result(value)
    real(dp), intent(in) :: x(:)
    value = mean_value(x)
  end function mean_measure

  real(dp) function sd_measure(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sample_sd(x)
  end function sd_measure

  real(dp) function semisd_measure(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mu
    mu = mean_value(x)
    value = sqrt(sum(merge((x - mu) ** 2, 0.0_dp, x <= mu)) / real(size(x), dp))
  end function semisd_measure

  real(dp) function var_measure(x, alpha) result(value)
    real(dp), intent(in) :: x(:), alpha
    value = -quantile_type7(x, alpha)
  end function var_measure

  real(dp) function es_measure(x, alpha, status) result(value)
    real(dp), intent(in) :: x(:), alpha
    integer, intent(out), optional :: status
    real(dp) :: q
    integer :: n_tail
    q = quantile_type7(x, alpha)
    n_tail = count(x <= q)
    if (n_tail == 0) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = -sum(pack(x, x <= q)) / real(n_tail, dp)
      if (present(status)) status = rpese_success
    end if
  end function es_measure

  real(dp) function sharpe_measure(x, risk_free, status) result(value)
    real(dp), intent(in) :: x(:), risk_free
    integer, intent(out), optional :: status
    real(dp) :: scale
    scale = sample_sd(x)
    if (scale <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = (mean_value(x) - risk_free) / scale
      if (present(status)) status = rpese_success
    end if
  end function sharpe_measure

  real(dp) function sortino_measure(x, risk_free, threshold, threshold_mode, status) result(value)
    real(dp), intent(in) :: x(:), risk_free, threshold
    character(len=*), intent(in) :: threshold_mode
    integer, intent(out), optional :: status
    real(dp) :: mu, denominator
    character(len=:), allocatable :: mode

    value = 0.0_dp
    mu = mean_value(x)
    mode = trim(lower_string(adjustl(threshold_mode)))
    if (mode == 'mean') then
      denominator = semisd_measure(x) * sqrt(2.0_dp)
      if (denominator > tiny(1.0_dp)) value = (mu - risk_free) / denominator
    else
      denominator = sqrt(lower_partial_moment(x, threshold, 2))
      if (denominator > tiny(1.0_dp)) value = (mu - threshold) / denominator
    end if
    if (denominator <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      if (present(status)) status = rpese_success
    end if
  end function sortino_measure

  real(dp) function downside_sharpe_measure(x, risk_free, source_compatibility, status) result(value)
    real(dp), intent(in) :: x(:), risk_free
    logical, intent(in) :: source_compatibility
    integer, intent(out), optional :: status
    real(dp) :: denominator, numerator
    denominator = semisd_measure(x) * sqrt(2.0_dp)
    numerator = mean_value(x)
    if (.not. source_compatibility) numerator = numerator - risk_free
    if (denominator <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = numerator / denominator
      if (present(status)) status = rpese_success
    end if
  end function downside_sharpe_measure

  real(dp) function es_ratio_measure(x, alpha, risk_free, status) result(value)
    real(dp), intent(in) :: x(:), alpha, risk_free
    integer, intent(out), optional :: status
    real(dp) :: tail
    integer :: stat
    tail = es_measure(x, alpha, stat)
    if (stat /= rpese_success .or. abs(tail) <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = (mean_value(x) - risk_free) / tail
      if (present(status)) status = rpese_success
    end if
  end function es_ratio_measure

  real(dp) function var_ratio_measure(x, alpha, risk_free, status) result(value)
    real(dp), intent(in) :: x(:), alpha, risk_free
    integer, intent(out), optional :: status
    real(dp) :: risk
    risk = var_measure(x, alpha)
    if (abs(risk) <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = (mean_value(x) - risk_free) / risk
      if (present(status)) status = rpese_success
    end if
  end function var_ratio_measure

  real(dp) function rachev_ratio_measure(x, alpha, beta, status) result(value)
    real(dp), intent(in) :: x(:), alpha, beta
    integer, intent(out), optional :: status
    real(dp) :: q_lower, q_upper, lower_tail, upper_tail
    integer :: n_lower, n_upper

    q_lower = quantile_type7(x, alpha)
    q_upper = quantile_type7(x, 1.0_dp - beta)
    n_lower = count(x <= q_lower)
    n_upper = count(x >= q_upper)
    if (n_lower == 0 .or. n_upper == 0) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
      return
    end if
    lower_tail = -sum(pack(x, x <= q_lower)) / real(n_lower, dp)
    upper_tail = sum(pack(x, x >= q_upper)) / real(n_upper, dp)
    if (abs(lower_tail) <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = upper_tail / lower_tail
      if (present(status)) status = rpese_success
    end if
  end function rachev_ratio_measure

  real(dp) function robust_mean_measure(x, family, efficiency, status) result(value)
    real(dp), intent(in) :: x(:), efficiency
    character(len=*), intent(in) :: family
    integer, intent(out), optional :: status
    real(dp) :: scale
    logical :: converged
    call robust_location_scale(x, value, scale, family, efficiency, converged)
    if (present(status)) status = merge(rpese_success, rpese_numerical_failure, converged)
  end function robust_mean_measure

  real(dp) function lpm_measure(x, threshold, order) result(value)
    real(dp), intent(in) :: x(:), threshold
    integer, intent(in) :: order
    value = lower_partial_moment(x, threshold, order)
  end function lpm_measure

  real(dp) function upm_measure(x, threshold, order, source_compatibility) result(value)
    real(dp), intent(in) :: x(:), threshold
    integer, intent(in) :: order
    logical, intent(in), optional :: source_compatibility
    logical :: source
    source = .false.
    if (present(source_compatibility)) source = source_compatibility
    if (source) then
      value = sum(merge((threshold - x) ** order, 0.0_dp, x >= threshold)) / real(size(x), dp)
    else
      value = upper_partial_moment(x, threshold, order)
    end if
  end function upm_measure

  real(dp) function omega_ratio_measure(x, threshold, status) result(value)
    real(dp), intent(in) :: x(:), threshold
    integer, intent(out), optional :: status
    real(dp) :: lpm1
    lpm1 = lower_partial_moment(x, threshold, 1)
    if (abs(lpm1) <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(status)) status = rpese_numerical_failure
    else
      value = 1.0_dp + (mean_value(x) - threshold) / lpm1
      if (present(status)) status = rpese_success
    end if
  end function omega_ratio_measure
end module rpese_measures
