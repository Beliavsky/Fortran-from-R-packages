! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_influence
  use rpeif_kinds, only : dp
  use rpeif_types, only : nuisance_parameters, rpeif_options, influence_result, &
    rpeif_success, rpeif_invalid_argument, rpeif_numerical_failure, rpeif_unknown_estimator
  use rpeif_nuisance, only : nuisance_parameters_fn
  use rpeif_stats, only : mean_value, sample_sd, quantile_type7, gaussian_kde_density, &
    lower_partial_moment, upper_partial_moment, lower_string
  use rpeif_robust, only : robust_clean, robust_mean_influence
  use rpeif_prewhiten, only : ar_prewhiten
  use robstattm_psi, only : rho_prime, tuning_for_efficiency
  implicit none
  private
  public :: influence_from_data, influence_from_nuisance, influence_series, evaluate_shape
  public :: supported_estimator
contains
  logical function supported_estimator(estimator) result(supported)
    character(len=*), intent(in) :: estimator
    character(len=:), allocatable :: name

    name = canonical_estimator(estimator)
    select case (name)
    case ('mean', 'sd', 'semisd', 'var', 'es', 'sr', 'sor', 'dsr', 'esratio', &
          'varratio', 'rachevratio', 'robmean', 'lpm', 'omegaratio')
      supported = .true.
    case default
      supported = .false.
    end select
  end function supported_estimator

  subroutine influence_from_data(estimator, x_eval, returns, values, options, status, message)
    character(len=*), intent(in) :: estimator
    real(dp), intent(in) :: x_eval(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    type(rpeif_options) :: opts
    character(len=:), allocatable :: name
    real(dp) :: alpha, beta, efficiency
    integer :: stat
    character(len=160) :: msg

    opts = rpeif_options()
    if (present(options)) opts = options
    name = canonical_estimator(estimator)
    call resolve_parameters(name, opts, alpha, beta, efficiency)
    stat = rpeif_success
    msg = ''
    allocate(values(size(x_eval)))
    values = 0.0_dp

    if (.not. supported_estimator(name)) then
      stat = rpeif_unknown_estimator
      msg = 'Unknown estimator.'
    else if (size(returns) < 2 .or. size(x_eval) == 0) then
      stat = rpeif_invalid_argument
      msg = 'At least two returns and one evaluation value are required.'
    else if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. beta <= 0.0_dp .or. beta >= 1.0_dp) then
      stat = rpeif_invalid_argument
      msg = 'Tail probabilities must lie strictly between zero and one.'
    else if (opts%lpm_order < 1 .or. opts%lpm_order > 2) then
      stat = rpeif_invalid_argument
      msg = 'The LPM order must be one or two.'
    else
      call compute_from_data(name, x_eval, returns, values, opts, alpha, beta, efficiency, stat, msg)
    end if

    if (present(status)) status = stat
    if (present(message)) message = trim(msg)
  end subroutine influence_from_data

  subroutine influence_from_nuisance(estimator, x_eval, pars, values, options, status, message)
    character(len=*), intent(in) :: estimator
    real(dp), intent(in) :: x_eval(:)
    type(nuisance_parameters), intent(in) :: pars
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    type(rpeif_options) :: opts
    character(len=:), allocatable :: name
    real(dp) :: alpha, beta, efficiency
    integer :: stat
    character(len=160) :: msg

    opts = rpeif_options()
    if (present(options)) opts = options
    name = canonical_estimator(estimator)
    call resolve_parameters(name, opts, alpha, beta, efficiency)
    allocate(values(size(x_eval)))
    values = 0.0_dp
    stat = rpeif_success
    msg = ''

    if (.not. supported_estimator(name)) then
      stat = rpeif_unknown_estimator
      msg = 'Unknown estimator.'
    else if (size(x_eval) == 0) then
      stat = rpeif_invalid_argument
      msg = 'At least one evaluation value is required.'
    else
      call compute_from_nuisance(name, x_eval, pars, values, opts, alpha, beta, efficiency, stat, msg)
    end if

    if (present(status)) status = stat
    if (present(message)) message = trim(msg)
  end subroutine influence_from_nuisance

  subroutine influence_series(estimator, returns, result, options)
    character(len=*), intent(in) :: estimator
    real(dp), intent(in) :: returns(:)
    type(influence_result), intent(out) :: result
    type(rpeif_options), intent(in), optional :: options
    type(rpeif_options) :: opts
    real(dp), allocatable :: working(:), raw_values(:), residuals(:)
    real(dp) :: efficiency
    integer :: stat, prewhite_status
    character(len=160) :: msg
    character(len=:), allocatable :: name

    opts = rpeif_options()
    if (present(options)) opts = options
    name = canonical_estimator(estimator)
    result%estimator = name
    if (size(returns) < 2) then
      result%status = rpeif_invalid_argument
      result%message = 'At least two returns are required.'
      allocate(result%x(0), result%values(0))
      return
    end if

    efficiency = opts%efficiency
    if (efficiency <= 0.0_dp) efficiency = 0.99_dp
    if (opts%clean_outliers) then
      call robust_clean(returns, working, efficiency, 'mopt', stat)
    else
      allocate(working(size(returns)))
      working = returns
    end if

    call influence_from_data(name, working, working, raw_values, opts, stat, msg)
    if (stat /= rpeif_success .and. stat /= rpeif_numerical_failure) then
      result%status = stat
      result%message = msg
      result%x = working
      result%values = raw_values
      return
    end if

    if (opts%prewhiten) then
      call ar_prewhiten(raw_values, opts%ar_order, residuals, status=prewhite_status)
      result%values = residuals
      if (prewhite_status /= 0 .and. stat == rpeif_success) then
        stat = rpeif_numerical_failure
        msg = 'Influence values computed; AR prewhitening used a fallback.'
      end if
    else
      result%values = raw_values
    end if
    result%x = working
    result%status = stat
    result%message = msg
  end subroutine influence_series

  subroutine evaluate_shape(estimator, result, options, nuisance, returns, k, step)
    character(len=*), intent(in) :: estimator
    type(influence_result), intent(out) :: result
    type(rpeif_options), intent(in), optional :: options
    type(nuisance_parameters), intent(in), optional :: nuisance
    real(dp), intent(in), optional :: returns(:)
    integer, intent(in), optional :: k
    real(dp), intent(in), optional :: step
    type(rpeif_options) :: opts
    type(nuisance_parameters) :: pars
    real(dp) :: lower, upper, increment, mu, sigma
    integer :: range_k, n, i, stat
    character(len=160) :: msg

    opts = rpeif_options()
    if (present(options)) opts = options
    range_k = 4
    if (present(k)) range_k = max(1, k)
    increment = 0.001_dp
    if (present(step)) increment = step
    if (increment <= 0.0_dp) then
      result%status = rpeif_invalid_argument
      result%message = 'Shape-grid step must be positive.'
      allocate(result%x(0), result%values(0))
      return
    end if

    if (present(returns)) then
      if (size(returns) < 2) then
        result%status = rpeif_invalid_argument
        result%message = 'At least two returns are required.'
        allocate(result%x(0), result%values(0))
        return
      end if
      mu = mean_value(returns)
      sigma = sample_sd(returns)
      lower = mu - real(range_k, dp) * sigma
      upper = mu + real(range_k, dp) * sigma
    else
      lower = 0.005_dp - real(range_k, dp) * 0.07_dp
      upper = 0.005_dp + real(range_k, dp) * 0.07_dp
    end if

    n = max(1, int(floor((upper - lower) / increment + 1.0e-10_dp)) + 1)
    allocate(result%x(n))
    do i = 1, n
      result%x(i) = lower + real(i - 1, dp) * increment
    end do

    if (present(returns)) then
      call influence_from_data(estimator, result%x, returns, result%values, opts, stat, msg)
    else
      if (present(nuisance)) then
        pars = nuisance
      else
        call nuisance_parameters_fn(pars, status=stat)
      end if
      call influence_from_nuisance(estimator, result%x, pars, result%values, opts, stat, msg)
    end if
    result%estimator = canonical_estimator(estimator)
    result%status = stat
    result%message = msg
  end subroutine evaluate_shape

  subroutine compute_from_data(name, x, returns, values, opts, alpha, beta, efficiency, status, message)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: x(:), returns(:), alpha, beta, efficiency
    real(dp), intent(out) :: values(:)
    type(rpeif_options), intent(in) :: opts
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp) :: mu, sd, semisd, semimean, q_alpha, es_alpha, q_beta, eg_beta
    real(dp) :: lpm1, lpm2, upm1, omega, ratio, fq_alpha, denominator, sr, sor, dsr
    real(dp), allocatable :: selected(:)
    integer :: i, n_selected, robust_status
    logical :: indicator

    status = rpeif_success
    message = ''
    mu = mean_value(returns)
    sd = sample_sd(returns)

    select case (name)
    case ('mean')
      values = x - mu

    case ('sd')
      if (sd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Standard deviation is zero.')
        return
      end if
      values = ((x - mu) ** 2 - sd ** 2) / (2.0_dp * sd)

    case ('semisd', 'dsr', 'sor')
      semisd = sqrt(sum(merge((returns - mu) ** 2, 0.0_dp, returns <= mu)) / real(size(returns), dp))
      semimean = sum(merge(returns - mu, 0.0_dp, returns <= mu)) / real(size(returns), dp)
      if (name == 'semisd') then
        if (semisd <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
          return
        end if
        do i = 1, size(x)
          indicator = x(i) <= mu
          values(i) = (merge((x(i) - mu) ** 2, 0.0_dp, indicator) - &
            2.0_dp * semimean * (x(i) - mu) - semisd ** 2) / (2.0_dp * semisd)
        end do
      else if (name == 'dsr') then
        if (semisd <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
          return
        end if
        dsr = (mu - opts%risk_free) / (semisd * sqrt(2.0_dp))
        do i = 1, size(x)
          indicator = x(i) <= mu
          values(i) = (-dsr * merge((x(i) - mu) ** 2, 0.0_dp, indicator) / (2.0_dp * semisd ** 2) + &
            (x(i) - mu) * (1.0_dp / semisd + dsr * semimean / semisd ** 2) + dsr / 2.0_dp) / sqrt(2.0_dp)
        end do
      else
        if (trim(lower_string(opts%sortino_threshold)) == 'mean') then
          if (semisd <= tiny(1.0_dp)) then
            call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
            return
          end if
          sor = (mu - opts%risk_free) / semisd
          do i = 1, size(x)
            indicator = x(i) <= mu
            values(i) = -sor * merge((x(i) - mu) ** 2, 0.0_dp, indicator) / (2.0_dp * semisd ** 2) + &
              (x(i) - mu) * (1.0_dp / semisd + sor * semimean / semisd ** 2) + sor / 2.0_dp
          end do
        else
          lpm2 = lower_partial_moment(returns, opts%threshold_constant, 2)
          if (lpm2 <= tiny(1.0_dp)) then
            call numerical_error(values, status, message, 'Second lower partial moment is zero.')
            return
          end if
          sor = (mu - opts%threshold_constant) / sqrt(lpm2)
          do i = 1, size(x)
            indicator = x(i) <= opts%threshold_constant
            values(i) = -sor * merge((x(i) - opts%threshold_constant) ** 2, 0.0_dp, indicator) / &
              (2.0_dp * lpm2) + (x(i) - mu) / sqrt(lpm2) + sor / 2.0_dp
          end do
        end if
      end if

    case ('var', 'es', 'esratio', 'varratio', 'rachevratio')
      q_alpha = quantile_type7(returns, alpha)
      select case (name)
      case ('var')
        fq_alpha = gaussian_kde_density(returns, q_alpha)
        if (fq_alpha <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Estimated density at the quantile is zero.')
          return
        end if
        do i = 1, size(x)
          values(i) = (merge(1.0_dp, 0.0_dp, x(i) <= q_alpha) - alpha) / fq_alpha
        end do
      case ('es', 'esratio', 'rachevratio')
        n_selected = count(returns <= q_alpha)
        if (n_selected == 0) then
          call numerical_error(values, status, message, 'No observations lie in the lower tail.')
          return
        end if
        allocate(selected(n_selected))
        selected = pack(returns, returns <= q_alpha)
        es_alpha = -mean_value(selected)
        deallocate(selected)
        if (name == 'es') then
          do i = 1, size(x)
            values(i) = merge((-x(i) + q_alpha) / alpha, 0.0_dp, x(i) <= q_alpha) - q_alpha - es_alpha
          end do
        else if (name == 'esratio') then
          if (abs(es_alpha) <= tiny(1.0_dp)) then
            call numerical_error(values, status, message, 'Expected shortfall is zero.')
            return
          end if
          ratio = (mu - opts%risk_free) / es_alpha
          do i = 1, size(x)
            values(i) = (x(i) - mu) / es_alpha - ratio / es_alpha * &
              (merge((-x(i) + q_alpha) / alpha, 0.0_dp, x(i) <= q_alpha) - q_alpha - es_alpha)
          end do
        else
          q_beta = quantile_type7(returns, 1.0_dp - beta)
          n_selected = count(returns >= q_beta)
          if (n_selected == 0 .or. abs(es_alpha) <= tiny(1.0_dp)) then
            call numerical_error(values, status, message, 'Tail expectation is undefined.')
            return
          end if
          allocate(selected(n_selected))
          selected = pack(returns, returns >= q_beta)
          eg_beta = mean_value(selected)
          deallocate(selected)
          ratio = eg_beta / es_alpha
          do i = 1, size(x)
            values(i) = (merge((x(i) - q_beta) / beta, 0.0_dp, x(i) >= q_beta) + q_beta - eg_beta) / es_alpha - &
              ratio / es_alpha * (merge(-(x(i) - q_alpha) / alpha, 0.0_dp, x(i) <= q_alpha) - q_alpha - es_alpha)
          end do
        end if
      case ('varratio')
        denominator = -q_alpha
        if (abs(denominator) <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Value at risk is zero.')
          return
        end if
        ratio = (mu - opts%risk_free) / denominator
        if (opts%source_compatibility) then
          fq_alpha = q_alpha
        else
          fq_alpha = gaussian_kde_density(returns, q_alpha)
        end if
        if (abs(fq_alpha) <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'VaR-ratio density denominator is zero.')
          return
        end if
        do i = 1, size(x)
          values(i) = (x(i) - mu) / denominator - ratio / denominator * &
            (merge(1.0_dp, 0.0_dp, x(i) <= q_alpha) - alpha) / fq_alpha
        end do
      end select

    case ('sr')
      if (sd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Standard deviation is zero.')
        return
      end if
      sr = mu
      if (.not. opts%source_compatibility) sr = mu - opts%risk_free
      values = (x - mu) / sd - 0.5_dp * sr / sd ** 3 * ((x - mu) ** 2 - sd ** 2)

    case ('lpm')
      lpm1 = lower_partial_moment(returns, opts%threshold_constant, 1)
      lpm2 = lower_partial_moment(returns, opts%threshold_constant, 2)
      do i = 1, size(x)
        indicator = x(i) <= opts%threshold_constant
        if (opts%lpm_order == 1) then
          values(i) = merge(opts%threshold_constant - x(i), 0.0_dp, indicator) - lpm1
        else
          values(i) = merge((opts%threshold_constant - x(i)) ** 2, 0.0_dp, indicator) - lpm2
        end if
      end do

    case ('omegaratio')
      lpm1 = lower_partial_moment(returns, opts%threshold_constant, 1)
      if (abs(lpm1) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'First lower partial moment is zero.')
        return
      end if
      upm1 = upper_partial_moment(returns, opts%threshold_constant, 1, opts%source_compatibility)
      omega = 1.0_dp + (mu - opts%threshold_constant) / lpm1
      do i = 1, size(x)
        values(i) = (merge(x(i) - opts%threshold_constant, 0.0_dp, x(i) >= opts%threshold_constant) - upm1) / lpm1 - &
          omega / lpm1 * (merge(opts%threshold_constant - x(i), 0.0_dp, x(i) <= opts%threshold_constant) - lpm1)
      end do

    case ('robmean')
      call robust_mean_influence(x, returns, selected, opts%robust_family, efficiency, robust_status)
      values = selected
      if (robust_status /= 0) then
        status = rpeif_numerical_failure
        message = 'Robust location calculation did not fully converge.'
      end if
    end select
  end subroutine compute_from_data

  subroutine compute_from_nuisance(name, x, pars, values, opts, alpha, beta, efficiency, status, message)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: x(:), alpha, beta, efficiency
    type(nuisance_parameters), intent(in) :: pars
    real(dp), intent(out) :: values(:)
    type(rpeif_options), intent(in) :: opts
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp) :: denominator, ratio, dsr, sr, tuning
    integer :: i
    logical :: indicator

    status = rpeif_success
    message = ''
    select case (name)
    case ('mean')
      values = x - pars%mu
    case ('sd')
      if (pars%sd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Standard deviation is zero.')
        return
      end if
      values = ((x - pars%mu) ** 2 - pars%sd ** 2) / (2.0_dp * pars%sd)
    case ('semisd')
      if (pars%semisd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
        return
      end if
      do i = 1, size(x)
        indicator = x(i) <= pars%mu
        values(i) = (merge((x(i) - pars%mu) ** 2, 0.0_dp, indicator) - &
          2.0_dp * pars%semimean * (x(i) - pars%mu) - pars%semisd ** 2) / (2.0_dp * pars%semisd)
      end do
    case ('dsr')
      if (pars%semisd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
        return
      end if
      if (opts%source_compatibility) then
        dsr = pars%dsr
      else
        dsr = (pars%mu - opts%risk_free) / (pars%semisd * sqrt(2.0_dp))
      end if
      do i = 1, size(x)
        indicator = x(i) <= pars%mu
        values(i) = (-dsr * merge((x(i) - pars%mu) ** 2, 0.0_dp, indicator) / (2.0_dp * pars%semisd ** 2) + &
          (x(i) - pars%mu) * (1.0_dp / pars%semisd + dsr * pars%semimean / pars%semisd ** 2) + dsr / 2.0_dp) / sqrt(2.0_dp)
      end do
    case ('es')
      do i = 1, size(x)
        values(i) = merge((-x(i) + pars%q_alpha) / alpha, 0.0_dp, x(i) <= pars%q_alpha) - &
          pars%q_alpha - pars%es_alpha
      end do
    case ('esratio')
      if (abs(pars%es_alpha) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Expected shortfall is zero.')
        return
      end if
      ratio = pars%es_ratio
      do i = 1, size(x)
        values(i) = (x(i) - pars%mu) / pars%es_alpha - ratio / pars%es_alpha * &
          (merge((-x(i) + pars%q_alpha) / alpha, 0.0_dp, x(i) <= pars%q_alpha) - pars%q_alpha - pars%es_alpha)
      end do
    case ('lpm')
      do i = 1, size(x)
        indicator = x(i) <= opts%threshold_constant
        if (opts%lpm_order == 1) then
          values(i) = merge(opts%threshold_constant - x(i), 0.0_dp, indicator) - pars%lpm1
        else
          values(i) = merge((opts%threshold_constant - x(i)) ** 2, 0.0_dp, indicator) - pars%lpm2
        end if
      end do
    case ('omegaratio')
      if (abs(pars%lpm1) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'First lower partial moment is zero.')
        return
      end if
      do i = 1, size(x)
        values(i) = (merge(x(i) - opts%threshold_constant, 0.0_dp, x(i) >= opts%threshold_constant) - pars%upm1) / pars%lpm1 - &
          pars%omega / pars%lpm1 * (merge(opts%threshold_constant - x(i), 0.0_dp, x(i) <= opts%threshold_constant) - pars%lpm1)
      end do
    case ('rachevratio')
      if (abs(pars%es_alpha) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Expected shortfall is zero.')
        return
      end if
      do i = 1, size(x)
        values(i) = (merge((x(i) - pars%q_beta) / beta, 0.0_dp, x(i) >= pars%q_beta) + &
          pars%q_beta - pars%eg_beta) / pars%es_alpha - &
          pars%rachev_ratio / pars%es_alpha * &
          (merge(-(x(i) - pars%q_alpha) / alpha, 0.0_dp, x(i) <= pars%q_alpha) - pars%q_alpha - pars%es_alpha)
      end do
    case ('robmean')
      if (pars%sd <= tiny(1.0_dp) .or. abs(pars%psi_prime) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Robust-mean nuisance scale or derivative is zero.')
        return
      end if
      tuning = tuning_for_efficiency(efficiency, opts%robust_family)
      do i = 1, size(x)
        values(i) = pars%sd * rho_prime((x(i) - pars%mu) / pars%sd, opts%robust_family, tuning) / pars%psi_prime
      end do
    case ('sor')
      if (trim(lower_string(opts%sortino_threshold)) == 'mean') then
        if (pars%semisd <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Semi-standard deviation is zero.')
          return
        end if
        do i = 1, size(x)
          indicator = x(i) <= pars%mu
          values(i) = -pars%sor_mu * merge((x(i) - pars%mu) ** 2, 0.0_dp, indicator) / (2.0_dp * pars%semisd ** 2) + &
            (x(i) - pars%mu) * (1.0_dp / pars%semisd + pars%sor_mu * pars%semimean / pars%semisd ** 2) + pars%sor_mu / 2.0_dp
        end do
      else
        if (pars%lpm2 <= tiny(1.0_dp)) then
          call numerical_error(values, status, message, 'Second lower partial moment is zero.')
          return
        end if
        do i = 1, size(x)
          values(i) = -pars%sor_c * merge((x(i) - opts%threshold_constant) ** 2, 0.0_dp, &
            x(i) <= opts%threshold_constant) / (2.0_dp * pars%lpm2) + &
            (x(i) - pars%mu) / sqrt(pars%lpm2) + pars%sor_c / 2.0_dp
        end do
      end if
    case ('sr')
      if (pars%sd <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Standard deviation is zero.')
        return
      end if
      sr = pars%mu
      if (.not. opts%source_compatibility) sr = pars%mu - opts%risk_free
      values = (x - pars%mu) / pars%sd - 0.5_dp * sr / pars%sd ** 3 * &
        ((x - pars%mu) ** 2 - pars%sd ** 2)
    case ('var')
      if (pars%fq_alpha <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'Density at the quantile is zero.')
        return
      end if
      do i = 1, size(x)
        values(i) = (merge(1.0_dp, 0.0_dp, x(i) <= pars%q_alpha) - alpha) / pars%fq_alpha
      end do
    case ('varratio')
      denominator = -pars%q_alpha
      if (abs(denominator) <= tiny(1.0_dp) .or. abs(pars%fq_alpha) <= tiny(1.0_dp)) then
        call numerical_error(values, status, message, 'VaR-ratio nuisance denominator is zero.')
        return
      end if
      do i = 1, size(x)
        values(i) = (x(i) - pars%mu) / denominator - pars%var_ratio / denominator * &
          (merge(1.0_dp, 0.0_dp, x(i) <= pars%q_alpha) - alpha) / pars%fq_alpha
      end do
    end select
  end subroutine compute_from_nuisance

  subroutine resolve_parameters(name, opts, alpha, beta, efficiency)
    character(len=*), intent(in) :: name
    type(rpeif_options), intent(in) :: opts
    real(dp), intent(out) :: alpha, beta, efficiency

    alpha = opts%alpha
    if (alpha <= 0.0_dp) then
      select case (name)
      case ('var', 'varratio', 'es')
        alpha = 0.05_dp
      case default
        alpha = 0.1_dp
      end select
    end if
    beta = opts%beta
    if (beta <= 0.0_dp) beta = 0.1_dp
    efficiency = opts%efficiency
    if (efficiency <= 0.0_dp) then
      if (name == 'robmean') then
        efficiency = 0.95_dp
      else
        efficiency = 0.99_dp
      end if
    end if
  end subroutine resolve_parameters

  function canonical_estimator(estimator) result(name)
    character(len=*), intent(in) :: estimator
    character(len=:), allocatable :: name
    character(len=:), allocatable :: raw

    raw = trim(lower_string(adjustl(estimator)))
    select case (raw)
    case ('omega', 'omegaratio', 'omega_ratio')
      name = 'omegaratio'
    case ('rachev', 'rachevratio', 'rachev_ratio', 'rachr')
      name = 'rachevratio'
    case ('robustmean', 'robmean', 'rob_mean')
      name = 'robmean'
    case ('es_ratio', 'esratio')
      name = 'esratio'
    case ('var_ratio', 'varratio')
      name = 'varratio'
    case ('semi_sd', 'semisd')
      name = 'semisd'
    case default
      name = raw
    end select
  end function canonical_estimator

  subroutine numerical_error(values, status, message, text)
    real(dp), intent(out) :: values(:)
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=*), intent(in) :: text

    values = 0.0_dp
    status = rpeif_numerical_failure
    message = text
  end subroutine numerical_error
end module rpeif_influence
