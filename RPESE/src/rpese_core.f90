! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_core
  use rpese_kinds, only : dp
  use rpese_types, only : rpese_options, se_result, matrix_se_result, periodogram_result, &
    rpese_success, rpese_invalid_argument, rpese_numerical_failure, rpese_dependency_failure, &
    se_if_iid, se_if_cor, se_if_cor_adapt, se_if_cor_pw, se_boot_iid, se_boot_cor, &
    fit_exponential, fit_gamma
  use rpese_measures, only : point_estimate, canonical_measure, supported_measure
  use rpese_timeseries, only : fit_periodogram, lag1_correlation, fit_ar1, polynomial_design
  use rpese_bootstrap, only : bootstrap_iid_se, bootstrap_block_se
  use rpeif, only : rpeif_options, influence_result, influence_from_data, influence_series, &
    robust_clean, rpeif_success
  use rpeglmen, only : enet_options, fit_result, glmnet_exp, fit_glm_gamma_net, rpe_success
  implicit none
  private
  public :: estimate_se, estimate_se_matrix, influence_values, spectral_variance
contains
  subroutine estimate_se(returns, estimator, method, result, options)
    real(dp), intent(in) :: returns(:)
    character(len=*), intent(in) :: estimator
    integer, intent(in) :: method
    type(se_result), intent(out) :: result
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    real(dp), allocatable :: raw_if(:), pw_if(:), cleaned_returns(:)
    real(dp) :: se_cor_value, se_pw_value, phi, intercept, weight
    integer :: stat, stat2
    character(len=160) :: msg

    opts = rpese_options()
    if (present(options)) opts = options
    result%method = method
    result%estimator = canonical_measure(estimator)
    result%status = rpese_success
    result%message = ''

    if (.not. supported_measure(estimator)) then
      result%status = rpese_invalid_argument
      result%message = 'Unknown risk or performance measure.'
      return
    end if
    if (size(returns) < 4) then
      result%status = rpese_invalid_argument
      result%message = 'At least four observations are required.'
      return
    end if

    call point_estimate(returns, estimator, result%estimate, opts, stat, msg)
    if (stat /= rpese_success) then
      result%status = stat
      result%message = msg
      return
    end if

    call influence_values(returns, estimator, raw_if, opts, clean=.false., status=stat, message=msg)
    if (stat /= rpese_success .and. stat /= rpese_numerical_failure) then
      result%status = rpese_dependency_failure
      result%message = msg
      return
    end if
    result%influence_correlation = lag1_correlation(raw_if)
    call fit_ar1(raw_if, phi, intercept, pw_if, stat2)
    result%prewhitened_influence_correlation = lag1_correlation(pw_if)

    if (opts%clean_outliers) then
      call robust_clean(returns, cleaned_returns, opts%robust_efficiency, 'mopt', stat2)
      result%return_correlation = lag1_correlation(cleaned_returns)
    else
      result%return_correlation = lag1_correlation(returns)
    end if

    select case (method)
    case (se_if_iid)
      result%standard_error = sqrt(sum(raw_if * raw_if) / real(size(raw_if) * size(raw_if), dp))
    case (se_if_cor)
      call correlated_if_se(returns, estimator, .false., result%standard_error, &
        result%coefficients, result%ar1_coefficient, opts, stat, msg)
    case (se_if_cor_pw)
      call correlated_if_se(returns, estimator, .true., result%standard_error, &
        result%coefficients, result%ar1_coefficient, opts, stat, msg)
    case (se_if_cor_adapt)
      call correlated_if_se(returns, estimator, .false., se_cor_value, result%coefficients, &
        result%ar1_coefficient, opts, stat, msg)
      if (stat == rpese_success) then
        call correlated_if_se(returns, estimator, .true., se_pw_value, phi=phi, &
          options=opts, status=stat2, message=msg)
        if (stat2 /= rpese_success) stat = stat2
      end if
      if (stat == rpese_success) then
        call fit_ar1(returns, phi, intercept, status=stat2)
        if (opts%source_compatibility) then
          if (phi >= 0.0_dp .and. phi < opts%adaptive_a) then
            weight = 0.0_dp
          else if (phi >= opts%adaptive_a .and. phi <= opts%adaptive_b) then
            weight = (phi - opts%adaptive_a) / (opts%adaptive_b - opts%adaptive_a)
          else
            weight = 1.0_dp
          end if
        else
          weight = (max(0.0_dp, phi) - opts%adaptive_a) / &
            max(tiny(1.0_dp), opts%adaptive_b - opts%adaptive_a)
          weight = max(0.0_dp, min(1.0_dp, weight))
        end if
        result%adaptive_weight = weight
        result%ar1_coefficient = phi
        result%standard_error = (1.0_dp - weight) * se_cor_value + weight * se_pw_value
      end if
    case (se_boot_iid)
      call bootstrap_iid_se(returns, estimator, result%standard_error, opts, stat, msg)
    case (se_boot_cor)
      call bootstrap_block_se(returns, estimator, result%standard_error, opts, stat, msg)
    case default
      stat = rpese_invalid_argument
      msg = 'Unknown standard-error method.'
    end select

    if (stat /= rpese_success) then
      result%status = stat
      result%message = msg
    else
      result%status = rpese_success
      result%message = 'completed'
    end if
  end subroutine estimate_se

  subroutine estimate_se_matrix(data, estimator, method, result, options)
    real(dp), intent(in) :: data(:, :)
    character(len=*), intent(in) :: estimator
    integer, intent(in) :: method
    type(matrix_se_result), intent(out) :: result
    type(rpese_options), intent(in), optional :: options
    integer :: j

    if (size(data, 1) < 4 .or. size(data, 2) < 1) then
      result%status = rpese_invalid_argument
      result%message = 'The data matrix must have at least four rows and one column.'
      allocate(result%column(0))
      return
    end if
    allocate(result%column(size(data, 2)))
    do j = 1, size(data, 2)
      call estimate_se(data(:, j), estimator, method, result%column(j), options)
      if (result%column(j)%status /= rpese_success) then
        result%status = result%column(j)%status
        write(result%message, '(a,i0,a,a)') 'Column ', j, ': ', trim(result%column(j)%message)
        return
      end if
    end do
    result%status = rpese_success
    result%message = 'completed'
  end subroutine estimate_se_matrix

  subroutine influence_values(returns, estimator, values, options, clean, status, message)
    real(dp), intent(in) :: returns(:)
    character(len=*), intent(in) :: estimator
    real(dp), allocatable, intent(out) :: values(:)
    type(rpese_options), intent(in) :: options
    logical, intent(in), optional :: clean
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    type(rpeif_options) :: if_options
    type(influence_result) :: if_result
    logical :: do_clean
    integer :: stat
    character(len=160) :: msg

    do_clean = options%clean_outliers
    if (present(clean)) do_clean = clean
    call map_if_options(options, if_options)
    if_options%clean_outliers = do_clean
    if_options%prewhiten = .false.
    call influence_series(canonical_measure(estimator), returns, if_result, if_options)
    values = if_result%values
    stat = if_result%status
    msg = if_result%message
    if (stat == rpeif_success) stat = rpese_success
    if (present(status)) status = stat
    if (present(message)) message = trim(msg)
  end subroutine influence_values

  subroutine correlated_if_se(returns, estimator, prewhiten, standard_error, coefficients, phi, &
      options, status, message)
    real(dp), intent(in) :: returns(:)
    character(len=*), intent(in) :: estimator
    logical, intent(in) :: prewhiten
    real(dp), intent(out) :: standard_error
    real(dp), allocatable, intent(out), optional :: coefficients(:)
    real(dp), intent(out), optional :: phi
    type(rpese_options), intent(in) :: options
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), allocatable :: if_values(:), residuals(:), fit_coefficients(:)
    real(dp) :: ar_coefficient, ar_intercept, variance
    integer :: stat
    character(len=160) :: msg

    call influence_values(returns, estimator, if_values, options, status=stat, message=msg)
    if (stat /= rpese_success .and. stat /= rpese_numerical_failure) then
      status = rpese_dependency_failure
      message = trim(msg)
      standard_error = 0.0_dp
      return
    end if
    ar_coefficient = 0.0_dp
    if (prewhiten) then
      call fit_ar1(if_values, ar_coefficient, ar_intercept, residuals, stat)
      call spectral_variance(residuals, variance, fit_coefficients, options, stat, msg)
      if (stat == rpese_success) then
        if (abs(1.0_dp - ar_coefficient) <= 1.0e-8_dp) then
          stat = rpese_numerical_failure
          msg = 'The prewhitening long-run variance correction is singular.'
        else
          variance = variance / (1.0_dp - ar_coefficient) ** 2
        end if
      end if
    else
      call spectral_variance(if_values, variance, fit_coefficients, options, stat, msg)
    end if
    if (present(phi)) phi = ar_coefficient
    if (present(coefficients)) coefficients = fit_coefficients
    if (stat == rpese_success .and. variance >= 0.0_dp) then
      standard_error = sqrt(variance)
      status = rpese_success
      message = 'completed'
    else
      standard_error = 0.0_dp
      status = stat
      message = trim(msg)
    end if
  end subroutine correlated_if_se

  subroutine spectral_variance(data, variance, coefficients, options, status, message)
    real(dp), intent(in) :: data(:)
    real(dp), intent(out) :: variance
    real(dp), allocatable, intent(out), optional :: coefficients(:)
    type(rpese_options), intent(in) :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    type(periodogram_result) :: periodogram
    type(enet_options) :: en_options
    type(fit_result) :: fit
    real(dp), allocatable :: design(:, :), means(:), scales(:)
    real(dp) :: log_spectrum_zero
    integer :: stat, j
    character(len=160) :: msg

    variance = 0.0_dp
    call fit_periodogram(data, periodogram, options)
    if (periodogram%status /= rpese_success) then
      stat = periodogram%status
      msg = periodogram%message
      if (present(coefficients)) allocate(coefficients(0))
      if (present(status)) status = stat
      if (present(message)) message = trim(msg)
      return
    end if
    call polynomial_design(periodogram%frequency, options%polynomial_degree, design, &
      options%standardize_design, means, scales, stat)
    if (stat /= rpese_success) then
      msg = 'Could not construct the polynomial frequency design.'
      if (present(coefficients)) allocate(coefficients(0))
      if (present(status)) status = stat
      if (present(message)) message = trim(msg)
      return
    end if

    en_options = enet_options()
    en_options%alpha = options%elastic_net_alpha
    en_options%num_lambda = max(1, options%num_lambda)
    en_options%max_iter = max(10, options%max_iterations)
    en_options%k_fold = max(2, options%cv_folds)
    en_options%k_fold_iter = max(1, options%cv_repeats)
    en_options%seed = options%seed
    en_options%has_intercept = .true.
    en_options%penalize_intercept = options%fitting_method == fit_gamma
    en_options%abs_tol = 1.0e-6_dp
    en_options%rel_tol = 1.0e-5_dp

    select case (options%fitting_method)
    case (fit_exponential)
      call glmnet_exp(design, periodogram%spectrum, fit, en_options)
    case (fit_gamma)
      call fit_glm_gamma_net(design, periodogram%spectrum, fit, en_options)
    case default
      stat = rpese_invalid_argument
      msg = 'Unknown periodogram fitting distribution.'
      if (present(coefficients)) allocate(coefficients(0))
      if (present(status)) status = stat
      if (present(message)) message = trim(msg)
      return
    end select

    if (.not. allocated(fit%coefficients) .or. fit%status /= rpe_success) then
      stat = rpese_dependency_failure
      msg = 'The penalized periodogram model did not converge: ' // trim(fit%message)
      if (present(coefficients)) allocate(coefficients(0))
    else
      if (options%standardize_design) then
        log_spectrum_zero = fit%coefficients(1)
        do j = 1, size(means)
          log_spectrum_zero = log_spectrum_zero - fit%coefficients(j + 1) * means(j) / scales(j)
        end do
      else
        log_spectrum_zero = fit%coefficients(1)
      end if
      if (log_spectrum_zero > log(huge(1.0_dp)) - 2.0_dp) then
        stat = rpese_numerical_failure
        msg = 'The fitted zero-frequency spectrum overflowed.'
      else
        variance = exp(log_spectrum_zero) / real(size(data), dp)
        stat = rpese_success
        msg = 'completed'
      end if
      if (present(coefficients)) coefficients = fit%coefficients
    end if
    if (present(status)) status = stat
    if (present(message)) message = trim(msg)
  end subroutine spectral_variance

  subroutine map_if_options(options, if_options)
    type(rpese_options), intent(in) :: options
    type(rpeif_options), intent(out) :: if_options
    if_options = rpeif_options()
    if_options%alpha = options%alpha
    if_options%beta = options%beta
    if_options%risk_free = options%risk_free
    if_options%threshold_constant = options%threshold_constant
    if_options%lpm_order = options%lpm_order
    if_options%sortino_threshold = options%sortino_threshold
    if_options%robust_family = options%robust_family
    if_options%efficiency = options%robust_efficiency
    if_options%clean_outliers = options%clean_outliers
    if_options%prewhiten = .false.
    if_options%source_compatibility = options%source_compatibility
  end subroutine map_if_options
end module rpese_core
