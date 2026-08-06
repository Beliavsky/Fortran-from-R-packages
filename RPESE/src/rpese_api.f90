! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_api
  use rpese_kinds, only : dp
  use rpese_types, only : rpese_options, se_result, se_if_iid
  use rpese_core, only : estimate_se
  implicit none
  private
  public :: mean_se, sd_se, semisd_se, var_se, es_se, sr_se, sor_se, dsr_se
  public :: esratio_se, varratio_se, rachevratio_se, robmean_se, lpm_se, omegaratio_se
contains
  subroutine mean_se(data, result, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    call dispatch(data, 'mean', result, method, options)
  end subroutine mean_se

  subroutine sd_se(data, result, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    call dispatch(data, 'sd', result, method, options)
  end subroutine sd_se

  subroutine semisd_se(data, result, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    call dispatch(data, 'semisd', result, method, options)
  end subroutine semisd_se

  subroutine var_se(data, result, confidence, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(confidence)) opts%alpha = 1.0_dp - confidence
    call estimate_se(data, 'var', optional_method(method), result, opts)
  end subroutine var_se

  subroutine es_se(data, result, confidence, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(confidence)) opts%alpha = 1.0_dp - confidence
    call estimate_se(data, 'es', optional_method(method), result, opts)
  end subroutine es_se

  subroutine sr_se(data, result, risk_free, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: risk_free
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(risk_free)) opts%risk_free = risk_free
    call estimate_se(data, 'sr', optional_method(method), result, opts)
  end subroutine sr_se

  subroutine sor_se(data, result, threshold, threshold_mode, risk_free, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: threshold, risk_free
    character(len=*), intent(in), optional :: threshold_mode
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(threshold)) opts%threshold_constant = threshold
    if (present(threshold_mode)) opts%sortino_threshold = threshold_mode
    if (present(risk_free)) opts%risk_free = risk_free
    call estimate_se(data, 'sor', optional_method(method), result, opts)
  end subroutine sor_se

  subroutine dsr_se(data, result, risk_free, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: risk_free
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(risk_free)) opts%risk_free = risk_free
    call estimate_se(data, 'dsr', optional_method(method), result, opts)
  end subroutine dsr_se

  subroutine esratio_se(data, result, alpha, risk_free, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, risk_free
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(alpha)) opts%alpha = alpha
    if (present(risk_free)) opts%risk_free = risk_free
    call estimate_se(data, 'esratio', optional_method(method), result, opts)
  end subroutine esratio_se

  subroutine varratio_se(data, result, alpha, risk_free, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, risk_free
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(alpha)) opts%alpha = alpha
    if (present(risk_free)) opts%risk_free = risk_free
    call estimate_se(data, 'varratio', optional_method(method), result, opts)
  end subroutine varratio_se

  subroutine rachevratio_se(data, result, alpha, beta, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(alpha)) opts%alpha = alpha
    if (present(beta)) opts%beta = beta
    call estimate_se(data, 'rachevratio', optional_method(method), result, opts)
  end subroutine rachevratio_se

  subroutine robmean_se(data, result, family, efficiency, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    character(len=*), intent(in), optional :: family
    real(dp), intent(in), optional :: efficiency
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(family)) opts%robust_family = family
    if (present(efficiency)) opts%robust_efficiency = efficiency
    opts%clean_outliers = .false.
    call estimate_se(data, 'robmean', optional_method(method), result, opts)
  end subroutine robmean_se

  subroutine lpm_se(data, result, threshold, order, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: threshold
    integer, intent(in), optional :: order, method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(threshold)) opts%threshold_constant = threshold
    if (present(order)) opts%lpm_order = order
    call estimate_se(data, 'lpm', optional_method(method), result, opts)
  end subroutine lpm_se

  subroutine omegaratio_se(data, result, threshold, method, options)
    real(dp), intent(in) :: data(:)
    type(se_result), intent(out) :: result
    real(dp), intent(in), optional :: threshold
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    if (present(threshold)) opts%threshold_constant = threshold
    call estimate_se(data, 'omegaratio', optional_method(method), result, opts)
  end subroutine omegaratio_se

  subroutine dispatch(data, estimator, result, method, options)
    real(dp), intent(in) :: data(:)
    character(len=*), intent(in) :: estimator
    type(se_result), intent(out) :: result
    integer, intent(in), optional :: method
    type(rpese_options), intent(in), optional :: options
    type(rpese_options) :: opts
    opts = rpese_options()
    if (present(options)) opts = options
    call estimate_se(data, estimator, optional_method(method), result, opts)
  end subroutine dispatch

  integer function optional_method(method) result(value)
    integer, intent(in), optional :: method
    value = se_if_iid
    if (present(method)) value = method
  end function optional_method
end module rpese_api
