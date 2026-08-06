! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_gogarch
  use ghyp_kinds, only : dp, i8
  use tsgarch, only : fit_options, garch_fit, estimate_garch, forecast_garch, garch_forecast
  use tsd_moments, only : dskewness, dkurtosis
  use tsmarch_types
  use tsmarch_linalg
  use tsmarch_ica, only : fastica, radical
  implicit none
  private
  public :: estimate_gogarch, forecast_gogarch, gogarch_covariance
  public :: gogarch_correlation, gogarch_coskewness, gogarch_cokurtosis
  public :: portfolio_variance, portfolio_skewness, portfolio_kurtosis

contains

  function estimate_gogarch(data, specification, options, seed) result(fit)
    real(dp), intent(in) :: data(:, :)
    type(gogarch_spec), intent(in) :: specification
    type(fit_options), intent(in), optional :: options
    integer(i8), intent(in), optional :: seed
    type(gogarch_fit) :: fit
    type(fit_options) :: opt
    real(dp), allocatable :: variances(:, :), skewness(:), kurtosis(:)
    integer :: j, t, n, m

    n = size(data, 1)
    m = size(data, 2)
    if (n < max(30, 3 * m) .or. m < 2) then
      fit%status = tsm_invalid_argument
      fit%message = 'GO-GARCH requires at least two series and enough observations for ICA'
      return
    end if
    fit%spec = specification
    if (trim(lower_string(specification%ica_method)) == 'radical') then
      if (present(seed)) then
        fit%ica = radical(data, seed=seed)
      else
        fit%ica = radical(data)
      end if
    else
      if (present(seed)) then
        fit%ica = fastica(data, seed=seed)
      else
        fit%ica = fastica(data)
      end if
    end if
    if (fit%ica%status /= tsm_success .and. fit%ica%status /= tsm_no_convergence) then
      fit%status = fit%ica%status
      fit%message = fit%ica%message
      return
    end if
    fit%factor_series = fit%ica%components
    allocate(fit%factors(m), variances(n, m), skewness(m), kurtosis(m))
    opt%compute_inference = .false.
    if (present(options)) opt = options
    fit%log_likelihood = 0.0_dp
    do j = 1, m
      fit%factors(j) = estimate_garch(fit%factor_series(:, j), specification%factor_spec, options=opt)
      if (fit%factors(j)%status /= 0) then
        fit%status = tsm_numerical_failure
        fit%message = 'a factor GARCH fit failed'
        return
      end if
      variances(:, j) = fit%factors(j)%filtered%variance
      skewness(j) = dskewness(fit%factors(j)%spec%distribution, fit%factors(j)%parameters%dist%skew, &
        fit%factors(j)%parameters%dist%shape, fit%factors(j)%parameters%dist%lambda)
      kurtosis(j) = dkurtosis(fit%factors(j)%spec%distribution, fit%factors(j)%parameters%dist%skew, &
        fit%factors(j)%parameters%dist%shape, fit%factors(j)%parameters%dist%lambda)
      if (.not. finite_value(skewness(j))) skewness(j) = 0.0_dp
      if (.not. finite_value(kurtosis(j))) kurtosis(j) = 0.0_dp
      fit%log_likelihood = fit%log_likelihood + fit%factors(j)%log_likelihood
    end do
    allocate(fit%covariance(m, m, n), fit%correlation(m, m, n))
    do t = 1, n
      fit%covariance(:, :, t) = gogarch_covariance(variances(t, :), fit%ica%mixing)
      fit%correlation(:, :, t) = covariance_to_correlation(fit%covariance(:, :, t))
    end do
    fit%coskewness = gogarch_coskewness(fit%ica%mixing, skewness, variances(n, :))
    fit%cokurtosis = gogarch_cokurtosis(fit%ica%mixing, kurtosis, variances(n, :))
    fit%status = tsm_success
    fit%message = 'ok'
  end function estimate_gogarch

  function forecast_gogarch(data, fit, horizon, paths, seed) result(out)
    real(dp), intent(in) :: data(:, :)
    type(gogarch_fit), intent(in) :: fit
    integer, intent(in) :: horizon
    integer, intent(in), optional :: paths
    integer(i8), intent(in), optional :: seed
    type(gogarch_forecast) :: out
    type(garch_forecast) :: factor_forecast
    real(dp), allocatable :: factor_mean(:, :), z(:, :), component_draw(:), asset_draw(:)
    integer :: m, h, j, p, npaths

    if (fit%status /= tsm_success .or. horizon < 1) then
      out%status = tsm_invalid_argument
      out%message = 'a successful GO-GARCH fit and positive horizon are required'
      return
    end if
    m = size(data, 2)
    npaths = 0
    if (present(paths)) npaths = max(0, paths)
    allocate(out%factor_variance(horizon, m), factor_mean(horizon, m))
    do j = 1, m
      factor_forecast = forecast_garch(fit%factor_series(:, j), fit%factors(j), horizon, &
        paths=max(100, npaths), seed=int(1009 + 17 * j, i8))
      if (factor_forecast%status /= 0) then
        out%status = tsm_numerical_failure
        out%message = 'a factor forecast failed'
        return
      end if
      out%factor_variance(:, j) = factor_forecast%variance
      factor_mean(:, j) = factor_forecast%mean
    end do
    allocate(out%covariance(m, m, horizon), out%correlation(m, m, horizon))
    do h = 1, horizon
      out%covariance(:, :, h) = gogarch_covariance(out%factor_variance(h, :), fit%ica%mixing)
      out%correlation(:, :, h) = covariance_to_correlation(out%covariance(:, :, h))
    end do
    if (npaths > 0) then
      if (present(seed)) call set_random_seed(seed)
      allocate(out%simulated(horizon, m, npaths), component_draw(m), asset_draw(m))
      do p = 1, npaths
        do h = 1, horizon
          z = random_normal_matrix(m, 1)
          component_draw = factor_mean(h, :) + sqrt(max(out%factor_variance(h, :), 0.0_dp)) * z(:, 1)
          asset_draw = fit%ica%mean + matmul(fit%ica%mixing, component_draw)
          out%simulated(h, :, p) = asset_draw
        end do
      end do
    end if
    out%status = tsm_success
    out%message = 'ok'
  end function forecast_gogarch

  function gogarch_covariance(factor_variance, mixing) result(covariance)
    real(dp), intent(in) :: factor_variance(:), mixing(:, :)
    real(dp), allocatable :: covariance(:, :), d(:, :)
    integer :: i, m
    m = size(mixing, 1)
    allocate(d(size(factor_variance), size(factor_variance)))
    d = 0.0_dp
    do i = 1, size(factor_variance)
      d(i, i) = max(factor_variance(i), 0.0_dp)
    end do
    covariance = matmul(mixing, matmul(d, transpose(mixing)))
  end function gogarch_covariance

  function gogarch_correlation(factor_variance, mixing) result(correlation)
    real(dp), intent(in) :: factor_variance(:), mixing(:, :)
    real(dp), allocatable :: correlation(:, :)
    correlation = covariance_to_correlation(gogarch_covariance(factor_variance, mixing))
  end function gogarch_correlation

  function gogarch_coskewness(mixing, factor_skewness, factor_variance) result(tensor)
    real(dp), intent(in) :: mixing(:, :), factor_skewness(:), factor_variance(:)
    real(dp), allocatable :: tensor(:, :, :)
    real(dp) :: moment3
    integer :: a, b, c, k, m
    m = size(mixing, 1)
    allocate(tensor(m, m, m))
    tensor = 0.0_dp
    do k = 1, size(mixing, 2)
      moment3 = factor_skewness(k) * max(factor_variance(k), 0.0_dp) ** 1.5_dp
      do c = 1, m
        do b = 1, m
          do a = 1, m
            tensor(a, b, c) = tensor(a, b, c) + mixing(a, k) * mixing(b, k) * mixing(c, k) * moment3
          end do
        end do
      end do
    end do
  end function gogarch_coskewness

  function gogarch_cokurtosis(mixing, factor_excess_kurtosis, factor_variance) result(tensor)
    real(dp), intent(in) :: mixing(:, :), factor_excess_kurtosis(:), factor_variance(:)
    real(dp), allocatable :: tensor(:, :, :, :), covariance(:, :)
    real(dp) :: cumulant4
    integer :: a, b, c, d, k, m
    m = size(mixing, 1)
    covariance = gogarch_covariance(factor_variance, mixing)
    allocate(tensor(m, m, m, m))
    do d = 1, m
      do c = 1, m
        do b = 1, m
          do a = 1, m
            tensor(a, b, c, d) = covariance(a, b) * covariance(c, d) + &
              covariance(a, c) * covariance(b, d) + covariance(a, d) * covariance(b, c)
          end do
        end do
      end do
    end do
    do k = 1, size(mixing, 2)
      cumulant4 = factor_excess_kurtosis(k) * max(factor_variance(k), 0.0_dp) ** 2
      do d = 1, m
        do c = 1, m
          do b = 1, m
            do a = 1, m
              tensor(a, b, c, d) = tensor(a, b, c, d) + mixing(a, k) * mixing(b, k) * &
                mixing(c, k) * mixing(d, k) * cumulant4
            end do
          end do
        end do
      end do
    end do
  end function gogarch_cokurtosis

  function portfolio_variance(covariance, weights) result(value)
    real(dp), intent(in) :: covariance(:, :), weights(:)
    real(dp) :: value
    value = dot_product(weights, matmul(covariance, weights))
  end function portfolio_variance

  function portfolio_skewness(coskewness, covariance, weights) result(value)
    real(dp), intent(in) :: coskewness(:, :, :), covariance(:, :), weights(:)
    real(dp) :: value, m3, variance
    integer :: a, b, c, m
    m = size(weights)
    m3 = 0.0_dp
    do c = 1, m
      do b = 1, m
        do a = 1, m
          m3 = m3 + weights(a) * weights(b) * weights(c) * coskewness(a, b, c)
        end do
      end do
    end do
    variance = portfolio_variance(covariance, weights)
    value = m3 / max(variance, tiny(1.0_dp)) ** 1.5_dp
  end function portfolio_skewness

  function portfolio_kurtosis(cokurtosis, covariance, weights) result(value)
    real(dp), intent(in) :: cokurtosis(:, :, :, :), covariance(:, :), weights(:)
    real(dp) :: value, m4, variance
    integer :: a, b, c, d, m
    m = size(weights)
    m4 = 0.0_dp
    do d = 1, m
      do c = 1, m
        do b = 1, m
          do a = 1, m
            m4 = m4 + weights(a) * weights(b) * weights(c) * weights(d) * cokurtosis(a, b, c, d)
          end do
        end do
      end do
    end do
    variance = portfolio_variance(covariance, weights)
    value = m4 / max(variance, tiny(1.0_dp)) ** 2
  end function portfolio_kurtosis

  logical function finite_value(x)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp), intent(in) :: x
    finite_value = ieee_is_finite(x)
  end function finite_value

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

end module tsmarch_gogarch
