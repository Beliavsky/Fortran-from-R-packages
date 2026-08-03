! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_inference
  use nlme_kinds, only : dp
  use nlme_linalg, only : solve_least_squares
  use nlme_status, only : NLME_SUCCESS
  use segmented_status
  use segmented_types
  use segmented_utils
  use segmented_fit, only : fit_segmented_lm, fit_segmented_glm
  implicit none
  private
  public :: predict_segmented, segment_slopes, segment_intercepts
  public :: breakpoint_confint, aapc, broken_line_values
  public :: davies_test, pscore_test, pwr_seg, select_breakpoints_bic
contains
  subroutine predict_segmented(fit, xnew, znew, prediction, status, linear_predictor)
    type(segmented_result), intent(in) :: fit
    real(dp), intent(in) :: xnew(:,:), znew(:,:)
    real(dp), allocatable, intent(out) :: prediction(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: linear_predictor
    real(dp), allocatable :: design(:,:)
    logical :: linear
    if (.not. allocated(fit%breakpoints) .or. .not. allocated(fit%coefficients)) then
      allocate(prediction(0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    if (size(xnew, 2) /= fit%n_base .or. size(znew, 2) /= fit%n_break .or. &
        size(xnew, 1) /= size(znew, 1)) then
      allocate(prediction(0))
      status = SEG_DIMENSION_ERROR
      return
    end if
    call make_design(xnew, znew, fit%breakpoints, fit%kind, design)
    allocate(prediction(size(xnew, 1)))
    prediction = matmul(design, fit%coefficients)
    linear = .false.
    if (present(linear_predictor)) linear = linear_predictor
    if (.not. linear) then
      select case (fit%family)
      case (FAMILY_BINOMIAL)
        prediction = 1.0_dp / (1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, prediction))))
      case (FAMILY_POISSON)
        prediction = exp(max(-35.0_dp, min(35.0_dp, prediction)))
      end select
    end if
    status = SEG_SUCCESS
  end subroutine predict_segmented

  subroutine segment_slopes(fit, base_coefficient, slopes, standard_errors, status)
    type(segmented_result), intent(in) :: fit
    integer, intent(in) :: base_coefficient
    real(dp), allocatable, intent(out) :: slopes(:), standard_errors(:)
    integer, intent(out) :: status
    real(dp), allocatable :: contrast(:)
    integer :: k, j, idx
    if (fit%kind /= SEGMENTED_CONTINUOUS .or. base_coefficient < 1 .or. &
        base_coefficient > fit%n_base .or. .not. allocated(fit%covariance)) then
      allocate(slopes(0), standard_errors(0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    k = fit%n_break
    allocate(slopes(k + 1), standard_errors(k + 1), contrast(size(fit%coefficients)))
    do j = 0, k
      contrast = 0.0_dp
      contrast(base_coefficient) = 1.0_dp
      if (j > 0) then
        do idx = 1, j
          contrast(fit%n_base + idx) = 1.0_dp
        end do
      end if
      slopes(j + 1) = dot_product(contrast, fit%coefficients)
      standard_errors(j + 1) = sqrt(max(0.0_dp, &
          dot_product(contrast, matmul(fit%covariance, contrast))))
    end do
    status = SEG_SUCCESS
  end subroutine segment_slopes

  subroutine segment_intercepts(fit, intercept_coefficient, base_slope_coefficient, &
      intercepts, status)
    type(segmented_result), intent(in) :: fit
    integer, intent(in) :: intercept_coefficient, base_slope_coefficient
    real(dp), allocatable, intent(out) :: intercepts(:)
    integer, intent(out) :: status
    integer :: j, k
    if (fit%kind /= SEGMENTED_CONTINUOUS .or. intercept_coefficient < 1 .or. &
        intercept_coefficient > fit%n_base .or. base_slope_coefficient < 1 .or. &
        base_slope_coefficient > fit%n_base) then
      allocate(intercepts(0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    k = fit%n_break
    allocate(intercepts(k + 1))
    intercepts(1) = fit%coefficients(intercept_coefficient)
    do j = 1, k
      intercepts(j + 1) = intercepts(j) - &
          fit%coefficients(fit%n_base + j) * fit%breakpoints(j)
    end do
    status = SEG_SUCCESS
  end subroutine segment_intercepts

  subroutine breakpoint_confint(fit, confidence, interval, status)
    type(segmented_result), intent(in) :: fit
    real(dp), intent(in) :: confidence
    real(dp), allocatable, intent(out) :: interval(:,:)
    integer, intent(out) :: status
    real(dp) :: zcrit
    if (confidence <= 0.0_dp .or. confidence >= 1.0_dp .or. &
        .not. allocated(fit%breakpoint_se)) then
      allocate(interval(0, 0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    zcrit = normal_quantile(0.5_dp + 0.5_dp * confidence)
    allocate(interval(size(fit%breakpoints), 2))
    interval(:, 1) = fit%breakpoints - zcrit * fit%breakpoint_se
    interval(:, 2) = fit%breakpoints + zcrit * fit%breakpoint_se
    status = SEG_SUCCESS
  end subroutine breakpoint_confint

  function aapc(fit, x_min, x_max, base_slope_coefficient, exponentiate, status) result(value)
    type(segmented_result), intent(in) :: fit
    real(dp), intent(in) :: x_min, x_max
    integer, intent(in) :: base_slope_coefficient
    logical, intent(in), optional :: exponentiate
    integer, intent(out), optional :: status
    real(dp) :: value, total, left, right
    real(dp), allocatable :: slopes(:), ses(:)
    integer :: st, j
    logical :: expit
    call segment_slopes(fit, base_slope_coefficient, slopes, ses, st)
    if (st /= SEG_SUCCESS .or. x_max <= x_min) then
      value = huge(1.0_dp)
      if (present(status)) status = SEG_INVALID_ARGUMENT
      return
    end if
    total = 0.0_dp
    left = x_min
    do j = 1, fit%n_break
      right = min(x_max, max(x_min, fit%breakpoints(j)))
      if (right > left) total = total + slopes(j) * (right - left)
      left = max(left, right)
    end do
    if (x_max > left) total = total + slopes(size(slopes)) * (x_max - left)
    value = total / (x_max - x_min)
    expit = .false.
    if (present(exponentiate)) expit = exponentiate
    if (expit) value = exp(value) - 1.0_dp
    if (present(status)) status = SEG_SUCCESS
  end function aapc

  pure function broken_line_values(x, intercept, base_slope, changes, breakpoints) result(y)
    real(dp), intent(in) :: x, intercept, base_slope
    real(dp), intent(in) :: changes(:), breakpoints(:)
    real(dp) :: y
    integer :: j
    y = intercept + base_slope * x
    do j = 1, min(size(changes), size(breakpoints))
      y = y + changes(j) * max(0.0_dp, x - breakpoints(j))
    end do
  end function broken_line_values

  subroutine davies_test(y, x, z, grid_points, result, weights, offset)
    real(dp), intent(in) :: y(:), x(:,:), z(:)
    integer, intent(in), optional :: grid_points
    type(test_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    integer :: ng, g, n, p, status
    real(dp), allocatable :: design(:,:), beta(:), cov(:,:), fitted(:), residuals(:)
    real(dp), allocatable :: zz(:,:), psi(:), w(:)
    real(dp) :: lo, hi, candidate, rss, sigma, loglik, tstat, best
    n = size(y)
    p = size(x, 2)
    if (size(x, 1) /= n .or. size(z) /= n .or. n <= p + 1) then
      result%status = SEG_DIMENSION_ERROR
      return
    end if
    ng = 25
    if (present(grid_points)) ng = max(5, grid_points)
    lo = quantile_value(z, 0.05_dp)
    hi = quantile_value(z, 0.95_dp)
    allocate(w(n), zz(n, 1), psi(1))
    w = 1.0_dp
    if (present(weights)) w = weights
    zz(:, 1) = z
    best = 0.0_dp
    result%breakpoint = 0.5_dp * (lo + hi)
    do g = 1, ng
      candidate = lo + real(g - 1, dp) / real(max(1, ng - 1), dp) * (hi - lo)
      psi(1) = candidate
      call make_design(x, zz, psi, SEGMENTED_CONTINUOUS, design)
      call weighted_lm_fit(y, design, beta, cov, fitted, residuals, rss, sigma, &
          loglik, status, w, offset)
      if (status /= SEG_SUCCESS) cycle
      tstat = abs(beta(p + 1)) / sqrt(max(tiny(1.0_dp), cov(p + 1, p + 1)))
      if (tstat > best) then
        best = tstat
        result%breakpoint = candidate
      end if
    end do
    result%statistic = best
    result%grid_evaluated = ng
    result%p_value = min(1.0_dp, 2.0_dp * real(ng, dp) * (1.0_dp - normal_cdf(best)))
    result%status = SEG_SUCCESS
  end subroutine davies_test

  subroutine pscore_test(y, x, z, breakpoint, result, weights, offset)
    real(dp), intent(in) :: y(:), x(:,:), z(:), breakpoint
    type(test_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    real(dp), allocatable :: design(:,:), beta(:), cov(:,:), fitted(:), residuals(:)
    real(dp), allocatable :: zz(:,:), psi(:), w(:)
    real(dp) :: rss, sigma, loglik
    integer :: n, p, status
    n = size(y)
    p = size(x, 2)
    if (size(x, 1) /= n .or. size(z) /= n .or. n <= p + 1) then
      result%status = SEG_DIMENSION_ERROR
      return
    end if
    allocate(w(n), zz(n, 1), psi(1))
    w = 1.0_dp
    if (present(weights)) w = weights
    zz(:, 1) = z
    psi(1) = breakpoint
    call make_design(x, zz, psi, SEGMENTED_CONTINUOUS, design)
    call weighted_lm_fit(y, design, beta, cov, fitted, residuals, rss, sigma, &
        loglik, status, w, offset)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    result%statistic = beta(p + 1) / sqrt(max(tiny(1.0_dp), cov(p + 1, p + 1)))
    result%p_value = 2.0_dp * (1.0_dp - normal_cdf(abs(result%statistic)))
    result%breakpoint = breakpoint
    result%grid_evaluated = 1
    result%status = SEG_SUCCESS
  end subroutine pscore_test

  function pwr_seg(z, x, breakpoint, slope_change, sigma, alpha, two_sided, status) result(power)
    real(dp), intent(in) :: z(:), x(:,:), breakpoint, slope_change, sigma
    real(dp), intent(in), optional :: alpha
    logical, intent(in), optional :: two_sided
    integer, intent(out), optional :: status
    real(dp) :: power, a, zcrit, ncp, rss
    real(dp), allocatable :: u(:), beta(:), cov(:,:), residuals(:)
    integer :: st, n
    logical :: two
    n = size(z)
    if (size(x, 1) /= n .or. sigma <= 0.0_dp .or. n <= size(x, 2)) then
      power = 0.0_dp
      if (present(status)) status = SEG_INVALID_ARGUMENT
      return
    end if
    allocate(u(n))
    u = max(0.0_dp, z - breakpoint)
    call solve_least_squares(x, u, beta, cov, rss, st)
    if (st /= NLME_SUCCESS) then
      power = 0.0_dp
      if (present(status)) status = SEG_SINGULAR
      return
    end if
    residuals = u - matmul(x, beta)
    ncp = abs(slope_change) * sqrt(sum(residuals**2)) / sigma
    a = 0.05_dp
    if (present(alpha)) a = alpha
    two = .true.
    if (present(two_sided)) two = two_sided
    if (two) then
      zcrit = normal_quantile(1.0_dp - 0.5_dp * a)
      power = normal_cdf(-zcrit - ncp) + 1.0_dp - normal_cdf(zcrit - ncp)
    else
      zcrit = normal_quantile(1.0_dp - a)
      power = 1.0_dp - normal_cdf(zcrit - ncp)
    end if
    if (present(status)) status = SEG_SUCCESS
  end function pwr_seg

  subroutine select_breakpoints_bic(y, x, z, max_breaks, family, best_fit, bic_values, &
      weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:)
    integer, intent(in) :: max_breaks, family
    type(segmented_result), intent(out) :: best_fit
    real(dp), allocatable, intent(out) :: bic_values(:)
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    type(segmented_result) :: fit
    real(dp), allocatable :: zz(:,:), psi(:), beta(:), cov(:,:)
    real(dp), allocatable :: fitted(:), residuals(:), w(:)
    real(dp) :: objective, sigma, loglik, best_bic
    integer :: k, j, status, n, p
    ctl = segmented_control()
    if (present(control)) ctl = control
    n = size(y)
    p = size(x, 2)
    if (size(x, 1) /= n .or. size(z) /= n .or. max_breaks < 1) then
      allocate(bic_values(0))
      best_fit%status = SEG_INVALID_ARGUMENT
      return
    end if
    allocate(bic_values(0:max_breaks), w(n))
    w = 1.0_dp
    if (present(weights)) w = weights
    if (family == FAMILY_GAUSSIAN) then
      call weighted_lm_fit(y, x, beta, cov, fitted, residuals, objective, sigma, &
          loglik, status, w, offset)
    else
      call glm_fit_fixed(y, x, family, beta, cov, fitted, residuals, objective, &
          loglik, status, w, offset, ctl%glm_max_iter, ctl%tolerance)
    end if
    if (status /= SEG_SUCCESS .and. status /= SEG_MAX_ITER) then
      best_fit%status = status
      return
    end if
    bic_values(0) = -2.0_dp * loglik + log(real(n, dp)) * real(p, dp)
    best_bic = bic_values(0)
    best_fit%coefficients = beta
    best_fit%covariance = cov
    best_fit%fitted = fitted
    best_fit%residuals = residuals
    best_fit%weights = w
    allocate(best_fit%breakpoints(0), best_fit%breakpoint_se(0))
    best_fit%objective = objective
    best_fit%log_likelihood = loglik
    best_fit%aic = -2.0_dp * loglik + 2.0_dp * real(p, dp)
    best_fit%bic = bic_values(0)
    best_fit%n_base = p
    best_fit%n_break = 0
    best_fit%family = family
    best_fit%kind = SEGMENTED_CONTINUOUS
    best_fit%status = SEG_SUCCESS
    best_fit%converged = .true.
    do k = 1, max_breaks
      allocate(zz(n, k), psi(k))
      do j = 1, k
        zz(:, j) = z
        psi(j) = quantile_value(z, real(j, dp) / real(k + 1, dp))
      end do
      if (family == FAMILY_GAUSSIAN) then
        call fit_segmented_lm(y, x, zz, psi, fit, weights, offset, ctl)
      else
        call fit_segmented_glm(y, x, zz, psi, family, fit, weights, offset, ctl)
      end if
      if (fit%status == SEG_SUCCESS .or. fit%status == SEG_MAX_ITER) then
        bic_values(k) = fit%bic
        if (fit%bic < best_bic) then
          best_bic = fit%bic
          best_fit = fit
        end if
      else
        bic_values(k) = huge(1.0_dp)
      end if
      deallocate(zz, psi)
    end do
  end subroutine select_breakpoints_bic
end module segmented_inference
