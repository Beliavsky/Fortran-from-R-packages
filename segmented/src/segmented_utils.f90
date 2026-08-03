! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp, pi_dp
  use nlme_linalg, only : solve_least_squares, inverse_spd
  use nlme_status, only : NLME_SUCCESS
  use segmented_status
  use segmented_types, only : FAMILY_GAUSSIAN, FAMILY_BINOMIAL, FAMILY_POISSON
  implicit none
  private
  public :: quantile_value, breakpoint_limits, hinge_matrix, step_matrix
  public :: make_design, weighted_lm_fit, glm_fit_fixed, normal_cdf
  public :: normal_quantile, sort_real, clamp_breakpoints, same_column
  public :: gaussian_loglik, glm_deviance
contains
  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_real

  function quantile_value(x, probability) result(value)
    real(dp), intent(in) :: x(:), probability
    real(dp) :: value, h, frac
    real(dp), allocatable :: work(:)
    integer :: lo, hi
    if (size(x) == 0) then
      value = 0.0_dp
      return
    end if
    allocate(work(size(x)))
    work = x
    call sort_real(work)
    h = 1.0_dp + real(size(x) - 1, dp) * min(1.0_dp, max(0.0_dp, probability))
    lo = max(1, min(size(x), int(floor(h))))
    hi = max(1, min(size(x), lo + 1))
    frac = h - real(lo, dp)
    value = (1.0_dp - frac) * work(lo) + frac * work(hi)
  end function quantile_value

  subroutine breakpoint_limits(z, lower_q, upper_q, lower, upper, status)
    real(dp), intent(in) :: z(:,:), lower_q, upper_q
    real(dp), allocatable, intent(out) :: lower(:), upper(:)
    integer, intent(out) :: status
    integer :: j
    if (size(z, 1) < 3 .or. size(z, 2) < 1 .or. lower_q < 0.0_dp .or. &
        upper_q > 1.0_dp .or. lower_q >= upper_q) then
      allocate(lower(0), upper(0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    allocate(lower(size(z, 2)), upper(size(z, 2)))
    do j = 1, size(z, 2)
      lower(j) = quantile_value(z(:, j), lower_q)
      upper(j) = quantile_value(z(:, j), upper_q)
      if (.not. ieee_is_finite(lower(j)) .or. .not. ieee_is_finite(upper(j)) .or. &
          lower(j) >= upper(j)) then
        status = SEG_INFEASIBLE_BREAKPOINT
        return
      end if
    end do
    status = SEG_SUCCESS
  end subroutine breakpoint_limits

  pure logical function same_column(a, b)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: scale
    if (size(a) /= size(b)) then
      same_column = .false.
      return
    end if
    scale = max(1.0_dp, maxval(abs(a)), maxval(abs(b)))
    same_column = maxval(abs(a - b)) <= 100.0_dp * epsilon(1.0_dp) * scale
  end function same_column

  subroutine clamp_breakpoints(z, lower, upper, psi)
    real(dp), intent(in) :: z(:,:), lower(:), upper(:)
    real(dp), intent(inout) :: psi(:)
    integer :: i, j, k
    real(dp) :: temp, gap
    psi = min(upper, max(lower, psi))
    do i = 1, size(psi)
      do j = i + 1, size(psi)
        if (same_column(z(:, i), z(:, j))) then
          if (psi(i) > psi(j)) then
            temp = psi(i)
            psi(i) = psi(j)
            psi(j) = temp
          end if
        end if
      end do
    end do
    do i = 1, size(psi)
      do j = i + 1, size(psi)
        if (same_column(z(:, i), z(:, j))) then
          gap = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, upper(i) - lower(i))
          if (psi(j) <= psi(i) + gap) psi(j) = min(upper(j), psi(i) + gap)
        end if
      end do
    end do
    do k = 1, size(psi)
      psi(k) = min(upper(k), max(lower(k), psi(k)))
    end do
  end subroutine clamp_breakpoints

  subroutine hinge_matrix(z, psi, u, v)
    real(dp), intent(in) :: z(:,:), psi(:)
    real(dp), allocatable, intent(out) :: u(:,:)
    real(dp), allocatable, intent(out), optional :: v(:,:)
    integer :: i, j
    allocate(u(size(z, 1), size(z, 2)))
    if (present(v)) allocate(v(size(z, 1), size(z, 2)))
    do j = 1, size(z, 2)
      do i = 1, size(z, 1)
        if (z(i, j) > psi(j)) then
          u(i, j) = z(i, j) - psi(j)
          if (present(v)) v(i, j) = -1.0_dp
        else
          u(i, j) = 0.0_dp
          if (present(v)) v(i, j) = 0.0_dp
        end if
      end do
    end do
  end subroutine hinge_matrix

  subroutine step_matrix(z, psi, h)
    real(dp), intent(in) :: z(:,:), psi(:)
    real(dp), allocatable, intent(out) :: h(:,:)
    integer :: i, j
    allocate(h(size(z, 1), size(z, 2)))
    do j = 1, size(z, 2)
      do i = 1, size(z, 1)
        h(i, j) = merge(1.0_dp, 0.0_dp, z(i, j) > psi(j))
      end do
    end do
  end subroutine step_matrix

  subroutine make_design(x, z, psi, kind, design, derivative)
    real(dp), intent(in) :: x(:,:), z(:,:), psi(:)
    integer, intent(in) :: kind
    real(dp), allocatable, intent(out) :: design(:,:)
    real(dp), allocatable, intent(out), optional :: derivative(:,:)
    real(dp), allocatable :: term(:,:), v(:,:)
    integer :: p, m
    p = size(x, 2)
    m = size(z, 2)
    if (kind == 1) then
      if (present(derivative)) then
        call hinge_matrix(z, psi, term, v)
      else
        call hinge_matrix(z, psi, term)
      end if
    else
      call step_matrix(z, psi, term)
      if (present(derivative)) allocate(v(size(z, 1), m), source = 0.0_dp)
    end if
    allocate(design(size(x, 1), p + m))
    design(:, :p) = x
    design(:, p + 1:) = term
    if (present(derivative)) call move_alloc(v, derivative)
  end subroutine make_design

  subroutine weighted_lm_fit(y, design, beta, covariance, fitted, residuals, rss, sigma, &
      loglik, status, weights, offset)
    real(dp), intent(in) :: y(:), design(:,:)
    real(dp), allocatable, intent(out) :: beta(:), covariance(:,:), fitted(:), residuals(:)
    real(dp), intent(out) :: rss, sigma, loglik
    integer, intent(out) :: status
    real(dp), intent(in), optional :: weights(:), offset(:)
    real(dp), allocatable :: yadj(:), cov0(:,:), w(:), off(:)
    real(dp) :: df, sumlogw
    integer :: n, p, st
    n = size(y)
    p = size(design, 2)
    if (size(design, 1) /= n .or. n <= p) then
      status = SEG_DIMENSION_ERROR
      return
    end if
    allocate(w(n), off(n), yadj(n))
    w = 1.0_dp
    off = 0.0_dp
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
        status = SEG_INVALID_ARGUMENT
        return
      end if
      w = weights
    end if
    if (present(offset)) then
      if (size(offset) /= n) then
        status = SEG_DIMENSION_ERROR
        return
      end if
      off = offset
    end if
    yadj = y - off
    call solve_least_squares(design, yadj, beta, cov0, rss, st, w)
    if (st /= NLME_SUCCESS) then
      status = SEG_SINGULAR
      return
    end if
    allocate(fitted(n), residuals(n), covariance(p, p))
    fitted = off + matmul(design, beta)
    residuals = y - fitted
    rss = sum(w * residuals**2)
    df = real(max(1, n - p), dp)
    sigma = sqrt(max(0.0_dp, rss / df))
    covariance = cov0 * sigma**2
    sumlogw = sum(log(w))
    loglik = gaussian_loglik(rss, n, sumlogw)
    status = SEG_SUCCESS
  end subroutine weighted_lm_fit

  pure function gaussian_loglik(rss, n, sumlogw) result(value)
    real(dp), intent(in) :: rss, sumlogw
    integer, intent(in) :: n
    real(dp) :: value, s2
    s2 = max(tiny(1.0_dp), rss / real(max(1, n), dp))
    value = -0.5_dp * (real(n, dp) * (log(2.0_dp * pi_dp * s2) + 1.0_dp) - sumlogw)
  end function gaussian_loglik

  subroutine glm_fit_fixed(y, design, family, beta, covariance, fitted, residuals, &
      objective, loglik, status, prior_weights, offset, max_iter, tolerance)
    real(dp), intent(in) :: y(:), design(:,:)
    integer, intent(in) :: family
    real(dp), allocatable, intent(out) :: beta(:), covariance(:,:), fitted(:), residuals(:)
    real(dp), intent(out) :: objective, loglik
    integer, intent(out) :: status
    real(dp), intent(in), optional :: prior_weights(:), offset(:)
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: b(:), bnew(:), cov0(:,:), eta(:), mu(:), varmu(:), dmu(:)
    real(dp), allocatable :: worky(:), workw(:), pw(:), off(:)
    real(dp) :: rss, sigma, tol, diff, scale
    integer :: n, p, iter, itmax, st
    n = size(y)
    p = size(design, 2)
    if (size(design, 1) /= n .or. n <= p) then
      status = SEG_DIMENSION_ERROR
      return
    end if
    if (family == FAMILY_GAUSSIAN) then
      call weighted_lm_fit(y, design, beta, covariance, fitted, residuals, rss, sigma, &
          loglik, status, prior_weights, offset)
      objective = rss
      return
    end if
    if (family /= FAMILY_BINOMIAL .and. family /= FAMILY_POISSON) then
      status = SEG_INVALID_ARGUMENT
      return
    end if
    if (family == FAMILY_BINOMIAL) then
      if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
        status = SEG_INVALID_ARGUMENT
        return
      end if
    else if (any(y < 0.0_dp)) then
      status = SEG_INVALID_ARGUMENT
      return
    end if
    allocate(pw(n), off(n), eta(n), mu(n), varmu(n), dmu(n), worky(n), workw(n))
    allocate(b(p), source = 0.0_dp)
    pw = 1.0_dp
    off = 0.0_dp
    if (present(prior_weights)) then
      if (size(prior_weights) /= n .or. any(prior_weights <= 0.0_dp)) then
        status = SEG_INVALID_ARGUMENT
        return
      end if
      pw = prior_weights
    end if
    if (present(offset)) then
      if (size(offset) /= n) then
        status = SEG_DIMENSION_ERROR
        return
      end if
      off = offset
    end if
    if (family == FAMILY_BINOMIAL) then
      mu = min(1.0_dp - 1.0e-4_dp, max(1.0e-4_dp, (y + 0.5_dp) / 2.0_dp))
      eta = log(mu / (1.0_dp - mu))
    else
      mu = max(0.1_dp, y + 0.1_dp)
      eta = log(mu)
    end if
    itmax = 100
    if (present(max_iter)) itmax = max_iter
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = tolerance
    do iter = 1, itmax
      if (family == FAMILY_BINOMIAL) then
        mu = 1.0_dp / (1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, eta))))
        varmu = max(1.0e-10_dp, mu * (1.0_dp - mu))
        dmu = varmu
      else
        mu = exp(max(-35.0_dp, min(35.0_dp, eta)))
        varmu = max(1.0e-10_dp, mu)
        dmu = mu
      end if
      workw = pw * dmu**2 / varmu
      worky = eta + (y - mu) / dmu
      call solve_least_squares(design, worky - off, bnew, cov0, rss, st, workw)
      if (st /= NLME_SUCCESS) then
        status = SEG_SINGULAR
        return
      end if
      scale = max(1.0_dp, maxval(abs(b)))
      diff = maxval(abs(bnew - b)) / scale
      b = bnew
      eta = off + matmul(design, b)
      if (diff <= tol) exit
    end do
    beta = b
    allocate(fitted(n), residuals(n), covariance(p, p))
    if (family == FAMILY_BINOMIAL) then
      fitted = 1.0_dp / (1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, eta))))
      varmu = max(1.0e-10_dp, fitted * (1.0_dp - fitted))
      dmu = varmu
    else
      fitted = exp(max(-35.0_dp, min(35.0_dp, eta)))
      varmu = max(1.0e-10_dp, fitted)
      dmu = fitted
    end if
    residuals = y - fitted
    workw = pw * dmu**2 / varmu
    call solve_least_squares(design, eta - off, bnew, cov0, rss, st, workw)
    if (st /= NLME_SUCCESS) then
      status = SEG_SINGULAR
      return
    end if
    covariance = cov0
    objective = glm_deviance(y, fitted, family, pw)
    if (family == FAMILY_BINOMIAL) then
      loglik = sum(pw * (y * log(max(tiny(1.0_dp), fitted)) + &
          (1.0_dp - y) * log(max(tiny(1.0_dp), 1.0_dp - fitted))))
    else
      loglik = sum(pw * (y * log(max(tiny(1.0_dp), fitted)) - fitted - log_gamma(y + 1.0_dp)))
    end if
    if (iter > itmax) then
      status = SEG_MAX_ITER
    else
      status = SEG_SUCCESS
    end if
  end subroutine glm_fit_fixed

  function glm_deviance(y, mu, family, weights) result(value)
    real(dp), intent(in) :: y(:), mu(:), weights(:)
    integer, intent(in) :: family
    real(dp) :: value, term
    integer :: i
    value = 0.0_dp
    select case (family)
    case (FAMILY_GAUSSIAN)
      value = sum(weights * (y - mu)**2)
    case (FAMILY_BINOMIAL)
      do i = 1, size(y)
        term = 0.0_dp
        if (y(i) > 0.0_dp) term = term + y(i) * log(y(i) / max(tiny(1.0_dp), mu(i)))
        if (y(i) < 1.0_dp) term = term + (1.0_dp - y(i)) * &
            log((1.0_dp - y(i)) / max(tiny(1.0_dp), 1.0_dp - mu(i)))
        value = value + 2.0_dp * weights(i) * term
      end do
    case (FAMILY_POISSON)
      do i = 1, size(y)
        if (y(i) > 0.0_dp) then
          term = y(i) * log(y(i) / max(tiny(1.0_dp), mu(i))) - (y(i) - mu(i))
        else
          term = mu(i)
        end if
        value = value + 2.0_dp * weights(i) * term
      end do
    end select
  end function glm_deviance

  pure elemental function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, &
      2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, &
      2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
          ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp * log(1.0_dp - p))
      x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
          ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    else
      q = p - 0.5_dp
      r = q * q
      x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
          (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
    end if
  end function normal_quantile
end module segmented_utils
