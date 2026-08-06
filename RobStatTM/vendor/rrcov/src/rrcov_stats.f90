! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_stats
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success, rrcov_invalid_argument, rrcov_dimension_error
  use rrcov_sort, only : sort_real
  use rrcov_linalg, only : mahalanobis_squared
  implicit none
  private
  public :: mean_vector, covariance_matrix, weighted_mean_covariance
  public :: median, quantile, mad_scale, qn_scale, sn_scale, tau_scale
  public :: spatial_median, standardize_columns, robust_standardize_columns
  public :: chi_square_cdf, chi_square_quantile, normal_cdf, f_cdf
  public :: regularized_gamma_p, regularized_beta, column_medians
contains
  pure function mean_vector(x) result(value)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: value(size(x, 2))
    if (size(x, 1) > 0) then
      value = sum(x, dim=1) / real(size(x, 1), dp)
    else
      value = 0.0_dp
    end if
  end function mean_vector

  function covariance_matrix(x, unbiased, center, status) result(value)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in), optional :: unbiased
    real(dp), intent(out), optional :: center(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :), centered(:, :), mu(:)
    logical :: use_unbiased
    real(dp) :: denominator
    integer :: n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(value(p, p), mu(p))
    value = 0.0_dp
    mu = mean_vector(x)
    if (present(center)) then
      if (size(center) == p) center = mu
    end if
    if (n < 1 .or. p < 1) then
      if (present(status)) status = rrcov_invalid_argument
      return
    end if
    allocate(centered(n, p))
    centered = x - spread(mu, 1, n)
    use_unbiased = .true.
    if (present(unbiased)) use_unbiased = unbiased
    if (use_unbiased .and. n > 1) then
      denominator = real(n - 1, dp)
    else
      denominator = real(n, dp)
    end if
    value = matmul(transpose(centered), centered) / max(denominator, 1.0_dp)
    value = 0.5_dp * (value + transpose(value))
    if (present(status)) status = rrcov_success
  end function covariance_matrix

  subroutine weighted_mean_covariance(x, weights, center, covariance, status, normalize)
    real(dp), intent(in) :: x(:, :), weights(:)
    real(dp), allocatable, intent(out) :: center(:), covariance(:, :)
    integer, intent(out) :: status
    logical, intent(in), optional :: normalize
    real(dp), allocatable :: delta(:)
    real(dp) :: sumw, denominator
    logical :: normalized
    integer :: i, n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(center(p), covariance(p, p), delta(p))
    center = 0.0_dp
    covariance = 0.0_dp
    if (size(weights) /= n .or. n < 1 .or. p < 1 .or. any(weights < 0.0_dp)) then
      status = rrcov_dimension_error
      return
    end if
    sumw = sum(weights)
    if (sumw <= tiny(1.0_dp)) then
      status = rrcov_invalid_argument
      return
    end if
    do i = 1, n
      center = center + weights(i) * x(i, :)
    end do
    center = center / sumw
    do i = 1, n
      delta = x(i, :) - center
      covariance = covariance + weights(i) * spread(delta, 2, p) * spread(delta, 1, p)
    end do
    normalized = .true.
    if (present(normalize)) normalized = normalize
    if (normalized) then
      denominator = sumw
    else
      denominator = max(sumw - 1.0_dp, 1.0_dp)
    end if
    covariance = covariance / denominator
    covariance = 0.5_dp * (covariance + transpose(covariance))
    status = rrcov_success
  end subroutine weighted_mean_covariance

  function median(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp), allocatable :: work(:)
    integer :: n
    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    work = x
    call sort_real(work)
    if (mod(n, 2) == 1) then
      value = work((n + 1) / 2)
    else
      value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
    end if
  end function median

  function quantile(x, probability) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: probability
    real(dp) :: value, position, fraction
    real(dp), allocatable :: work(:)
    integer :: lower, n
    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    work = x
    call sort_real(work)
    position = 1.0_dp + min(1.0_dp, max(0.0_dp, probability)) * real(n - 1, dp)
    lower = floor(position)
    fraction = position - real(lower, dp)
    if (lower >= n) then
      value = work(n)
    else
      value = (1.0_dp - fraction) * work(lower) + fraction * work(lower + 1)
    end if
  end function quantile

  function mad_scale(x, center, consistency) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out), optional :: center
    logical, intent(in), optional :: consistency
    real(dp) :: value, mu
    logical :: apply_consistency
    mu = median(x)
    value = median(abs(x - mu))
    apply_consistency = .true.
    if (present(consistency)) apply_consistency = consistency
    if (apply_consistency) value = 1.482602218505602_dp * value
    if (present(center)) center = mu
  end function mad_scale

  function qn_scale(x, consistency) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: consistency
    real(dp) :: value, correction
    real(dp), allocatable :: differences(:)
    logical :: apply_consistency
    integer :: n, i, j, k, m, h, order
    n = size(x)
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    m = n * (n - 1) / 2
    allocate(differences(m))
    k = 0
    do j = 2, n
      do i = 1, j - 1
        k = k + 1
        differences(k) = abs(x(j) - x(i))
      end do
    end do
    call sort_real(differences)
    h = n / 2 + 1
    order = h * (h - 1) / 2
    order = min(m, max(1, order))
    value = differences(order)
    apply_consistency = .true.
    if (present(consistency)) apply_consistency = consistency
    if (apply_consistency) then
      correction = qn_finite_correction(n)
      value = 2.2219_dp * correction * value
    end if
  end function qn_scale

  pure function qn_finite_correction(n) result(value)
    integer, intent(in) :: n
    real(dp) :: value
    select case (n)
    case (2); value = 0.399_dp
    case (3); value = 0.994_dp
    case (4); value = 0.512_dp
    case (5); value = 0.844_dp
    case (6); value = 0.611_dp
    case (7); value = 0.857_dp
    case (8); value = 0.669_dp
    case (9); value = 0.872_dp
    case default
      if (mod(n, 2) == 1) then
        value = real(n, dp) / (real(n, dp) + 1.4_dp)
      else
        value = real(n, dp) / (real(n, dp) + 3.8_dp)
      end if
    end select
  end function qn_finite_correction

  function sn_scale(x, consistency) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: consistency
    real(dp) :: value, correction
    real(dp), allocatable :: row_medians(:), differences(:)
    logical :: apply_consistency
    integer :: i, n
    n = size(x)
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    allocate(row_medians(n), differences(n))
    do i = 1, n
      differences = abs(x(i) - x)
      row_medians(i) = median(differences)
    end do
    value = median(row_medians)
    apply_consistency = .true.
    if (present(consistency)) apply_consistency = consistency
    if (apply_consistency) then
      if (n <= 9) then
        select case (n)
        case (2); correction = 0.743_dp
        case (3); correction = 1.851_dp
        case (4); correction = 0.954_dp
        case (5); correction = 1.351_dp
        case (6); correction = 0.993_dp
        case (7); correction = 1.198_dp
        case (8); correction = 1.005_dp
        case (9); correction = 1.131_dp
        end select
      else if (mod(n, 2) == 1) then
        correction = real(n, dp) / (real(n, dp) - 0.9_dp)
      else
        correction = 1.0_dp
      end if
      value = 1.1926_dp * correction * value
    end if
  end function sn_scale

  function tau_scale(x, location) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out), optional :: location
    real(dp), allocatable :: u(:), weights(:)
    real(dp) :: value, med, sigma0, mu, denominator
    integer :: n
    real(dp), parameter :: c1 = 4.5_dp
    real(dp), parameter :: c2_squared = 9.0_dp
    real(dp), parameter :: expectation = 0.9247153921761315_dp
    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      if (present(location)) location = 0.0_dp
      return
    end if
    med = median(x)
    sigma0 = median(abs(x - med))
    if (sigma0 <= tiny(1.0_dp)) then
      value = 0.0_dp
      if (present(location)) location = med
      return
    end if
    allocate(u(n), weights(n))
    u = abs(x - med) / (c1 * sigma0)
    weights = max(0.0_dp, 1.0_dp - u * u) ** 2
    denominator = sum(weights)
    if (denominator <= tiny(1.0_dp)) then
      mu = med
    else
      mu = sum(x * weights) / denominator
    end if
    u = min(((x - mu) / sigma0) ** 2, c2_squared)
    value = sigma0 * sqrt(sum(u) / (real(n, dp) * expectation))
    if (present(location)) location = mu
  end function tau_scale

  function column_medians(x) result(value)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: value(size(x, 2))
    integer :: j
    do j = 1, size(x, 2)
      value(j) = median(x(:, j))
    end do
  end function column_medians

  subroutine spatial_median(x, center, status, tolerance, max_iterations)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: center(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: next(:), delta(:)
    real(dp) :: distance, denominator, tol, error
    integer :: i, iteration, maxit, n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(center(p), next(p), delta(p))
    if (n < 1 .or. p < 1) then
      center = 0.0_dp
      status = rrcov_invalid_argument
      return
    end if
    center = column_medians(x)
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 500
    if (present(max_iterations)) maxit = max(1, max_iterations)
    do iteration = 1, maxit
      next = 0.0_dp
      denominator = 0.0_dp
      do i = 1, n
        delta = x(i, :) - center
        distance = sqrt(sum(delta * delta))
        if (distance <= tol) cycle
        next = next + x(i, :) / distance
        denominator = denominator + 1.0_dp / distance
      end do
      if (denominator <= tiny(1.0_dp)) exit
      next = next / denominator
      error = sqrt(sum((next - center) ** 2))
      center = next
      if (error <= tol * max(1.0_dp, sqrt(sum(center * center)))) exit
    end do
    if (iteration > maxit) then
      status = 4
    else
      status = rrcov_success
    end if
  end subroutine spatial_median

  subroutine standardize_columns(x, z, center, scale, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: z(:, :), center(:), scale(:)
    integer, intent(out) :: status
    real(dp), allocatable :: covariance(:, :)
    integer :: j, n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(z(n, p), center(p), scale(p))
    covariance = covariance_matrix(x, center=center, status=status)
    do j = 1, p
      scale(j) = sqrt(max(covariance(j, j), 0.0_dp))
      if (scale(j) <= tiny(1.0_dp)) scale(j) = 1.0_dp
      z(:, j) = (x(:, j) - center(j)) / scale(j)
    end do
  end subroutine standardize_columns

  subroutine robust_standardize_columns(x, z, center, scale, status, use_tau)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: z(:, :), center(:), scale(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: use_tau
    logical :: tau
    integer :: j, n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(z(n, p), center(p), scale(p))
    tau = .false.
    if (present(use_tau)) tau = use_tau
    do j = 1, p
      if (tau) then
        scale(j) = tau_scale(x(:, j), center(j))
      else
        center(j) = median(x(:, j))
        scale(j) = mad_scale(x(:, j), consistency=.true.)
      end if
      if (scale(j) <= tiny(1.0_dp)) scale(j) = 1.0_dp
      z(:, j) = (x(:, j) - center(j)) / scale(j)
    end do
    status = rrcov_success
  end subroutine robust_standardize_columns

  pure function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure function regularized_gamma_p(a, x) result(value)
    real(dp), intent(in) :: a, x
    real(dp) :: value
    real(dp) :: sum_value, delta, ap, b, c, d, h, an, log_factor
    integer :: i
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 1.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (x <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    log_factor = -x + a * log(x) - log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sum_value = 1.0_dp / a
      delta = sum_value
      do i = 1, maxit
        ap = ap + 1.0_dp
        delta = delta * x / ap
        sum_value = sum_value + delta
        if (abs(delta) <= abs(sum_value) * eps) exit
      end do
      value = sum_value * exp(log_factor)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do i = 1, maxit
        an = -real(i, dp) * (real(i, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp / d
        delta = d * c
        h = h * delta
        if (abs(delta - 1.0_dp) <= eps) exit
      end do
      value = 1.0_dp - exp(log_factor) * h
    end if
    value = min(1.0_dp, max(0.0_dp, value))
  end function regularized_gamma_p

  pure function chi_square_cdf(x, degrees_freedom) result(value)
    real(dp), intent(in) :: x, degrees_freedom
    real(dp) :: value
    if (x <= 0.0_dp) then
      value = 0.0_dp
    else
      value = regularized_gamma_p(0.5_dp * degrees_freedom, 0.5_dp * x)
    end if
  end function chi_square_cdf

  function chi_square_quantile(probability, degrees_freedom) result(value)
    real(dp), intent(in) :: probability, degrees_freedom
    real(dp) :: value, lower, upper, midpoint, p
    integer :: i
    p = min(1.0_dp - 1.0e-14_dp, max(1.0e-14_dp, probability))
    lower = 0.0_dp
    upper = max(1.0_dp, degrees_freedom)
    do while (chi_square_cdf(upper, degrees_freedom) < p)
      upper = 2.0_dp * upper
      if (upper > 1.0e12_dp) exit
    end do
    do i = 1, 200
      midpoint = 0.5_dp * (lower + upper)
      if (chi_square_cdf(midpoint, degrees_freedom) < p) then
        lower = midpoint
      else
        upper = midpoint
      end if
      if (upper - lower <= 1.0e-12_dp * max(1.0_dp, midpoint)) exit
    end do
    value = 0.5_dp * (lower + upper)
  end function chi_square_quantile

  pure function regularized_beta(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: value, bt
    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      value = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      value = bt * beta_continued_fraction(x, a, b) / a
    else
      value = 1.0_dp - bt * beta_continued_fraction(1.0_dp - x, b, a) / b
    end if
    value = min(1.0_dp, max(0.0_dp, value))
  end function regularized_beta

  pure function beta_continued_fraction(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: value
    real(dp) :: qab, qap, qam, c, d, h, aa, delta
    integer :: m, m2
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 1.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, maxit
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
        ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
        ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      delta = d * c
      h = h * delta
      if (abs(delta - 1.0_dp) <= eps) exit
    end do
    value = h
  end function beta_continued_fraction

  pure function f_cdf(x, df1, df2) result(value)
    real(dp), intent(in) :: x, df1, df2
    real(dp) :: value, z
    if (x <= 0.0_dp) then
      value = 0.0_dp
    else
      z = df1 * x / (df1 * x + df2)
      value = regularized_beta(z, 0.5_dp * df1, 0.5_dp * df2)
    end if
  end function f_cdf
end module rrcov_stats
