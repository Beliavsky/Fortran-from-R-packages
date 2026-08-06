! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_utils
  use ghyp_kinds, only : dp, i8
  use ghyp_linalg, only : inverse_spd, logdet_spd, cholesky_lower
  use tsd_math, only : regularized_gamma_p
  use tsmarch_types, only : escc_result, tsm_success, tsm_invalid_argument, tsm_numerical_failure
  use tsmarch_linalg
  implicit none
  private
  public :: ewma_covariance, lw_covariance, cor2cov, make_psd
  public :: multivariate_normal_density, multivariate_student_density
  public :: rmvnorm, rmvt, combn_fast, escc_test
  public :: aggregate_mean, aggregate_sigma, lag_matrix
  public :: lower_triangle, triangle_to_symmetric

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  function ewma_covariance(x, lambda, demean) result(covariance)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: lambda
    logical, intent(in), optional :: demean
    real(dp), allocatable :: covariance(:, :, :), z(:, :), mu(:)
    logical :: center
    integer :: t, n, m
    n = size(x, 1)
    m = size(x, 2)
    allocate(covariance(m, m, n), z(n, m))
    center = .true.
    if (present(demean)) center = demean
    if (center) then
      mu = column_mean(x)
      z = x - spread(mu, 1, n)
    else
      z = x
    end if
    covariance(:, :, 1) = sample_covariance(z)
    do t = 2, n
      covariance(:, :, t) = lambda * covariance(:, :, t - 1) + &
        (1.0_dp - lambda) * outer_product(z(t - 1, :), z(t - 1, :))
    end do
  end function ewma_covariance

  function lw_covariance(x) result(shrunk)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable :: shrunk(:, :), sample(:, :), target(:, :), centered(:, :), mu(:)
    real(dp) :: muvar, phi, rho, gamma, shrinkage
    integer :: n, m, t, i
    n = size(x, 1)
    m = size(x, 2)
    sample = sample_covariance(x)
    allocate(target(m, m), centered(n, m))
    target = 0.0_dp
    muvar = sum([(sample(i, i), i = 1, m)]) / real(m, dp)
    do i = 1, m
      target(i, i) = muvar
    end do
    mu = column_mean(x)
    centered = x - spread(mu, 1, n)
    phi = 0.0_dp
    do t = 1, n
      phi = phi + sum((outer_product(centered(t, :), centered(t, :)) - sample) ** 2)
    end do
    phi = phi / real(n, dp)
    gamma = sum((sample - target) ** 2)
    rho = 0.0_dp
    if (gamma > tiny(1.0_dp)) then
      shrinkage = min(1.0_dp, max(0.0_dp, (phi - rho) / (real(n, dp) * gamma)))
    else
      shrinkage = 1.0_dp
    end if
    shrunk = shrinkage * target + (1.0_dp - shrinkage) * sample
  end function lw_covariance

  function cor2cov(correlation, sigma) result(covariance)
    real(dp), intent(in) :: correlation(:, :), sigma(:)
    real(dp), allocatable :: covariance(:, :)
    integer :: i, j, m
    m = size(correlation, 1)
    allocate(covariance(m, m))
    do j = 1, m
      do i = 1, m
        covariance(i, j) = correlation(i, j) * sigma(i) * sigma(j)
      end do
    end do
  end function cor2cov

  function make_psd(a, floor_value) result(out)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: floor_value
    real(dp), allocatable :: out(:, :)
    real(dp) :: floorv
    floorv = 1.0e-8_dp
    if (present(floor_value)) floorv = floor_value
    out = nearest_correlation(a, floorv)
  end function make_psd

  function multivariate_normal_density(x, mean, covariance, log_density) result(value)
    real(dp), intent(in) :: x(:), mean(:), covariance(:, :)
    logical, intent(in), optional :: log_density
    real(dp) :: value, q, ld
    real(dp), allocatable :: inv(:, :), d(:)
    logical :: ok, logd
    integer :: m
    m = size(x)
    logd = .false.
    if (present(log_density)) logd = log_density
    if (size(mean) /= m .or. size(covariance, 1) /= m .or. size(covariance, 2) /= m) then
      value = -huge(1.0_dp)
      if (.not. logd) value = 0.0_dp
      return
    end if
    call inverse_spd(covariance, inv, ok)
    if (.not. ok) then
      value = -huge(1.0_dp)
      if (.not. logd) value = 0.0_dp
      return
    end if
    ld = logdet_spd(covariance, ok)
    d = x - mean
    q = dot_product(d, matmul(inv, d))
    value = -0.5_dp * (real(m, dp) * log(2.0_dp * pi) + ld + q)
    if (.not. logd) value = exp(value)
  end function multivariate_normal_density

  function multivariate_student_density(x, mean, covariance, shape, log_density) result(value)
    real(dp), intent(in) :: x(:), mean(:), covariance(:, :), shape
    logical, intent(in), optional :: log_density
    real(dp) :: value, q, ld
    real(dp), allocatable :: inv(:, :), d(:)
    logical :: ok, logd
    integer :: m
    m = size(x)
    logd = .false.
    if (present(log_density)) logd = log_density
    if (shape <= 2.0_dp .or. size(mean) /= m .or. size(covariance, 1) /= m .or. size(covariance, 2) /= m) then
      value = -huge(1.0_dp)
      if (.not. logd) value = 0.0_dp
      return
    end if
    call inverse_spd(covariance, inv, ok)
    if (.not. ok) then
      value = -huge(1.0_dp)
      if (.not. logd) value = 0.0_dp
      return
    end if
    ld = logdet_spd(covariance, ok)
    d = x - mean
    q = dot_product(d, matmul(inv, d))
    value = log_gamma(0.5_dp * (shape + real(m, dp))) - log_gamma(0.5_dp * shape) - &
      0.5_dp * (ld + real(m, dp) * log(pi * (shape - 2.0_dp))) - &
      0.5_dp * (shape + real(m, dp)) * log(1.0_dp + q / (shape - 2.0_dp))
    if (.not. logd) value = exp(value)
  end function multivariate_student_density

  function rmvnorm(n, mean, covariance, seed) result(draws)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:), covariance(:, :)
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: draws(:, :), z(:, :), l(:, :)
    logical :: ok
    integer :: i
    allocate(z(n, size(mean)))
    if (present(seed)) then
      z = random_normal_matrix(n, size(mean), seed)
    else
      z = random_normal_matrix(n, size(mean))
    end if
    call cholesky_lower(covariance, l, ok)
    allocate(draws(n, size(mean)))
    if (.not. ok) then
      draws = huge(1.0_dp)
      return
    end if
    do i = 1, n
      draws(i, :) = mean + matmul(l, z(i, :))
    end do
  end function rmvnorm

  function rmvt(n, mean, covariance, shape, seed) result(draws)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:), covariance(:, :), shape
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: draws(:, :), z(:, :), l(:, :)
    real(dp) :: chisq
    logical :: ok
    integer :: i
    allocate(z(n, size(mean)))
    if (present(seed)) then
      z = random_normal_matrix(n, size(mean), seed)
    else
      z = random_normal_matrix(n, size(mean))
    end if
    call cholesky_lower(covariance, l, ok)
    allocate(draws(n, size(mean)))
    if (.not. ok .or. shape <= 2.0_dp) then
      draws = huge(1.0_dp)
      return
    end if
    do i = 1, n
      chisq = chi_square_draw(shape)
      draws(i, :) = mean + matmul(l, z(i, :)) * sqrt((shape - 2.0_dp) / max(chisq, tiny(1.0_dp)))
    end do
  end function rmvt

  function combn_fast(n, m) result(combinations)
    integer, intent(in) :: n, m
    integer, allocatable :: combinations(:, :)
    integer, allocatable :: current(:)
    integer :: count, position
    if (m < 0 .or. m > n) then
      allocate(combinations(0, 0))
      return
    end if
    count = binomial_integer(n, m)
    allocate(combinations(m, count), current(m))
    position = 0
    call generate_combinations(1, 1)
  contains
    recursive subroutine generate_combinations(start_value, depth)
      integer, intent(in) :: start_value, depth
      integer :: value
      if (depth > m) then
        position = position + 1
        combinations(:, position) = current
        return
      end if
      do value = start_value, n - m + depth
        current(depth) = value
        call generate_combinations(value + 1, depth + 1)
      end do
    end subroutine generate_combinations
  end function combn_fast

  function escc_test(whitened_residuals, lags) result(out)
    real(dp), intent(in) :: whitened_residuals(:, :)
    integer, intent(in), optional :: lags
    type(escc_result) :: out
    real(dp), allocatable :: products(:, :), regressors(:, :), response(:), xtx(:, :), xty(:), beta(:), inv(:, :), noise(:)
    real(dp) :: sigma2
    integer :: lag_count, n, m, pairs, i, j, k, t, row
    logical :: ok
    n = size(whitened_residuals, 1)
    m = size(whitened_residuals, 2)
    lag_count = 1
    if (present(lags)) lag_count = lags
    if (n <= lag_count + 2 .or. m < 2 .or. lag_count < 1) then
      out%status = tsm_invalid_argument
      out%message = 'ESCC requires two or more series and enough observations for the requested lags'
      return
    end if
    pairs = m * (m - 1) / 2
    allocate(products(n, pairs))
    k = 0
    do i = 1, m - 1
      do j = i + 1, m
        k = k + 1
        products(:, k) = whitened_residuals(:, i) * whitened_residuals(:, j)
      end do
    end do
    allocate(regressors(pairs * (n - lag_count), lag_count + 1), response(pairs * (n - lag_count)))
    row = 0
    do k = 1, pairs
      do t = lag_count + 1, n
        row = row + 1
        regressors(row, 1) = 1.0_dp
        do j = 1, lag_count
          regressors(row, j + 1) = products(t - j, k)
        end do
        response(row) = products(t, k)
      end do
    end do
    xtx = matmul(transpose(regressors), regressors)
    xty = matmul(transpose(regressors), response)
    call matrix_inverse_general(xtx, inv, ok)
    if (.not. ok) then
      out%status = tsm_numerical_failure
      out%message = 'ESCC artificial regression is singular'
      return
    end if
    beta = matmul(inv, xty)
    noise = response - matmul(regressors, beta)
    sigma2 = sum(noise ** 2) / real(max(1, size(response) - lag_count - 1), dp)
    out%statistic = dot_product(beta, matmul(xtx, beta)) / sqrt(max(sigma2, tiny(1.0_dp)))
    out%degrees_freedom = lag_count + 1
    out%p_value = 1.0_dp - regularized_gamma_p(0.5_dp * real(out%degrees_freedom, dp), 0.5_dp * out%statistic)
    out%p_value = min(max(out%p_value, 0.0_dp), 1.0_dp)
    out%status = tsm_success
    out%message = 'ok'
  end function escc_test

  function aggregate_mean(mean_matrix, weights) result(value)
    real(dp), intent(in) :: mean_matrix(:, :), weights(:)
    real(dp), allocatable :: value(:)
    value = matmul(mean_matrix, weights)
  end function aggregate_mean

  function aggregate_sigma(covariance, weights) result(value)
    real(dp), intent(in) :: covariance(:, :, :), weights(:)
    real(dp), allocatable :: value(:)
    integer :: t
    allocate(value(size(covariance, 3)))
    do t = 1, size(covariance, 3)
      value(t) = sqrt(max(dot_product(weights, matmul(covariance(:, :, t), weights)), 0.0_dp))
    end do
  end function aggregate_sigma

  subroutine lag_matrix(x, lags, include_constant, y, regressors)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lags
    logical, intent(in), optional :: include_constant
    real(dp), allocatable, intent(out) :: y(:), regressors(:, :)
    logical :: constant
    integer :: t, j, offset
    constant = .false.
    if (present(include_constant)) constant = include_constant
    offset = merge(1, 0, constant)
    allocate(y(size(x) - lags), regressors(size(x) - lags, lags + offset))
    do t = lags + 1, size(x)
      y(t - lags) = x(t)
      if (constant) regressors(t - lags, 1) = 1.0_dp
      do j = 1, lags
        regressors(t - lags, offset + j) = x(t - j)
      end do
    end do
  end subroutine lag_matrix

  function lower_triangle(a, include_diagonal) result(values)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in), optional :: include_diagonal
    real(dp), allocatable :: values(:)
    logical :: diag
    integer :: n, count, i, j, k
    n = size(a, 1)
    diag = .false.
    if (present(include_diagonal)) diag = include_diagonal
    count = merge(n * (n + 1) / 2, n * (n - 1) / 2, diag)
    allocate(values(count))
    k = 0
    do j = 1, n
      do i = j + merge(0, 1, diag), n
        k = k + 1
        values(k) = a(i, j)
      end do
    end do
  end function lower_triangle

  function triangle_to_symmetric(values, n, include_diagonal) result(a)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: n
    logical, intent(in), optional :: include_diagonal
    real(dp), allocatable :: a(:, :)
    logical :: diag
    integer :: i, j, k
    diag = .false.
    if (present(include_diagonal)) diag = include_diagonal
    allocate(a(n, n))
    a = 0.0_dp
    if (.not. diag) then
      do i = 1, n
        a(i, i) = 1.0_dp
      end do
    end if
    k = 0
    do j = 1, n
      do i = j + merge(0, 1, diag), n
        k = k + 1
        if (k <= size(values)) then
          a(i, j) = values(k)
          a(j, i) = values(k)
        end if
      end do
    end do
  end function triangle_to_symmetric

  integer function binomial_integer(n, k) result(value)
    integer, intent(in) :: n, k
    integer :: i, kk
    kk = min(k, n - k)
    value = 1
    do i = 1, kk
      value = value * (n - kk + i) / i
    end do
  end function binomial_integer

  function chi_square_draw(df) result(value)
    real(dp), intent(in) :: df
    real(dp) :: value, u
    real(dp), allocatable :: z(:, :)
    integer :: k
    k = max(1, nint(df))
    z = random_normal_matrix(k, 1)
    value = sum(z(:, 1) ** 2)
    if (abs(real(k, dp) - df) > 0.25_dp) then
      call random_number(u)
      value = value * df / real(k, dp) * (0.75_dp + 0.5_dp * u)
    end if
    value = max(value, tiny(1.0_dp))
  end function chi_square_draw

end module tsmarch_utils
