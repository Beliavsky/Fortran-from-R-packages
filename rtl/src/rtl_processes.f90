! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_processes
  use rtl_kinds, only: dp
  use rtl_types, only: path_result, ou_fit_result, multivariate_result
  use rtl_stats, only: fill_normal, seed_random, random_poisson, random_lognormal
  use rtl_stats, only: mean_value, sample_sd, kendall_correlation, nearest_psd
  use rtl_stats, only: cholesky_factor, linear_interpolate
  implicit none
  private

  public :: sim_gbm, sim_ou, sim_ou_time, sim_ou_jump, fit_ou
  public :: sim_multivariates
  public :: simGBM, simOU, simOUt, simOUJ, fitOU, simMultivariates

  interface simGBM
    module procedure sim_gbm
  end interface simGBM

  interface simOU
    module procedure sim_ou
  end interface simOU

  interface simOUt
    module procedure sim_ou_time
  end interface simOUt

  interface simOUJ
    module procedure sim_ou_jump
  end interface simOUJ

  interface fitOU
    module procedure fit_ou
  end interface fitOU

  interface simMultivariates
    module procedure sim_multivariates
  end interface simMultivariates

contains

  function sim_gbm(nsims, s0, drift, sigma, maturity, dt, normal_increments, seed) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: s0, drift, sigma, maturity, dt
    real(dp), intent(in), optional :: normal_increments(:, :)
    integer, intent(in), optional :: seed
    type(path_result) :: output
    real(dp), allocatable :: diffusion(:, :)
    integer :: periods, i, j

    call validate_path_inputs(nsims, maturity, dt, output)
    if (.not. output%status%ok) return
    if (s0 <= 0.0_dp .or. sigma < 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "GBM requires positive s0 and nonnegative sigma"
      return
    end if
    periods = nint(maturity / dt)
    if (present(seed)) call seed_random(seed)
    allocate(diffusion(periods, nsims))
    if (present(normal_increments)) then
      if (size(normal_increments, 1) /= periods .or. size(normal_increments, 2) /= nsims) then
        output%status%ok = .false.
        output%status%message = "normal increment dimensions do not match periods and nsims"
        return
      end if
      diffusion = normal_increments
    else
      call fill_normal(diffusion, sqrt(dt))
    end if
    allocate(output%time(0:periods), output%values(0:periods, nsims))
    do i = 0, periods
      output%time(i) = real(i, dp) * dt
    end do
    output%values(0, :) = s0
    do j = 1, nsims
      do i = 1, periods
        output%values(i, j) = output%values(i - 1, j) * exp((drift - 0.5_dp * sigma**2) * dt + &
          sigma * diffusion(i, j))
      end do
    end do
  end function sim_gbm

  function sim_ou(nsims, s0, mu, theta, sigma, maturity, dt, epsilon, seed) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: s0, mu, theta, sigma, maturity, dt
    real(dp), intent(in), optional :: epsilon(:, :)
    integer, intent(in), optional :: seed
    type(path_result) :: output
    real(dp), allocatable :: diffusion(:, :)
    integer :: periods, i

    call validate_path_inputs(nsims, maturity, dt, output)
    if (.not. output%status%ok) return
    periods = nint(maturity / dt)
    if (present(seed)) call seed_random(seed)
    allocate(diffusion(periods, nsims))
    if (present(epsilon)) then
      if (size(epsilon, 1) /= periods .or. size(epsilon, 2) /= nsims) then
        output%status%ok = .false.
        output%status%message = "epsilon dimensions do not match periods and nsims"
        return
      end if
      diffusion = epsilon
    else
      call fill_normal(diffusion, sqrt(dt))
    end if
    allocate(output%time(0:periods), output%values(0:periods, nsims))
    output%values(0, :) = s0
    do i = 0, periods
      output%time(i) = real(i, dp) * dt
    end do
    do i = 1, periods
      output%values(i, :) = output%values(i - 1, :) + &
        theta * (mu - output%values(i - 1, :)) * dt + sigma * diffusion(i, :)
    end do
  end function sim_ou

  function sim_ou_time(nsims, s0, mu_time, mu_value, theta, sigma, maturity, dt, &
      epsilon, seed) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: s0, mu_time(:), mu_value(:), theta, sigma, maturity, dt
    real(dp), intent(in), optional :: epsilon(:, :)
    integer, intent(in), optional :: seed
    type(path_result) :: output
    real(dp), allocatable :: diffusion(:, :), mean_path(:)
    integer :: periods, i

    call validate_path_inputs(nsims, maturity, dt, output)
    if (.not. output%status%ok) return
    if (size(mu_time) /= size(mu_value) .or. size(mu_time) < 2) then
      output%status%ok = .false.
      output%status%message = "mu_time and mu_value must have equal length of at least two"
      return
    end if
    periods = nint(maturity / dt)
    allocate(diffusion(periods, nsims), mean_path(0:periods))
    if (present(seed)) call seed_random(seed)
    if (present(epsilon)) then
      if (size(epsilon, 1) /= periods .or. size(epsilon, 2) /= nsims) then
        output%status%ok = .false.
        output%status%message = "epsilon dimensions do not match periods and nsims"
        return
      end if
      diffusion = epsilon
    else
      call fill_normal(diffusion, sqrt(dt))
    end if
    allocate(output%time(0:periods), output%values(0:periods, nsims))
    output%values(0, :) = s0
    do i = 0, periods
      output%time(i) = real(i, dp) * dt
      mean_path(i) = linear_interpolate(mu_time, mu_value, output%time(i))
    end do
    do i = 1, periods
      output%values(i, :) = output%values(i - 1, :) + &
        theta * (mean_path(i) - output%values(i - 1, :)) * dt + sigma * diffusion(i, :)
    end do
  end function sim_ou_time

  function sim_ou_jump(nsims, s0, mu, theta, sigma, jump_probability, jump_average_size, &
      jump_log_sd, maturity, dt, seed, legacy_recycled_jump) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: s0, mu, theta, sigma, jump_probability
    real(dp), intent(in) :: jump_average_size, jump_log_sd, maturity, dt
    integer, intent(in), optional :: seed
    logical, intent(in), optional :: legacy_recycled_jump
    type(path_result) :: output
    real(dp), allocatable :: diffusion(:, :), jumps(:, :)
    real(dp) :: recycled_size = 0.0_dp
    integer :: periods, i, j, count
    logical :: legacy

    call validate_path_inputs(nsims, maturity, dt, output)
    if (.not. output%status%ok) return
    if (jump_probability < 0.0_dp .or. jump_average_size <= 0.0_dp .or. jump_log_sd < 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "invalid OU-jump parameters"
      return
    end if
    periods = nint(maturity / dt)
    legacy = .false.
    if (present(legacy_recycled_jump)) legacy = legacy_recycled_jump
    if (present(seed)) call seed_random(seed)
    allocate(diffusion(periods, nsims), jumps(periods, nsims))
    call fill_normal(diffusion, sqrt(dt))
    jumps = 0.0_dp
    if (legacy) recycled_size = random_lognormal(log(jump_average_size), jump_log_sd)
    do j = 1, nsims
      do i = 1, periods
        count = random_poisson(jump_probability * dt)
        if (count > 0) then
          if (legacy) then
            jumps(i, j) = real(count, dp) * recycled_size
          else
            jumps(i, j) = sum_independent_jumps(count, jump_average_size, jump_log_sd)
          end if
        end if
      end do
    end do
    allocate(output%time(0:periods), output%values(0:periods, nsims))
    output%values(0, :) = s0
    do i = 0, periods
      output%time(i) = real(i, dp) * dt
    end do
    do i = 1, periods
      output%values(i, :) = output%values(i - 1, :) + theta * &
        (mu - jump_probability * jump_average_size - output%values(i - 1, :)) * dt + &
        sigma * diffusion(i, :) + jumps(i, :)
    end do
  end function sim_ou_jump

  function fit_ou(spread, dt) result(output)
    real(dp), intent(in) :: spread(:), dt
    type(ou_fit_result) :: output
    integer :: n
    real(dp) :: sx, sy, sxx, syy, sxy, denominator, ratio, a, sigma_h2

    n = size(spread)
    if (n < 3 .or. dt <= 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "fit_ou requires at least three values and positive dt"
      return
    end if
    sx = sum(spread(1:n - 1))
    sy = sum(spread(2:n))
    sxx = sum(spread(1:n - 1)**2)
    syy = sum(spread(2:n)**2)
    sxy = sum(spread(1:n - 1) * spread(2:n))
    denominator = real(n - 1, dp) * (sxx - sxy) - (sx**2 - sx * sy)
    if (abs(denominator) <= 1.0e-14_dp) then
      output%status%ok = .false.
      output%status%message = "OU mean estimate is singular"
      return
    end if
    output%mu = (sy * sxx - sx * sxy) / denominator
    denominator = sxx - 2.0_dp * output%mu * sx + real(n - 1, dp) * output%mu**2
    ratio = (sxy - output%mu * sx - output%mu * sy + real(n - 1, dp) * output%mu**2) / denominator
    if (ratio <= 0.0_dp .or. ratio >= 1.0_dp) then
      output%status%ok = .false.
      output%status%message = "OU autoregressive coefficient must lie in (0,1)"
      return
    end if
    output%theta = -log(ratio) / dt
    a = ratio
    sigma_h2 = (syy - 2.0_dp * a * sxy + a**2 * sxx - &
      2.0_dp * output%mu * (1.0_dp - a) * (sy - a * sx) + &
      real(n - 1, dp) * output%mu**2 * (1.0_dp - a)**2) / real(n - 1, dp)
    output%sigma = sqrt(max(0.0_dp, sigma_h2 * 2.0_dp * (-log(a)) / (1.0_dp - a**2) / dt))
    output%half_life_periods = log(2.0_dp) / output%theta / dt
  end function fit_ou

  function sim_multivariates(nsims, prices, s0, seed, use_last_start) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: prices(:, :)
    real(dp), intent(in), optional :: s0(:)
    integer, intent(in), optional :: seed
    logical, intent(in), optional :: use_last_start
    type(multivariate_result) :: output
    real(dp), allocatable :: changes(:, :), start(:), lower(:, :), z(:, :)
    integer :: n, p, i
    logical :: ok, legacy_start

    n = size(prices, 1)
    p = size(prices, 2)
    if (n < 3 .or. p < 1 .or. nsims < 1) then
      output%status%ok = .false.
      output%status%message = "sim_multivariates requires at least three rows and one asset"
      return
    end if
    allocate(changes(n - 1, p), output%mean(p), output%sd(p), output%correlation(p, p))
    allocate(output%covariance(p, p), output%simulations(nsims, p), start(p), lower(p, p), z(nsims, p))
    changes = prices(2:n, :) - prices(1:n - 1, :)
    do i = 1, p
      output%mean(i) = mean_value(changes(:, i))
      output%sd(i) = sample_sd(changes(:, i))
    end do
    call kendall_correlation(changes, output%correlation)
    do i = 1, p
      output%covariance(i, :) = output%sd(i) * output%sd * output%correlation(i, :)
    end do
    call nearest_psd(output%covariance, output%covariance_adjusted)
    call cholesky_factor(output%covariance, lower, ok)
    if (.not. ok) then
      output%status%ok = .false.
      output%status%message = "could not factor adjusted covariance"
      return
    end if
    legacy_start = .false.
    if (present(use_last_start)) legacy_start = use_last_start
    if (present(s0)) then
      if (size(s0) /= p) then
        output%status%ok = .false.
        output%status%message = "s0 length does not match number of assets"
        return
      end if
      start = s0
    else
      start = 0.0_dp
    end if
    if (legacy_start) start = prices(n, :)
    if (present(seed)) call seed_random(seed)
    call fill_normal(z)
    output%simulations = spread(start, 1, nsims) + matmul(z, transpose(lower))
  end function sim_multivariates

  subroutine validate_path_inputs(nsims, maturity, dt, output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: maturity, dt
    type(path_result), intent(inout) :: output
    real(dp) :: periods
    if (nsims < 1 .or. maturity <= 0.0_dp .or. dt <= 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "nsims, maturity, and dt must be positive"
      return
    end if
    periods = maturity / dt
    if (abs(periods - real(nint(periods), dp)) > 1.0e-10_dp) then
      output%status%ok = .false.
      output%status%message = "maturity/dt must be an integer"
    end if
  end subroutine validate_path_inputs

  real(dp) function sum_independent_jumps(count, average_size, log_sd) result(value)
    integer, intent(in) :: count
    real(dp), intent(in) :: average_size, log_sd
    integer :: i
    value = 0.0_dp
    do i = 1, count
      value = value + random_lognormal(log(average_size), log_sd)
    end do
  end function sum_independent_jumps

end module rtl_processes
