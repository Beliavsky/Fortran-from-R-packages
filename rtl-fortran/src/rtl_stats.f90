! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rtl_kinds, only: dp, pi, sqrt_two
  implicit none
  private

  public :: normal_pdf, normal_cdf, mean_value, sample_sd, covariance_value
  public :: kendall_correlation, covariance_matrix, nearest_psd
  public :: seed_random, fill_normal, random_normal, random_poisson
  public :: random_lognormal, cholesky_factor, solve_linear, linear_interpolate
  public :: natural_cubic_interpolate, beta_value, cumulative_product

contains

  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(-x / sqrt_two)
  end function normal_cdf

  pure real(dp) function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: center
    if (size(x) < 2) then
      value = 0.0_dp
      return
    end if
    center = mean_value(x)
    value = sqrt(max(0.0_dp, sum((x - center)**2) / real(size(x) - 1, dp)))
  end function sample_sd

  pure real(dp) function covariance_value(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: mx, my
    integer :: n
    n = min(size(x), size(y))
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    mx = sum(x(1:n)) / real(n, dp)
    my = sum(y(1:n)) / real(n, dp)
    value = sum((x(1:n) - mx) * (y(1:n) - my)) / real(n - 1, dp)
  end function covariance_value

  pure real(dp) function beta_value(asset, benchmark, mask) result(value)
    real(dp), intent(in) :: asset(:), benchmark(:)
    logical, intent(in), optional :: mask(:)
    real(dp), allocatable :: xa(:), xb(:)
    integer :: i, n, k

    n = min(size(asset), size(benchmark))
    if (present(mask)) n = min(n, size(mask))
    allocate(xa(n), xb(n))
    k = 0
    do i = 1, n
      if (present(mask)) then
        if (.not. mask(i)) cycle
      end if
      if (.not. ieee_is_finite(asset(i)) .or. .not. ieee_is_finite(benchmark(i))) cycle
      k = k + 1
      xa(k) = asset(i)
      xb(k) = benchmark(i)
    end do
    if (k < 2 .or. covariance_value(xb(1:k), xb(1:k)) <= 0.0_dp) then
      value = 0.0_dp
    else
      value = covariance_value(xa(1:k), xb(1:k)) / covariance_value(xb(1:k), xb(1:k))
    end if
  end function beta_value

  pure subroutine covariance_matrix(x, covariance, means, sds)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: covariance(:, :)
    real(dp), intent(out), optional :: means(:), sds(:)
    integer :: i, j, p

    p = size(x, 2)
    covariance = 0.0_dp
    do i = 1, p
      if (present(means)) means(i) = mean_value(x(:, i))
      if (present(sds)) sds(i) = sample_sd(x(:, i))
      do j = i, p
        covariance(i, j) = covariance_value(x(:, i), x(:, j))
        covariance(j, i) = covariance(i, j)
      end do
    end do
  end subroutine covariance_matrix

  pure subroutine kendall_correlation(x, correlation)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: correlation(:, :)
    integer :: i, j, a, b, p, n
    integer :: concordant, discordant
    real(dp) :: dx, dy

    n = size(x, 1)
    p = size(x, 2)
    correlation = 0.0_dp
    do i = 1, p
      correlation(i, i) = 1.0_dp
      do j = i + 1, p
        concordant = 0
        discordant = 0
        do a = 1, n - 1
          do b = a + 1, n
            dx = x(b, i) - x(a, i)
            dy = x(b, j) - x(a, j)
            if (dx * dy > 0.0_dp) concordant = concordant + 1
            if (dx * dy < 0.0_dp) discordant = discordant + 1
          end do
        end do
        if (concordant + discordant > 0) then
          correlation(i, j) = real(concordant - discordant, dp) / &
            real(concordant + discordant, dp)
        else
          correlation(i, j) = 0.0_dp
        end if
        correlation(j, i) = correlation(i, j)
      end do
    end do
  end subroutine kendall_correlation

  subroutine nearest_psd(a, adjusted)
    real(dp), intent(inout) :: a(:, :)
    logical, intent(out), optional :: adjusted
    real(dp), allocatable :: v(:, :), d(:), work(:, :)
    real(dp) :: max_off, theta, c, s, app, aqq, apq, floor_value
    integer :: n, i, j, p, q, iter
    logical :: changed

    n = size(a, 1)
    allocate(v(n, n), d(n), work(n, n))
    work = 0.5_dp * (a + transpose(a))
    v = 0.0_dp
    do i = 1, n
      v(i, i) = 1.0_dp
    end do

    do iter = 1, 100 * max(1, n * n)
      max_off = 0.0_dp
      p = 1
      q = 1
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(work(i, j)) > max_off) then
            max_off = abs(work(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (max_off < 1.0e-13_dp) exit
      app = work(p, p)
      aqq = work(q, q)
      apq = work(p, q)
      theta = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
      c = cos(theta)
      s = sin(theta)
      call jacobi_rotate(work, p, q, c, s)
      call rotate_columns(v, p, q, c, s)
    end do

    do i = 1, n
      d(i) = work(i, i)
    end do
    floor_value = max(1.0e-12_dp, maxval(abs(d)) * 1.0e-12_dp)
    changed = any(d < floor_value)
    d = max(d, floor_value)
    a = matmul(v, matmul(diagonal_matrix(d), transpose(v)))
    a = 0.5_dp * (a + transpose(a))
    if (present(adjusted)) adjusted = changed
  contains
    subroutine jacobi_rotate(matrix, ip, iq, cc, ss)
      real(dp), intent(inout) :: matrix(:, :)
      integer, intent(in) :: ip, iq
      real(dp), intent(in) :: cc, ss
      real(dp) :: mpk, mqk
      integer :: k
      do k = 1, size(matrix, 1)
        if (k /= ip .and. k /= iq) then
          mpk = matrix(ip, k)
          mqk = matrix(iq, k)
          matrix(ip, k) = cc * mpk - ss * mqk
          matrix(k, ip) = matrix(ip, k)
          matrix(iq, k) = ss * mpk + cc * mqk
          matrix(k, iq) = matrix(iq, k)
        end if
      end do
      mpk = matrix(ip, ip)
      mqk = matrix(iq, iq)
      matrix(ip, ip) = cc * cc * mpk - 2.0_dp * ss * cc * matrix(ip, iq) + ss * ss * mqk
      matrix(iq, iq) = ss * ss * mpk + 2.0_dp * ss * cc * matrix(ip, iq) + cc * cc * mqk
      matrix(ip, iq) = 0.0_dp
      matrix(iq, ip) = 0.0_dp
    end subroutine jacobi_rotate

    subroutine rotate_columns(matrix, ip, iq, cc, ss)
      real(dp), intent(inout) :: matrix(:, :)
      integer, intent(in) :: ip, iq
      real(dp), intent(in) :: cc, ss
      real(dp) :: vp, vq
      integer :: k
      do k = 1, size(matrix, 1)
        vp = matrix(k, ip)
        vq = matrix(k, iq)
        matrix(k, ip) = cc * vp - ss * vq
        matrix(k, iq) = ss * vp + cc * vq
      end do
    end subroutine rotate_columns

    pure function diagonal_matrix(values) result(matrix)
      real(dp), intent(in) :: values(:)
      real(dp) :: matrix(size(values), size(values))
      integer :: k
      matrix = 0.0_dp
      do k = 1, size(values)
        matrix(k, k) = values(k)
      end do
    end function diagonal_matrix
  end subroutine nearest_psd

  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer, allocatable :: values(:)
    integer :: n, i
    call random_seed(size=n)
    allocate(values(n))
    do i = 1, n
      values(i) = modulo(seed + 104729 * i, huge(1) - 1)
      if (values(i) == 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine seed_random

  real(dp) function random_normal() result(value)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, 1.0e-15_dp)
    value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
  end function random_normal

  subroutine fill_normal(x, sd)
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(in), optional :: sd
    real(dp) :: scale
    integer :: i, j
    scale = 1.0_dp
    if (present(sd)) scale = sd
    do j = 1, size(x, 2)
      do i = 1, size(x, 1)
        x(i, j) = scale * random_normal()
      end do
    end do
  end subroutine fill_normal

  integer function random_poisson(lambda) result(value)
    real(dp), intent(in) :: lambda
    real(dp) :: limit_value, product_value
    integer :: k
    if (lambda <= 0.0_dp) then
      value = 0
      return
    end if
    if (lambda < 30.0_dp) then
      limit_value = exp(-lambda)
      product_value = 1.0_dp
      k = 0
      do
        k = k + 1
        call random_number(product_value)
        exit
      end do
      product_value = 1.0_dp
      k = 0
      do while (product_value > limit_value)
        k = k + 1
        call multiply_uniform(product_value)
      end do
      value = k - 1
    else
      value = max(0, nint(lambda + sqrt(lambda) * random_normal()))
    end if
  contains
    subroutine multiply_uniform(product_input)
      real(dp), intent(inout) :: product_input
      real(dp) :: u
      call random_number(u)
      product_input = product_input * u
    end subroutine multiply_uniform
  end function random_poisson

  real(dp) function random_lognormal(log_mean, log_sd) result(value)
    real(dp), intent(in) :: log_mean, log_sd
    value = exp(log_mean + log_sd * random_normal())
  end function random_lognormal

  subroutine cholesky_factor(a, lower, ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: lower(:, :)
    logical, intent(out) :: ok
    real(dp) :: subtotal
    integer :: i, j, k, n

    n = size(a, 1)
    lower = 0.0_dp
    ok = .true.
    do i = 1, n
      do j = 1, i
        subtotal = a(i, j)
        do k = 1, j - 1
          subtotal = subtotal - lower(i, k) * lower(j, k)
        end do
        if (i == j) then
          if (subtotal <= 0.0_dp) then
            ok = .false.
            return
          end if
          lower(i, j) = sqrt(subtotal)
        else
          lower(i, j) = subtotal / lower(j, j)
        end if
      end do
    end do
  end subroutine cholesky_factor

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:, :), row_tmp(:)
    real(dp) :: factor, pivot_value
    integer :: i, j, k, pivot, n

    n = size(b)
    allocate(aug(n, n + 1), row_tmp(n + 1))
    aug(:, 1:n) = a
    aug(:, n + 1) = b
    ok = .true.
    do k = 1, n
      pivot = k
      pivot_value = abs(aug(k, k))
      do i = k + 1, n
        if (abs(aug(i, k)) > pivot_value) then
          pivot = i
          pivot_value = abs(aug(i, k))
        end if
      end do
      if (pivot_value <= 1.0e-14_dp) then
        ok = .false.
        x = 0.0_dp
        return
      end if
      if (pivot /= k) then
        row_tmp = aug(k, :)
        aug(k, :) = aug(pivot, :)
        aug(pivot, :) = row_tmp
      end if
      do i = k + 1, n
        factor = aug(i, k) / aug(k, k)
        aug(i, k:n + 1) = aug(i, k:n + 1) - factor * aug(k, k:n + 1)
      end do
    end do
    do i = n, 1, -1
      x(i) = aug(i, n + 1)
      do j = i + 1, n
        x(i) = x(i) - aug(i, j) * x(j)
      end do
      x(i) = x(i) / aug(i, i)
    end do
  end subroutine solve_linear

  pure real(dp) function linear_interpolate(x, y, xout) result(value)
    real(dp), intent(in) :: x(:), y(:), xout
    integer :: i, n
    n = size(x)
    if (n == 1) then
      value = y(1)
      return
    end if
    if (xout <= x(1)) then
      value = y(1) + (xout - x(1)) * (y(2) - y(1)) / (x(2) - x(1))
      return
    end if
    if (xout >= x(n)) then
      value = y(n - 1) + (xout - x(n - 1)) * (y(n) - y(n - 1)) / (x(n) - x(n - 1))
      return
    end if
    do i = 1, n - 1
      if (xout >= x(i) .and. xout <= x(i + 1)) then
        value = y(i) + (xout - x(i)) * (y(i + 1) - y(i)) / (x(i + 1) - x(i))
        return
      end if
    end do
    value = y(n)
  end function linear_interpolate

  function natural_cubic_interpolate(x, y, xout) result(value)
    real(dp), intent(in) :: x(:), y(:), xout
    real(dp) :: value
    real(dp), allocatable :: y2(:), u(:)
    real(dp) :: sig, p, h, a, b, xx
    integer :: n, i, k_low, k_high, k

    n = size(x)
    if (n < 3) then
      value = linear_interpolate(x, y, xout)
      return
    end if
    allocate(y2(n), u(n - 1))
    y2(1) = 0.0_dp
    u(1) = 0.0_dp
    do i = 2, n - 1
      sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
      p = sig * y2(i - 1) + 2.0_dp
      y2(i) = (sig - 1.0_dp) / p
      u(i) = (6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
        (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do
    y2(n) = 0.0_dp
    do k = n - 1, 1, -1
      y2(k) = y2(k) * y2(k + 1) + u(k)
    end do

    xx = min(max(xout, x(1)), x(n))
    k_low = 1
    k_high = n
    do while (k_high - k_low > 1)
      k = (k_high + k_low) / 2
      if (x(k) > xx) then
        k_high = k
      else
        k_low = k
      end if
    end do
    h = x(k_high) - x(k_low)
    a = (x(k_high) - xx) / h
    b = (xx - x(k_low)) / h
    value = a * y(k_low) + b * y(k_high) + &
      ((a**3 - a) * y2(k_low) + (b**3 - b) * y2(k_high)) * h * h / 6.0_dp
  end function natural_cubic_interpolate

  pure subroutine cumulative_product(x, result)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: result(:)
    integer :: i
    if (size(x) == 0) return
    result(1) = x(1)
    do i = 2, size(x)
      result(i) = result(i - 1) * x(i)
    end do
  end subroutine cumulative_product

end module rtl_stats
