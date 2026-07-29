! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_math
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use acdm_kinds, only : dp, pi, tiny_pos, ACDM_SUCCESS, ACDM_BAD_INPUT, &
                         ACDM_NUMERIC_FAILURE, ACDM_SINGULAR
  implicit none
  private

  type, public :: rng_state
    integer(int64) :: state = 88172645463393265_int64
  end type rng_state

  public :: seed_rng, random_uniform, random_normal, random_gamma
  public :: normal_cdf, normal_quantile, gamma_p, beta_i
  public :: gamma_quantile, beta_quantile, chi_square_cdf
  public :: solve_linear, invert_matrix, least_squares
  public :: sample_mean, sample_variance, empirical_quantile
  public :: natural_spline_second, natural_spline_eval
  public :: autocorrelation

contains

  subroutine seed_rng(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: seed
    integer(int64) :: s

    s = int(seed, int64)
    if (s == 0_int64) s = 88172645463393265_int64
    rng%state = ieor(s, int(z'9E3779B97F4A7C15', int64))
    call warm_rng(rng)
  end subroutine seed_rng

  subroutine warm_rng(rng)
    type(rng_state), intent(inout) :: rng
    integer :: i
    real(dp) :: dummy

    do i = 1, 16
      dummy = random_uniform(rng)
    end do
  end subroutine warm_rng

  function next_uint64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(int64) :: x

    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
  end function next_uint64

  function random_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(int64) :: x, mask

    mask = int(z'001FFFFFFFFFFFFF', int64)
    x = iand(next_uint64(rng), mask)
    u = (real(x, dp) + 0.5_dp) / 9007199254740992.0_dp
    u = max(tiny_pos, min(1.0_dp - tiny_pos, u))
  end function random_uniform

  function random_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z
    real(dp) :: u1, u2

    u1 = random_uniform(rng)
    u2 = random_uniform(rng)
    z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
  end function random_normal

  recursive function random_gamma(rng, shape) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp) :: x
    real(dp) :: d, c, z, v, u

    if (shape <= 0.0_dp) then
      x = ieee_value(x, ieee_quiet_nan)
      return
    end if

    if (shape < 1.0_dp) then
      x = random_gamma(rng, shape + 1.0_dp) * &
          random_uniform(rng)**(1.0_dp / shape)
      return
    end if

    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      do
        z = random_normal(rng)
        v = 1.0_dp + c * z
        if (v > 0.0_dp) exit
      end do
      v = v**3
      u = random_uniform(rng)
      if (u < 1.0_dp - 0.0331_dp * z**4) exit
      if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
    end do
    x = d * v
  end function random_gamma

  pure elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p

    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x
    real(dp) :: q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp - plow

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + &
            c(5)) * q + c(6)) / &
          ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + &
            a(5)) * r + a(6)) * q / &
          (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + &
            b(5)) * r + 1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - p))
      x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + &
             c(5)) * q + c(6)) / &
           ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    end if
  end function normal_quantile

  pure function gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p
    integer, parameter :: itmax = 400
    real(dp), parameter :: eps = 2.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    integer :: n
    real(dp) :: sumv, del, ap, b, c, d, h, an, gln

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(p, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if

    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do n = 1, itmax
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv) * eps) exit
      end do
      p = sumv * exp(-x + a * log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do n = 1, itmax
        an = -real(n, dp) * (real(n, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) <= eps) exit
      end do
      p = 1.0_dp - exp(-x + a * log(x) - gln) * h
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function gamma_p

  pure function beta_i(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: p
    real(dp) :: bt

    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
      p = ieee_value(p, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if

    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
             a * log(x) + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      p = bt * beta_cf(a, b, x) / a
    else
      p = 1.0_dp - bt * beta_cf(b, a, 1.0_dp - x) / b
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function beta_i

  pure function beta_cf(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    real(dp) :: h
    integer, parameter :: maxit = 400
    real(dp), parameter :: eps = 2.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, qab, qam, qap

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
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
  end function beta_cf

  function gamma_quantile(p, shape) result(x)
    real(dp), intent(in) :: p, shape
    real(dp) :: x
    real(dp) :: lo, hi, mid
    integer :: i

    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    else if (shape <= 0.0_dp) then
      x = ieee_value(x, ieee_quiet_nan)
      return
    end if

    lo = 0.0_dp
    hi = max(1.0_dp, shape)
    do while (gamma_p(shape, hi) < p .and. hi < 1.0e300_dp)
      hi = hi * 2.0_dp
    end do
    do i = 1, 160
      mid = 0.5_dp * (lo + hi)
      if (gamma_p(shape, mid) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function gamma_quantile

  function beta_quantile(p, a, b) result(x)
    real(dp), intent(in) :: p, a, b
    real(dp) :: x
    real(dp) :: lo, hi, mid
    integer :: i

    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = 1.0_dp
      return
    else if (a <= 0.0_dp .or. b <= 0.0_dp) then
      x = ieee_value(x, ieee_quiet_nan)
      return
    end if

    lo = 0.0_dp
    hi = 1.0_dp
    do i = 1, 160
      mid = 0.5_dp * (lo + hi)
      if (beta_i(mid, a, b) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function beta_quantile

  pure function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = gamma_p(0.5_dp * df, 0.5_dp * x)
    end if
  end function chi_square_cdf

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:, :), rowtmp(:)
    real(dp) :: factor, pivot
    integer :: n, i, j, k, ip

    n = size(b)
    status = ACDM_SUCCESS
    if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) then
      status = ACDM_BAD_INPUT
      return
    end if

    allocate(aug(n, n + 1), rowtmp(n + 1))
    aug(:, 1:n) = a
    aug(:, n + 1) = b

    do k = 1, n
      ip = k
      pivot = abs(aug(k, k))
      do i = k + 1, n
        if (abs(aug(i, k)) > pivot) then
          pivot = abs(aug(i, k))
          ip = i
        end if
      end do
      if (pivot <= 1.0e-14_dp) then
        status = ACDM_SINGULAR
        x = 0.0_dp
        return
      end if
      if (ip /= k) then
        rowtmp = aug(k, :)
        aug(k, :) = aug(ip, :)
        aug(ip, :) = rowtmp
      end if
      do i = k + 1, n
        factor = aug(i, k) / aug(k, k)
        aug(i, k:n + 1) = aug(i, k:n + 1) - &
                          factor * aug(k, k:n + 1)
      end do
    end do

    x = 0.0_dp
    do i = n, 1, -1
      x(i) = aug(i, n + 1)
      do j = i + 1, n
        x(i) = x(i) - aug(i, j) * x(j)
      end do
      x(i) = x(i) / aug(i, i)
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j, st

    n = size(a, 1)
    if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. &
        size(ainv, 2) /= n) then
      status = ACDM_BAD_INPUT
      return
    end if

    allocate(e(n), x(n))
    ainv = 0.0_dp
    status = ACDM_SUCCESS
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, x, st)
      if (st /= ACDM_SUCCESS) then
        status = st
        return
      end if
      ainv(:, j) = x
    end do
  end subroutine invert_matrix

  subroutine least_squares(xmat, y, beta, residual, status, ridge)
    real(dp), intent(in) :: xmat(:, :), y(:)
    real(dp), intent(out) :: beta(:), residual(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: ridge
    real(dp), allocatable :: xtx(:, :), xty(:)
    real(dp) :: r
    integer :: p, i

    p = size(xmat, 2)
    if (size(xmat, 1) /= size(y) .or. size(beta) /= p .or. &
        size(residual) /= size(y)) then
      status = ACDM_BAD_INPUT
      return
    end if

    allocate(xtx(p, p), xty(p))
    xtx = matmul(transpose(xmat), xmat)
    xty = matmul(transpose(xmat), y)
    r = 1.0e-10_dp
    if (present(ridge)) r = max(0.0_dp, ridge)
    do i = 1, p
      xtx(i, i) = xtx(i, i) + r
    end do
    call solve_linear(xtx, xty, beta, status)
    if (status == ACDM_SUCCESS) residual = y - matmul(xmat, beta)
  end subroutine least_squares

  pure function sample_mean(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m

    if (size(x) == 0) then
      m = ieee_value(m, ieee_quiet_nan)
    else
      m = sum(x) / real(size(x), dp)
    end if
  end function sample_mean

  pure function sample_variance(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v, m

    if (size(x) < 2) then
      v = ieee_value(v, ieee_quiet_nan)
    else
      m = sample_mean(x)
      v = sum((x - m)**2) / real(size(x) - 1, dp)
    end if
  end function sample_variance

  function empirical_quantile(x, prob) result(q)
    real(dp), intent(in) :: x(:), prob
    real(dp) :: q
    real(dp), allocatable :: y(:)
    real(dp) :: h, frac
    integer :: n, lo, hi

    n = size(x)
    if (n == 0 .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      q = ieee_value(q, ieee_quiet_nan)
      return
    end if
    allocate(y(n))
    y = x
    call quicksort(y, 1, n)
    if (n == 1) then
      q = y(1)
      return
    end if
    h = 1.0_dp + real(n - 1, dp) * prob
    lo = floor(h)
    hi = ceiling(h)
    frac = h - real(lo, dp)
    q = (1.0_dp - frac) * y(lo) + frac * y(hi)
  end function empirical_quantile

  recursive subroutine quicksort(a, left, right)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(dp) :: pivot, temp

    if (left >= right) return
    i = left
    j = right
    pivot = a((left + right) / 2)
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp = a(i)
        a(i) = a(j)
        a(j) = temp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
  end subroutine quicksort

  subroutine natural_spline_second(x, y, y2, status)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: y2(:)
    integer, intent(out) :: status
    real(dp), allocatable :: u(:)
    real(dp) :: sig, p
    integer :: n, i, k

    n = size(x)
    status = ACDM_SUCCESS
    if (n < 2 .or. size(y) /= n .or. size(y2) /= n) then
      status = ACDM_BAD_INPUT
      return
    end if
    if (any(x(2:n) <= x(1:n - 1))) then
      status = ACDM_BAD_INPUT
      return
    end if

    allocate(u(n))
    y2(1) = 0.0_dp
    u(1) = 0.0_dp
    do i = 2, n - 1
      sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
      p = sig * y2(i - 1) + 2.0_dp
      y2(i) = (sig - 1.0_dp) / p
      u(i) = (6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
              (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / &
              (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do
    y2(n) = 0.0_dp
    do k = n - 1, 1, -1
      y2(k) = y2(k) * y2(k + 1) + u(k)
    end do
  end subroutine natural_spline_second

  pure function natural_spline_eval(xa, ya, y2a, x) result(y)
    real(dp), intent(in) :: xa(:), ya(:), y2a(:), x
    real(dp) :: y
    integer :: klo, khi, k, n
    real(dp) :: h, a, b

    n = size(xa)
    klo = 1
    khi = n
    if (x <= xa(1)) then
      klo = 1
      khi = 2
    else if (x >= xa(n)) then
      klo = n - 1
      khi = n
    else
      do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (xa(k) > x) then
          khi = k
        else
          klo = k
        end if
      end do
    end if
    h = xa(khi) - xa(klo)
    if (h <= 0.0_dp) then
      y = ieee_value(y, ieee_quiet_nan)
      return
    end if
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * h**2 / 6.0_dp
  end function natural_spline_eval

  function autocorrelation(x, max_lag) result(acf)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: max_lag
    real(dp), allocatable :: acf(:)
    real(dp) :: m, denom
    integer :: lag, n

    n = size(x)
    allocate(acf(0:max_lag))
    if (n == 0) then
      acf = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    m = sample_mean(x)
    denom = sum((x - m)**2)
    if (denom <= 0.0_dp) then
      acf = 0.0_dp
      acf(0) = 1.0_dp
      return
    end if
    acf(0) = 1.0_dp
    do lag = 1, max_lag
      if (lag >= n) then
        acf(lag) = 0.0_dp
      else
        acf(lag) = sum((x(1:n - lag) - m) * &
                       (x(1 + lag:n) - m)) / denom
      end if
    end do
  end function autocorrelation

end module acdm_math
