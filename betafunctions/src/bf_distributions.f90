! SPDX-License-Identifier: CC0-1.0
module bf_distributions
  use bf_kinds, only: dp
  use bf_special, only: log_beta_fn, reg_incomplete_beta, inv_reg_incomplete_beta, &
                        log_choose, binomial_pmf, binomial_cdf_le
  use bf_quadrature, only: integrate_gk
  use bf_random, only: random_beta, random_discrete
  implicit none
  private

  public :: beta4_pdf, beta4_cdf, beta4_quantile, beta4_random
  public :: beta_ms_pdf, beta_ms_cdf, beta_ms_quantile, beta_ms_random
  public :: ams, bms, labmsu, uabmsl, beta_mode, beta_median
  public :: mla, mlb, mlm
  public :: beta_binomial_pmf, beta_binomial_cdf, beta_binomial_random
  public :: gamma_binomial_pdf, gamma_binomial_cdf, gamma_binomial_quantile, gamma_binomial_random
  public :: compound_binomial_pmf, compound_binomial_cdf, compound_binomial_random
  public :: beta_compound_binomial_pmf, beta_compound_binomial_random
  public :: beta_times_binomial_tail, beta_times_gamma_binomial_tail, beta_times_beta_tail

  type :: beta_binom_ctx
    integer :: x = 0
    integer :: n = 0
    real(dp) :: l = 0.0_dp, u = 1.0_dp, alpha = 1.0_dp, beta = 1.0_dp
  end type beta_binom_ctx

  type :: gamma_binom_ctx
    real(dp) :: size = 1.0_dp, prob = 0.5_dp
  end type gamma_binom_ctx

  type :: beta_cbinom_ctx
    integer :: x = 0, n = 0
    real(dp) :: k = 0.0_dp
    real(dp) :: l = 0.0_dp, u = 1.0_dp, alpha = 1.0_dp, beta = 1.0_dp
  end type beta_cbinom_ctx

contains

  pure real(dp) function beta4_pdf(x, l, u, alpha, beta) result(v)
    real(dp), intent(in) :: x, l, u, alpha, beta
    real(dp) :: z, logv

    if (alpha <= 0.0_dp .or. beta <= 0.0_dp .or. u <= l) then
      v = 0.0_dp
      return
    end if
    if (x < l .or. x > u) then
      v = 0.0_dp
      return
    end if
    z = (x - l) / (u - l)
    if (z <= 0.0_dp) then
      if (alpha < 1.0_dp) then
        v = huge(1.0_dp)
      else if (abs(alpha - 1.0_dp) <= epsilon(1.0_dp)) then
        v = exp(-log(u - l) - log_beta_fn(alpha, beta))
      else
        v = 0.0_dp
      end if
      return
    else if (z >= 1.0_dp) then
      if (beta < 1.0_dp) then
        v = huge(1.0_dp)
      else if (abs(beta - 1.0_dp) <= epsilon(1.0_dp)) then
        v = exp(-log(u - l) - log_beta_fn(alpha, beta))
      else
        v = 0.0_dp
      end if
      return
    end if
    logv = (alpha - 1.0_dp) * log(z) + (beta - 1.0_dp) * log(1.0_dp - z) - &
           log_beta_fn(alpha, beta) - log(u - l)
    v = exp(logv)
  end function beta4_pdf

  pure real(dp) function beta4_cdf(q, l, u, alpha, beta, lower_tail) result(v)
    real(dp), intent(in) :: q, l, u, alpha, beta
    logical, intent(in), optional :: lower_tail
    logical :: lower
    real(dp) :: z

    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    if (q <= l) then
      v = 0.0_dp
    else if (q >= u) then
      v = 1.0_dp
    else
      z = (q - l) / (u - l)
      v = reg_incomplete_beta(z, alpha, beta)
    end if
    if (.not. lower) v = 1.0_dp - v
  end function beta4_cdf

  real(dp) function beta4_quantile(p, l, u, alpha, beta, lower_tail) result(x)
    real(dp), intent(in) :: p, l, u, alpha, beta
    logical, intent(in), optional :: lower_tail
    logical :: lower
    real(dp) :: z

    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    z = inv_reg_incomplete_beta(p, alpha, beta)
    if (lower) then
      x = l + (u - l) * z
    else
      ! Compatibility with qBeta.4P(): the upstream code uses 1-qbeta(p),
      ! rather than qbeta(1-p), for the upper-tail branch.
      x = l + (u - l) * (1.0_dp - z)
    end if
  end function beta4_quantile

  subroutine beta4_random(x, l, u, alpha, beta)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: l, u, alpha, beta
    integer :: i
    do i = 1, size(x)
      x(i) = l + (u - l) * random_beta(alpha, beta)
    end do
  end subroutine beta4_random

  pure real(dp) function ams(mean, variance, l, u) result(alpha)
    real(dp), intent(in) :: mean, variance
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    alpha = ((ll - mean) * (ll * (mean - uu) - mean**2 + mean * uu - variance)) / &
            (variance * (ll - uu))
  end function ams

  pure real(dp) function bms(mean, variance, l, u) result(beta)
    real(dp), intent(in) :: mean, variance
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    beta = ((mean - uu) * (ll * (uu - mean) + mean**2 - mean * uu + variance)) / &
           (variance * (uu - ll))
  end function bms

  pure real(dp) function labmsu(alpha, beta, mean, variance, u, have_variance, have_u) result(l)
    real(dp), intent(in) :: alpha, beta, mean, variance, u
    logical, intent(in) :: have_variance, have_u
    if (have_u .and. .not. have_variance) then
      l = (alpha * mean - alpha * u + beta * mean) / beta
    else if (.not. have_u .and. have_variance) then
      l = mean - alpha * sqrt(variance * (alpha + beta + 1.0_dp)) / sqrt(alpha * beta)
    else
      l = (alpha * (mean - u)**3 + beta * variance * (beta * u - mean + u)) / &
          (beta**2 * variance)
    end if
  end function labmsu

  pure real(dp) function uabmsl(alpha, beta, mean, variance, l, have_variance, have_l) result(u)
    real(dp), intent(in) :: alpha, beta, mean, variance, l
    logical, intent(in) :: have_variance, have_l
    if (have_l .and. .not. have_variance) then
      u = beta * (mean - l) / alpha + mean
    else if (.not. have_l .and. have_variance) then
      u = mean + beta * sqrt(variance * (alpha + beta + 1.0_dp)) / sqrt(alpha * beta)
    else
      u = (alpha * variance * (alpha * l + l - mean) - beta * (l - mean)**3) / &
          (alpha**2 * variance)
    end if
  end function uabmsl

  pure real(dp) function beta_mode(alpha, beta, l, u) result(x)
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    x = ((alpha - 1.0_dp) / (alpha + beta - 2.0_dp)) * (uu - ll) + ll
  end function beta_mode

  pure real(dp) function beta_median(alpha, beta, l, u) result(x)
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    x = (alpha - 1.0_dp / 3.0_dp) / (alpha + beta - 2.0_dp / 3.0_dp) * (uu - ll) + ll
  end function beta_median

  pure real(dp) function mla(alpha, beta, x, n, use_xn) result(a)
    real(dp), intent(in) :: alpha, beta, x, n
    logical, intent(in) :: use_xn
    real(dp) :: nn
    if (use_xn) then
      a = x * (n - 2.0_dp) + 1.0_dp
    else
      nn = alpha + beta
      a = (alpha * (nn - 2.0_dp) + nn) / nn
    end if
  end function mla

  pure real(dp) function mlb(alpha, beta, x, n, use_xn) result(b)
    real(dp), intent(in) :: alpha, beta, x, n
    logical, intent(in) :: use_xn
    real(dp) :: nn, a
    if (use_xn) then
      b = n - (x * (n - 2.0_dp) + 1.0_dp)
    else
      nn = alpha + beta
      a = (alpha * (nn - 2.0_dp) + nn) / nn
      b = nn - a
    end if
  end function mlb

  pure real(dp) function mlm(alpha, beta, x, n, use_xn) result(m)
    real(dp), intent(in) :: alpha, beta, x, n
    logical, intent(in) :: use_xn
    real(dp) :: nn, a, b
    if (use_xn) then
      m = (x * (n - 2.0_dp) + 1.0_dp) / n
    else
      nn = alpha + beta
      a = (alpha * (nn - 2.0_dp) + nn) / nn
      b = nn - a
      m = a / (a + b)
    end if
  end function mlm

  pure real(dp) function beta_ms_pdf(x, mean, variance, l, u) result(v)
    real(dp), intent(in) :: x, mean, variance
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    v = beta4_pdf(x, ll, uu, ams(mean, variance, ll, uu), bms(mean, variance, ll, uu))
  end function beta_ms_pdf

  pure real(dp) function beta_ms_cdf(q, mean, variance, l, u, lower_tail) result(v)
    real(dp), intent(in) :: q, mean, variance
    real(dp), intent(in), optional :: l, u
    logical, intent(in), optional :: lower_tail
    real(dp) :: ll, uu
    logical :: lower
    ll = 0.0_dp
    uu = 1.0_dp
    lower = .true.
    if (present(l)) ll = l
    if (present(u)) uu = u
    if (present(lower_tail)) lower = lower_tail
    v = beta4_cdf(q, ll, uu, ams(mean, variance, ll, uu), bms(mean, variance, ll, uu), lower)
  end function beta_ms_cdf

  real(dp) function beta_ms_quantile(p, mean, variance, l, u, lower_tail) result(x)
    real(dp), intent(in) :: p, mean, variance
    real(dp), intent(in), optional :: l, u
    logical, intent(in), optional :: lower_tail
    real(dp) :: ll, uu
    logical :: lower
    ll = 0.0_dp
    uu = 1.0_dp
    lower = .true.
    if (present(l)) ll = l
    if (present(u)) uu = u
    if (present(lower_tail)) lower = lower_tail
    x = beta4_quantile(p, ll, uu, ams(mean, variance, ll, uu), bms(mean, variance, ll, uu), lower)
  end function beta_ms_quantile

  subroutine beta_ms_random(x, mean, variance, l, u)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: mean, variance
    real(dp), intent(in), optional :: l, u
    real(dp) :: ll, uu
    ll = 0.0_dp
    uu = 1.0_dp
    if (present(l)) ll = l
    if (present(u)) uu = u
    call beta4_random(x, ll, uu, ams(mean, variance, ll, uu), bms(mean, variance, ll, uu))
  end subroutine beta_ms_random

  real(dp) function beta_binomial_integrand(y, ctx) result(v)
    real(dp), intent(in) :: y
    class(*), intent(in) :: ctx
    select type (c => ctx)
    type is (beta_binom_ctx)
      v = binomial_pmf(c%x, c%n, y) * beta4_pdf(y, c%l, c%u, c%alpha, c%beta)
    class default
      v = 0.0_dp
    end select
  end function beta_binomial_integrand

  real(dp) function beta_binomial_pmf(x, n, l, u, alpha, beta) result(v)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: l, u, alpha, beta
    type(beta_binom_ctx) :: ctx
    if (x < 0 .or. x > n) then
      v = 0.0_dp
      return
    end if
    ctx = beta_binom_ctx(x, n, l, u, alpha, beta)
    v = integrate_gk(beta_binomial_integrand, ctx, 0.0_dp, 1.0_dp, 1.0e-11_dp, 1.0e-9_dp)
  end function beta_binomial_pmf

  real(dp) function beta_binomial_cdf(q, n, l, u, alpha, beta, lower_tail) result(v)
    integer, intent(in) :: q, n
    real(dp), intent(in) :: l, u, alpha, beta
    logical, intent(in), optional :: lower_tail
    logical :: lower
    integer :: j
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    if (q <= 0) then
      if (lower) then
        v = 0.0_dp
      else
        v = 1.0_dp
      end if
      return
    end if
    v = 0.0_dp
    do j = q, n
      v = v + beta_binomial_pmf(j, n, l, u, alpha, beta)
    end do
    ! Upstream pBetaBinom uses P[X < q] for lower.tail=TRUE.
    if (lower) v = 1.0_dp - v
    v = max(0.0_dp, min(1.0_dp, v))
  end function beta_binomial_cdf

  subroutine beta_binomial_random(x, n, l, u, alpha, beta)
    integer, intent(out) :: x(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: l, u, alpha, beta
    real(dp), allocatable :: w(:)
    integer :: i, j
    allocate(w(0:n))
    do j = 0, n
      w(j) = beta_binomial_pmf(j, n, l, u, alpha, beta)
    end do
    do i = 1, size(x)
      x(i) = random_discrete(w) - 1
    end do
  end subroutine beta_binomial_random

  pure real(dp) function gamma_binomial_kernel(x, size, prob) result(v)
    real(dp), intent(in) :: x, size, prob
    real(dp) :: lv
    if (x < 0.0_dp .or. x > size .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      v = 0.0_dp
      return
    end if
    if (prob <= 0.0_dp) then
      if (x <= 0.0_dp) then
        v = 1.0_dp
      else
        v = 0.0_dp
      end if
      return
    else if (prob >= 1.0_dp) then
      if (x >= size) then
        v = 1.0_dp
      else
        v = 0.0_dp
      end if
      return
    end if
    lv = log_choose(size, x) + x * log(prob) + (size - x) * log(1.0_dp - prob)
    v = exp(lv)
  end function gamma_binomial_kernel

  real(dp) function gamma_binomial_integrand(y, ctx) result(v)
    real(dp), intent(in) :: y
    class(*), intent(in) :: ctx
    select type (c => ctx)
    type is (gamma_binom_ctx)
      v = gamma_binomial_kernel(y, c%size, c%prob)
    class default
      v = 0.0_dp
    end select
  end function gamma_binomial_integrand

  real(dp) function gamma_binomial_pdf(x, size, prob, normalized) result(v)
    real(dp), intent(in) :: x, size, prob
    logical, intent(in), optional :: normalized
    logical :: norm
    type(gamma_binom_ctx) :: ctx
    real(dp) :: den
    norm = .false.
    if (present(normalized)) norm = normalized
    v = gamma_binomial_kernel(x, size, prob)
    if (norm .and. v > 0.0_dp) then
      ctx = gamma_binom_ctx(size, prob)
      den = integrate_gk(gamma_binomial_integrand, ctx, 0.0_dp, size, 1.0e-11_dp, 1.0e-9_dp)
      if (den > 0.0_dp) v = v / den
    end if
  end function gamma_binomial_pdf

  real(dp) function gamma_binomial_cdf(q, size, prob, lower_tail) result(v)
    real(dp), intent(in) :: q, size, prob
    logical, intent(in), optional :: lower_tail
    logical :: lower
    type(gamma_binom_ctx) :: ctx
    real(dp) :: num, den, qq
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    qq = max(0.0_dp, min(size, q))
    ctx = gamma_binom_ctx(size, prob)
    den = integrate_gk(gamma_binomial_integrand, ctx, 0.0_dp, size, 1.0e-11_dp, 1.0e-9_dp)
    num = integrate_gk(gamma_binomial_integrand, ctx, 0.0_dp, qq, 1.0e-11_dp, 1.0e-9_dp)
    if (den <= 0.0_dp) then
      v = 0.0_dp
    else
      v = num / den
    end if
    if (.not. lower) v = 1.0_dp - v
    v = max(0.0_dp, min(1.0_dp, v))
  end function gamma_binomial_cdf

  real(dp) function gamma_binomial_quantile(p, size, prob, lower_tail, precision) result(x)
    real(dp), intent(in) :: p, size, prob
    logical, intent(in), optional :: lower_tail
    real(dp), intent(in), optional :: precision
    logical :: lower
    real(dp) :: lo, hi, mid, y, tol
    integer :: iter
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    tol = 1.0e-7_dp
    if (present(precision)) tol = precision
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = size
      return
    end if
    lo = 0.0_dp
    hi = size
    do iter = 1, 200
      mid = 0.5_dp * (lo + hi)
      y = gamma_binomial_cdf(mid, size, prob, lower)
      if (abs(y - p) <= tol) exit
      if (lower) then
        if (y < p) then
          lo = mid
        else
          hi = mid
        end if
      else
        if (y < p) then
          hi = mid
        else
          lo = mid
        end if
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function gamma_binomial_quantile

  subroutine gamma_binomial_random(x, size_par, prob, precision)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: size_par, prob
    real(dp), intent(in), optional :: precision
    real(dp) :: u, tol
    integer :: i
    tol = 1.0e-4_dp
    if (present(precision)) tol = precision
    do i = 1, size(x)
      call random_number(u)
      x(i) = gamma_binomial_quantile(u, size_par, prob, .true., tol)
    end do
  end subroutine gamma_binomial_random

  pure real(dp) function compound_binomial_pmf(x, n, k, p) result(v)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: k, p
    v = binomial_pmf(x, n, p) - k * p * (1.0_dp - p) * &
        (binomial_pmf(x, n - 2, p) - 2.0_dp * binomial_pmf(x - 1, n - 2, p) + &
         binomial_pmf(x - 2, n - 2, p))
  end function compound_binomial_pmf

  pure real(dp) function compound_binomial_cdf(q, n, k, p, lower_tail) result(v)
    integer, intent(in) :: q, n
    real(dp), intent(in) :: k, p
    logical, intent(in), optional :: lower_tail
    logical :: lower
    integer :: j
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    v = 0.0_dp
    do j = max(q, 0), n
      v = v + compound_binomial_pmf(j, n, k, p)
    end do
    if (lower) v = 1.0_dp - v
  end function compound_binomial_cdf

  subroutine compound_binomial_random(x, n, k, p)
    integer, intent(out) :: x(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: k
    real(dp), intent(in) :: p(:)
    real(dp), allocatable :: w(:)
    real(dp) :: shift
    integer :: i, j, ip
    allocate(w(0:n))
    do i = 1, size(x)
      ip = min(i, size(p))
      do j = 0, n
        w(j) = compound_binomial_pmf(j, n, k, p(ip))
      end do
      if (minval(w) < 0.0_dp) then
        shift = abs(minval(w))
        w = w + shift
      end if
      if (sum(w) > 0.0_dp) w = w / sum(w)
      x(i) = random_discrete(w) - 1
    end do
  end subroutine compound_binomial_random

  real(dp) function beta_cbinom_integrand(y, ctx) result(v)
    real(dp), intent(in) :: y
    class(*), intent(in) :: ctx
    select type (c => ctx)
    type is (beta_cbinom_ctx)
      v = compound_binomial_pmf(c%x, c%n, c%k, y) * beta4_pdf(y, c%l, c%u, c%alpha, c%beta)
    class default
      v = 0.0_dp
    end select
  end function beta_cbinom_integrand

  real(dp) function beta_compound_binomial_pmf(x, n, k, l, u, alpha, beta) result(v)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: k, l, u, alpha, beta
    type(beta_cbinom_ctx) :: ctx
    if (x < 0 .or. x > n) then
      v = 0.0_dp
      return
    end if
    ctx = beta_cbinom_ctx(x, n, k, l, u, alpha, beta)
    v = integrate_gk(beta_cbinom_integrand, ctx, 0.0_dp, 1.0_dp, 1.0e-11_dp, 1.0e-9_dp)
  end function beta_compound_binomial_pmf

  subroutine beta_compound_binomial_random(x, n, k, l, u, alpha, beta)
    integer, intent(out) :: x(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: k, l, u, alpha, beta
    real(dp), allocatable :: w(:)
    real(dp) :: shift
    integer :: i, j
    allocate(w(0:n))
    do j = 0, n
      w(j) = beta_compound_binomial_pmf(j, n, k, l, u, alpha, beta)
    end do
    if (minval(w) < 0.0_dp) then
      shift = abs(minval(w))
      w = w + shift
    end if
    if (sum(w) > 0.0_dp) w = w / sum(w)
    do i = 1, size(x)
      x(i) = random_discrete(w) - 1
    end do
  end subroutine beta_compound_binomial_random

  real(dp) function beta_times_binomial_tail(x, l, u, alpha, beta, n, c, lower_tail) result(v)
    real(dp), intent(in) :: x, l, u, alpha, beta, c
    integer, intent(in) :: n
    logical, intent(in), optional :: lower_tail
    logical :: lower
    real(dp) :: pc
    integer :: kcut
    lower = .false.
    if (present(lower_tail)) lower = lower_tail
    kcut = nint(real(n, dp) * c) - 1
    pc = binomial_cdf_le(kcut, n, x)
    if (lower) then
      v = beta4_pdf(x, l, u, alpha, beta) * pc
    else
      v = beta4_pdf(x, l, u, alpha, beta) * (1.0_dp - pc)
    end if
  end function beta_times_binomial_tail

  real(dp) function beta_times_gamma_binomial_tail(x, l, u, alpha, beta, n, c, lower_tail) result(v)
    real(dp), intent(in) :: x, l, u, alpha, beta, n, c
    logical, intent(in), optional :: lower_tail
    logical :: lower
    real(dp) :: pc
    lower = .false.
    if (present(lower_tail)) lower = lower_tail
    pc = gamma_binomial_cdf(n * c, n, x, lower)
    v = beta4_pdf(x, l, u, alpha, beta) * pc
  end function beta_times_gamma_binomial_tail

  real(dp) function beta_times_beta_tail(x, l, u, alpha, beta, n, c, lower_tail) result(v)
    real(dp), intent(in) :: x, l, u, alpha, beta, n, c
    logical, intent(in), optional :: lower_tail
    logical :: lower
    real(dp) :: pc
    lower = .false.
    if (present(lower_tail)) lower = lower_tail
    if (x <= 0.0_dp .or. x >= 1.0_dp) then
      pc = merge(0.0_dp, 1.0_dp, x <= 0.0_dp)
    else
      pc = reg_incomplete_beta(c, x * n, (1.0_dp - x) * n)
      if (.not. lower) pc = 1.0_dp - pc
    end if
    if (lower) then
      v = beta4_pdf(x, l, u, alpha, beta) * pc
    else
      v = beta4_pdf(x, l, u, alpha, beta) * pc
    end if
  end function beta_times_beta_tail

end module bf_distributions
