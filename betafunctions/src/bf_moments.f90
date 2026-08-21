! SPDX-License-Identifier: CC0-1.0
module bf_moments
  use bf_kinds, only: dp
  use bf_distributions, only: beta_binomial_pmf, ams, bms
  use bf_special, only: binomial_pmf
  implicit none
  private
  public :: moment_set, beta_params
  public :: beta_moments, binomial_moments, beta_binomial_moments, observed_moments
  public :: beta4_fit, beta2_fit, descending_factorial, ascending_factorial
  public :: tsm, hb_tsm, beta_true_score_fit, beta_true_score_moments, hb_beta_true_score_fit

  type :: moment_set
    real(dp), allocatable :: raw(:)
    real(dp), allocatable :: central(:)
    real(dp), allocatable :: standardized(:)
  end type moment_set

  type :: beta_params
    real(dp) :: alpha = 0.0_dp
    real(dp) :: beta = 0.0_dp
    real(dp) :: l = 0.0_dp
    real(dp) :: u = 1.0_dp
    real(dp) :: etl = 0.0_dp
    real(dp) :: k = 0.0_dp
    real(dp) :: n = 0.0_dp
    logical :: used_failsafe = .false.
  end type beta_params

contains

  pure real(dp) function int_choose(n, k) result(v)
    integer, intent(in) :: n, k
    integer :: j, kk
    v = 1.0_dp
    if (k < 0 .or. k > n) then
      v = 0.0_dp
      return
    end if
    kk = min(k, n - k)
    do j = 1, kk
      v = v * real(n - kk + j, dp) / real(j, dp)
    end do
  end function int_choose

  pure real(dp) function rising_ratio(a, b, j) result(v)
    real(dp), intent(in) :: a, b
    integer, intent(in) :: j
    integer :: i
    v = 1.0_dp
    do i = 0, j - 1
      v = v * (a + real(i, dp)) / (b + real(i, dp))
    end do
  end function rising_ratio

  subroutine raw_to_all(raw, out)
    real(dp), intent(in) :: raw(:)
    type(moment_set), intent(out) :: out
    integer :: r, j, m
    real(dp) :: mu, var

    m = size(raw)
    allocate(out%raw(m), out%central(m), out%standardized(m))
    out%raw = raw
    mu = raw(1)
    do r = 1, m
      out%central(r) = (-mu)**r
      do j = 1, r
        out%central(r) = out%central(r) + int_choose(r, j) * (-mu)**(r - j) * raw(j)
      end do
    end do
    var = 0.0_dp
    if (m >= 2) var = out%central(2)
    do r = 1, m
      if (var > 0.0_dp) then
        out%standardized(r) = out%central(r) / sqrt(var)**r
      else
        out%standardized(r) = 0.0_dp
      end if
    end do
  end subroutine raw_to_all

  subroutine beta_moments(alpha, beta, l, u, orders, out)
    real(dp), intent(in) :: alpha, beta, l, u
    integer, intent(in) :: orders
    type(moment_set), intent(out) :: out
    real(dp), allocatable :: raw(:)
    integer :: r, j

    allocate(raw(orders))
    do r = 1, orders
      raw(r) = 0.0_dp
      do j = 0, r
        raw(r) = raw(r) + int_choose(r, j) * l**(r - j) * (u - l)**j * &
                 rising_ratio(alpha, alpha + beta, j)
      end do
    end do
    call raw_to_all(raw, out)
  end subroutine beta_moments

  subroutine binomial_moments(n, p, orders, out)
    integer, intent(in) :: n, orders
    real(dp), intent(in) :: p
    type(moment_set), intent(out) :: out
    real(dp), allocatable :: raw(:)
    integer :: r, k

    allocate(raw(orders))
    raw = 0.0_dp
    do r = 1, orders
      do k = 0, n
        raw(r) = raw(r) + real(k, dp)**r * binomial_pmf(k, n, p)
      end do
    end do
    call raw_to_all(raw, out)
  end subroutine binomial_moments

  subroutine beta_binomial_moments(n, l, u, alpha, beta, orders, out)
    integer, intent(in) :: n, orders
    real(dp), intent(in) :: l, u, alpha, beta
    type(moment_set), intent(out) :: out
    real(dp), allocatable :: raw(:), weights(:)
    integer :: r, k

    allocate(raw(orders), weights(0:n))
    do k = 0, n
      weights(k) = beta_binomial_pmf(k, n, l, u, alpha, beta)
    end do
    raw = 0.0_dp
    do r = 1, orders
      do k = 0, n
        raw(r) = raw(r) + real(k, dp)**r * weights(k)
      end do
    end do
    call raw_to_all(raw, out)
  end subroutine beta_binomial_moments

  subroutine observed_moments(x, orders, correct, out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: orders
    logical, intent(in), optional :: correct
    type(moment_set), intent(out) :: out
    logical :: corr
    real(dp) :: mu, svar, denom
    integer :: r, n

    corr = .true.
    if (present(correct)) corr = correct
    n = size(x)
    allocate(out%raw(orders), out%central(orders), out%standardized(orders))
    mu = sum(x) / real(n, dp)
    if (n > 1) then
      svar = sum((x - mu)**2) / real(n - 1, dp)
    else
      svar = 0.0_dp
    end if
    do r = 1, orders
      out%raw(r) = sum(x**r) / real(n, dp)
      denom = real(n, dp)
      if (corr) denom = real(max(n - 1, 1), dp)
      out%central(r) = sum((x - mu)**r) / denom
      if (svar > 0.0_dp) then
        out%standardized(r) = sum((x - mu)**r / sqrt(svar)**r) / real(n, dp)
      else
        out%standardized(r) = 0.0_dp
      end if
    end do
  end subroutine observed_moments

  function beta4_fit(scores, mean_value, variance, skewness, kurtosis, use_moments) result(par)
    real(dp), intent(in) :: scores(:)
    real(dp), intent(in) :: mean_value, variance, skewness, kurtosis
    logical, intent(in) :: use_moments
    type(beta_params) :: par
    type(moment_set) :: m
    real(dp) :: m1, s2, g3, g4, r, disc

    if (use_moments) then
      m1 = mean_value
      s2 = variance
      g3 = skewness
      g4 = kurtosis
    else
      call observed_moments(scores, 4, .true., m)
      m1 = m%raw(1)
      s2 = m%central(2)
      g3 = m%standardized(3)
      g4 = m%standardized(4)
    end if

    r = 6.0_dp * (g4 - g3**2 - 1.0_dp) / (6.0_dp + 3.0_dp * g3**2 - 2.0_dp * g4)
    disc = 1.0_dp - 24.0_dp * (r + 1.0_dp) / &
           ((r + 2.0_dp) * (r + 3.0_dp) * g4 - 3.0_dp * (r - 6.0_dp) * (r + 1.0_dp))
    if (g3 < 0.0_dp) then
      par%alpha = 0.5_dp * r * (1.0_dp + sqrt(disc))
      par%beta = 0.5_dp * r * (1.0_dp - sqrt(disc))
    else
      par%beta = 0.5_dp * r * (1.0_dp + sqrt(disc))
      par%alpha = 0.5_dp * r * (1.0_dp - sqrt(disc))
    end if
    par%l = m1 - par%alpha * sqrt(s2 * (par%alpha + par%beta + 1.0_dp)) / &
            sqrt(par%alpha * par%beta)
    par%u = m1 + par%beta * sqrt(s2 * (par%alpha + par%beta + 1.0_dp)) / &
            sqrt(par%alpha * par%beta)
  end function beta4_fit

  function beta2_fit(scores, mean_value, variance, l, u, use_scores) result(par)
    real(dp), intent(in) :: scores(:)
    real(dp), intent(in) :: mean_value, variance, l, u
    logical, intent(in) :: use_scores
    type(beta_params) :: par
    real(dp) :: m1, s2
    integer :: n

    if (use_scores) then
      n = size(scores)
      m1 = sum(scores) / real(n, dp)
      s2 = sum((scores - m1)**2) / real(max(n - 1, 1), dp)
    else
      m1 = mean_value
      s2 = variance
    end if
    par%alpha = ams(m1, s2, l, u)
    par%beta = bms(m1, s2, l, u)
    par%l = l
    par%u = u
  end function beta2_fit

  pure real(dp) function descending_factorial(x, r) result(v)
    real(dp), intent(in) :: x
    integer, intent(in) :: r
    integer :: i
    if (x < real(r, dp)) then
      v = 0.0_dp
      return
    end if
    if (r <= 1) then
      v = x**r
      return
    end if
    v = 1.0_dp
    do i = 0, r - 1
      v = v * (x - real(i, dp))
    end do
  end function descending_factorial

  pure real(dp) function ascending_factorial(x, r) result(v)
    real(dp), intent(in) :: x
    integer, intent(in) :: r
    integer :: i
    if (r <= 1) then
      v = x**r
      return
    end if
    v = 1.0_dp
    do i = 0, r - 1
      v = v * (x + real(i, dp))
    end do
  end function ascending_factorial

  real(dp) function tsm(x, r, n) result(v)
    real(dp), intent(in) :: x(:), n
    integer, intent(in) :: r
    real(dp) :: sx
    integer :: i
    if (r == 1) then
      v = sum(x) / real(size(x), dp) / n
      return
    end if
    sx = 0.0_dp
    do i = 1, size(x)
      sx = sx + descending_factorial(x(i), r)
    end do
    sx = sx / real(size(x), dp)
    v = (sx / descending_factorial(n - 2.0_dp, r - 2)) / descending_factorial(n, 2)
  end function tsm

  subroutine hb_tsm(x, r, n, k, m)
    real(dp), intent(in) :: x(:), n, k
    integer, intent(in) :: r
    real(dp), intent(out) :: m(r)
    real(dp) :: sx, d2
    integer :: i, j

    if (r <= 0) return
    m(1) = sum(x) / real(size(x), dp) / n
    do i = 2, r
      sx = 0.0_dp
      do j = 1, size(x)
        sx = sx + descending_factorial(x(j), i)
      end do
      sx = sx / real(size(x), dp)
      d2 = descending_factorial(real(i, dp), 2)
      m(i) = ((sx / descending_factorial(n - 2.0_dp, i - 2)) + k * d2 * m(i - 1)) / &
             (descending_factorial(n, 2) + k * d2)
    end do
  end subroutine hb_tsm

  subroutine beta_true_score_moments(x, min_value, max_value, etl, reliability, have_etl, out, n_eff)
    real(dp), intent(in) :: x(:), min_value, max_value, etl, reliability
    logical, intent(in) :: have_etl
    type(moment_set), intent(out) :: out
    real(dp), intent(out) :: n_eff
    real(dp), allocatable :: xs(:)
    real(dp) :: meanx, varx, m(4), s2
    integer :: nobs, i

    nobs = size(x)
    meanx = sum(x) / real(nobs, dp)
    varx = sum((x - meanx)**2) / real(max(nobs - 1, 1), dp)
    if (have_etl) then
      n_eff = etl
    else
      n_eff = ((meanx - min_value) * (max_value - meanx) - reliability * varx) / &
              (varx * (1.0_dp - reliability))
    end if
    allocate(xs(nobs))
    xs = (x - min_value) / (max_value - min_value) * n_eff
    call hb_tsm(xs, 4, n_eff, 0.0_dp, m)
    allocate(out%raw(4), out%central(4), out%standardized(4))
    out%raw = m
    out%central(1) = 0.0_dp
    out%central(2) = m(2) - m(1)**2
    out%central(3) = m(3) - 3.0_dp * m(1) * m(2) + 2.0_dp * m(1)**3
    out%central(4) = m(4) - 4.0_dp * m(1) * m(3) + 6.0_dp * m(1)**2 * m(2) - 3.0_dp * m(1)**4
    s2 = out%central(2)
    out%standardized = 0.0_dp
    if (s2 > 0.0_dp) then
      do i = 2, 4
        out%standardized(i) = out%central(i) / sqrt(s2)**i
      end do
    end if
  end subroutine beta_true_score_moments

  function beta_true_score_fit(x, min_value, max_value, etl, reliability, have_etl, &
                               four_parameter, failsafe, l, u) result(par)
    real(dp), intent(in) :: x(:), min_value, max_value, etl, reliability, l, u
    logical, intent(in) :: have_etl, four_parameter, failsafe
    type(beta_params) :: par
    real(dp), allocatable :: xs(:)
    real(dp) :: n_eff, meanx, varx, s2, g3, g4
    real(dp) :: m(4)
    type(beta_params) :: p4
    integer :: nobs

    nobs = size(x)
    meanx = sum(x) / real(nobs, dp)
    varx = sum((x - meanx)**2) / real(max(nobs - 1, 1), dp)
    if (have_etl) then
      n_eff = etl
    else
      n_eff = ((meanx - min_value) * (max_value - meanx) - reliability * varx) / &
              (varx * (1.0_dp - reliability))
    end if
    allocate(xs(nobs))
    xs = (x - min_value) / (max_value - min_value) * n_eff
    call hb_tsm(xs, 4, n_eff, 0.0_dp, m)
    s2 = m(2) - m(1)**2
    if (s2 <= 0.0_dp) then
      par = beta2_fit(xs(1:0), m(1), max(s2, tiny(1.0_dp)), l, u, .false.)
      par%etl = n_eff
      return
    end if
    g3 = (m(3) - 3.0_dp * m(1) * m(2) + 2.0_dp * m(1)**3) / sqrt(s2)**3
    g4 = (m(4) - 4.0_dp * m(1) * m(3) + 6.0_dp * m(1)**2 * m(2) - 3.0_dp * m(1)**4) / sqrt(s2)**4
    if (four_parameter) then
      p4 = beta4_fit(xs(1:0), m(1), s2, g3, g4, .true.)
      if (failsafe .and. (p4%l < 0.0_dp .or. p4%u > 1.0_dp)) then
        par = beta2_fit(xs(1:0), m(1), s2, l, u, .false.)
        par%used_failsafe = .true.
      else
        par = p4
      end if
    else
      par = beta2_fit(xs(1:0), m(1), s2, l, u, .false.)
    end if
    par%etl = n_eff
  end function beta_true_score_fit

  function hb_beta_true_score_fit(x, n, k, four_parameter, failsafe, l, u) result(par)
    real(dp), intent(in) :: x(:), n, k, l, u
    logical, intent(in) :: four_parameter, failsafe
    type(beta_params) :: par, p4
    real(dp) :: m(4), s2, g3, g4

    call hb_tsm(x, 4, n, k, m)
    s2 = m(2) - m(1)**2
    g3 = (m(3) - 3.0_dp * m(1) * m(2) + 2.0_dp * m(1)**3) / sqrt(s2)**3
    g4 = (m(4) - 4.0_dp * m(1) * m(3) + 6.0_dp * m(1)**2 * m(2) - 3.0_dp * m(1)**4) / sqrt(s2)**4
    if (four_parameter) then
      p4 = beta4_fit(x(1:0), m(1), s2, g3, g4, .true.)
      if (failsafe .and. (p4%l < 0.0_dp .or. p4%u > 1.0_dp)) then
        par = beta2_fit(x(1:0), m(1), s2, l, u, .false.)
        par%used_failsafe = .true.
      else
        par = p4
      end if
    else
      par = beta2_fit(x(1:0), m(1), s2, l, u, .false.)
    end if
    par%k = k
    par%n = n
  end function hb_beta_true_score_fit

end module bf_moments
