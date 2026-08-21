! SPDX-License-Identifier: CC0-1.0
module bf_special
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bf_kinds, only: dp
  implicit none
  private
  public :: log_beta_fn, beta_fn, reg_incomplete_beta, inv_reg_incomplete_beta
  public :: log_choose, choose_general, gchoose, binomial_pmf, binomial_cdf_le, binomial_prob_range
  public :: reg_lower_gamma, reg_upper_gamma, chi_square_sf

contains

  pure real(dp) function log_beta_fn(a, b) result(v)
    real(dp), intent(in) :: a, b
    v = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
  end function log_beta_fn

  pure real(dp) function beta_fn(a, b) result(v)
    real(dp), intent(in) :: a, b
    v = exp(log_beta_fn(a, b))
  end function beta_fn

  pure real(dp) function beta_cf(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
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

  pure real(dp) function reg_incomplete_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt

    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
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
  end function reg_incomplete_beta

  real(dp) function inv_reg_incomplete_beta(p, a, b, tol) result(x)
    real(dp), intent(in) :: p, a, b
    real(dp), intent(in), optional :: tol
    real(dp) :: lo, hi, mid, target, eps
    integer :: iter

    target = max(0.0_dp, min(1.0_dp, p))
    if (target <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (target >= 1.0_dp) then
      x = 1.0_dp
      return
    end if
    eps = 1.0e-12_dp
    if (present(tol)) eps = max(tol, 10.0_dp * epsilon(1.0_dp))

    lo = 0.0_dp
    hi = 1.0_dp
    do iter = 1, 200
      mid = 0.5_dp * (lo + hi)
      if (reg_incomplete_beta(mid, a, b) < target) then
        lo = mid
      else
        hi = mid
      end if
      if (hi - lo <= eps * max(1.0_dp, abs(mid))) exit
    end do
    x = 0.5_dp * (lo + hi)
  end function inv_reg_incomplete_beta

  pure real(dp) function log_choose(n, k) result(v)
    real(dp), intent(in) :: n, k
    if (k < 0.0_dp .or. k > n) then
      v = -huge(1.0_dp)
    else
      v = log_gamma(n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(n - k + 1.0_dp)
    end if
  end function log_choose

  pure real(dp) function choose_general(n, k) result(v)
    real(dp), intent(in) :: n, k
    real(dp) :: lv
    lv = log_choose(n, k)
    if (lv <= -0.5_dp * huge(1.0_dp)) then
      v = 0.0_dp
    else
      v = exp(lv)
    end if
  end function choose_general

  pure real(dp) function gchoose(n, k) result(v)
    real(dp), intent(in) :: n, k
    v = choose_general(n, k)
  end function gchoose

  pure real(dp) function binomial_pmf(k, n, p) result(v)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p
    real(dp) :: lp

    if (k < 0 .or. k > n .or. p < 0.0_dp .or. p > 1.0_dp) then
      v = 0.0_dp
      return
    end if
    if (p <= 0.0_dp) then
      v = merge(1.0_dp, 0.0_dp, k == 0)
      return
    else if (p >= 1.0_dp) then
      v = merge(1.0_dp, 0.0_dp, k == n)
      return
    end if
    lp = log_choose(real(n, dp), real(k, dp)) + real(k, dp) * log(p) + &
         real(n - k, dp) * log(1.0_dp - p)
    v = exp(lp)
  end function binomial_pmf

  pure real(dp) function binomial_cdf_le(k, n, p) result(v)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p
    integer :: j

    if (k < 0) then
      v = 0.0_dp
      return
    else if (k >= n) then
      v = 1.0_dp
      return
    end if
    if (p <= 0.0_dp) then
      v = 1.0_dp
      return
    else if (p >= 1.0_dp) then
      v = 0.0_dp
      return
    end if

    ! Direct summation is robust for the moderate test lengths targeted here.
    v = 0.0_dp
    do j = 0, k
      v = v + binomial_pmf(j, n, p)
    end do
    v = max(0.0_dp, min(1.0_dp, v))
  end function binomial_cdf_le

  pure real(dp) function binomial_prob_range(lo, hi, n, p) result(v)
    integer, intent(in) :: lo, hi, n
    real(dp), intent(in) :: p
    if (hi < lo) then
      v = 0.0_dp
    else
      v = binomial_cdf_le(hi, n, p) - binomial_cdf_le(lo - 1, n, p)
      v = max(0.0_dp, min(1.0_dp, v))
    end if
  end function binomial_prob_range

  pure real(dp) function reg_lower_gamma(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit = 500
    real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an

    if (a <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if

    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do n = 1, maxit
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) < abs(sumv) * eps) exit
      end do
      p = sumv * exp(-x + a * log(x) - log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / tiny(1.0_dp)
      d = 1.0_dp / b
      h = d
      do n = 1, maxit
        an = -real(n, dp) * (real(n, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
        c = b + an / c
        if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) < eps) exit
      end do
      p = 1.0_dp - exp(-x + a * log(x) - log_gamma(a)) * h
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_lower_gamma

  pure real(dp) function reg_upper_gamma(a, x) result(q)
    real(dp), intent(in) :: a, x
    q = 1.0_dp - reg_lower_gamma(a, x)
    q = max(0.0_dp, min(1.0_dp, q))
  end function reg_upper_gamma

  pure real(dp) function chi_square_sf(x, df) result(q)
    real(dp), intent(in) :: x
    integer, intent(in) :: df
    if (x <= 0.0_dp) then
      q = 1.0_dp
    else if (df <= 0) then
      q = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      q = reg_upper_gamma(0.5_dp * real(df, dp), 0.5_dp * x)
    end if
  end function chi_square_sf

end module bf_special
