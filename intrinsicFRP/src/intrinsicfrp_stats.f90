! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_stats
  use intrinsicfrp_kinds, only: dp, i8
  implicit none
  private
  public :: normal_cdf, normal_quantile, chi_square_cdf, chi_square_quantile
  public :: gamma_p, rng_state, rng_seed, rng_uniform, rng_normal

  type, public :: rng_state
    integer(i8) :: state = 88172645463393265_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

contains

  pure function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376d1, 2.209460984245205d2, -2.759285104469687d2, &
       1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0 ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406d1, 1.615858368580409d2, -1.556989798598866d2, &
       6.680131188771972d1, -1.328068155288572d1 ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293d-3, -3.223964580411365d-1, -2.400758277161838d0, &
      -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0 ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462d-3, 3.224671290700398d-1, 2.445134137142996d0, &
       3.754408661907416d0 ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
          ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
          (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - p))
      x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
           ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
    end if
  end function normal_quantile

  function gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p
    real(dp) :: sumv, del, ap, b, c, d, h, an, gln
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x <= tiny(1.0_dp)) then
      p = 0.0_dp
      return
    end if
    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do n = 1, 10000
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv) * 10.0_dp * epsilon(1.0_dp)) exit
      end do
      p = sumv * exp(-x + a * log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / tiny(1.0_dp)
      d = 1.0_dp / b
      h = d
      do n = 1, 10000
        an = -real(n, dp) * (real(n, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
        c = b + an / c
        if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) exit
      end do
      p = 1.0_dp - exp(-x + a * log(x) - gln) * h
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function gamma_p

  function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: p
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = gamma_p(0.5_dp * df, 0.5_dp * x)
    end if
  end function chi_square_cdf

  function chi_square_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: x, lo, hi, z
    integer :: iter
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    z = normal_quantile(p)
    x = df * max(0.05_dp, (1.0_dp - 2.0_dp / (9.0_dp * df) + &
      z * sqrt(2.0_dp / (9.0_dp * df))) ** 3)
    lo = 0.0_dp
    hi = max(1.0_dp, 2.0_dp * x)
    do while (chi_square_cdf(hi, df) < p)
      hi = 2.0_dp * hi
      if (hi > huge(1.0_dp) / 4.0_dp) exit
    end do
    do iter = 1, 100
      x = 0.5_dp * (lo + hi)
      if (chi_square_cdf(x, df) < p) then
        lo = x
      else
        hi = x
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function chi_square_quantile

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    if (seed == 0_i8) then
      rng%state = 88172645463393265_i8
    else
      rng%state = seed
    end if
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(i8) :: x
    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
    u = real(iand(x, int(z'001FFFFFFFFFFFFF', i8)), dp) / &
      real(int(z'0020000000000000', i8), dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z
    real(dp) :: u1, u2, r
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp * acos(-1.0_dp) * u2)
    rng%spare = r * sin(2.0_dp * acos(-1.0_dp) * u2)
    rng%has_spare = .true.
  end function rng_normal

end module intrinsicfrp_stats
