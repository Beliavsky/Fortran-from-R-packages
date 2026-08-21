! SPDX-License-Identifier: CC0-1.0
module bf_random
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bf_kinds, only: dp
  implicit none
  private
  public :: random_normal, random_gamma, random_beta, random_binomial, random_discrete

contains

  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
  end function random_normal

  recursive real(dp) function random_gamma(shape) result(g)
    real(dp), intent(in) :: shape
    real(dp) :: d, c, x, v, u

    if (shape <= 0.0_dp) then
      g = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      g = random_gamma(shape + 1.0_dp) * u**(1.0_dp / shape)
      return
    end if

    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      x = random_normal()
      v = 1.0_dp + c * x
      if (v <= 0.0_dp) cycle
      v = v**3
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp * x**4) exit
      if (log(u) < 0.5_dp * x*x + d * (1.0_dp - v + log(v))) exit
    end do
    g = d * v
  end function random_gamma

  real(dp) function random_beta(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: g1, g2
    g1 = random_gamma(a)
    g2 = random_gamma(b)
    x = g1 / (g1 + g2)
  end function random_beta

  integer function random_binomial(n, p) result(k)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: i
    real(dp) :: u
    k = 0
    do i = 1, n
      call random_number(u)
      if (u < p) k = k + 1
    end do
  end function random_binomial

  integer function random_discrete(weights) result(idx)
    real(dp), intent(in) :: weights(:)
    real(dp) :: u, s, total
    integer :: i

    total = sum(max(weights, 0.0_dp))
    if (total <= 0.0_dp) then
      idx = 1
      return
    end if
    call random_number(u)
    u = u * total
    s = 0.0_dp
    do i = 1, size(weights)
      s = s + max(weights(i), 0.0_dp)
      if (u <= s) then
        idx = i
        return
      end if
    end do
    idx = size(weights)
  end function random_discrete

end module bf_random
