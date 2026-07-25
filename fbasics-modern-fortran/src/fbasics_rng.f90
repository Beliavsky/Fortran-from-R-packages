! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_rng
  use fbasics_kinds, only: dp, pi
  implicit none
  private
  integer(kind=8), save :: lcg_state = 123456789_8
  logical, save :: have_spare = .false.
  real(dp), save :: spare_normal = 0.0_dp
  public :: set_lcg_seed, get_lcg_seed, runif_lcg, rnorm_lcg, rt_lcg
  public :: random_gamma, random_chisq, random_inverse_gaussian
contains
  subroutine set_lcg_seed(seed)
    integer(kind=8), intent(in) :: seed
    lcg_state = modulo(max(abs(seed),1_8)-1_8, 2147483646_8) + 1_8
    have_spare = .false.
  end subroutine set_lcg_seed

  integer(kind=8) function get_lcg_seed() result(seed)
    seed = lcg_state
  end function get_lcg_seed

  real(dp) function runif_lcg() result(u)
    integer(kind=8), parameter :: a = 16807_8, m = 2147483647_8
    lcg_state = modulo(a * lcg_state, m)
    u = real(lcg_state, dp) / real(m, dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function runif_lcg

  real(dp) function rnorm_lcg() result(z)
    real(dp) :: u1, u2, r
    if (have_spare) then
      z = spare_normal
      have_spare = .false.
      return
    end if
    u1 = runif_lcg()
    u2 = runif_lcg()
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp*pi*u2)
    spare_normal = r * sin(2.0_dp*pi*u2)
    have_spare = .true.
  end function rnorm_lcg

  recursive real(dp) function random_gamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = random_gamma(shape + 1.0_dp, scale) * runif_lcg()**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rnorm_lcg()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v**3
      u = runif_lcg()
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp - v + log(v))) exit
    end do
    x = scale*d*v
  end function random_gamma

  real(dp) function random_chisq(df) result(x)
    real(dp), intent(in) :: df
    x = random_gamma(0.5_dp*df, 2.0_dp)
  end function random_chisq

  real(dp) function rt_lcg(df) result(x)
    real(dp), intent(in) :: df
    x = rnorm_lcg()/sqrt(random_chisq(df)/df)
  end function rt_lcg

  real(dp) function random_inverse_gaussian(mu, lambda) result(x)
    real(dp), intent(in) :: mu, lambda
    real(dp) :: v, y, candidate, u
    if (mu <= 0.0_dp .or. lambda <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    v = rnorm_lcg()
    y = v*v
    candidate = mu + (mu*mu*y)/(2.0_dp*lambda) - &
      (mu/(2.0_dp*lambda))*sqrt(4.0_dp*mu*lambda*y + mu*mu*y*y)
    u = runif_lcg()
    if (u <= mu/(mu + candidate)) then
      x = candidate
    else
      x = mu*mu/candidate
    end if
  end function random_inverse_gaussian
end module fbasics_rng
