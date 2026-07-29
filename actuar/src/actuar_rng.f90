! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_rng
  use actuar_kinds, only : dp, i8, pi
  use actuar_special, only : nan_dp
  implicit none
  private
  public :: seed_rng, runif, rnorm, rexp, rgamma_shape, rbeta_ab
  public :: rpois, rbinom, rnbinom, rinvgauss_rng

  integer(i8), save :: state = 88172645463393265_i8
  logical, save :: have_spare = .false.
  real(dp), save :: spare = 0.0_dp
contains

  subroutine seed_rng(seed)
    integer(i8), intent(in) :: seed
    state = merge(seed, 88172645463393265_i8, seed /= 0_i8)
    have_spare = .false.
  end subroutine seed_rng

  function next_u64() result(x)
    integer(i8) :: x
    x = state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    state = x
  end function next_u64

  function runif() result(u)
    real(dp) :: u
    integer(i8) :: x
    x = iand(next_u64(), int(z'7FFFFFFFFFFFFFFF', i8))
    u = (real(x, dp) + 0.5_dp) / (real(huge(0_i8), dp) + 1.0_dp)
    u = max(tiny(1.0_dp), min(1.0_dp-tiny(1.0_dp), u))
  end function runif

  function rnorm() result(z)
    real(dp) :: z, u, v, s
    if (have_spare) then
      z = spare
      have_spare = .false.
      return
    end if
    do
      u = 2.0_dp*runif()-1.0_dp
      v = 2.0_dp*runif()-1.0_dp
      s = u*u+v*v
      if (s > 0.0_dp .and. s < 1.0_dp) exit
    end do
    s = sqrt(-2.0_dp*log(s)/s)
    z = u*s
    spare = v*s
    have_spare = .true.
  end function rnorm

  function rexp(rate) result(x)
    real(dp), intent(in), optional :: rate
    real(dp) :: x, r
    r = 1.0_dp
    if (present(rate)) r = rate
    if (r <= 0.0_dp) then
      x = nan_dp()
    else
      x = -log(runif())/r
    end if
  end function rexp

  recursive function rgamma_shape(shape, scale) result(x)
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: x, d, c, z, u, sc
    sc = 1.0_dp
    if (present(scale)) sc = scale
    if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
      x = nan_dp(); return
    end if
    if (shape < 1.0_dp) then
      x = rgamma_shape(shape+1.0_dp)*runif()**(1.0_dp/shape)*sc
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rnorm()
        if (1.0_dp+c*z > 0.0_dp) exit
      end do
      x = (1.0_dp+c*z)**3
      u = runif()
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-x+log(x))) exit
    end do
    x = d*x*sc
  end function rgamma_shape

  function rbeta_ab(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: x, y
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      x = nan_dp(); return
    end if
    x = rgamma_shape(a)
    y = rgamma_shape(b)
    x = x/(x+y)
  end function rbeta_ab

  function rpois(lambda) result(k)
    real(dp), intent(in) :: lambda
    integer :: k
    real(dp) :: l, p, z, x
    if (lambda < 0.0_dp) then
      k = -1; return
    end if
    if (lambda == 0.0_dp) then
      k = 0; return
    end if
    if (lambda < 30.0_dp) then
      l = exp(-lambda); p = 1.0_dp; k = 0
      do
        k = k+1
        p = p*runif()
        if (p <= l) exit
      end do
      k = k-1
    else
      do
        z = rnorm()
        x = lambda + sqrt(lambda)*z
        if (x >= 0.0_dp) exit
      end do
      k = nint(x)
    end if
  end function rpois

  function rbinom(n, p) result(k)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: k, i
    if (n < 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      k = -1; return
    end if
    k = 0
    do i = 1, n
      if (runif() < p) k = k+1
    end do
  end function rbinom

  function rnbinom(size, prob) result(k)
    real(dp), intent(in) :: size, prob
    integer :: k
    real(dp) :: lambda
    if (size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
      k = -1; return
    end if
    lambda = rgamma_shape(size, (1.0_dp-prob)/prob)
    k = rpois(lambda)
  end function rnbinom

  function rinvgauss_rng(mu, phi) result(x)
    real(dp), intent(in) :: mu, phi
    real(dp) :: x, y, z, candidate
    if (mu <= 0.0_dp .or. phi <= 0.0_dp) then
      x = nan_dp(); return
    end if
    z = rnorm(); y = z*z
    candidate = mu + 0.5_dp*mu*mu*y/phi - &
      0.5_dp*mu/phi*sqrt(4.0_dp*mu*phi*y + mu*mu*y*y)
    if (runif() <= mu/(mu+candidate)) then
      x = candidate
    else
      x = mu*mu/candidate
    end if
  end function rinvgauss_rng

end module actuar_rng
