! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_rng
  use iso_fortran_env, only : int64
  use sn_kinds, only : dp, pi, tiny_dp
  use sn_math, only : normal_cdf, normal_quantile
  implicit none
  private

  type, public :: sn_rng_state
    integer(int64) :: s(4) = [ &
      int(z'123456789abcdef0',int64), int(z'0fedcba987654321',int64), &
      int(z'9e3779b97f4a7c15',int64), int(z'6a09e667f3bcc909',int64)]
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
    procedure :: gamma => rng_gamma
    procedure :: chi_square => rng_chi_square
    procedure :: truncated_normal_lower => rng_truncated_normal_lower
  end type sn_rng_state

  public :: splitmix64

contains

  pure integer(int64) function rotl(x, k) result(y)
    integer(int64), intent(in) :: x
    integer, intent(in) :: k
    y = ior(shiftl(x,k),shiftr(x,64-k))
  end function rotl

  integer(int64) function splitmix64(x) result(z)
    integer(int64), intent(inout) :: x
    x = x + int(z'9e3779b97f4a7c15',int64)
    z = x
    z = ieor(z,shiftr(z,30))*int(z'bf58476d1ce4e5b9',int64)
    z = ieor(z,shiftr(z,27))*int(z'94d049bb133111eb',int64)
    z = ieor(z,shiftr(z,31))
  end function splitmix64

  subroutine rng_seed(self, seed_value)
    class(sn_rng_state), intent(inout) :: self
    integer(int64), intent(in) :: seed_value
    integer(int64) :: x
    integer :: i
    x = seed_value
    do i=1,4
      self%s(i) = splitmix64(x)
    end do
    if (all(self%s == 0_int64)) self%s(1) = 1_int64
    self%has_spare = .false.
    self%spare = 0.0_dp
  end subroutine rng_seed

  integer(int64) function next_u64(self) result(value)
    class(sn_rng_state), intent(inout) :: self
    integer(int64) :: t
    value = rotl(self%s(2)*5_int64,7)*9_int64
    t = shiftl(self%s(2),17)
    self%s(3) = ieor(self%s(3),self%s(1))
    self%s(4) = ieor(self%s(4),self%s(2))
    self%s(2) = ieor(self%s(2),self%s(3))
    self%s(1) = ieor(self%s(1),self%s(4))
    self%s(3) = ieor(self%s(3),t)
    self%s(4) = rotl(self%s(4),45)
  end function next_u64

  real(dp) function rng_uniform(self) result(u)
    class(sn_rng_state), intent(inout) :: self
    integer(int64) :: z
    z = next_u64(self)
    u = real(shiftr(z,11),dp)*1.11022302462515654042e-16_dp
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp-epsilon(1.0_dp)
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(sn_rng_state), intent(inout) :: self
    real(dp) :: u1, u2, r
    if (self%has_spare) then
      z = self%spare
      self%has_spare = .false.
      return
    end if
    u1 = max(self%uniform(),tiny_dp)
    u2 = self%uniform()
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*pi*u2)
    self%spare = r*sin(2.0_dp*pi*u2)
    self%has_spare = .true.
  end function rng_normal

  recursive real(dp) function rng_gamma(self, shape, scale) result(x)
    class(sn_rng_state), intent(inout) :: self
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: d, c, z, v, u, scl

    scl = 1.0_dp
    if (present(scale)) scl = scale
    if (shape <= 0.0_dp .or. scl <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = self%gamma(shape+1.0_dp,1.0_dp)*self%uniform()**(1.0_dp/shape)*scl
      return
    end if

    d = shape-1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = self%normal()
      v = 1.0_dp+c*z
      if (v <= 0.0_dp) cycle
      v = v*v*v
      u = self%uniform()
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v*scl
  end function rng_gamma

  real(dp) function rng_chi_square(self, nu) result(x)
    class(sn_rng_state), intent(inout) :: self
    real(dp), intent(in) :: nu
    x = self%gamma(0.5_dp*nu,2.0_dp)
  end function rng_chi_square

  real(dp) function rng_truncated_normal_lower(self, lower) result(x)
    class(sn_rng_state), intent(inout) :: self
    real(dp), intent(in) :: lower
    real(dp) :: p0, u
    p0 = normal_cdf(lower)
    u = p0+(1.0_dp-p0)*self%uniform()
    x = normal_quantile(u)
  end function rng_truncated_normal_lower

end module sn_rng
