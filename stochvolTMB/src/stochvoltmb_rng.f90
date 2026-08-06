! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_rng
  use stochvoltmb_kinds, only : dp, pi, tiny_dp
  implicit none
  private

  type, public :: sv_rng_state
    integer(kind=8) :: state = 88172645463393265_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
    procedure :: gamma => rng_gamma
    procedure :: student_t => rng_student_t
  end type sv_rng_state

contains

  subroutine rng_seed(self, seed)
    class(sv_rng_state), intent(inout) :: self
    integer(kind=8), intent(in) :: seed
    self%state = seed
    if (self%state == 0_8) self%state = 88172645463393265_8
    self%has_spare = .false.
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(u)
    class(sv_rng_state), intent(inout) :: self
    integer(kind=8) :: x
    x = self%state
    x = ieor(x, shiftl(x,13))
    x = ieor(x, shiftr(x,7))
    x = ieor(x, shiftl(x,17))
    self%state = x
    u = real(iand(x, int(z'7FFFFFFFFFFFFFFF',kind=8)),dp) / real(huge(0_8),dp)
    u = min(1.0_dp-tiny_dp,max(tiny_dp,u))
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(sv_rng_state), intent(inout) :: self
    real(dp) :: u1, u2, r
    if (self%has_spare) then
      z = self%spare
      self%has_spare = .false.
      return
    end if
    u1 = self%uniform()
    u2 = self%uniform()
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*pi*u2)
    self%spare = r*sin(2.0_dp*pi*u2)
    self%has_spare = .true.
  end function rng_normal

  recursive real(dp) function rng_gamma(self, shape) result(x)
    class(sv_rng_state), intent(inout) :: self
    real(dp), intent(in) :: shape
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = self%gamma(shape+1.0_dp)*self%uniform()**(1.0_dp/shape)
      return
    end if
    d = shape-1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = self%normal()
        v = 1.0_dp+c*z
        if (v > 0.0_dp) exit
      end do
      v = v*v*v
      u = self%uniform()
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v
  end function rng_gamma

  real(dp) function rng_student_t(self, df) result(x)
    class(sv_rng_state), intent(inout) :: self
    real(dp), intent(in) :: df
    x = self%normal()/sqrt(2.0_dp*self%gamma(0.5_dp*df)/df)
  end function rng_student_t

end module stochvoltmb_rng
