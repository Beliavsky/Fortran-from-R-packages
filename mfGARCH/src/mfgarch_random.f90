! SPDX-License-Identifier: MIT
module mfgarch_random
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch_kinds, only : dp
  use mfgarch_math, only : pi_dp
  implicit none
  private

  type, public :: mfgarch_rng
    private
    integer(int64) :: state = 88172645463393265_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  contains
    procedure, public :: seed => rng_seed
    procedure, public :: uniform => rng_uniform
    procedure, public :: normal => rng_normal
    procedure, public :: gamma => rng_gamma
    procedure, public :: student_t => rng_student_t
  end type mfgarch_rng

contains

  subroutine rng_seed(self, seed)
    class(mfgarch_rng), intent(inout) :: self
    integer(int64), intent(in) :: seed

    self%state = seed
    if (self%state == 0_int64) self%state = 88172645463393265_int64
    self%has_spare = .false.
  end subroutine rng_seed

  function rng_uniform(self) result(value)
    class(mfgarch_rng), intent(inout) :: self
    real(dp) :: value
    integer(int64) :: x, positive

    x = self%state
    x = ieor(x, shiftl(x,13))
    x = ieor(x, shiftr(x,7))
    x = ieor(x, shiftl(x,17))
    self%state = x
    positive = iand(x, huge(1_int64))
    value = (real(positive, dp) + 0.5_dp) / (real(huge(1_int64), dp) + 1.0_dp)
    value = min(1.0_dp - epsilon(1.0_dp), max(tiny(1.0_dp), value))
  end function rng_uniform

  function rng_normal(self) result(value)
    class(mfgarch_rng), intent(inout) :: self
    real(dp) :: value
    real(dp) :: radius, angle

    if (self%has_spare) then
      value = self%spare
      self%has_spare = .false.
      return
    end if
    radius = sqrt(-2.0_dp * log(self%uniform()))
    angle = 2.0_dp * pi_dp * self%uniform()
    value = radius * cos(angle)
    self%spare = radius * sin(angle)
    self%has_spare = .true.
  end function rng_normal

  recursive function rng_gamma(self, shape) result(value)
    class(mfgarch_rng), intent(inout) :: self
    real(dp), intent(in) :: shape
    real(dp) :: value
    real(dp) :: d, c, x, v, u

    if (shape <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      value = self%gamma(shape + 1.0_dp) * self%uniform()**(1.0_dp / shape)
      return
    end if
    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      do
        x = self%normal()
        v = 1.0_dp + c * x
        if (v > 0.0_dp) exit
      end do
      v = v**3
      u = self%uniform()
      if (u < 1.0_dp - 0.0331_dp * x**4) exit
      if (log(u) < 0.5_dp * x**2 + d * (1.0_dp - v + log(v))) exit
    end do
    value = d * v
  end function rng_gamma

  function rng_student_t(self, degrees_of_freedom) result(value)
    class(mfgarch_rng), intent(inout) :: self
    real(dp), intent(in) :: degrees_of_freedom
    real(dp) :: value, chi_square

    if (degrees_of_freedom <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    chi_square = 2.0_dp * self%gamma(0.5_dp * degrees_of_freedom)
    value = self%normal() / sqrt(chi_square / degrees_of_freedom)
  end function rng_student_t

end module mfgarch_random
