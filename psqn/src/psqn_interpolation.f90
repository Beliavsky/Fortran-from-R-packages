! SPDX-License-Identifier: Apache-2.0
module psqn_interpolation
  use psqn_types, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: intrapolate_type

  type :: intrapolate_type
    real(dp) :: f0 = 0.0_dp
    real(dp) :: d0 = 0.0_dp
    real(dp) :: xold = 0.0_dp
    real(dp) :: fold = 0.0_dp
    real(dp) :: xnew = 0.0_dp
    real(dp) :: fnew = 0.0_dp
    logical :: has_two_values = .false.
  contains
    procedure :: init => intrapolate_init
    procedure :: get_value => intrapolate_get_value
    procedure :: update => intrapolate_update
  end type intrapolate_type

contains

  subroutine intrapolate_init(self, f0, d0, x, f)
    class(intrapolate_type), intent(inout) :: self
    real(dp), intent(in) :: f0, d0, x, f
    self%f0 = f0
    self%d0 = d0
    self%xnew = x
    self%fnew = f
    self%has_two_values = .false.
  end subroutine intrapolate_init

  real(dp) function intrapolate_get_value(self, v1, v2) result(val)
    class(intrapolate_type), intent(in) :: self
    real(dp), intent(in) :: v1, v2
    real(dp) :: a, b, small, fac, f1, f2, vb, va, deter

    a = min(v1, v2)
    b = max(v1, v2)
    small = 0.01_dp * (b - a)

    if (.not. self%has_two_values) then
      if (self%fnew - self%f0 - self%d0 * self%xnew == 0.0_dp) then
        val = a + 0.5_dp * (b - a)
      else
        val = -self%d0 * self%xnew * self%xnew / &
          (2.0_dp * (self%fnew - self%f0 - self%d0 * self%xnew))
      end if
    else
      fac = self%xnew**2 * self%xold**2 * (self%xnew - self%xold)
      if (fac == 0.0_dp) fac = 1.0_dp
      f1 = self%fnew - self%f0 - self%d0 * self%xnew
      f2 = self%fold - self%f0 - self%d0 * self%xold
      vb = (-self%xold**3 * f1 + self%xnew**3 * f2) / fac
      va = (self%xold**2 * f1 - self%xnew**2 * f2) / fac
      deter = vb * vb - 3.0_dp * va * self%d0
      if (deter < 0.0_dp .or. va == 0.0_dp) then
        val = a + 0.5_dp * (b - a)
      else
        val = (-vb + sqrt(deter)) / (3.0_dp * va)
      end if
    end if

    if (val < a + small .or. val > b - small .or. .not. ieee_is_finite(val)) then
      val = a + 0.5_dp * (b - a)
    end if
  end function intrapolate_get_value

  subroutine intrapolate_update(self, x, f)
    class(intrapolate_type), intent(inout) :: self
    real(dp), intent(in) :: x, f
    self%xold = self%xnew
    self%fold = self%fnew
    self%xnew = x
    self%fnew = f
    self%has_two_values = .true.
  end subroutine intrapolate_update

end module psqn_interpolation
