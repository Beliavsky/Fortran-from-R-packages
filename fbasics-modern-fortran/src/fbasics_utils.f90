! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_utils
  use fbasics_kinds, only: dp
  implicit none
  private
  public :: heaviside, sign_function, delta_function, boxcar, ramp
contains
  pure elemental real(dp) function heaviside(x) result(v)
    real(dp),intent(in)::x
    if(x<0.0_dp)then;v=0.0_dp;else if(x>0.0_dp)then;v=1.0_dp;else;v=0.5_dp;end if
  end function
  pure elemental real(dp) function sign_function(x) result(v)
    real(dp),intent(in)::x
    if(x<0.0_dp)then;v=-1.0_dp;else if(x>0.0_dp)then;v=1.0_dp;else;v=0.0_dp;end if
  end function
  pure elemental real(dp) function delta_function(x,tolerance) result(v)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::tolerance
    real(dp)::tol
    tol=sqrt(epsilon(1.0_dp));if(present(tolerance))tol=tolerance
    if(abs(x)<=tol)then;v=1.0_dp/tol;else;v=0.0_dp;end if
  end function
  pure elemental real(dp) function boxcar(x,left,right) result(v)
    real(dp),intent(in)::x,left,right
    v=heaviside(x-left)-heaviside(x-right)
  end function
  pure elemental real(dp) function ramp(x) result(v)
    real(dp),intent(in)::x
    v=max(x,0.0_dp)
  end function
end module fbasics_utils
