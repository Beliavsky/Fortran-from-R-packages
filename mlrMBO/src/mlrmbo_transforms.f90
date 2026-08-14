module mlrmbo_transforms
  use mlrmbo_kinds, only : dp
  implicit none
  private
  public :: trafo_log, trafo_sqrt, inverse_trafo_log, inverse_trafo_sqrt
contains
  elemental real(dp) function trafo_log(x) result(y)
    real(dp), intent(in) :: x
    if(x<=0.0_dp) error stop 'trafo_log: x must be positive'
    y=log(x)
  end function trafo_log
  elemental real(dp) function trafo_sqrt(x) result(y)
    real(dp), intent(in) :: x
    if(x<0.0_dp) error stop 'trafo_sqrt: x must be nonnegative'
    y=sqrt(x)
  end function trafo_sqrt
  elemental real(dp) function inverse_trafo_log(y) result(x)
    real(dp), intent(in) :: y
    x=exp(y)
  end function inverse_trafo_log
  elemental real(dp) function inverse_trafo_sqrt(y) result(x)
    real(dp), intent(in) :: y
    x=y*y
  end function inverse_trafo_sqrt
end module mlrmbo_transforms
