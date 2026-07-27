! SPDX-License-Identifier: GPL-3.0-only
module test_custom_callbacks
  use svdnf_kinds, only : dp
  use svdnf_types, only : svm_dynamics
  implicit none
  private
  public :: custom_mu_y, custom_sigma_y, custom_mu_x, custom_sigma_x, custom_setter
contains
  function custom_mu_y(x,p) result(value)
    real(dp), intent(in) :: x,p(:)
    real(dp) :: value
    value=p(1)+0.1_dp*x
  end function custom_mu_y
  function custom_sigma_y(x,p) result(value)
    real(dp), intent(in) :: x,p(:)
    real(dp) :: value
    value=p(1)+0.0_dp*x
  end function custom_sigma_y
  function custom_mu_x(x,p) result(value)
    real(dp), intent(in) :: x,p(:)
    real(dp) :: value
    value=p(1)*x
  end function custom_mu_x
  function custom_sigma_x(x,p) result(value)
    real(dp), intent(in) :: x,p(:)
    real(dp) :: value
    value=p(1)+0.0_dp*x
  end function custom_sigma_x
  subroutine custom_setter(dynamics,parameters,ok)
    type(svm_dynamics), intent(inout) :: dynamics
    real(dp), intent(in) :: parameters(:)
    logical, intent(out) :: ok
    ok=size(parameters)==1
    if (ok) dynamics%mu_y_parameters=[parameters(1)]
  end subroutine custom_setter
end module test_custom_callbacks
