! SPDX-License-Identifier: GPL-2.0-only
module calibrar_test_functions
  use calibrar_kinds, only : dp
  use calibrar_random, only : rand_normal
  implicit none
  private
  public :: sphere_n, quadratic_objective, shifted_quadratic, quadratic_gradient
  public :: two_component_objective, rosenbrock_objective
contains
  function sphere_n(x, noise_sd, aggregate) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: noise_sd
    logical, intent(in), optional :: aggregate
    real(dp) :: v
    real(dp) :: sd
    integer :: i
    sd=0.1_dp
    if(present(noise_sd)) sd=noise_sd
    v=0.0_dp
    do i=1,size(x)
      v=v+(x(i)+sd*rand_normal())**2
    end do
    if(present(aggregate)) then
      if(.not.aggregate) error stop "sphere_n: use component-wise code for aggregate=.false."
    end if
  end function sphere_n

  function quadratic_objective(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v=sum(x*x)+10.0_dp
  end function quadratic_objective

  function shifted_quadratic(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v=sum((x-1.5_dp)**2)
  end function shifted_quadratic

  subroutine quadratic_gradient(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g=2.0_dp*x
  end subroutine quadratic_gradient

  subroutine two_component_objective(x,f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    if(size(f)<2) error stop "two_component_objective: need two outputs"
    f(1)=sum(x*x)
    f(2)=sum((x-0.5_dp)**2)
  end subroutine two_component_objective

  function rosenbrock_objective(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    integer :: i
    v=0.0_dp
    do i=1,size(x)-1
      v=v+100.0_dp*(x(i+1)-x(i)*x(i))**2+(1.0_dp-x(i))**2
    end do
  end function rosenbrock_objective
end module calibrar_test_functions
