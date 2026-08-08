! SPDX-License-Identifier: Apache-2.0
program test_hessian
  use psqn, only : dp, psqn_element_spec, psqn_generic_hess
  implicit none
  type(psqn_element_spec) :: specs(2)
  real(dp) :: x(3), h(3,3), expected(3,3)

  specs(1)%idx = [1,2]
  specs(2)%idx = [2,3]
  x = [0.3_dp,-0.7_dp,2.1_dp]
  call psqn_generic_hess(x, specs, elem, h)
  expected = 0.0_dp
  expected(1,1) = 1.0_dp
  expected(2,2) = 2.0_dp
  expected(3,3) = 1.0_dp
  if (maxval(abs(h-expected)) > 1.0e-7_dp) error stop "numerical Hessian mismatch"
  print '(a)', 'test_hessian: PASS'
contains
  subroutine elem(i, x, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    logical, intent(in) :: comp_grad
    select case(i)
    case(1)
      f = 0.5_dp*(x(1)-1.0_dp)**2 + 0.5_dp*(x(2)+2.0_dp)**2
      if (comp_grad) g = [x(1)-1.0_dp, x(2)+2.0_dp]
    case(2)
      f = 0.5_dp*(x(1)+2.0_dp)**2 + 0.5_dp*(x(2)-3.0_dp)**2
      if (comp_grad) g = [x(1)+2.0_dp, x(2)-3.0_dp]
    end select
  end subroutine elem
end program test_hessian
