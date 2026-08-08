! SPDX-License-Identifier: Apache-2.0
program test_generic
  use psqn, only : dp, psqn_element_spec, psqn_options, psqn_info, psqn_optimize_generic, psqn_converged
  implicit none
  type(psqn_element_spec) :: specs(2)
  type(psqn_options) :: opt
  type(psqn_info) :: res
  real(dp) :: x(3), h(3,3)

  specs(1)%idx = [1,2]
  specs(2)%idx = [2,3]
  x = 0.0_dp
  opt%max_it = 200
  opt%gr_tol = 1.0e-9_dp
  call psqn_optimize_generic(x, specs, elem, res, opt, hess=h)
  if (res%info /= psqn_converged) error stop "generic PSQN did not converge"
  if (maxval(abs(x - [1.0_dp,-2.0_dp,3.0_dp])) > 1.0e-6_dp) error stop "generic PSQN wrong solution"
  print '(a)', 'test_generic: PASS'
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
end program test_generic
