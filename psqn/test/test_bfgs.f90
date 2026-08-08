! SPDX-License-Identifier: Apache-2.0
program test_bfgs
  use psqn, only : dp, psqn_info, psqn_bfgs_options, psqn_bfgs_optimize, psqn_converged
  implicit none
  real(dp) :: x(2)
  type(psqn_info) :: res
  type(psqn_bfgs_options) :: opt

  x = [-1.2_dp, 1.0_dp]
  opt%max_it = 300
  opt%gr_tol = 1.0e-8_dp
  call psqn_bfgs_optimize(x, rosenbrock, res, opt)
  if (res%info /= psqn_converged) error stop "BFGS did not converge"
  if (maxval(abs(x - [1.0_dp, 1.0_dp])) > 1.0e-5_dp) error stop "BFGS wrong solution"
  if (abs(res%value) > 1.0e-10_dp) error stop "BFGS wrong value"
  print '(a)', 'test_bfgs: PASS'
contains
  subroutine rosenbrock(x, f, g, comp_grad)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    logical, intent(in) :: comp_grad
    f = 100.0_dp * (x(2) - x(1)*x(1))**2 + (1.0_dp - x(1))**2
    if (comp_grad) then
      g(1) = -400.0_dp*x(1)*(x(2)-x(1)*x(1)) - 2.0_dp*(1.0_dp-x(1))
      g(2) = 200.0_dp*(x(2)-x(1)*x(1))
    end if
  end subroutine rosenbrock
end program test_bfgs
