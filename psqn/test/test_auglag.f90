! SPDX-License-Identifier: Apache-2.0
program test_auglag
  use psqn, only : dp, psqn_element_spec, psqn_options, psqn_auglag_options, psqn_auglag_info, &
                   psqn_aug_lagrang_generic, psqn_converged
  implicit none
  type(psqn_element_spec) :: specs(1), cspecs(1)
  type(psqn_options) :: opt
  type(psqn_auglag_options) :: aopt
  type(psqn_auglag_info) :: res
  real(dp) :: x(2), multipliers(1)

  specs(1)%idx = [1,2]
  cspecs(1)%idx = [1,2]
  x = [1.0_dp,2.0_dp]
  multipliers = 0.0_dp
  opt%max_it = 100
  opt%gr_tol = 1.0e-10_dp
  aopt%max_it_outer = 80
  aopt%violations_norm_thresh = 1.0e-7_dp
  aopt%tau = 2.0_dp
  call psqn_aug_lagrang_generic(x, specs, objective, cspecs, constraint, multipliers, res, opt, aopt)
  if (res%info /= psqn_converged) error stop "auglag did not converge"
  if (maxval(abs(x - [0.0_dp,1.0_dp])) > 2.0e-5_dp) error stop "auglag wrong solution"
  print '(a)', 'test_auglag: PASS'
contains
  subroutine objective(i, x, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    logical, intent(in) :: comp_grad
    f = 0.5_dp*(x(1)-1.0_dp)**2 + 0.5_dp*(x(2)-2.0_dp)**2
    if (comp_grad) g = [x(1)-1.0_dp, x(2)-2.0_dp]
  end subroutine objective
  subroutine constraint(i, x, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    logical, intent(in) :: comp_grad
    f = x(1) + x(2) - 1.0_dp
    if (comp_grad) g = 1.0_dp
  end subroutine constraint
end program test_auglag
