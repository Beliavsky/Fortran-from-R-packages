! SPDX-License-Identifier: Apache-2.0
program test_structured
  use psqn, only : dp, psqn_options, psqn_info, psqn_optimize_structured, psqn_converged, psqn_pre_block
  implicit none
  integer :: private_dims(2)
  type(psqn_options) :: opt
  type(psqn_info) :: res
  real(dp) :: x(3)

  private_dims = 1
  x = 0.0_dp
  opt%max_it = 200
  opt%pre_method = psqn_pre_block
  call psqn_optimize_structured(x, 1, private_dims, elem, res, opt)
  if (res%info /= psqn_converged) error stop "structured PSQN did not converge"
  if (maxval(abs(x - [1.0_dp,3.0_dp,0.0_dp])) > 1.0e-6_dp) error stop "structured PSQN wrong solution"
  print '(a)', 'test_structured: PASS'
contains
  subroutine elem(i, x, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    logical, intent(in) :: comp_grad
    real(dp) :: a
    if (i == 1) then
      a = 2.0_dp
    else
      a = -1.0_dp
    end if
    f = 0.5_dp*(x(2)-x(1)-a)**2 + 0.25_dp*(x(1)-1.0_dp)**2
    if (comp_grad) then
      g(1) = -(x(2)-x(1)-a) + 0.5_dp*(x(1)-1.0_dp)
      g(2) = x(2)-x(1)-a
    end if
  end subroutine elem
end program test_structured
