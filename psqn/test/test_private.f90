! SPDX-License-Identifier: Apache-2.0
program test_private
  use psqn, only : dp, psqn_optimize_private_structured
  implicit none
  integer :: private_dims(2)
  real(dp) :: x(3), value

  private_dims = 1
  x = [1.0_dp, 0.0_dp, 0.0_dp]
  call psqn_optimize_private_structured(x, 1, private_dims, elem, value, gr_tol=1.0e-10_dp)
  if (maxval(abs(x - [1.0_dp,3.0_dp,0.0_dp])) > 1.0e-6_dp) error stop "private optimize wrong solution"
  if (abs(value) > 1.0e-10_dp) error stop "private optimize wrong value"
  print '(a)', 'test_private: PASS'
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
    f = 0.5_dp*(x(2)-x(1)-a)**2
    if (comp_grad) then
      g(1) = -(x(2)-x(1)-a)
      g(2) = x(2)-x(1)-a
    end if
  end subroutine elem
end program test_private
