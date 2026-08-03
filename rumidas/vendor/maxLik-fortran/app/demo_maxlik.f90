! SPDX-License-Identifier: GPL-2.0-or-later
program demo_maxlik
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik
  implicit none

  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  character(len=16), parameter :: methods(5) = [character(len=16) :: 'nr', 'bfgs', 'bfgsr', 'cg', 'nm']
  real(dp) :: start(2)
  integer :: i

  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  problem%hessian => hessian
  control%iterlim = 400
  start = [-3.0_dp, 4.0_dp]

  print '(a)', 'maxLik modern Fortran demonstration'
  print '(a)', 'method       maximum          x(1)          x(2)   code'
  do i = 1, size(methods)
    call max_lik(problem, start, result, trim(methods(i)), control)
    print '(a10,3f14.6,i7)', trim(methods(i)), result%maximum, result%estimate, result%code
  end do

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -0.5_dp * ((x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = [-(x(1) - 1.0_dp), -2.0_dp * (x(2) + 2.0_dp)]
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    h = 0.0_dp
    h(1, 1) = -1.0_dp + 0.0_dp * x(1)
    h(2, 2) = -2.0_dp
    status = 0
  end subroutine hessian

end program demo_maxlik
