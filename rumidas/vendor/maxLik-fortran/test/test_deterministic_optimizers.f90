! SPDX-License-Identifier: GPL-2.0-or-later
program test_deterministic_optimizers
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik
  implicit none

  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp) :: start(2)
  character(len=16), parameter :: methods(5) = [character(len=16) :: 'nr', 'bfgs', 'bfgsr', 'cg', 'nm']
  integer :: i

  start = [-3.0_dp, 4.0_dp]
  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  problem%hessian => hessian
  control%iterlim = 400
  control%gradtol = 1.0e-7_dp
  control%tol = 1.0e-9_dp

  do i = 1, size(methods)
    call max_lik(problem, start, result, trim(methods(i)), control)
    call assert_true(result%converged, trim(methods(i)) // ' converged')
    call assert_close(result%estimate(1), 1.0_dp, 2.0e-4_dp, trim(methods(i)) // ' x1')
    call assert_close(result%estimate(2), -2.0_dp, 2.0e-4_dp, trim(methods(i)) // ' x2')
    call assert_close(result%maximum, 0.0_dp, 1.0e-7_dp, trim(methods(i)) // ' maximum')
  end do
  print '(a)', 'test_deterministic_optimizers: PASS'

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

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      if (allocated(result%estimate)) print '(a,2f14.6)', 'estimate:', result%estimate
      print '(a,i0,2a)', 'code=', result%code, ' message=', trim(result%message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_deterministic_optimizers
