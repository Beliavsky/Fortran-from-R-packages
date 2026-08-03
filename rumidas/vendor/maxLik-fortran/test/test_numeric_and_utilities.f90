! SPDX-License-Identifier: GPL-2.0-or-later
program test_numeric_and_utilities
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, max_lik, &
    numeric_jacobian, pack_active_parameters, unpack_active_parameters, progressive_condition_numbers, set_bounds
  implicit none

  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp), allocatable :: jacobian(:, :), subset(:), full(:), condition(:)
  real(dp) :: matrix(2, 2)
  integer :: status

  call initialize_problem(problem, 2, objective)
  control%iterlim = 300
  control%gradtol = 1.0e-6_dp
  call max_lik(problem, [-2.0_dp, 3.0_dp], result, 'nr', control)
  call assert_true(result%converged, 'numeric Newton converged')
  call assert_close(result%estimate(1), 1.0_dp, 2.0e-4_dp, 'numeric Newton x1')
  call assert_close(result%estimate(2), -2.0_dp, 2.0e-4_dp, 'numeric Newton x2')

  call set_bounds(problem, [-10.0_dp, -10.0_dp], [0.0_dp, 10.0_dp], status)
  call assert_true(status == 0, 'set bounds')
  call max_lik(problem, [-2.0_dp, 3.0_dp], result, 'bfgs', control)
  call assert_true(result%converged, 'bounded BFGS converged')
  call assert_close(result%estimate(1), 0.0_dp, 2.0e-6_dp, 'bounded x1')
  call assert_close(result%estimate(2), -2.0_dp, 2.0e-5_dp, 'bounded x2')

  call numeric_jacobian(vector_function, [2.0_dp, 3.0_dp], 2, jacobian, status)
  call assert_true(status == 0, 'numeric Jacobian status')
  call assert_close(jacobian(1, 1), 4.0_dp, 1.0e-5_dp, 'Jacobian 11')
  call assert_close(jacobian(1, 2), 1.0_dp, 1.0e-5_dp, 'Jacobian 12')
  call assert_close(jacobian(2, 1), 3.0_dp, 1.0e-5_dp, 'Jacobian 21')
  call assert_close(jacobian(2, 2), 2.0_dp, 1.0e-5_dp, 'Jacobian 22')

  call pack_active_parameters([1.0_dp, 2.0_dp, 3.0_dp], [.true., .false., .true.], subset, status)
  call assert_true(status == 0 .and. all(abs(subset - [1.0_dp, 3.0_dp]) < 1.0e-14_dp), 'pack active')
  call unpack_active_parameters([4.0_dp, 5.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], &
    [.true., .false., .true.], full, status)
  call assert_true(status == 0 .and. all(abs(full - [4.0_dp, 2.0_dp, 5.0_dp]) < 1.0e-14_dp), 'unpack active')

  matrix = reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])
  call progressive_condition_numbers(matrix, condition, status)
  call assert_true(status == 0, 'condition number status')
  call assert_close(condition(1), 1.0_dp, 1.0e-12_dp, 'condition one column')
  call assert_close(condition(2), 2.0_dp, 1.0e-12_dp, 'condition two columns')
  print '(a)', 'test_numeric_and_utilities: PASS'

contains

  subroutine objective(x, value, callback_status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: callback_status
    value = -0.5_dp * ((x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2)
    callback_status = 0
  end subroutine objective

  subroutine vector_function(x, values, callback_status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
    integer, intent(out) :: callback_status
    values = [x(1)**2 + x(2), x(1) * x(2)]
    callback_status = 0
  end subroutine vector_function

  subroutine assert_true(condition_value, label)
    logical, intent(in) :: condition_value
    character(len=*), intent(in) :: label
    if (.not. condition_value) then
      print '(a)', 'FAIL: ' // label
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_numeric_and_utilities
