program test_errors
   use sparse_index_tracking, only : dp, sparse_index_fit, fit_sparse_index_tracking, &
                                     sit_dimension_error, sit_infeasible_bounds, &
                                     sit_invalid_argument, sit_iteration_limit
   implicit none

   real(dp) :: x(8, 3), target(8), bad_initial(3)
   type(sparse_index_fit) :: fit
   integer :: i

   do i = 1, 8
      x(i, :) = [real(i, dp), real(i * i, dp), sin(real(i, dp))] * 1.0e-3_dp
   end do
   target = 0.5_dp * x(:, 1) + 0.5_dp * x(:, 2)

   call fit_sparse_index_tracking(x(:, 1:1), target, 1.0e-7_dp, fit)
   call assert_true(fit%info == sit_dimension_error, 'univariate x was accepted')

   call fit_sparse_index_tracking(x, target, 0.0_dp, fit)
   call assert_true(fit%info == sit_invalid_argument, 'zero lambda was accepted')

   call fit_sparse_index_tracking(x, target, 1.0e-7_dp, fit, upper_bound=0.2_dp)
   call assert_true(fit%info == sit_infeasible_bounds, 'infeasible cap was accepted')

   call fit_sparse_index_tracking(x, target, 1.0e-7_dp, fit, measure='hete')
   call assert_true(fit%info == sit_invalid_argument, 'missing Huber parameter was accepted')

   call fit_sparse_index_tracking(x, target, 1.0e-7_dp, fit, measure='unknown')
   call assert_true(fit%info == sit_invalid_argument, 'unknown measure was accepted')

   bad_initial = [0.8_dp, 0.3_dp, -0.1_dp]
   call fit_sparse_index_tracking(x, target, 1.0e-7_dp, fit, initial_weights=bad_initial)
   call assert_true(fit%info == sit_invalid_argument, 'invalid initial weights were accepted')

   call fit_sparse_index_tracking(x, target, 1.0e-7_dp, fit, max_iterations=1)
   call assert_true(fit%info == sit_iteration_limit, 'iteration limit was not reported')
   call assert_true(abs(sum(fit%weights) - 1.0_dp) < 1.0e-11_dp, 'limited fit did not return feasible weights')

   print '(a)', 'test_errors: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_errors
