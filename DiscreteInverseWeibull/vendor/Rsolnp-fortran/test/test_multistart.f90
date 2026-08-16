! SPDX-License-Identifier: GPL-2.0-only
program test_multistart
   use rsolnp, only : dp, solnp_problem, multistart_result, solnp_control, &
      solnp_problem_suite, startpars, csolnp_ms, solnp_success
   implicit none

   type(solnp_problem) :: problem
   type(multistart_result) :: result
   type(solnp_control) :: control
   real(dp), allocatable :: starts1(:, :), starts2(:, :)
   integer :: status
   character(len=160) :: message

   call solnp_problem_suite('hs', 5, problem, status, message)
   call check(status == solnp_success, 'retrieve HS05')
   call startpars(problem, 8, starts1, seed=17, status=status, message=message)
   call startpars(problem, 8, starts2, seed=17, status=status, message=message)
   call check(maxval(abs(starts1 - starts2)) <= epsilon(1.0_dp), 'deterministic starts')
   call check(all(starts1 >= spread(problem%lower, 1, 8)), 'starts above lower bounds')
   call check(all(starts1 <= spread(problem%upper, 1, 8)), 'starts below upper bounds')

   control%max_iter = 80
   control%min_iter = 300
   control%tol = 1.0e-8_dp
   call csolnp_ms(problem, 8, result, control, seed=17)
   call check(result%best_index > 0, 'multistart best index')
   call check(result%best%objective < -1.9_dp, 'multistart HS05 objective')
   call check(result%best%convergence == solnp_success, 'multistart convergence')

   print '(a)', 'test_multistart: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check
end program test_multistart
