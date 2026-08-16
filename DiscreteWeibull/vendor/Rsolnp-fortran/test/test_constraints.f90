! SPDX-License-Identifier: GPL-2.0-only
program test_constraints
   use rsolnp, only : dp, solnp_problem, solnp_result, solnp_control, csolnp, &
      solnp_problem_suite, solnp_success
   implicit none

   type(solnp_problem) :: problem
   type(solnp_result) :: result
   type(solnp_control) :: control
   integer :: status
   character(len=160) :: message

   control%max_iter = 200
   control%min_iter = 500
   control%tol = 1.0e-8_dp
   control%restoration_iter = 150

   call solnp_problem_suite('hs', 6, problem, status, message)
   call check(status == solnp_success, 'retrieve HS06')
   call csolnp(problem, result, control)
   call check(result%convergence == solnp_success, 'HS06 convergence')
   call check(result%objective < 1.0e-10_dp, 'HS06 objective')
   call check(result%kkt%eq_violation < 1.0e-7_dp, 'HS06 equality')

   call solnp_problem_suite('Other', 2, problem, status, message)
   call check(status == solnp_success, 'retrieve box')
   call csolnp(problem, result, control)
   call check(result%convergence == solnp_success, 'box convergence')
   call check(abs(result%objective + 48.11252243_dp) < 1.0e-4_dp, 'box objective')
   call check(result%kkt%eq_violation < 2.0e-7_dp, 'box equality')

   call solnp_problem_suite('hs', 11, problem, status, message)
   call check(status == solnp_success, 'retrieve HS11')
   call csolnp(problem, result, control)
   call check(result%kkt%ineq_violation < 1.0e-5_dp, 'HS11 inequality')
   call check(result%objective < -8.0_dp, 'HS11 objective improvement')

   print '(a)', 'test_constraints: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check
end program test_constraints
