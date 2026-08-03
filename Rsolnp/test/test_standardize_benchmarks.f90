! SPDX-License-Identifier: GPL-2.0-only
program test_standardize_benchmarks
   use rsolnp, only : dp, solnp_problem, problem_table_entry, solnp_result, &
      solnp_problem_suite, solnp_problems_table, solnp_standardize_problem, &
      csolnp, kkt_diagnose, kkt_diagnostics, solnp_control, solnp_success
   implicit none

   type(solnp_problem) :: problem, standardized
   type(problem_table_entry), allocatable :: table(:)
   type(solnp_result) :: result
   type(kkt_diagnostics) :: diagnostics
   type(solnp_control) :: control
   integer :: status
   character(len=160) :: message

   call solnp_problems_table(table)
   call check(size(table) == 77, 'benchmark table size')
   call check(count(table%implemented) == 18, 'translated benchmark count')

   call solnp_problem_suite('Other', 9, problem, status, message)
   call check(status == solnp_success, 'retrieve Wright9')
   call solnp_standardize_problem(problem, standardized, status, message)
   call check(status == solnp_success, 'standardize Wright9')
   call check(standardized%standard_form, 'standard-form flag')
   call check(standardized%n_ineq == 6, 'two-sided inequalities expanded')
   call check(maxval(abs(standardized%ineq_upper)) <= epsilon(1.0_dp), 'standard upper bounds')

   call solnp_problem_suite('Other', 6, problem, status, message)
   control%max_iter = 250
   control%min_iter = 500
   control%tol = 1.0e-8_dp
   call csolnp(problem, result, control)
   call kkt_diagnose(problem, result, diagnostics, 1.0e-6_dp)
   call check(abs(result%objective - problem%best_fn) < 2.0e-5_dp, 'Powell objective')
   call check(diagnostics%eq_violation < 1.0e-6_dp, 'Powell feasibility')

   print '(a)', 'test_standardize_benchmarks: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check
end program test_standardize_benchmarks
