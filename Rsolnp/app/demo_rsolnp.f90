! SPDX-License-Identifier: GPL-2.0-only
program demo_rsolnp
   use rsolnp, only : dp, solnp_problem, solnp_result, solnp_control, &
      solnp_problem_suite, csolnp
   implicit none

   type(solnp_problem) :: problem
   type(solnp_result) :: result
   type(solnp_control) :: control
   integer :: status
   character(len=160) :: message

   call solnp_problem_suite('Other', 6, problem, status, message)
   if (status /= 0) error stop trim(message)
   control%max_iter = 250
   control%min_iter = 500
   control%tol = 1.0e-8_dp
   control%trace = 1
   call csolnp(problem, result, control)

   write(*, '(a)') ''
   write(*, '(a)') 'Powell exponential benchmark'
   write(*, '(a,es16.8)') 'objective:      ', result%objective
   write(*, '(a,5f13.7)') 'parameters:     ', result%pars
   write(*, '(a,es12.4)') 'eq violation:   ', result%kkt%eq_violation
   write(*, '(a,i0,2a)') 'convergence:    ', result%convergence, ' - ', trim(result%message)
end program demo_rsolnp
