! SPDX-License-Identifier: GPL-2.0-only
program multistart_example
   use rsolnp, only : solnp_problem, solnp_control, multistart_result, &
      solnp_problem_suite, csolnp_ms
   implicit none

   type(solnp_problem) :: problem
   type(solnp_control) :: control
   type(multistart_result) :: result
   integer :: status
   character(len=160) :: message

   call solnp_problem_suite('hs', 5, problem, status, message)
   if (status /= 0) error stop trim(message)
   control%max_iter = 100
   control%min_iter = 400
   call csolnp_ms(problem, 12, result, control, seed=2026)

   write(*, '(a,i0)') 'best start:  ', result%best_index
   write(*, '(a,es16.8)') 'objective:   ', result%best%objective
   write(*, '(a,2f14.8)') 'parameters:  ', result%best%pars
end program multistart_example
