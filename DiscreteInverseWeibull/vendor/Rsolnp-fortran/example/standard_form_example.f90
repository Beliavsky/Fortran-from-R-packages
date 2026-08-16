! SPDX-License-Identifier: GPL-2.0-only
program standard_form_example
   use rsolnp, only : solnp_problem, solnp_problem_suite, solnp_standardize_problem
   implicit none

   type(solnp_problem) :: original, standardized
   integer :: status
   character(len=160) :: message

   call solnp_problem_suite('Other', 9, original, status, message)
   if (status /= 0) error stop trim(message)
   call solnp_standardize_problem(original, standardized, status, message)
   if (status /= 0) error stop trim(message)

   write(*, '(a,i0)') 'original two-sided inequalities: ', original%n_ineq
   write(*, '(a,i0)') 'standard one-sided inequalities: ', standardized%n_ineq
   write(*, '(a,l1)') 'standard-form flag: ', standardized%standard_form
end program standard_form_example
