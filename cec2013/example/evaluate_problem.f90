program evaluate_problem
   use cec2013
   implicit none
   type(cec2013_context) :: ctx
   real(dp) :: x(10), f
   integer :: status

   x = 0.0_dp
   call ctx%init(10, 'data', status)
   if (status /= CEC2013_OK) error stop 'could not load CEC2013 data'
   f = cec2013_evaluate(ctx, 1, x, status)
   print '(a,es18.10)', 'CEC2013 F1 at x=0: ', f
end program evaluate_problem
