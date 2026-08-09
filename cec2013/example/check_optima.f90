program check_optima
   use cec2013
   implicit none
   type(cec2013_context) :: ctx
   real(dp) :: x(10), f
   integer :: p, status

   call ctx%init(10, 'data', status)
   if (status /= CEC2013_OK) error stop 'could not load CEC2013 data'
   x = ctx%shift(1:10)
   do p = 1, 28
      f = cec2013_evaluate(ctx, p, x, status)
      print '(i3,2x,f10.3,2x,f10.3)', p, f, cec2013_optimum_value(p)
   end do
end program check_optima
