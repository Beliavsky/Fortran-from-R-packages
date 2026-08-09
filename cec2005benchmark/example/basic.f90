program basic
   use cec2005benchmark
   implicit none
   real(dp) :: x(10), f
   integer :: ios
   x = 0.0_dp
   f = cec2005_eval(1,x,'data',.false.,ios)
   if (ios /= 0) error stop 'CEC2005 initialization failed'
   print '(a,es16.8)', 'F1(0) = ', f
end program basic
