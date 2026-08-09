program reuse_context
   use cec2005benchmark
   implicit none
   type(cec2005_context) :: ctx
   real(dp) :: x(3,10), f(3)
   integer :: i, ios
   do i = 1, 10
      x(:,i) = [0.1_dp*real(i,dp), -0.05_dp*real(i,dp), 0.02_dp*real(i,dp)]
   end do
   call ctx%init(21,10,'data',.false.,ios)
   if (ios /= 0) error stop 'CEC2005 initialization failed'
   call ctx%evaluate_batch(x,f,ios)
   if (ios /= 0) error stop 'CEC2005 batch evaluation failed'
   print '(a,3(1x,es14.6))', 'F21 values:', f
end program reuse_context
