program test_context_batch_noise
   use cec2005benchmark
   implicit none
   type(cec2005_context) :: ctx
   real(dp) :: x(3,10), f1(3), f2(3), clean, noisy
   integer :: i, ios

   do i = 1, 10
      x(1,i) = -3.0_dp + 0.2_dp*real(i,dp)
      x(2,i) =  1.5_dp - 0.1_dp*real(i,dp)
      x(3,i) =  0.05_dp*real(i*i,dp)
   end do
   call ctx%init(10,10,'data',.false.,ios)
   if (ios /= 0) error stop 'context initialization failed'
   call ctx%evaluate_batch(x,f1,ios)
   if (ios /= 0) error stop 'context batch failed'
   call cec2005_eval_batch(10,x,f2,'data',.false.,ios)
   if (ios /= 0) error stop 'convenience batch failed'
   if (maxval(abs(f1-f2)) > 1.0e-12_dp*max(1.0_dp,maxval(abs(f1)))) then
      error stop 'context and convenience batch differ'
   end if

   call ctx%init(4,10,'data',.false.,ios)
   clean = ctx%evaluate(x(1,:))
   call cec2005_seed(12345)
   call ctx%set_noise(.true.)
   noisy = ctx%evaluate(x(1,:))
   if (noisy < clean) error stop 'F4 noise factor should be nonnegative'

   print *, 'PASS test_context_batch_noise'
end program test_context_batch_noise
