program test_batch_status
   use cec2013
   implicit none
   type(cec2013_context) :: ctx
   real(dp) :: x(5,3), f(3), fsingle
   integer :: i, status

   call ctx%init(5, 'data', status)
   if (status /= CEC2013_OK) error stop 'context initialization failed'
   x(:,1) = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp]
   x(:,2) = ctx%shift(1:5)
   x(:,3) = [-7.0_dp, 6.0_dp, -5.0_dp, 4.0_dp, -3.0_dp]
   call cec2013_evaluate_batch(ctx, 8, x, f, status)
   if (status /= CEC2013_OK) error stop 'batch evaluation failed'
   do i = 1, 3
      fsingle = cec2013_evaluate(ctx, 8, x(:,i), status)
      if (abs(f(i)-fsingle) > 1.0e-11_dp*max(1.0_dp,abs(fsingle))) then
         error stop 'batch/single mismatch'
      end if
   end do
   fsingle = cec2013_evaluate(ctx, 29, x(:,1), status)
   if (status /= CEC2013_BAD_PROBLEM) error stop 'invalid problem accepted'
   if (fsingle < huge(1.0_dp)/2.0_dp) error stop 'invalid problem did not return sentinel'
   print *, 'PASS test_batch_status'
end program test_batch_status
