program test_optima
   use cec2013
   implicit none
   type(cec2013_context) :: ctx
   real(dp) :: x(10), got, want
   integer :: i, status

   call ctx%init(10, 'data', status)
   if (status /= CEC2013_OK) error stop 'context initialization failed'
   x = ctx%shift(1:10)
   do i = 1, 28
      got = cec2013_evaluate(ctx, i, x, status)
      want = cec2013_optimum_value(i)
      if (status /= CEC2013_OK) error stop 'evaluation failed'
      if (abs(got-want) > 1.0e-9_dp*max(1.0_dp,abs(want))) then
         write(*,'(a,i0,2(1x,es24.16))') 'optimum mismatch for problem', i, got, want
         error stop 1
      end if
   end do
   print *, 'PASS test_optima'
end program test_optima
