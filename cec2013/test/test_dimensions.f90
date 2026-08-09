program test_dimensions
   use cec2013
   implicit none
   integer, parameter :: dims(12) = [2,5,10,20,30,40,50,60,70,80,90,100]
   type(cec2013_context) :: ctx
   real(dp), allocatable :: x(:)
   real(dp) :: got
   integer :: i, status

   do i = 1, size(dims)
      call ctx%init(dims(i), 'data', status)
      if (status /= CEC2013_OK) error stop 'valid dimension did not initialize'
      allocate(x(dims(i)))
      x = ctx%shift(1:dims(i))
      got = cec2013_evaluate(ctx, 1, x, status)
      if (status /= CEC2013_OK .or. abs(got+1400.0_dp) > 1.0e-10_dp) then
         error stop 'dimension regression failed'
      end if
      deallocate(x)
   end do

   call ctx%init(3, 'data', status)
   if (status /= CEC2013_BAD_DIMENSION) error stop 'invalid dimension accepted'
   print *, 'PASS test_dimensions'
end program test_dimensions
