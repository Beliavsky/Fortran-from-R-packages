program test_validation
   use cec2005benchmark
   implicit none
   type(cec2005_context) :: ctx
   real(dp) :: x3(3), f
   integer :: ios

   x3 = 0.0_dp
   call ctx%init(0,10,'data',.false.,ios)
   if (ios /= 1) error stop 'invalid function id not rejected'
   call ctx%init(1,3,'data',.false.,ios)
   if (ios /= 2) error stop 'invalid dimension not rejected'
   call ctx%init(1,2,'data',.false.,ios)
   if (ios /= 0) error stop 'valid context rejected'
   f = ctx%evaluate(x3,ios)
   if (ios /= 1) error stop 'wrong vector length not rejected'
   if (f <= 0.0_dp) then
      ! The sentinel is huge(); this branch only silences unused-result concerns.
      error stop 'unexpected invalid-evaluation sentinel'
   end if
   print *, 'PASS test_validation'
end program test_validation
