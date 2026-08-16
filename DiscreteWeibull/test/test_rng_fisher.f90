program test_rng_fisher
   use discrete_weibull
   implicit none
   integer, parameter :: n=50000
   integer(i64), allocatable :: x(:)
   real(dp) :: emp,theo
   type(fisher_result) :: fr
   integer :: fails

   fails = 0
   allocate(x(n))
   call seed_rng(13579)
   call rdweibull(x,0.45_dp,1.0_dp,.false.)
   emp = sum(real(x,dp))/real(n,dp)
   theo = 1.0_dp/(1.0_dp-0.45_dp)
   if (abs(emp-theo) > 0.04_dp) fails=fails+1

   fr = varFisher(x(:1000),.false.)
   if (fr%status /= 0) fails=fails+1
   if (fr%information(1,1) <= 0.0_dp .or. fr%information(2,2) <= 0.0_dp) &
      fails=fails+1
   if (fr%inverse(1,1) <= 0.0_dp .or. fr%inverse(2,2) <= 0.0_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_rng_fisher: FAIL", fails, emp, fr%mle
      error stop 1
   end if
   print *, "test_rng_fisher: PASS"

contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: ns,i
      integer, allocatable :: s(:)
      call random_seed(size=ns)
      allocate(s(ns))
      do i=1,ns
         s(i)=mod(seed+65537*i,2147483646)+1
      end do
      call random_seed(put=s)
   end subroutine seed_rng
end program test_rng_fisher
