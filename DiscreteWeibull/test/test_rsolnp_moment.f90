program test_rsolnp_moment
   use discrete_weibull
   implicit none
   integer(i64), allocatable :: x(:)
   type(dweibull_fit_result) :: fit
   integer :: fails

   fails = 0
   allocate(x(40))
   call seed_rng(24680)
   call rdweibull(x,0.55_dp,1.3_dp,.false.)
   fit = estdweibull(x,"M",.false.,1.0e-5_dp,5000)
   if (fit%status /= 0 .and. fit%status /= 1) fails=fails+1
   if (fit%pars(1) < 0.0_dp .or. fit%pars(1) >= 1.0_dp) fails=fails+1
   if (fit%pars(2) <= 0.0_dp .or. fit%pars(2) > 100.0_dp) fails=fails+1
   if (fit%objective > 0.5_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_rsolnp_moment: FAIL", fails, fit%status, fit%pars, fit%objective
      error stop 1
   end if
   print *, "test_rsolnp_moment: PASS"

contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: s(:)
      call random_seed(size=n)
      allocate(s(n))
      do i=1,n
         s(i)=mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=s)
   end subroutine seed_rng
end program test_rsolnp_moment
