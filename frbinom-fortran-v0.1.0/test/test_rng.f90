program test_rng
   use frbinom
   implicit none
   integer, parameter :: n=120000
   integer, allocatable :: x(:)
   real(dp), allocatable :: pmf(:)
   real(dp) :: em,tm,ev,tv
   integer :: k,fails,status

   fails = 0
   allocate(x(n))

   call set_frbinom_seed(12345)
   call rfrbinom(x,30,0.6_dp,0.7_dp,0.2_dp)
   call frbinom_pmf_table(30,0.6_dp,0.7_dp,0.2_dp,.false.,pmf,status)
   tm = 0.0_dp
   tv = 0.0_dp
   do k = 0, 30
      tm = tm+real(k,dp)*pmf(k)
   end do
   do k = 0, 30
      tv = tv+(real(k,dp)-tm)**2*pmf(k)
   end do
   em = sum(real(x,dp))/real(n,dp)
   ev = sum((real(x,dp)-em)**2)/real(n-1,dp)
   if (abs(em-tm) > 0.08_dp) fails=fails+1
   if (abs(ev-tv) > 0.70_dp) fails=fails+1

   call set_frbinom_seed(54321)
   call rfrbinom2(x,30,0.8_dp,0.2_dp,0.1_dp)
   call frbinom2_pmf_table(30,0.8_dp,0.2_dp,0.1_dp,.false.,pmf,status)
   tm = 0.0_dp
   tv = 0.0_dp
   do k = 0, 30
      tm = tm+real(k,dp)*pmf(k)
   end do
   do k = 0, 30
      tv = tv+(real(k,dp)-tm)**2*pmf(k)
   end do
   em = sum(real(x,dp))/real(n,dp)
   ev = sum((real(x,dp)-em)**2)/real(n-1,dp)
   if (abs(em-tm) > 0.025_dp) fails=fails+1
   if (abs(ev-tv) > 0.08_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_rng: FAIL", fails, em,tm,ev,tv
      error stop 1
   end if
   print *, "test_rng: PASS"
end program test_rng
