program test_observed_info
   use discrete_laplace
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer, parameter :: n=1200
   integer(int64) :: xr(n)
   integer :: x(n),i,fails,st
   real(dp) :: c(2,2)
   fails=0
   call set_rng_seed(97531)
   call rdlaplace2(xr,0.35_dp,0.55_dp)
   do i=1,n; x(i)=int(xr(i)); end do
   c=iofi2(x,st)
   if(st/=0 .or. (.not. all(ieee_is_finite(c)))) fails=fails+1
   if(c(1,1)<=0.0_dp .or. c(2,2)<=0.0_dp) fails=fails+1
   if(abs(c(1,2)-c(2,1))>1.0e-12_dp) fails=fails+1
   if(fails/=0) then
      print *, 'test_observed_info: FAIL',fails,c
      error stop 1
   end if
   print *, 'test_observed_info: PASS'
end program
