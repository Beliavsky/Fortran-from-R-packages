program test_estimation
   use discrete_laplace
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer, parameter :: n=4000
   integer(int64) :: xr(n)
   integer :: x(n),i,fails,st
   real(dp) :: par(2),pp(2),pm(2),pmm(2)
   type(estdlaplace_result) :: f1
   fails=0
   call set_rng_seed(12345)
   call rdlaplace(xr,0.45_dp,0.25_dp)
   do i=1,n; x(i)=int(xr(i)); end do
   f1=estdlaplace(x)
   if(abs(f1%hatp-0.45_dp)>0.04_dp .or. abs(f1%hatq-0.25_dp)>0.04_dp) fails=fails+1

   call set_rng_seed(24680)
   call rdlaplace2(xr,0.45_dp,0.25_dp)
   do i=1,n; x(i)=int(xr(i)); end do
   par=estdlaplace2(x,'ML',status=st)
   if(st/=0 .or. abs(par(1)-0.45_dp)>0.05_dp .or. abs(par(2)-0.25_dp)>0.05_dp) fails=fails+1
   pm=estdlaplace2(x,'M',status=st)
   if(st/=0 .or. abs(pm(1)-0.45_dp)>0.07_dp .or. abs(pm(2)-0.25_dp)>0.07_dp) fails=fails+1
   pp=estdlaplace2(x,'P',status=st)
   if(st/=0 .or. any(pp<=0.0_dp) .or. any(pp>=1.0_dp)) fails=fails+1
   pmm=estdlaplace2(x,'MM',status=st)
   if(st/=0 .or. any(pmm<=0.0_dp) .or. any(pmm>=1.0_dp)) fails=fails+1
   if(fails/=0) then
      print *, 'test_estimation: FAIL',fails,par,pm,pp,pmm
      error stop 1
   end if
   print *, 'test_estimation: PASS'
end program
