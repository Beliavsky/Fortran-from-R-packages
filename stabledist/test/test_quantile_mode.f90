program test_quantile_mode
   use r_compat, only: dp, qnorm
   use stabledist
   implicit none
   real(dp) :: q,p,u,qlev,q0,m,d0,dl,dr
   integer :: i
   real(dp), parameter :: probs(5)=[0.05_dp,0.2_dp,0.5_dp,0.8_dp,0.95_dp]

   do i=1,size(probs)
      q=qstable(probs(i),1.5_dp,0.4_dp,tol=2e-7_dp,integ_tol=2e-9_dp,subdivisions=1200)
      p=pstable(q,1.5_dp,0.4_dp,tol=2e-9_dp,subdivisions=1200)
      if(abs(p-probs(i))>2e-6_dp)then
         print *,'p-q inversion',probs(i),q,p
         error stop 1
      end if
   end do

   ! Exact Levy quantile: F(x)=erfc(sqrt(1/(2x))).  S0 has a -1 shift here.
   p=0.6_dp
   u=qnorm(1.0_dp-p/2.0_dp)
   qlev=1.0_dp/(u*u)
   q0=qstable(p,0.5_dp,1.0_dp,tol=1e-10_dp,integ_tol=1e-10_dp,subdivisions=2000,maxiter=300)
   if(abs(q0-(qlev-1.0_dp))>2e-7_dp)then
      print *,'Levy/S0 quantile',q0,qlev-1.0_dp
      error stop 1
   end if

   m=stable_mode(1.2_dp,0.8_dp,tol=1e-9_dp)
   if(abs(m-(-0.2697162360190157_dp))>2e-6_dp)then
      print *,'mode ref',m
      error stop 1
   end if

   ! S2 is mode-centered by construction.
   d0=dstable(0.4_dp,1.2_dp,0.8_dp,gamma=1.7_dp,delta=0.4_dp,pm=2,tol=1e-9_dp)
   dl=dstable(0.39_dp,1.2_dp,0.8_dp,gamma=1.7_dp,delta=0.4_dp,pm=2,tol=1e-9_dp)
   dr=dstable(0.41_dp,1.2_dp,0.8_dp,gamma=1.7_dp,delta=0.4_dp,pm=2,tol=1e-9_dp)
   if(d0<dl .or. d0<dr)error stop 'S2 not mode-centered'

   print *, 'test_quantile_mode: PASS'
end program test_quantile_mode
