program fit_example
   use fattailsr, only : dp, kiener_parameters, make_k4, qkiener, fit_kiener_k4
   implicit none
   integer, parameter :: n=199
   real(dp) :: x(n), p
   type(kiener_parameters) :: truth, fitted
   integer :: i
   truth=make_k4(0.2_dp,1.1_dp,5.0_dp,0.2_dp)
   do i=1,n
      p=real(i,dp)/real(n+1,dp)
      x(i)=qkiener(p,truth)
   end do
   fitted=fit_kiener_k4(x,maxk=20.0_dp,mink=1.53_dp,maxe=0.5_dp)
   print '(a,4(1x,f10.5))','truth  m g k e:',truth%m,truth%g,truth%k,truth%e
   print '(a,4(1x,f10.5))','fitted m g k e:',fitted%m,fitted%g,fitted%k,fitted%e
end program fit_example
