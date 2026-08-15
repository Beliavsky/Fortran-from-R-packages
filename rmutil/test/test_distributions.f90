program test_distributions
   use rmutil
   implicit none
   real(dp) :: p, x, sm
   integer :: k

   call check_close(pinvgauss(1.4_dp,2.0_dp,0.3_dp),0.4438189599194026_dp,2.0e-13_dp,"pinvgauss")
   call check_close(dinvgauss(1.4_dp,2.0_dp,0.3_dp),0.3950261383958223_dp,2.0e-13_dp,"dinvgauss")
   call check_close(plaplace(-0.3_dp,0.2_dp,1.3_dp),0.3403561991616927_dp,2.0e-14_dp,"plaplace")
   call check_close(ppareto(3.2_dp,2.1_dp,4.5_dp),0.8033667499402304_dp,2.0e-14_dp,"ppareto")
   call check_close(pggamma(2.2_dp,1.7_dp,3.0_dp,1.2_dp),0.47247434923591164_dp,2.0e-13_dp,"pggamma")
   call check_close(ppowexp(1.3_dp,0.4_dp,2.2_dp,1.0_dp),0.7280014976064876_dp,2.0e-13_dp,"ppowexp normal")
   call check_close(dbetabinom(3,10,0.3_dp,4.0_dp),0.13505265664000032_dp,3.0e-14_dp,"dbetabinom")
   call check_close(pbetabinom(3,10,0.3_dp,4.0_dp),0.6257311866880011_dp,3.0e-14_dp,"pbetabinom")
   call check_close(dgammacount(2,3.2_dp,1.4_dp),0.23563415559132306_dp,3.0e-13_dp,"dgammacount")
   call check_close(pgammacount(2,3.2_dp,1.4_dp),0.38349461299768356_dp,3.0e-13_dp,"pgammacount")

   p=0.73_dp; x=qinvgauss(p,2.0_dp,0.3_dp)
   call check_close(pinvgauss(x,2.0_dp,0.3_dp),p,2.0e-10_dp,"inverse Gaussian quantile")
   x=qsimplex(0.67_dp,0.4_dp,0.2_dp)
   call check_close(psimplex(x,0.4_dp,0.2_dp),0.67_dp,4.0e-7_dp,"simplex quantile")
   x=qpowexp(0.83_dp,-0.2_dp,1.7_dp,0.7_dp)
   call check_close(ppowexp(x,-0.2_dp,1.7_dp,0.7_dp),0.83_dp,2.0e-11_dp,"power exponential quantile")

   sm=0.0_dp
   do k=0,12
      sm=sm+ddoublebinom(k,12,0.35_dp,1.4_dp)
   end do
   call check_close(sm,1.0_dp,2.0e-13_dp,"double binomial normalization")

   sm=0.0_dp
   do k=0,60
      sm=sm+ddoublepois(k,7.0_dp,0.8_dp)
   end do
   call check(abs(sm-1.0_dp)<3.0e-8_dp,"double Poisson normalization")

   sm=0.0_dp
   do k=0,60
      sm=sm+dmultpois(k,4.0_dp,0.85_dp)
   end do
   call check(abs(sm-1.0_dp)<3.0e-8_dp,"multiplicative Poisson normalization")

   call check_close(ppvfpois(0,3.0_dp,1.2_dp,0.4_dp),dpvfpois(0,3.0_dp,1.2_dp,0.4_dp), &
      2.0e-14_dp,"PVF CDF at zero")
   p=0.8_dp; k=qpvfpois(p,4.0_dp,1.1_dp,0.35_dp)
   call check(ppvfpois(k,4.0_dp,1.1_dp,0.35_dp)>=p,"PVF quantile")

   x=qginvgauss(0.45_dp,1.3_dp,0.7_dp,0.8_dp)
   call check_close(pginvgauss(x,1.3_dp,0.7_dp,0.8_dp),0.45_dp,3.0e-5_dp,"GIG quantile")

   print *, "test_distributions: PASS"
contains
   subroutine check_close(x,y,tol,msg)
      real(dp),intent(in)::x,y,tol; character(*),intent(in)::msg
      call check(abs(x-y)<=tol,msg)
   end subroutine check_close
   subroutine check(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok)then
         print *, "FAIL: ", trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_distributions
