program test_roundtrip
   use rmutil
   implicit none
   real(dp), parameter :: p0=0.37_dp
   real(dp) :: x

   x=qinvgauss(p0,2.0_dp,0.3_dp); call checkp(pinvgauss(x,2.0_dp,0.3_dp),"invgauss",2e-10_dp)
   x=qlaplace(p0,0.2_dp,1.3_dp); call checkp(plaplace(x,0.2_dp,1.3_dp),"laplace",2e-12_dp)
   x=qlevy(p0,0.2_dp,1.2_dp); call checkp(plevy(x,0.2_dp,1.2_dp),"levy",2e-11_dp)
   x=qpareto(p0,2.1_dp,4.5_dp); call checkp(ppareto(x,2.1_dp,4.5_dp),"pareto",2e-12_dp)
   x=qsimplex(p0,0.35_dp,0.4_dp); call checkp(psimplex(x,0.35_dp,0.4_dp),"simplex",3e-7_dp)
   x=qtwosidedpower(p0,0.35_dp,2.2_dp); call checkp(ptwosidedpower(x,0.35_dp,2.2_dp),"two-sided power",2e-12_dp)
   x=qboxcox(p0,1.2_dp,0.6_dp,0.7_dp); call checkp(pboxcox(x,1.2_dp,0.6_dp,0.7_dp),"Box-Cox",3e-10_dp)
   x=qburr(p0,1.7_dp,2.1_dp,1.4_dp); call checkp(pburr(x,1.7_dp,2.1_dp,1.4_dp),"Burr",2e-12_dp)
   x=qgextval(p0,2.0_dp,1.5_dp,0.8_dp); call checkp(pgextval(x,2.0_dp,1.5_dp,0.8_dp),"generalized extreme value",3e-10_dp)
   x=qggamma(p0,1.7_dp,3.0_dp,1.2_dp); call checkp(pggamma(x,1.7_dp,3.0_dp,1.2_dp),"generalized gamma",3e-9_dp)
   x=qginvgauss(p0,1.4_dp,0.8_dp,0.7_dp); call checkp(pginvgauss(x,1.4_dp,0.8_dp,0.7_dp),"generalized inverse Gaussian",5e-6_dp)
   x=qglogis(p0,0.2_dp,1.3_dp,1.7_dp); call checkp(pglogis(x,0.2_dp,1.3_dp,1.7_dp),"generalized logistic",2e-12_dp)
   x=qgweibull(p0,1.8_dp,2.0_dp,1.3_dp); call checkp(pgweibull(x,1.8_dp,2.0_dp,1.3_dp),"generalized Weibull",2e-12_dp)
   x=qhjorth(p0,1.2_dp,0.4_dp,0.9_dp); call checkp(phjorth(x,1.2_dp,0.4_dp,0.9_dp),"Hjorth",3e-10_dp)
   x=qpowexp(p0,0.2_dp,1.4_dp,1.3_dp); call checkp(ppowexp(x,0.2_dp,1.4_dp,1.3_dp),"power exponential",3e-9_dp)
   ! The upstream qskewlaplace branch point is 0.5; use f=1 where this is exact.
   x=qskewlaplace(p0,0.2_dp,1.3_dp,1.0_dp); call checkp(pskewlaplace(x,0.2_dp,1.3_dp,1.0_dp),"skew Laplace",2e-12_dp)

   call check(qbetabinom(0.5_dp,10,0.3_dp,4.0_dp)>=0,"beta-binomial quantile")
   call check(qdoublebinom(0.5_dp,10,0.3_dp,1.4_dp)>=0,"double-binomial quantile")
   call check(qmultbinom(0.5_dp,10,0.3_dp,0.7_dp)>=0,"multiplicative-binomial quantile")
   call check(qdoublepois(0.5_dp,3.0_dp,1.3_dp)>=0,"double-Poisson quantile")
   call check(qmultpois(0.5_dp,3.0_dp,0.8_dp)>=0,"multiplicative-Poisson quantile")
   call check(qpvfpois(0.5_dp,3.0_dp,0.8_dp,0.4_dp)>=0,"PVF-Poisson quantile")
   call check(qgammacount(0.5_dp,3.0_dp,1.3_dp)>=0,"gamma-count quantile")
   call check(qconsul(0.5_dp,3.0_dp,0.2_dp)>=0,"Consul quantile")

   print *, "test_roundtrip: PASS"
contains
   subroutine checkp(p,name,tol)
      real(dp),intent(in)::p,tol
      character(*),intent(in)::name
      call check(abs(p-p0)<tol,name)
   end subroutine checkp
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *,"FAIL: ",trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_roundtrip
