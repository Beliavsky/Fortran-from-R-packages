program test_core
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_ph
   use matrixdist_dph
   use matrixdist_iph
   implicit none
   real(dp) :: a(1),s(1,1),x,p,q,tol,z(1)
   integer :: k(1)
   a=1.0_dp
   s(1,1)=-2.0_dp
   tol=2.0e-10_dp
   x=0.7_dp
   call chk(ph_density(x,a,s),2.0_dp*exp(-1.4_dp),tol,'ph density')
   call chk(ph_cdf(x,a,s),1.0_dp-exp(-1.4_dp),tol,'ph cdf')
   call chk(ph_survival(x,a,s),exp(-1.4_dp),tol,'ph survival')
   call chk(ph_mean(a,s),0.5_dp,tol,'ph mean')
   call chk(ph_variance(a,s),0.25_dp,tol,'ph variance')
   call chk(ph_moment(3,a,s),0.75_dp,tol,'ph moment3')
   call chk(ph_laplace(0.3_dp,a,s),2.0_dp/2.3_dp,tol,'ph laplace')
   p=0.8_dp
   q=ph_quantile(p,a,s)
   call chk(q,-log(0.2_dp)/2.0_dp,2e-8_dp,'ph quantile')
   call chk(iph_density(x,a,s,'weibull',[2.0_dp]),4.0_dp*x*exp(-2.0_dp*x*x),2e-9_dp,'weibull density')
   call chk(iph_cdf(x,a,s,'weibull',[2.0_dp]),1.0_dp-exp(-2.0_dp*x*x),2e-9_dp,'weibull cdf')
   s(1,1)=-1.0_dp
   call chk(mgev_density(x,a,s,[0.0_dp,1.0_dp,0.0_dp]),exp(-x-exp(-x)),2e-9_dp,'gev density')
   call chk(mgev_cdf(x,a,s,[0.0_dp,1.0_dp,0.0_dp]),exp(-exp(-x)),2e-9_dp,'gev cdf')


   s(1,1)=0.6_dp
   k=3
   z=0.8_dp
   call chk(dph_density(k(1),a,s),0.4_dp*0.6_dp**2,tol,'dph density')
   call chk(dph_cdf(k(1),a,s),1.0_dp-0.6_dp**3,tol,'dph cdf')
   call chk(dph_mean(a,s),2.5_dp,tol,'dph mean')
   call chk(dph_variance(a,s),3.75_dp,tol,'dph variance')
   call chk(dph_pgf(z(1),a,s),0.4_dp*z(1)/(1.0_dp-0.6_dp*z(1)),tol,'dph pgf')
   print *, 'test_core: PASS'
contains
   subroutine chk(got,want,eps,label)
      real(dp),intent(in)::got,want,eps
      character(len=*),intent(in)::label
      if(abs(got-want)>eps*max(1.0_dp,abs(want))) then
         print *,trim(label),got,want
         error stop 1
      end if
   end subroutine
end program
