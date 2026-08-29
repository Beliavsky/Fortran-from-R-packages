program test_multivariate
   use r_compat, only: dp
   use matrixdist_multivariate
   use matrixdist_miph
   implicit none
   real(dp)::a(1),s11(1,1),s12(1,1),s22(1,1),f,want,c(2,2),m(2),tol
   real(dp)::s3(1,1,2),rew(1,2),cm(2,2),beta(3,2),x(2)
   integer::kk(2)
   character(len=12)::kinds(2)
   a=1.0_dp
   s11=-2.0_dp
   s12=2.0_dp
   s22=-3.0_dp
   tol=3e-9_dp
   f=bivph_density(0.4_dp,0.7_dp,a,s11,s12,s22)
   want=6.0_dp*exp(-2.0_dp*0.4_dp-3.0_dp*0.7_dp)
   call chk(f,want,tol,'bivph density')
   m=bivph_mean(a,s11,s12,s22)
   call chk(m(1),0.5_dp,tol,'bivph mean1')
   call chk(m(2),1.0_dp/3.0_dp,tol,'bivph mean2')
   c=bivph_cov(a,s11,s12,s22)
   call chk(c(1,2),0.0_dp,2e-8_dp,'bivph cov')

   s11=0.6_dp
   s12=0.4_dp
   s22=0.7_dp
   f=bivdph_density(3,2,a,s11,s12,s22)
   want=0.4_dp*0.6_dp**2*0.3_dp*0.7_dp
   call chk(f,want,tol,'bivdph density')
   m=bivdph_mean(a,s11,s12,s22)
   call chk(m(1),2.5_dp,tol,'bivdph mean1')
   call chk(m(2),1.0_dp/0.3_dp,tol,'bivdph mean2')
   c=bivdph_cov(a,s11,s12,s22)
   call chk(c(1,2),0.0_dp,2e-8_dp,'bivdph cov')

   s3(1,1,1)=-2.0_dp
   s3(1,1,2)=-3.0_dp
   x=[0.4_dp,0.7_dp]
   f=mph_density_point(x,a,s3)
   want=(2.0_dp*exp(-0.8_dp))*(3.0_dp*exp(-2.1_dp))
   call chk(f,want,tol,'mph density')
   m=mph_mean(a,s3)
   call chk(m(1),0.5_dp,tol,'mph mean1')
   call chk(m(2),1.0_dp/3.0_dp,tol,'mph mean2')
   c=mph_cov(a,s3)
   call chk(c(1,2),0.0_dp,2e-8_dp,'mph cov')

   s3(1,1,1)=0.6_dp
   s3(1,1,2)=0.7_dp
   kk=[3,2]
   f=mdph_density_point(kk,a,s3)
   want=(0.4_dp*0.6_dp**2)*(0.3_dp*0.7_dp)
   call chk(f,want,tol,'mdph density')
   m=mdph_mean(a,s3)
   call chk(m(1),2.5_dp,tol,'mdph mean1')
   call chk(m(2),1.0_dp/0.3_dp,tol,'mdph mean2')

   s11=-2.0_dp
   s12=2.0_dp
   s22=-3.0_dp
   kinds=['weibull     ','weibull     ']
   beta=0.0_dp
   beta(1,:)=[2.0_dp,1.5_dp]
   f=biviph_density(x,a,s11,s12,s22,kinds,beta)
   want=bivph_density(x(1)**2,x(2)**1.5_dp,a,s11,s12,s22)*(2*x(1))*(1.5_dp*sqrt(x(2)))
   call chk(f,want,2e-8_dp,'biviph density')

   s11=-2.0_dp
   rew(1,:)=[1.0_dp,2.0_dp]
   m=mphstar_mean(a,s11,rew)
   call chk(m(1),0.5_dp,tol,'mphstar mean1')
   call chk(m(2),1.0_dp,tol,'mphstar mean2')
   cm=mphstar_cov(a,s11,rew)
   call chk(cm(1,2),0.5_dp,2e-8_dp,'mphstar covariance')
   print *, 'test_multivariate: PASS'
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
