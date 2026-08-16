module rfast_distributions
   use, intrinsic :: iso_fortran_env, only : int64
   use rfast_special, only : dp, pi, log1p_r, nan_r
   use rfast_linalg, only : cholesky_lower, inverse_matrix, logdet_spd, mahalanobis_sq
   use zigg, only : ziggurat_rng
   implicit none
   private
   type(ziggurat_rng), save :: rng
   public :: set_seed, rnorm_vector, rexp_vector, runif_vector
   public :: dmvnorm, dmvt, rmvnorm, rmvt, rmvlaplace
   public :: gamma_random, chisq_random, vonmises_random

contains

   subroutine set_seed(seed)
      integer(int64),intent(in)::seed
      call rng%set_seed(seed)
   end subroutine set_seed

   function rnorm_vector(n,mean,sd) result(x)
      integer,intent(in)::n;real(dp),intent(in),optional::mean,sd;real(dp),allocatable::x(:);real(dp)::m,s
      integer::i
      m=0.0_dp;if(present(mean))m=mean;s=1.0_dp;if(present(sd))s=sd;allocate(x(max(0,n)))
      do i=1,n;x(i)=m+s*rng%rnorm();end do
   end function rnorm_vector

   function rexp_vector(n,rate) result(x)
      integer,intent(in)::n;real(dp),intent(in),optional::rate;real(dp),allocatable::x(:);real(dp)::r;integer::i
      r=1.0_dp;if(present(rate))r=rate;allocate(x(max(0,n)));do i=1,n;x(i)=rng%rexp()/r;end do
   end function rexp_vector

   function runif_vector(n,a,b) result(x)
      integer,intent(in)::n;real(dp),intent(in),optional::a,b;real(dp),allocatable::x(:);real(dp)::lo,hi;integer::i
      lo=0.0_dp;if(present(a))lo=a;hi=1.0_dp;if(present(b))hi=b;allocate(x(max(0,n)))
      do i=1,n;x(i)=lo+(hi-lo)*rng%runi();end do
   end function runif_vector

   recursive real(dp) function gamma_random(shape,scale) result(x)
      real(dp),intent(in)::shape;real(dp),intent(in),optional::scale
      real(dp)::d,c,z,u,s
      s=1.0_dp;if(present(scale))s=scale
      if(shape<=0.0_dp)then;x=nan_r();return;end if
      if(shape<1.0_dp)then
         u=max(tiny(1.0_dp),rng%runi());x=gamma_random(shape+1.0_dp)*u**(1.0_dp/shape)*s;return
      end if
      d=shape-1.0_dp/3.0_dp;c=1.0_dp/sqrt(9.0_dp*d)
      do
         z=rng%rnorm();u=rng%runi()
         if(1.0_dp+c*z<=0.0_dp)cycle
         x=(1.0_dp+c*z)**3
         if(u<1.0_dp-0.0331_dp*z**4)exit
         if(log(u)<0.5_dp*z*z+d*(1.0_dp-x+log(x)))exit
      end do
      x=d*x*s
   end function gamma_random

   real(dp) function chisq_random(df) result(x)
      real(dp),intent(in)::df;x=gamma_random(0.5_dp*df,2.0_dp)
   end function chisq_random

   function dmvnorm(x,mu,sigma,log_density) result(d)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:);logical,intent(in),optional::log_density
      real(dp)::d(size(x,1)),inv(size(sigma,1),size(sigma,2)),z(size(mu)),ld
      logical::lg;integer::i,info,p
      lg=.false.;if(present(log_density))lg=log_density;p=size(mu)
      call inverse_matrix(sigma,inv,info);ld=logdet_spd(sigma,info)
      if(info/=0)then;d=nan_r();return;end if
      do i=1,size(x,1)
         z=x(i,:)-mu;d(i)=-0.5_dp*(real(p,dp)*log(2.0_dp*pi)+ld+dot_product(z,matmul(inv,z)))
         if(.not.lg)d(i)=exp(d(i))
      end do
   end function dmvnorm

   function dmvt(x,mu,sigma,df,log_density) result(d)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),df;logical,intent(in),optional::log_density
      real(dp)::d(size(x,1)),inv(size(sigma,1),size(sigma,2)),z(size(mu)),ld,q
      logical::lg;integer::i,info,p
      lg=.false.;if(present(log_density))lg=log_density;p=size(mu)
      call inverse_matrix(sigma,inv,info);ld=logdet_spd(sigma,info)
      if(info/=0)then;d=nan_r();return;end if
      do i=1,size(x,1)
         z=x(i,:)-mu;q=dot_product(z,matmul(inv,z))
         d(i)=log_gamma(0.5_dp*(df+p))-log_gamma(0.5_dp*df)-0.5_dp*(real(p,dp)*log(df*pi)+ld) &
              -0.5_dp*(df+p)*log1p_r(q/df)
         if(.not.lg)d(i)=exp(d(i))
      end do
   end function dmvt

   function rmvnorm(n,mu,sigma) result(x)
      integer,intent(in)::n;real(dp),intent(in)::mu(:),sigma(:,:);real(dp)::x(n,size(mu)),l(size(mu),size(mu)),z(size(mu))
      integer::i,j,info
      call cholesky_lower(sigma,l,info)
      if(info/=0)then;x=nan_r();return;end if
      do i=1,n
         do j=1,size(mu);z(j)=rng%rnorm();end do
         x(i,:)=mu+matmul(l,z)
      end do
   end function rmvnorm

   function rmvt(n,mu,sigma,df) result(x)
      integer,intent(in)::n;real(dp),intent(in)::mu(:),sigma(:,:),df;real(dp)::x(n,size(mu)),l(size(mu),size(mu)),z(size(mu)),w
      integer::i,j,info
      call cholesky_lower(sigma,l,info);if(info/=0)then;x=nan_r();return;end if
      do i=1,n
         do j=1,size(mu);z(j)=rng%rnorm();end do;w=sqrt(chisq_random(df)/df);x(i,:)=mu+matmul(l,z)/w
      end do
   end function rmvt

   function rmvlaplace(n,mu,sigma) result(x)
      integer,intent(in)::n;real(dp),intent(in)::mu(:),sigma(:,:);real(dp)::x(n,size(mu)),l(size(mu),size(mu)),z(size(mu)),w
      integer::i,j,info
      call cholesky_lower(sigma,l,info);if(info/=0)then;x=nan_r();return;end if
      do i=1,n
         do j=1,size(mu);z(j)=rng%rnorm();end do;w=sqrt(2.0_dp*rng%rexp());x(i,:)=mu+w*matmul(l,z)
      end do
   end function rmvlaplace

   real(dp) function vonmises_random(mu,kappa) result(theta)
      real(dp),intent(in)::mu,kappa;real(dp)::a,b,r,u1,u2,u3,z,f,c
      if(kappa<1.0e-8_dp)then;theta=mu+2.0_dp*pi*(rng%runi()-0.5_dp);return;end if
      a=1.0_dp+sqrt(1.0_dp+4.0_dp*kappa*kappa);b=(a-sqrt(2.0_dp*a))/(2.0_dp*kappa);r=(1.0_dp+b*b)/(2.0_dp*b)
      do
         u1=rng%runi();z=cos(pi*u1);f=(1.0_dp+r*z)/(r+z);c=kappa*(r-f);u2=rng%runi()
         if(c*(2.0_dp-c)-u2>0.0_dp)exit
         if(log(c/u2)+1.0_dp-c>=0.0_dp)exit
      end do
      u3=rng%runi();if(u3>0.5_dp)then;theta=mu+acos(f);else;theta=mu-acos(f);end if
      theta=modulo(theta+pi,2.0_dp*pi)-pi
   end function vonmises_random

end module rfast_distributions
