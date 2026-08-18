! SPDX-License-Identifier: GPL-3.0-only
module mcmcpack_rng
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_math, only : normal_cdf,normal_quantile
   use mcmcpack_linalg, only : chol_lower
   implicit none
   private
   public :: set_seed, runif, rnorm, rgamma_mt, rchisq, rbeta, rinvgamma_rng
   public :: rdirichlet_rng, rmvnorm, rtruncnorm, rinvgauss
contains
   subroutine set_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: put(:)
      call random_seed(size=n); allocate(put(n))
      do i=1,n
         put(i)=modulo(seed+104729*i+37*i*i,huge(1)-1)+1
      end do
      call random_seed(put=put)
   end subroutine set_seed

   real(dp) function runif() result(u)
      call random_number(u)
      if (u <= 0.0_dp) u=tiny(1.0_dp)
      if (u >= 1.0_dp) u=1.0_dp-epsilon(1.0_dp)
   end function runif

   real(dp) function rnorm(mu,sd) result(x)
      real(dp),intent(in),optional::mu,sd
      real(dp)::u1,u2,m,s
      m=0.0_dp; s=1.0_dp; if(present(mu))m=mu; if(present(sd))s=sd
      u1=runif(); u2=runif(); x=m+s*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function rnorm

   recursive real(dp) function rgamma_mt(shape,scale) result(x)
      real(dp),intent(in)::shape
      real(dp),intent(in),optional::scale
      real(dp)::d,c,z,u,sc
      sc=1.0_dp; if(present(scale))sc=scale
      if(shape<=0.0_dp .or. sc<=0.0_dp) then; x=0.0_dp; return; end if
      if(shape<1.0_dp) then
         x=rgamma_mt(shape+1.0_dp,sc)*runif()**(1.0_dp/shape); return
      end if
      d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
      do
         do; z=rnorm(); if(1.0_dp+c*z>0.0_dp)exit; end do
         x=(1.0_dp+c*z)**3; u=runif()
         if(u<1.0_dp-0.0331_dp*z**4)exit
         if(log(u)<0.5_dp*z*z+d*(1.0_dp-x+log(x)))exit
      end do
      x=sc*d*x
   end function rgamma_mt

   real(dp) function rchisq(df) result(x)
      real(dp),intent(in)::df
      x=rgamma_mt(0.5_dp*df,2.0_dp)
   end function rchisq

   real(dp) function rbeta(a,b) result(x)
      real(dp),intent(in)::a,b
      real(dp)::g1,g2
      g1=rgamma_mt(a); g2=rgamma_mt(b); x=g1/(g1+g2)
   end function rbeta

   real(dp) function rinvgamma_rng(shape,scale) result(x)
      real(dp),intent(in)::shape,scale
      x=1.0_dp/rgamma_mt(shape,1.0_dp/scale)
   end function rinvgamma_rng

   subroutine rdirichlet_rng(alpha,x)
      real(dp),intent(in)::alpha(:)
      real(dp),intent(out)::x(size(alpha))
      integer::i
      do i=1,size(alpha); x(i)=rgamma_mt(alpha(i)); end do
      x=x/sum(x)
   end subroutine rdirichlet_rng

   subroutine rmvnorm(mean,cov,x,info)
      real(dp),intent(in)::mean(:),cov(:,:)
      real(dp),intent(out)::x(size(mean))
      integer,intent(out)::info
      real(dp),allocatable::l(:,:),z(:)
      integer::i,n
      n=size(mean); allocate(l(n,n),z(n)); call chol_lower(cov,l,info); if(info/=0)return
      do i=1,n; z(i)=rnorm(); end do
      x=mean+matmul(l,z)
   end subroutine rmvnorm

   real(dp) function rtruncnorm(mu,sd,lower,upper) result(x)
      real(dp),intent(in)::mu,sd,lower,upper
      real(dp)::a,b,u
      if(sd<=0.0_dp .or. lower>=upper) then; x=mu; return; end if
      a=normal_cdf((lower-mu)/sd); b=normal_cdf((upper-mu)/sd)
      if(b-a > 1.0e-14_dp) then
         u=a+(b-a)*runif(); u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
         x=mu+sd*normal_quantile(u)
      else
         ! Tail fallback: rejection from the untruncated normal.
         do
            x=rnorm(mu,sd); if(x>=lower .and. x<=upper)exit
         end do
      end if
   end function rtruncnorm

   real(dp) function rinvgauss(mu,lambda) result(x)
      real(dp),intent(in)::mu,lambda
      real(dp)::v,y,cand
      v=rnorm(); y=v*v
      cand=mu+(mu*mu*y)/(2.0_dp*lambda)-mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*y+mu*mu*y*y)
      if(runif()<=mu/(mu+cand)) then; x=cand; else; x=mu*mu/cand; end if
   end function rinvgauss
end module mcmcpack_rng
