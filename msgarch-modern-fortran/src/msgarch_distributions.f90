! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use msgarch_kinds, only : dp, pi
   use msgarch_rng, only : random_uniform, random_normal, random_gamma, random_student_t
   use msgarch_special, only : normal_cdf, normal_quantile, student_t_cdf, student_t_quantile, &
      regularized_gamma_p, gamma_quantile
   implicit none
   private
   public :: distribution_valid, innovation_logpdf, innovation_pdf, innovation_cdf
   public :: innovation_quantile, random_innovation, distribution_moments
contains
   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower=text
      do i=1,len(text)
         code=iachar(text(i:i))
         if(code>=iachar('A').and.code<=iachar('Z'))lower(i:i)=achar(code+32)
      end do
   end function lower_string

   pure function distribution_valid(distribution,shape,skew) result(ok)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: shape,skew
      logical :: ok
      character(len=:), allocatable :: dist
      dist=trim(adjustl(lower_string(distribution)))
      select case(dist)
      case('norm'); ok=.true.
      case('std'); ok=shape>2.0_dp
      case('ged'); ok=shape>0.0_dp
      case('snorm'); ok=skew>0.0_dp
      case('sstd'); ok=shape>2.0_dp.and.skew>0.0_dp
      case('sged'); ok=shape>0.0_dp.and.skew>0.0_dp
      case default; ok=.false.
      end select
   end function distribution_valid

   pure function base_abs_moment(base,order,shape) result(moment)
      character(len=*), intent(in) :: base
      real(dp), intent(in) :: order,shape
      real(dp) :: moment,lambda
      select case(trim(base))
      case('norm')
         moment=2.0_dp**(0.5_dp*order)*exp(log_gamma(0.5_dp*(order+1.0_dp))-0.5_dp*log(pi))
      case('std')
         if(shape<=max(2.0_dp,order))then
            moment=huge(1.0_dp)
         else
            moment=(shape-2.0_dp)**(0.5_dp*order)*exp(log_gamma(0.5_dp*(order+1.0_dp))+ &
               log_gamma(0.5_dp*(shape-order))-0.5_dp*log(pi)-log_gamma(0.5_dp*shape))
         end if
      case('ged')
         lambda=sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         moment=lambda**order*exp(log_gamma((order+1.0_dp)/shape)-log_gamma(1.0_dp/shape))
      case default
         moment=huge(1.0_dp)
      end select
   end function base_abs_moment

   pure subroutine fs_location_scale(base,shape,skew,mu,sigma)
      character(len=*), intent(in) :: base
      real(dp), intent(in) :: shape,skew
      real(dp), intent(out) :: mu,sigma
      real(dp) :: m1,second
      m1=base_abs_moment(base,1.0_dp,shape)
      mu=m1*(skew-1.0_dp/skew)
      second=(skew**3+skew**(-3))/(skew+1.0_dp/skew)
      sigma=sqrt(max(second-mu*mu,1.0e-14_dp))
   end subroutine fs_location_scale

   pure function base_logpdf(x,base,shape) result(value)
      real(dp), intent(in) :: x,shape
      character(len=*), intent(in) :: base
      real(dp) :: value,lambda
      select case(trim(base))
      case('norm')
         value=-0.5_dp*(log(2.0_dp*pi)+x*x)
      case('std')
         value=log_gamma(0.5_dp*(shape+1.0_dp))-log_gamma(0.5_dp*shape)- &
            0.5_dp*log(pi*(shape-2.0_dp))-0.5_dp*(shape+1.0_dp)*log(1.0_dp+x*x/(shape-2.0_dp))
      case('ged')
         lambda=2.0_dp**(-1.0_dp/shape)*sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         value=log(shape)-log(lambda)-(1.0_dp+1.0_dp/shape)*log(2.0_dp)- &
            log_gamma(1.0_dp/shape)-0.5_dp*(abs(x)/lambda)**shape
      case default
         value=-huge(1.0_dp)
      end select
   end function base_logpdf

   function base_cdf(x,base,shape) result(p)
      real(dp), intent(in) :: x,shape
      character(len=*), intent(in) :: base
      real(dp) :: p,lambda,q
      select case(trim(base))
      case('norm')
         p=normal_cdf(x)
      case('std')
         p=student_t_cdf(x*sqrt(shape/(shape-2.0_dp)),shape)
      case('ged')
         lambda=2.0_dp**(-1.0_dp/shape)*sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         q=regularized_gamma_p(1.0_dp/shape,0.5_dp*(abs(x)/lambda)**shape)
         if(x<0.0_dp)then;p=0.5_dp*(1.0_dp-q);else;p=0.5_dp*(1.0_dp+q);end if
      case default
         p=0.0_dp
      end select
      p=max(0.0_dp,min(1.0_dp,p))
   end function base_cdf

   function base_quantile(p,base,shape) result(x)
      real(dp), intent(in) :: p,shape
      character(len=*), intent(in) :: base
      real(dp) :: x,lambda,q
      select case(trim(base))
      case('norm')
         x=normal_quantile(p)
      case('std')
         x=student_t_quantile(p,shape)*sqrt((shape-2.0_dp)/shape)
      case('ged')
         lambda=2.0_dp**(-1.0_dp/shape)*sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         if(p<0.5_dp)then
            q=gamma_quantile(1.0_dp-2.0_dp*p,1.0_dp/shape)
            x=-lambda*(2.0_dp*q)**(1.0_dp/shape)
         else
            q=gamma_quantile(2.0_dp*p-1.0_dp,1.0_dp/shape)
            x=lambda*(2.0_dp*q)**(1.0_dp/shape)
         end if
      case default
         x=0.0_dp
      end select
   end function base_quantile

   pure function innovation_logpdf(z,distribution,shape,skew) result(value)
      real(dp), intent(in) :: z,shape,skew
      character(len=*), intent(in) :: distribution
      real(dp) :: value,mu,sigma,raw,transformed
      character(len=:), allocatable :: dist,base
      dist=trim(adjustl(lower_string(distribution)))
      if(.not.distribution_valid(dist,shape,skew))then;value=-huge(1.0_dp);return;end if
      select case(dist)
      case('norm','std','ged')
         value=base_logpdf(z,dist,shape)
      case('snorm','sstd','sged')
         base=dist(2:)
         call fs_location_scale(base,shape,skew,mu,sigma)
         raw=mu+sigma*z
         if(raw>=0.0_dp)then;transformed=raw/skew;else;transformed=raw*skew;end if
         value=log(sigma)+log(2.0_dp)-log(skew+1.0_dp/skew)+base_logpdf(transformed,base,shape)
      case default
         value=-huge(1.0_dp)
      end select
   end function innovation_logpdf

   pure function innovation_pdf(z,distribution,shape,skew) result(value)
      real(dp), intent(in) :: z,shape,skew
      character(len=*), intent(in) :: distribution
      real(dp) :: value,lp
      lp=innovation_logpdf(z,distribution,shape,skew)
      if(lp<log(tiny(1.0_dp)))then;value=0.0_dp;else;value=exp(lp);end if
   end function innovation_pdf

   function innovation_cdf(z,distribution,shape,skew) result(p)
      real(dp), intent(in) :: z,shape,skew
      character(len=*), intent(in) :: distribution
      real(dp) :: p,mu,sigma,raw,num,cutoff
      character(len=:), allocatable :: dist,base
      dist=trim(adjustl(lower_string(distribution)))
      select case(dist)
      case('norm','std','ged')
         p=base_cdf(z,dist,shape)
      case('snorm','sstd','sged')
         base=dist(2:)
         call fs_location_scale(base,shape,skew,mu,sigma)
         raw=sigma*z+mu
         num=1.0_dp/(skew+1.0_dp/skew)
         cutoff=-mu/sigma
         if(z<cutoff)then
            p=2.0_dp/skew*num*base_cdf(raw*skew,base,shape)
         else
            p=2.0_dp*num*(skew*base_cdf(raw/skew,base,shape)+1.0_dp/skew)-1.0_dp
         end if
      case default
         p=0.0_dp
      end select
      p=max(0.0_dp,min(1.0_dp,p))
   end function innovation_cdf

   function innovation_quantile(p,distribution,shape,skew) result(z)
      real(dp), intent(in) :: p,shape,skew
      character(len=*), intent(in) :: distribution
      real(dp) :: z,mu,sigma,pcut,u,xi2
      character(len=:), allocatable :: dist,base
      dist=trim(adjustl(lower_string(distribution)))
      select case(dist)
      case('norm','std','ged')
         z=base_quantile(p,dist,shape)
      case('snorm','sstd','sged')
         base=dist(2:)
         call fs_location_scale(base,shape,skew,mu,sigma)
         xi2=skew*skew
         pcut=1.0_dp/(1.0_dp+xi2)
         if(p<pcut)then
            u=0.5_dp*p*(xi2+1.0_dp)
            z=(base_quantile(u,base,shape)/skew-mu)/sigma
         else
            u=0.5_dp*p*(1.0_dp+1.0_dp/xi2)-0.5_dp/xi2+0.5_dp
            z=(base_quantile(u,base,shape)*skew-mu)/sigma
         end if
      case default
         z=0.0_dp
      end select
   end function innovation_quantile

   function random_innovation(distribution,shape,skew) result(z)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: shape,skew
      real(dp) :: z,lambda,magnitude,raw,mu,sigma,pplus
      character(len=:), allocatable :: dist,base
      magnitude=0.0_dp
      dist=trim(adjustl(lower_string(distribution)))
      select case(dist)
      case('norm');z=random_normal()
      case('std');z=sqrt((shape-2.0_dp)/shape)*random_student_t(shape)
      case('ged')
         lambda=2.0_dp**(-1.0_dp/shape)*sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         magnitude=lambda*(2.0_dp*random_gamma(1.0_dp/shape))**(1.0_dp/shape)
         if(random_uniform()<0.5_dp)magnitude=-magnitude
         z=magnitude
      case('snorm','sstd','sged')
         base=dist(2:)
         select case(base)
         case('norm');magnitude=abs(random_normal())
         case('std');magnitude=abs(sqrt((shape-2.0_dp)/shape)*random_student_t(shape))
         case('ged')
            lambda=2.0_dp**(-1.0_dp/shape)*sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
            magnitude=lambda*(2.0_dp*random_gamma(1.0_dp/shape))**(1.0_dp/shape)
         end select
         pplus=skew*skew/(1.0_dp+skew*skew)
         if(random_uniform()<pplus)then;raw=skew*magnitude;else;raw=-magnitude/skew;end if
         call fs_location_scale(base,shape,skew,mu,sigma)
         z=(raw-mu)/sigma
      case default
         z=0.0_dp
      end select
   end function random_innovation

   subroutine distribution_moments(distribution,shape,skew,eabs,ezineg,ez2ineg)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: shape,skew
      real(dp), intent(out) :: eabs,ezineg,ez2ineg
      integer, parameter :: ngrid=8000
      real(dp) :: bound,h,z,w,pdf
      integer :: i
      character(len=:), allocatable :: dist
      dist=trim(adjustl(lower_string(distribution)))
      if(dist=='norm'.or.dist=='std'.or.dist=='ged')then
         eabs=base_abs_moment(dist,1.0_dp,shape)
         ezineg=-0.5_dp*eabs
         ez2ineg=0.5_dp
         return
      end if
      bound=15.0_dp
      if(index(dist,'std')>0)bound=60.0_dp
      h=2.0_dp*bound/real(ngrid,dp)
      eabs=0.0_dp;ezineg=0.0_dp;ez2ineg=0.0_dp
      do i=0,ngrid
         z=-bound+h*real(i,dp)
         if(i==0.or.i==ngrid)then;w=1.0_dp;else if(mod(i,2)==0)then;w=2.0_dp;else;w=4.0_dp;end if
         pdf=innovation_pdf(z,dist,shape,skew)
         eabs=eabs+w*abs(z)*pdf
         if(z<0.0_dp)then
            ezineg=ezineg+w*z*pdf
            ez2ineg=ez2ineg+w*z*z*pdf
         end if
      end do
      eabs=eabs*h/3.0_dp;ezineg=ezineg*h/3.0_dp;ez2ineg=ez2ineg*h/3.0_dp
      if(.not.ieee_is_finite(eabs))eabs=huge(1.0_dp)
   end subroutine distribution_moments
end module msgarch_distributions
