! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_scale
   use robustbase_kinds, only: dp
   use robustbase_sort, only: median, quantile_type7, sort_real, weighted_high_median
   implicit none
   private
   public :: mad_scale, iqr_scale, qn_scale, sn_scale, huber_location, huberize_vector, &
             tau_huber, trimmed_mean_abs_dev, scale_tau2
contains
   function mad_scale(x, center, constant) result(s)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: center, constant
      real(dp) :: s, c, k
      real(dp), allocatable :: d(:)
      c=median(x); if(present(center)) c=center
      k=1.482602218505602_dp; if(present(constant)) k=constant
      d=abs(x-c); s=k*median(d)
   end function mad_scale

   function iqr_scale(x) result(s)
      real(dp), intent(in) :: x(:)
      real(dp) :: s
      s=(quantile_type7(x,0.75_dp)-quantile_type7(x,0.25_dp))/1.348979500392164_dp
   end function iqr_scale

   function qn_finite(n) result(c)
      integer,intent(in)::n
      real(dp)::c
      real(dp),parameter::small(11)=[0.399356_dp,0.99365_dp,0.51321_dp,0.84401_dp,0.61220_dp,0.85877_dp,0.66993_dp,0.87344_dp,0.72014_dp,0.88906_dp,0.75743_dp]
      if(n<=12) then
         c=small(n-1)
      else if(mod(n,2)==1) then
         c=1.0_dp/((1.60188_dp+(-2.1284_dp-5.172_dp/real(n,dp))/real(n,dp))/real(n,dp)+1.0_dp)
      else
         c=1.0_dp/((3.67561_dp+(1.9654_dp+(6.987_dp-77.0_dp/real(n,dp))/real(n,dp))/real(n,dp))/real(n,dp)+1.0_dp)
      end if
   end function qn_finite

   function qn_scale(x, finite_correction, constant) result(s)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: finite_correction
      real(dp), intent(in), optional :: constant
      real(dp) :: s, cc
      real(dp), allocatable :: d(:)
      integer :: n,i,j,k,h,kk
      logical :: fc
      n=size(x)
      if(n<=1) then; s=0.0_dp; return; end if
      allocate(d(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; d(k)=abs(x(i)-x(j))
         end do
      end do
      call sort_real(d)
      h=n/2+1; kk=h*(h-1)/2
      cc=2.21914_dp; if(present(constant)) cc=constant
      s=cc*d(kk)
      fc=.true.; if(present(finite_correction)) fc=finite_correction
      if(fc) s=s*qn_finite(n)
   end function qn_scale

   function sn_scale(x, finite_correction, constant) result(s)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: finite_correction
      real(dp), intent(in), optional :: constant
      real(dp) :: s,cc,fn
      real(dp), allocatable :: a(:),d(:)
      integer :: n,i
      logical :: fc
      real(dp),parameter::small(8)=[0.743_dp,1.851_dp,0.954_dp,1.351_dp,0.993_dp,1.198_dp,1.005_dp,1.131_dp]
      n=size(x)
      if(n<=1) then; s=0.0_dp; return; end if
      allocate(a(n),d(n))
      do i=1,n
         d=abs(x-x(i)); a(i)=median(d)
      end do
      cc=1.1926_dp; if(present(constant)) cc=constant
      s=cc*median(a)
      fc=.true.; if(present(finite_correction)) fc=finite_correction
      if(fc) then
         if(n<=9) then; fn=small(n-1)
         else if(mod(n,2)==1) then; fn=real(n,dp)/(real(n,dp)-0.9_dp)
         else; fn=1.0_dp
         end if
         s=s*fn
      end if
   end function sn_scale

   subroutine huber_location(x, mu, scale, iterations, k, weights, tol, standard_error)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::mu,scale
      integer,intent(out)::iterations
      real(dp),intent(in),optional::k,weights(:),tol
      real(dp),intent(out),optional::standard_error
      real(dp)::kk,tt,mu1,sumw
      real(dp),allocatable::w(:),y(:)
      integer::n
      n=size(x); kk=1.5_dp; if(present(k))kk=k; tt=1.0e-6_dp;if(present(tol))tt=tol
      allocate(w(n),y(n)); w=1.0_dp; if(present(weights)) then
         if(size(weights)/=n .or. any(weights<0.0_dp)) error stop "huber_location: invalid weights"
         w=weights
      end if
      sumw=sum(w)
      if(sumw<=0.0_dp) then; mu=0.0_dp;scale=0.0_dp;iterations=0;if(present(standard_error))standard_error=0.0_dp;return;end if
      mu=weighted_high_median(x,w); scale=weighted_high_median(abs(x-mu),w)*1.482602218505602_dp
      iterations=0
      if(scale>0.0_dp) then
         do
            iterations=iterations+1
            y=min(max(x,mu-kk*scale),mu+kk*scale)
            mu1=sum(w*y)/sumw
            if(abs(mu1-mu)<tt*scale .or. iterations>=1000) exit
            mu=mu1
         end do
         mu=mu1
      end if
      if(present(standard_error)) then
         if(scale>0.0_dp) then
            standard_error=scale*sqrt(tau_huber(x,mu,scale,kk)/real(n,dp))
         else
            standard_error=0.0_dp
         end if
      end if
   end subroutine huber_location

   function tau_huber(x,mu,scale,k) result(tau)
      real(dp),intent(in)::x(:),mu,scale,k
      real(dp)::tau
      real(dp),allocatable::r(:),psi(:)
      integer::np
      if(scale<=0.0_dp) then;tau=0.0_dp;return;end if
      r=(x-mu)/scale; psi=max(-k,min(k,r)); np=count(abs(r)<=k)
      if(np==0) then;tau=huge(1.0_dp);else;tau=real(size(x),dp)*sum(psi*psi)/real(np*np,dp);end if
   end function tau_huber

   function trimmed_mean_abs_dev(x,center,trim) result(s)
      real(dp),intent(in)::x(:),center,trim
      real(dp)::s
      real(dp),allocatable::d(:)
      integer::n,k
      if(trim<0.0_dp .or. trim>0.5_dp) error stop "trimmed_mean_abs_dev: invalid trim"
      d=abs(x-center);call sort_real(d);n=size(d);k=int(trim*real(n,dp))
      if(n-2*k<=0) then;s=0.0_dp;else;s=sum(d(k+1:n-k))/real(n-2*k,dp);end if
   end function trimmed_mean_abs_dev

   subroutine huberize_vector(x,y,center,cutoff,k,scale_out)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      real(dp),intent(in),optional::center,cutoff,k
      real(dp),intent(out),optional::scale_out
      real(dp)::m,c,kk,s
      integer::it
      if(size(y)/=size(x)) error stop "huberize_vector: size mismatch"
      kk=1.5_dp;if(present(k))kk=k
      if(present(center)) then;m=center;else;call huber_location(x,m,s,it,k=kk);end if
      s=qn_scale(x)
      c=kk;if(present(cutoff))c=cutoff
      y=min(max(x,m-c*s),m+c*s)
      if(present(scale_out))scale_out=s
   end subroutine huberize_vector

   function scale_tau2(x, center) result(s)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::center
      real(dp)::s,m,s0,u,w
      integer::i,n
      n=size(x);m=median(x);if(present(center))m=center
      s0=mad_scale(x,m)
      if(s0<=0.0_dp) then;s=0.0_dp;return;end if
      s=0.0_dp
      do i=1,n
         u=(x(i)-m)/(4.5_dp*s0)
         if(abs(u)<1.0_dp) then
            w=(1.0_dp-u*u)**2
            s=s+(x(i)-m)**2*w*w
         end if
      end do
      s=sqrt(s/real(max(1,n),dp))/sqrt(0.437_dp)
   end function scale_tau2
end module robustbase_scale
