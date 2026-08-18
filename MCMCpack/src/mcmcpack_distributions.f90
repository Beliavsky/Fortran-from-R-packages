! SPDX-License-Identifier: GPL-3.0-only
! Distribution utilities translated from MCMCpack R/distn.R.
module mcmcpack_distributions
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_math, only : log_choose
   use mcmcpack_rng, only : rchisq,rnorm,rinvgamma_rng,rdirichlet_rng,runif
   use mcmcpack_linalg, only : chol_lower,inv_spd,logdet_spd,trace_mat
   implicit none
   private
   public :: dinvgamma, rinvgamma, ddirichlet, rdirichlet
   public :: dwish, rwish, diwish, riwish
   public :: noncenhypergeom_pmf, dnoncenhypergeom, rnoncenhypergeom
contains
   elemental real(dp) function dinvgamma(x,shape,scale) result(v)
      real(dp),intent(in)::x,shape
      real(dp),intent(in),optional::scale
      real(dp)::b
      b=1.0_dp; if(present(scale))b=scale
      if(x<=0.0_dp .or. shape<=0.0_dp .or. b<=0.0_dp) then; v=0.0_dp; return; end if
      v=exp(shape*log(b)-log_gamma(shape)-(shape+1.0_dp)*log(x)-b/x)
   end function dinvgamma

   real(dp) function rinvgamma(shape,scale) result(x)
      real(dp),intent(in)::shape
      real(dp),intent(in),optional::scale
      real(dp)::b
      b=1.0_dp; if(present(scale))b=scale
      x=rinvgamma_rng(shape,b)
   end function rinvgamma

   real(dp) function ddirichlet(x,alpha) result(v)
      real(dp),intent(in)::x(:),alpha(:)
      if(size(x)/=size(alpha) .or. any(alpha<=0.0_dp) .or. any(x<=0.0_dp) .or. any(x>1.0_dp) .or. &
         abs(sum(x)-1.0_dp)>100.0_dp*epsilon(1.0_dp)) then
         v=0.0_dp; return
      end if
      v=exp(log_gamma(sum(alpha))-sum(log_gamma(alpha))+sum((alpha-1.0_dp)*log(x)))
   end function ddirichlet

   subroutine rdirichlet(alpha,x)
      real(dp),intent(in)::alpha(:)
      real(dp),intent(out)::x(size(alpha))
      call rdirichlet_rng(alpha,x)
   end subroutine rdirichlet

   pure real(dp) function log_multigamma(a,p) result(v)
      real(dp),intent(in)::a
      integer,intent(in)::p
      integer::j
      v=real(p*(p-1),dp)*0.25_dp*log(pi)
      do j=1,p; v=v+log_gamma(a+0.5_dp*real(1-j,dp)); end do
   end function log_multigamma

   real(dp) function dwish(w,vdf,s) result(d)
      real(dp),intent(in)::w(:,:),s(:,:)
      real(dp),intent(in)::vdf
      real(dp),allocatable::sinv(:,:)
      real(dp)::ldw,lds,lval
      integer::p,info
      p=size(w,1); d=0.0_dp
      if(size(w,2)/=p .or. any(shape(s)/=[p,p]) .or. vdf<real(p,dp))return
      allocate(sinv(p,p)); call logdet_spd(w,ldw,info); if(info/=0)return
      call logdet_spd(s,lds,info); if(info/=0)return
      call inv_spd(s,sinv,info); if(info/=0)return
      lval=0.5_dp*(vdf-real(p+1,dp))*ldw-0.5_dp*trace_mat(matmul(sinv,w)) &
           -0.5_dp*vdf*real(p,dp)*log(2.0_dp)-0.5_dp*vdf*lds-log_multigamma(0.5_dp*vdf,p)
      d=exp(lval)
   end function dwish

   subroutine rwish(vdf,s,w,info)
      real(dp),intent(in)::vdf,s(:,:)
      real(dp),intent(out)::w(size(s,1),size(s,2))
      integer,intent(out)::info
      integer::p,i,j
      real(dp),allocatable::l(:,:),a(:,:),b(:,:)
      p=size(s,1); w=0.0_dp
      if(size(s,2)/=p .or. vdf<real(p,dp))then; info=-1; return; end if
      allocate(l(p,p),a(p,p),b(p,p)); call chol_lower(s,l,info); if(info/=0)return
      a=0.0_dp
      do i=1,p
         a(i,i)=sqrt(rchisq(vdf-real(i-1,dp)))
         do j=1,i-1; a(i,j)=rnorm(); end do
      end do
      b=matmul(a,transpose(l))
      w=matmul(transpose(b),b)
   end subroutine rwish

   real(dp) function diwish(w,vdf,s) result(d)
      real(dp),intent(in)::w(:,:),s(:,:)
      real(dp),intent(in)::vdf
      real(dp),allocatable::winv(:,:)
      real(dp)::ldw,lds,lval
      integer::p,info
      p=size(w,1); d=0.0_dp
      if(size(w,2)/=p .or. any(shape(s)/=[p,p]) .or. vdf<real(p,dp))return
      allocate(winv(p,p)); call logdet_spd(w,ldw,info); if(info/=0)return
      call logdet_spd(s,lds,info); if(info/=0)return
      call inv_spd(w,winv,info); if(info/=0)return
      lval=0.5_dp*vdf*lds-0.5_dp*(vdf+real(p+1,dp))*ldw-0.5_dp*trace_mat(matmul(s,winv)) &
           -0.5_dp*vdf*real(p,dp)*log(2.0_dp)-log_multigamma(0.5_dp*vdf,p)
      d=exp(lval)
   end function diwish

   subroutine riwish(vdf,s,w,info)
      real(dp),intent(in)::vdf,s(:,:)
      real(dp),intent(out)::w(size(s,1),size(s,2))
      integer,intent(out)::info
      real(dp),allocatable::sinv(:,:),q(:,:)
      integer::p
      p=size(s,1); allocate(sinv(p,p),q(p,p))
      call inv_spd(s,sinv,info); if(info/=0)return
      call rwish(vdf,sinv,q,info); if(info/=0)return
      call inv_spd(q,w,info)
   end subroutine riwish

   subroutine noncenhypergeom_pmf(n1,n2,m1,psi,values,prob,info)
      integer,intent(in)::n1,n2,m1
      real(dp),intent(in)::psi
      integer,allocatable,intent(out)::values(:)
      real(dp),allocatable,intent(out)::prob(:)
      integer,intent(out)::info
      integer::ll,uu,i,k
      real(dp),allocatable::lp(:)
      real(dp)::mx
      info=0
      if(n1<0 .or. n2<0 .or. m1<0 .or. m1>n1+n2 .or. psi<=0.0_dp)then; info=1; allocate(values(0),prob(0)); return; end if
      ll=max(0,m1-n2); uu=min(n1,m1); allocate(values(uu-ll+1),prob(uu-ll+1),lp(uu-ll+1))
      do i=ll,uu
         k=i-ll+1; values(k)=i
         lp(k)=log_choose(n1,i)+log_choose(n2,m1-i)+real(i,dp)*log(psi)
      end do
      mx=maxval(lp); prob=exp(lp-mx); prob=prob/sum(prob)
   end subroutine noncenhypergeom_pmf

   real(dp) function dnoncenhypergeom(x,n1,n2,m1,psi) result(v)
      integer,intent(in)::x,n1,n2,m1
      real(dp),intent(in)::psi
      integer,allocatable::vals(:)
      real(dp),allocatable::p(:)
      integer::info,k
      call noncenhypergeom_pmf(n1,n2,m1,psi,vals,p,info); v=0.0_dp
      if(info/=0)return
      do k=1,size(vals); if(vals(k)==x)then; v=p(k); return; end if; end do
   end function dnoncenhypergeom

   integer function rnoncenhypergeom(n1,n2,m1,psi) result(x)
      integer,intent(in)::n1,n2,m1
      real(dp),intent(in)::psi
      integer,allocatable::vals(:)
      real(dp),allocatable::p(:)
      integer::info,k
      real(dp)::u,c
      call noncenhypergeom_pmf(n1,n2,m1,psi,vals,p,info)
      if(info/=0)then; x=0; return; end if
      u=runif(); c=0.0_dp
      do k=1,size(vals); c=c+p(k); if(u<=c)then; x=vals(k); return; end if; end do
      x=vals(size(vals))
   end function rnoncenhypergeom
end module mcmcpack_distributions
