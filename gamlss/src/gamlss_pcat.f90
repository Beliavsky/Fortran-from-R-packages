! Penalized categorical level fusion corresponding to upstream pcat().
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_pcat
   use gamlss_kinds, only : dp
   use gamlss_linalg, only : penalized_weighted_least_squares, invert_matrix
   use gamlss_smoothers, only : all_pair_difference_matrix
   implicit none
   private
   public :: pcat_result_t, pcat_design, fit_pcat, pcat_fused_groups

   type,public :: pcat_result_t
      real(dp),allocatable :: coefficients(:),fitted(:),difference_weights(:)
      real(dp) :: lambda=1.0_dp
      real(dp) :: edf=0.0_dp
      integer :: iterations=0
      integer :: status=0
      logical :: converged=.false.
   end type pcat_result_t
contains

   subroutine pcat_design(category,nlevels,x,status)
      integer,intent(in)::category(:)
      integer,intent(in),optional::nlevels
      real(dp),allocatable,intent(out)::x(:,:)
      integer,intent(out),optional::status
      integer::lev,i
      lev=maxval(category)
      if(present(nlevels))lev=nlevels
      if(size(category)==0.or.lev<1.or.minval(category)<1.or.maxval(category)>lev)then
         allocate(x(0,0));if(present(status))status=1;return
      end if
      allocate(x(size(category),lev));x=0.0_dp
      do i=1,size(category);x(i,category(i))=1.0_dp;end do
      if(present(status))status=0
   end subroutine pcat_design

   subroutine fit_pcat(x,y,w,result,lambda,estimate_lambda,lp,kappa,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      type(pcat_result_t),intent(out)::result
      real(dp),intent(in),optional::lambda,lp,kappa,tol
      logical,intent(in),optional::estimate_lambda
      integer,intent(in),optional::max_iter
      real(dp),allocatable::d(:,:),omega(:),pen(:,:),beta(:),cov(:,:),u(:),oldbeta(:)
      real(dp)::lam,pnorm,kap,crit,sig2,tau2,newlam,delta
      integer::p,n,it,nit,istat,j
      logical::auto
      n=size(x,1);p=size(x,2)
      if(size(y)/=n.or.size(w)/=n.or.p<1.or.any(w<0.0_dp))then;result%status=1;return;end if
      lam=1.0e-3_dp;if(present(lambda))lam=max(1.0e-7_dp,lambda)
      auto=.true.;if(present(estimate_lambda))auto=estimate_lambda
      pnorm=0.0_dp;if(present(lp))pnorm=max(0.0_dp,min(2.0_dp,lp))
      kap=1.0e-5_dp;if(present(kappa))kap=max(1.0e-12_dp,kappa)
      nit=100;if(present(max_iter))nit=max_iter
      crit=1.0e-4_dp;if(present(tol))crit=tol
      call all_pair_difference_matrix(p,d);allocate(omega(size(d,1)));omega=1.0_dp
      allocate(beta(p));beta=0.0_dp
      do it=1,nit
         oldbeta=beta
         call adaptive_penalty(d,omega,lam,pen)
         do j=1,p;pen(j,j)=pen(j,j)+1.0e-6_dp;end do
         call penalized_weighted_least_squares(x,y,w,pen,beta,cov,istat)
         if(istat/=0)then;result%status=2;return;end if
         u=matmul(d,beta)
         omega=1.0_dp/(abs(u)**(2.0_dp-pnorm)+kap*kap)
         result%edf=pcat_edf(x,w,pen)
         if(auto)then
            sig2=sum(w*(y-matmul(x,beta))**2)/max(1.0_dp,sum(merge(1.0_dp,0.0_dp,w>0.0_dp))-result%edf)
            tau2=sum(omega*u*u)/max(1.0e-6_dp,result%edf)
            tau2=max(tau2,1.0e-7_dp)
            newlam=min(1.0e7_dp,max(1.0e-7_dp,sig2/tau2))
            lam=sqrt(lam*newlam)
         end if
         delta=maxval(abs(beta-oldbeta));result%iterations=it
         if(delta<crit)then;result%converged=.true.;exit;end if
      end do
      result%coefficients=beta;allocate(result%fitted(n));result%fitted=matmul(x,beta)
      result%difference_weights=omega;result%lambda=lam;result%status=0
   end subroutine fit_pcat

   subroutine adaptive_penalty(d,omega,lambda,pen)
      real(dp),intent(in)::d(:,:),omega(:),lambda
      real(dp),allocatable,intent(out)::pen(:,:)
      integer::i,j,k,p
      p=size(d,2);allocate(pen(p,p));pen=0.0_dp
      do i=1,size(d,1)
         do j=1,p;do k=1,p
            pen(j,k)=pen(j,k)+lambda*omega(i)*d(i,j)*d(i,k)
         end do;end do
      end do
   end subroutine adaptive_penalty

   real(dp) function pcat_edf(x,w,pen) result(edf)
      real(dp),intent(in)::x(:,:),w(:),pen(:,:)
      real(dp),allocatable::xtwx(:,:),a(:,:),ainv(:,:),h(:,:)
      integer::i,j,k,p,istat
      p=size(x,2);allocate(xtwx(p,p));xtwx=0.0_dp
      do i=1,size(x,1);do j=1,p;do k=1,p
         xtwx(j,k)=xtwx(j,k)+w(i)*x(i,j)*x(i,k)
      end do;end do;end do
      a=xtwx+pen;call invert_matrix(a,ainv,istat)
      if(istat/=0)then;edf=real(p,dp);return;end if
      h=matmul(ainv,xtwx);edf=0.0_dp;do i=1,p;edf=edf+h(i,i);end do
   end function pcat_edf

   subroutine pcat_fused_groups(result,tolerance,groups,ngroups)
      type(pcat_result_t),intent(in)::result
      real(dp),intent(in),optional::tolerance
      integer,allocatable,intent(out)::groups(:)
      integer,intent(out)::ngroups
      real(dp)::tol
      integer::i,j,g,p
      logical::matched
      tol=1.0e-4_dp;if(present(tolerance))tol=max(0.0_dp,tolerance)
      if(.not.allocated(result%coefficients))then;allocate(groups(0));ngroups=0;return;end if
      p=size(result%coefficients);allocate(groups(p));groups=0;g=0
      do i=1,p
         matched=.false.
         do j=1,i-1
            if(abs(result%coefficients(i)-result%coefficients(j))<=tol)then
               groups(i)=groups(j);matched=.true.;exit
            end if
         end do
         if(.not.matched)then;g=g+1;groups(i)=g;end if
      end do
      ngroups=g
   end subroutine pcat_fused_groups

end module gamlss_pcat
