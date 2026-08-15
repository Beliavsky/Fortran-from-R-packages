! LMS/centile computational layer using BCCG/BCT/BCPE and P-splines.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_lms
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_BCCG,GAMLSS_BCT,GAMLSS_BCPE
   use gamlss_boxcox, only : qBCCG,qBCT,qBCPE
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model,predict_gamlss_parameters
   use gamlss_smoothers, only : p_spline_spec_t,fit_p_spline_basis,predict_p_spline_basis
   implicit none
   private
   public :: lms_result_t, fit_lms, predict_lms, lms_centiles

   type, public :: lms_result_t
      type(gamlss_result_t) :: fit
      type(p_spline_spec_t) :: mu_spline,sigma_spline,nu_spline
      integer :: family=GAMLSS_BCCG
   end type lms_result_t

contains

   subroutine fit_lms(y,x,result,family,df_mu,df_sigma,df_nu,lambda_mu,lambda_sigma,lambda_nu,control,status)
      real(dp),intent(in)::y(:),x(:)
      type(lms_result_t),intent(out)::result
      integer,intent(in),optional::family,df_mu,df_sigma,df_nu
      real(dp),intent(in),optional::lambda_mu,lambda_sigma,lambda_nu
      type(gamlss_control_t),intent(in),optional::control
      integer,intent(out),optional::status
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),pt(:,:)
      integer::fam,dm,ds,dn,istat
      real(dp)::lm,ls,ln
      fam=GAMLSS_BCCG;if(present(family))fam=family
      if(fam/=GAMLSS_BCCG.and.fam/=GAMLSS_BCT.and.fam/=GAMLSS_BCPE)then
         if(present(status))status=1;return
      end if
      dm=10;if(present(df_mu))dm=df_mu;ds=6;if(present(df_sigma))ds=df_sigma;dn=6;if(present(df_nu))dn=df_nu
      lm=10.0_dp;if(present(lambda_mu))lm=lambda_mu
      ls=10.0_dp;if(present(lambda_sigma))ls=lambda_sigma
      ln=10.0_dp;if(present(lambda_nu))ln=lambda_nu
      call fit_p_spline_basis(x,result%mu_spline,xm,df=dm,status=istat);if(istat/=0)goto 900
      call fit_p_spline_basis(x,result%sigma_spline,xs,df=ds,status=istat);if(istat/=0)goto 900
      call fit_p_spline_basis(x,result%nu_spline,xn,df=dn,status=istat);if(istat/=0)goto 900
      result%family=fam
      if(fam==GAMLSS_BCCG)then
         call fit_gamlss_model(y,xm,fam,result%fit,GAMLSS_METHOD_RS,x_sigma=xs,x_nu=xn, &
            penalty_mu=result%mu_spline%penalty,penalty_sigma=result%sigma_spline%penalty, &
            penalty_nu=result%nu_spline%penalty,lambda_mu=lm,lambda_sigma=ls,lambda_nu=ln,control=control)
      else
         allocate(xt(size(y),1),pt(1,1));xt=1.0_dp;pt=0.0_dp
         call fit_gamlss_model(y,xm,fam,result%fit,GAMLSS_METHOD_RS,x_sigma=xs,x_nu=xn,x_tau=xt, &
            penalty_mu=result%mu_spline%penalty,penalty_sigma=result%sigma_spline%penalty, &
            penalty_nu=result%nu_spline%penalty,penalty_tau=pt,lambda_mu=lm,lambda_sigma=ls, &
            lambda_nu=ln,lambda_tau=0.0_dp,control=control)
      end if
      istat=result%fit%status
900   if(present(status))status=istat
   end subroutine fit_lms

   subroutine predict_lms(model,x,mu,sigma,nu,tau,status)
      type(lms_result_t),intent(in)::model
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::mu(:),sigma(:),nu(:)
      real(dp),allocatable,intent(out),optional::tau(:)
      integer,intent(out),optional::status
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:)
      integer::istat
      call predict_p_spline_basis(x,model%mu_spline,xm,istat);if(istat/=0)goto 900
      call predict_p_spline_basis(x,model%sigma_spline,xs,istat);if(istat/=0)goto 900
      call predict_p_spline_basis(x,model%nu_spline,xn,istat);if(istat/=0)goto 900
      if(model%family==GAMLSS_BCCG)then
         call predict_gamlss_parameters(model%family,model%fit,xm,mu,xs,sigma,xn,nu,status=istat)
      else
         allocate(xt(size(x),1));xt=1.0_dp
         if(present(tau))then
            call predict_gamlss_parameters(model%family,model%fit,xm,mu,xs,sigma,xn,nu,xt,tau,status=istat)
         else
            block
               real(dp),allocatable::tt(:)
               call predict_gamlss_parameters(model%family,model%fit,xm,mu,xs,sigma,xn,nu,xt,tt,status=istat)
            end block
         end if
      end if
900   if(present(status))status=istat
   end subroutine predict_lms

   subroutine lms_centiles(model,x,prob,centiles,status)
      type(lms_result_t),intent(in)::model
      real(dp),intent(in)::x(:),prob(:)
      real(dp),allocatable,intent(out)::centiles(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::mu(:),sigma(:),nu(:),tau(:)
      integer::i,j,istat
      call predict_lms(model,x,mu,sigma,nu,tau,istat)
      if(istat/=0)then;if(present(status))status=istat;allocate(centiles(0,0));return;end if
      allocate(centiles(size(x),size(prob)))
      do i=1,size(x);do j=1,size(prob)
         select case(model%family)
         case(GAMLSS_BCCG);centiles(i,j)=qBCCG(prob(j),mu(i),sigma(i),nu(i))
         case(GAMLSS_BCT);centiles(i,j)=qBCT(prob(j),mu(i),sigma(i),nu(i),tau(i))
         case(GAMLSS_BCPE);centiles(i,j)=qBCPE(prob(j),mu(i),sigma(i),nu(i),tau(i))
         end select
      end do;end do
      if(present(status))status=0
   end subroutine lms_centiles

end module gamlss_lms
