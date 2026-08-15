! Full-family quantile residuals and matrix-first cross-validation for GAMLSS.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_validation_v05
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_types, only : gamlss_result_t,gamlss_control_t
   use gamlss_core, only : fit_gamlss_model
   use gamlss_family_support, only : family_cdf,family_is_discrete
   use gamlss_special, only : normal_quantile
   implicit none
   private
   public :: gamlss_cv_result_t,cross_validate_gamlss,randomized_quantile_residuals_all

   type,public :: gamlss_cv_result_t
      real(dp),allocatable :: fold_loglik(:),point_loglik(:)
      integer,allocatable :: fold_n(:),fold_status(:)
      real(dp),allocatable :: fitted_mu(:),fitted_sigma(:),fitted_nu(:),fitted_tau(:)
      real(dp) :: total_loglik=0.0_dp
      real(dp) :: mean_log_score=huge(1.0_dp)
      integer :: successful_folds=0
      integer :: status=0
   end type gamlss_cv_result_t
contains

   subroutine randomized_quantile_residuals_all(y,result,residuals,status,randomize)
      real(dp),intent(in)::y(:)
      type(gamlss_result_t),intent(in)::result
      real(dp),allocatable,intent(out)::residuals(:)
      integer,intent(out),optional::status
      logical,intent(in),optional::randomize
      real(dp)::a,b,c,d,flo,fhi,p,u
      integer::i,istat
      logical::disc,rnd
      istat=0;rnd=.true.;if(present(randomize))rnd=randomize
      if(.not.allocated(result%mu%fitted).or.size(result%mu%fitted)/=size(y))then
         allocate(residuals(0));if(present(status))status=1;return
      end if
      disc=family_is_discrete(result%family);allocate(residuals(size(y)))
      do i=1,size(y)
         a=result%mu%fitted(i);b=1.0_dp;c=1.0_dp;d=1.0_dp
         if(allocated(result%sigma%fitted))b=result%sigma%fitted(i)
         if(allocated(result%nu%fitted))c=result%nu%fitted(i)
         if(allocated(result%tau%fitted))d=result%tau%fitted(i)
         if(disc)then
            flo=family_cdf(result%family,y(i)-1.0_dp,a,b,c,d)
            fhi=family_cdf(result%family,y(i),a,b,c,d)
            if(flo<0.0_dp.or.fhi<0.0_dp)then;istat=2;exit;end if
            if(rnd)then;call random_number(u);p=flo+u*max(0.0_dp,fhi-flo)
            else;p=0.5_dp*(flo+fhi);end if
         else
            p=family_cdf(result%family,y(i),a,b,c,d)
            if(p<0.0_dp)then;istat=2;exit;end if
         end if
         p=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p));residuals(i)=normal_quantile(p)
      end do
      if(istat/=0)then
         if(allocated(result%residuals).and.size(result%residuals)==size(y))then
            residuals=result%residuals
         else
            residuals=0.0_dp
         end if
      end if
      if(present(status))status=istat
   end subroutine randomized_quantile_residuals_all

   subroutine cross_validate_gamlss(y,x_mu,fold,family,result,x_sigma,x_nu,x_tau,weights, &
      offset_mu,offset_sigma,offset_nu,offset_tau,control)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::fold(:),family
      type(gamlss_cv_result_t),intent(out)::result
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional::control
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),w(:),om(:),os(:),on(:),ot(:)
      real(dp),allocatable::xmt(:,:),xst(:,:),xnt(:,:),xtt(:,:),yt(:),wt(:),omt(:),ost(:),ont(:),ott(:)
      logical,allocatable::train(:),test(:)
      type(gamlss_result_t)::fit
      integer::n,np,nfold,f,i,nt,istat
      real(dp)::a,b,c,d,e1,e2,e3,e4,lp,sw

      n=size(y);np=family_npar(family)
      if(n<4.or.size(x_mu,1)/=n.or.size(fold)/=n.or.np<1.or.np>4.or.any(fold<1))then
         result%status=1;return
      end if
      nfold=maxval(fold);if(nfold<2)then;result%status=2;return;end if
      call full_designs(n,np,x_mu,x_sigma,x_nu,x_tau,xm,xs,xn,xt,istat)
      if(istat/=0)then;result%status=3;return;end if
      allocate(w(n),om(n),os(n),on(n),ot(n));w=1.0_dp;om=0.0_dp;os=0.0_dp;on=0.0_dp;ot=0.0_dp
      if(present(weights))then;if(size(weights)/=n)then;result%status=4;return;end if;w=weights;end if
      if(present(offset_mu))then;if(size(offset_mu)/=n)then;result%status=5;return;end if;om=offset_mu;end if
      if(present(offset_sigma))then;if(size(offset_sigma)/=n)then;result%status=5;return;end if;os=offset_sigma;end if
      if(present(offset_nu))then;if(size(offset_nu)/=n)then;result%status=5;return;end if;on=offset_nu;end if
      if(present(offset_tau))then;if(size(offset_tau)/=n)then;result%status=5;return;end if;ot=offset_tau;end if
      allocate(result%fold_loglik(nfold),result%fold_n(nfold),result%fold_status(nfold),result%point_loglik(n))
      allocate(result%fitted_mu(n),result%fitted_sigma(n),result%fitted_nu(n),result%fitted_tau(n))
      result%fold_loglik=0.0_dp;result%fold_n=0;result%fold_status=1;result%point_loglik=0.0_dp
      result%fitted_mu=0.0_dp;result%fitted_sigma=1.0_dp;result%fitted_nu=1.0_dp;result%fitted_tau=1.0_dp
      result%successful_folds=0;result%total_loglik=0.0_dp;sw=0.0_dp
      allocate(train(n),test(n))
      do f=1,nfold
         test=(fold==f);train=.not.test;nt=count(test)
         result%fold_n(f)=nt
         if(nt==0.or.count(train)<=size(xm,2))then;result%fold_status(f)=10;cycle;end if
         call subset_vector(y,train,yt);call subset_vector(w,train,wt)
         call subset_vector(om,train,omt);call subset_vector(os,train,ost)
         call subset_vector(on,train,ont);call subset_vector(ot,train,ott)
         call subset_matrix(xm,train,xmt);call subset_matrix(xs,train,xst)
         call subset_matrix(xn,train,xnt);call subset_matrix(xt,train,xtt)
         call fit_gamlss_model(yt,xmt,family,fit,x_sigma=xst,x_nu=xnt,x_tau=xtt,weights=wt, &
            offset_mu=omt,offset_sigma=ost,offset_nu=ont,offset_tau=ott,control=control)
         result%fold_status(f)=fit%status
         if(fit%status/=0)cycle
         result%successful_folds=result%successful_folds+1
         do i=1,n
            if(.not.test(i))cycle
            e1=dot_product(xm(i,:),fit%mu%coefficients)+om(i);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
            if(np>=2)e2=dot_product(xs(i,:),fit%sigma%coefficients)+os(i)
            if(np>=3)e3=dot_product(xn(i,:),fit%nu%coefficients)+on(i)
            if(np>=4)e4=dot_product(xt(i,:),fit%tau%coefficients)+ot(i)
            call map_parameters(family,e1,e2,e3,e4,a,b,c,d)
            lp=family_logpdf(family,y(i),a,b,c,d)
            result%point_loglik(i)=w(i)*lp;result%fold_loglik(f)=result%fold_loglik(f)+w(i)*lp
            result%fitted_mu(i)=a;result%fitted_sigma(i)=b;result%fitted_nu(i)=c;result%fitted_tau(i)=d
            sw=sw+w(i)
         end do
         result%total_loglik=result%total_loglik+result%fold_loglik(f)
      end do
      result%status=merge(0,6,result%successful_folds==nfold)
      if(sw>0.0_dp)result%mean_log_score=-result%total_loglik/sw
   end subroutine cross_validate_gamlss

   subroutine full_designs(n,np,xmu,xs_in,xn_in,xt_in,xm,xs,xn,xt,status)
      integer,intent(in)::n,np
      real(dp),intent(in)::xmu(:,:)
      real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
      real(dp),allocatable,intent(out)::xm(:,:),xs(:,:),xn(:,:),xt(:,:)
      integer,intent(out)::status
      status=0;xm=xmu
      if(present(xs_in))then;if(size(xs_in,1)/=n)then;status=1;return;end if;xs=xs_in
      else;allocate(xs(n,1));xs=1.0_dp;end if
      if(present(xn_in))then;if(size(xn_in,1)/=n)then;status=2;return;end if;xn=xn_in
      else;allocate(xn(n,1));xn=1.0_dp;end if
      if(present(xt_in))then;if(size(xt_in,1)/=n)then;status=3;return;end if;xt=xt_in
      else;allocate(xt(n,1));xt=1.0_dp;end if
   end subroutine full_designs

   subroutine subset_vector(a,mask,b)
      real(dp),intent(in)::a(:)
      logical,intent(in)::mask(:)
      real(dp),allocatable,intent(out)::b(:)
      integer::i,k
      allocate(b(count(mask)));k=0
      do i=1,size(a);if(mask(i))then;k=k+1;b(k)=a(i);end if;end do
   end subroutine subset_vector

   subroutine subset_matrix(a,mask,b)
      real(dp),intent(in)::a(:,:)
      logical,intent(in)::mask(:)
      real(dp),allocatable,intent(out)::b(:,:)
      integer::i,k
      allocate(b(count(mask),size(a,2)));k=0
      do i=1,size(a,1);if(mask(i))then;k=k+1;b(k,:)=a(i,:);end if;end do
   end subroutine subset_matrix
end module gamlss_validation_v05
