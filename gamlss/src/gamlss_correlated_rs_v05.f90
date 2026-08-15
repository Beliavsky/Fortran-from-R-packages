! Correlated working-response RS iterations for non-Gaussian GAMLSS families.
! The local Fisher/IRLS working covariance is D(1/sqrt(w))*V*R*V*D(1/sqrt(w)),
! where R and V are supplied by the nlme correlation/variance structures.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_correlated_rs_v05
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_types, only : gamlss_result_t,gamlss_parameter_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use nlme_types, only : correlation_spec,variance_spec,gls_result,nlme_control,VAR_FIXED
   use nlme_gls, only : fit_gls,initialize_correlation,initialize_variance
   use nlme_variance, only : variance_sd
   implicit none
   private
   public :: correlated_rs_result_t,fit_gamlss_correlated_rs

   type :: rs_block_t
      real(dp),allocatable :: x(:,:),offset(:),beta(:),eta(:),covariance(:,:)
   end type rs_block_t

   type,public :: correlated_rs_result_t
      type(gamlss_result_t) :: model
      type(correlation_spec) :: correlation
      type(variance_spec) :: variance
      real(dp),allocatable :: correlation_parameters(:)
      real(dp),allocatable :: variance_sd(:)
      integer :: iterations=0
      integer :: status=0
      logical :: converged=.false.
   end type correlated_rs_result_t
contains

   subroutine fit_gamlss_correlated_rs(y,x_mu,family,result,correlation,variance,x_sigma,x_nu,x_tau, &
      weights,offset_mu,offset_sigma,offset_nu,offset_tau,time,group,var_covariate,var_group,coordinates, &
      control,nlme_control_in,max_outer,tolerance)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::family
      type(correlated_rs_result_t),intent(out)::result
      type(correlation_spec),intent(in),optional::correlation
      type(variance_spec),intent(in),optional::variance
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      real(dp),intent(in),optional::time(:),var_covariate(:),coordinates(:,:)
      integer,intent(in),optional::group(:),var_group(:)
      type(gamlss_control_t),intent(in),optional::control
      type(nlme_control),intent(in),optional::nlme_control_in
      integer,intent(in),optional::max_outer
      real(dp),intent(in),optional::tolerance

      type(rs_block_t)::b(4)
      type(gamlss_result_t)::initial
      type(gamlss_control_t)::gctl
      type(nlme_control)::nctl
      type(correlation_spec)::corr,corr_use
      type(variance_spec)::vbase,vwork
      type(gls_result)::gr
      real(dp),allocatable::w(:),tv(:),vc(:),base_sd(:),ww(:),z(:),var_work(:)
      real(dp),allocatable::oldbeta(:),oldeta(:)
      integer,allocatable::gv(:),vg(:)
      integer::n,np,j,it,nouter,istat,k
      real(dp)::crit,olddev,newdev,step,outer_change
      logical::outer_converged

      n=size(y);np=family_npar(family)
      if(n<=0.or.np<1.or.np>4.or.size(x_mu,1)/=n)then;result%status=1;return;end if
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,b,istat)
      if(istat/=0)then;result%status=2;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,b,istat)
      if(istat/=0)then;result%status=3;return;end if
      allocate(w(n));w=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;result%status=4;return;end if
         w=weights
      end if
      allocate(tv(n),vc(n),gv(n),vg(n))
      do k=1,n;tv(k)=real(k,dp);end do
      vc=abs(y);gv=1;vg=1
      if(present(time))then;if(size(time)/=n)then;result%status=5;return;end if;tv=time;end if
      if(present(var_covariate))then
         if(size(var_covariate)/=n)then;result%status=6;return;end if;vc=var_covariate
      end if
      if(present(group))then;if(size(group)/=n)then;result%status=7;return;end if;gv=group;end if
      if(present(var_group))then;if(size(var_group)/=n)then;result%status=8;return;end if;vg=var_group;end if
      if(present(coordinates))then
         if(size(coordinates,1)/=n)then;result%status=9;return;end if
      end if

      corr=correlation_spec();if(present(correlation))corr=correlation
      if(present(coordinates))then
         call initialize_correlation(corr,tv,gv,coordinates)
      else
         call initialize_correlation(corr,tv,gv)
      end if
      vbase=variance_spec();if(present(variance))vbase=variance
      call initialize_variance(vbase,vg)
      call variance_sd(vbase,vc,vg,base_sd,istat)
      if(istat/=0)then;result%status=10;return;end if
      result%variance_sd=base_sd

      gctl=gamlss_control_t();if(present(control))gctl=control
      ! Obtain stable family-specific starting values from ordinary RS.
      call fit_gamlss_model(y,x_mu,family,initial,method=GAMLSS_METHOD_RS,x_sigma=x_sigma,x_nu=x_nu,x_tau=x_tau, &
         weights=weights,offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
         control=gctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      call copy_initial(initial,np,b)

      nctl=nlme_control();if(present(nlme_control_in))nctl=nlme_control_in
      nctl%reml=.false.
      nouter=max(3,gctl%n_cyc);if(present(max_outer))nouter=max(1,max_outer)
      crit=max(1.0e-8_dp,gctl%c_crit);if(present(tolerance))crit=max(0.0_dp,tolerance)
      olddev=independent_deviance(y,w,family,b,np)
      outer_converged=.false.
      allocate(ww(n),z(n),var_work(n))

      do it=1,nouter
         outer_change=0.0_dp
         do j=1,np
            call working_values(j,y,w,family,b,np,z,ww)
            var_work=(base_sd*base_sd)/max(ww,1.0e-12_dp)
            vwork=variance_spec();vwork%kind=VAR_FIXED;vwork%fixed=.true.
            corr_use=corr
            ! Estimate shared correlation on the location working response only.
            if(j>1)corr_use%fixed=.true.
            if(present(coordinates))then
               call fit_gls(z-b(j)%offset,b(j)%x,gr,correlation=corr_use,variance=vwork,time=tv,group=gv, &
                  var_covariate=var_work,var_group=vg,coordinates=coordinates,control=nctl)
            else
               call fit_gls(z-b(j)%offset,b(j)%x,gr,correlation=corr_use,variance=vwork,time=tv,group=gv, &
                  var_covariate=var_work,var_group=vg,control=nctl)
            end if
            if(gr%status/=0)then;result%status=30+gr%status;return;end if
            oldbeta=b(j)%beta;oldeta=b(j)%eta
            step=max(0.0_dp,min(1.0_dp,parameter_step(j,gctl)))
            b(j)%beta=oldbeta+step*(gr%beta-oldbeta)
            b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
            newdev=independent_deviance(y,w,family,b,np)
            if(gctl%autostep.and.newdev>olddev)then
               do k=1,10
                  b(j)%beta=0.5_dp*(b(j)%beta+oldbeta)
                  b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
                  newdev=independent_deviance(y,w,family,b,np)
                  if(newdev<=olddev)exit
               end do
            end if
            outer_change=max(outer_change,maxval(abs(b(j)%beta-oldbeta)/(1.0_dp+abs(oldbeta))))
            b(j)%covariance=gr%beta_cov
            olddev=newdev
            if(j==1.and.allocated(gr%correlation_parameters).and..not.corr%fixed)then
               if(allocated(corr%par))deallocate(corr%par)
               allocate(corr%par(size(gr%correlation_parameters)));corr%par=gr%correlation_parameters
            end if
         end do
         newdev=independent_deviance(y,w,family,b,np)
         if(outer_change<crit)then
            outer_converged=.true.;exit
         end if
         olddev=newdev
      end do
      call fill_result(y,w,family,b,np,result%model)
      if(allocated(corr%par).and..not.corr%fixed)then
         result%model%df_fit=result%model%df_fit+real(size(corr%par),dp)
         result%model%df_residual=sum(w)-result%model%df_fit
         result%model%aic=result%model%global_deviance+2.0_dp*result%model%df_fit
         result%model%sbc=result%model%global_deviance+log(max(1.0_dp,sum(w)))*result%model%df_fit
      end if
      result%model%iterations=it;result%model%converged=outer_converged
      result%iterations=it;result%converged=outer_converged;result%status=0
      result%correlation=corr;result%variance=vbase
      if(allocated(corr%par))then
         allocate(result%correlation_parameters(size(corr%par)));result%correlation_parameters=corr%par
      else
         allocate(result%correlation_parameters(0))
      end if
   end subroutine fit_gamlss_correlated_rs

   subroutine make_designs(n,np,xmu,xs_in,xn_in,xt_in,b,status)
      integer,intent(in)::n,np
      real(dp),intent(in)::xmu(:,:)
      real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
      type(rs_block_t),intent(out)::b(4)
      integer,intent(out)::status
      status=0;b(1)%x=xmu
      if(np>=2)then
         if(present(xs_in))then
            if(size(xs_in,1)/=n)then;status=1;return;end if;b(2)%x=xs_in
         else;allocate(b(2)%x(n,1));b(2)%x=1.0_dp;end if
      end if
      if(np>=3)then
         if(present(xn_in))then
            if(size(xn_in,1)/=n)then;status=2;return;end if;b(3)%x=xn_in
         else;allocate(b(3)%x(n,1));b(3)%x=1.0_dp;end if
      end if
      if(np>=4)then
         if(present(xt_in))then
            if(size(xt_in,1)/=n)then;status=3;return;end if;b(4)%x=xt_in
         else;allocate(b(4)%x(n,1));b(4)%x=1.0_dp;end if
      end if
   end subroutine make_designs

   subroutine make_offsets(n,np,om,os,on,ot,b,status)
      integer,intent(in)::n,np
      real(dp),intent(in),optional::om(:),os(:),on(:),ot(:)
      type(rs_block_t),intent(inout)::b(4)
      integer,intent(out)::status
      integer::j
      status=0
      do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
      if(present(om))then;if(size(om)/=n)then;status=1;return;end if;b(1)%offset=om;end if
      if(np>=2.and.present(os))then;if(size(os)/=n)then;status=2;return;end if;b(2)%offset=os;end if
      if(np>=3.and.present(on))then;if(size(on)/=n)then;status=3;return;end if;b(3)%offset=on;end if
      if(np>=4.and.present(ot))then;if(size(ot)/=n)then;status=4;return;end if;b(4)%offset=ot;end if
   end subroutine make_offsets

   subroutine copy_initial(model,np,b)
      type(gamlss_result_t),intent(in)::model
      integer,intent(in)::np
      type(rs_block_t),intent(inout)::b(4)
      b(1)%beta=model%mu%coefficients;b(1)%eta=matmul(b(1)%x,b(1)%beta)+b(1)%offset
      if(np>=2)then;b(2)%beta=model%sigma%coefficients;b(2)%eta=matmul(b(2)%x,b(2)%beta)+b(2)%offset;end if
      if(np>=3)then;b(3)%beta=model%nu%coefficients;b(3)%eta=matmul(b(3)%x,b(3)%beta)+b(3)%offset;end if
      if(np>=4)then;b(4)%beta=model%tau%coefficients;b(4)%eta=matmul(b(4)%x,b(4)%beta)+b(4)%offset;end if
      b(1)%covariance=model%mu%covariance
      if(np>=2)b(2)%covariance=model%sigma%covariance
      if(np>=3)b(3)%covariance=model%nu%covariance
      if(np>=4)b(4)%covariance=model%tau%covariance
   end subroutine copy_initial

   subroutine working_values(j,y,w,family,b,np,z,ww)
      integer,intent(in)::j,family,np
      real(dp),intent(in)::y(:),w(:)
      type(rs_block_t),intent(in)::b(4)
      real(dp),intent(out)::z(:),ww(:)
      real(dp)::score,hess
      integer::i
      do i=1,size(y)
         call eta_derivative(i,j,y(i),family,b,np,score,hess)
         hess=min(hess,-1.0e-10_dp)
         ww(i)=w(i)*min(1.0e10_dp,max(1.0e-10_dp,-hess))
         z(i)=b(j)%eta(i)+score/max(1.0e-10_dp,-hess)
         if(.not.ieee_is_finite(z(i)))z(i)=b(j)%eta(i)
      end do
   end subroutine working_values

   subroutine eta_derivative(i,j,yi,family,b,np,score,hess)
      integer,intent(in)::i,j,family,np
      real(dp),intent(in)::yi
      type(rs_block_t),intent(in)::b(4)
      real(dp),intent(out)::score,hess
      real(dp)::e(4),h,l0,lp,lm
      integer::k
      e=0.0_dp;do k=1,np;e(k)=b(k)%eta(i);end do
      h=1.0e-4_dp*(1.0_dp+abs(e(j)))
      l0=loglik_eta(yi,family,e)
      e(j)=e(j)+h;lp=loglik_eta(yi,family,e)
      e(j)=e(j)-2.0_dp*h;lm=loglik_eta(yi,family,e)
      if(.not.(ieee_is_finite(l0).and.ieee_is_finite(lp).and.ieee_is_finite(lm)))then
         score=0.0_dp;hess=-1.0e-8_dp
      else
         score=(lp-lm)/(2.0_dp*h);hess=(lp-2.0_dp*l0+lm)/(h*h)
      end if
   end subroutine eta_derivative

   real(dp) function loglik_eta(yi,family,e) result(lp)
      real(dp),intent(in)::yi,e(4)
      integer,intent(in)::family
      real(dp)::a,s,c,d
      call map_parameters(family,e(1),e(2),e(3),e(4),a,s,c,d)
      lp=family_logpdf(family,yi,a,s,c,d)
   end function loglik_eta

   real(dp) function independent_deviance(y,w,family,b,np) result(dev)
      real(dp),intent(in)::y(:),w(:)
      integer,intent(in)::family,np
      type(rs_block_t),intent(in)::b(4)
      real(dp)::e(4),lp
      integer::i,j
      dev=0.0_dp
      do i=1,size(y)
         e=0.0_dp;do j=1,np;e(j)=b(j)%eta(i);end do
         lp=loglik_eta(y(i),family,e)
         if(.not.ieee_is_finite(lp))then;dev=huge(1.0_dp)/100.0_dp;return;end if
         dev=dev-2.0_dp*w(i)*lp
      end do
   end function independent_deviance

   real(dp) function max_block_score_change(y,w,family,b,np) result(v)
      real(dp),intent(in)::y(:),w(:)
      integer,intent(in)::family,np
      type(rs_block_t),intent(in)::b(4)
      real(dp)::score,hess
      integer::i,j
      v=0.0_dp
      do j=1,np
         do i=1,size(y)
            call eta_derivative(i,j,y(i),family,b,np,score,hess)
            v=max(v,abs(w(i)*score)/(1.0_dp+abs(hess)))
         end do
      end do
   end function max_block_score_change

   real(dp) function parameter_step(j,ctl) result(s)
      integer,intent(in)::j
      type(gamlss_control_t),intent(in)::ctl
      select case(j)
      case(1);s=ctl%mu_step
      case(2);s=ctl%sigma_step
      case(3);s=ctl%nu_step
      case(4);s=ctl%tau_step
      end select
   end function parameter_step

   subroutine fill_result(y,w,family,b,np,model)
      real(dp),intent(in)::y(:),w(:)
      integer,intent(in)::family,np
      type(rs_block_t),intent(in)::b(4)
      type(gamlss_result_t),intent(out)::model
      real(dp)::e(4),a,s,c,d,lp,score,hess
      integer::i,j,p
      model%family=family;model%method=GAMLSS_METHOD_RS;model%status=0
      model%global_deviance=independent_deviance(y,w,family,b,np)
      model%penalized_deviance=model%global_deviance;model%df_fit=0.0_dp
      do j=1,np;model%df_fit=model%df_fit+real(size(b(j)%beta),dp);end do
      model%df_residual=sum(w)-model%df_fit
      model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(max(1.0_dp,sum(w)))*model%df_fit
      allocate(model%residuals(size(y)),model%case_deviance(size(y)))
      do i=1,size(y)
         e=0.0_dp;do j=1,np;e(j)=b(j)%eta(i);end do
         lp=loglik_eta(y(i),family,e);model%case_deviance(i)=-2.0_dp*lp
         call eta_derivative(i,1,y(i),family,b,np,score,hess)
         model%residuals(i)=score/sqrt(max(1.0e-12_dp,-hess))
      end do
      call copy_block(b(1),model%mu)
      if(np>=2)call copy_block(b(2),model%sigma)
      if(np>=3)call copy_block(b(3),model%nu)
      if(np>=4)call copy_block(b(4),model%tau)
      do i=1,size(y)
         e=0.0_dp;do j=1,np;e(j)=b(j)%eta(i);end do
         call map_parameters(family,e(1),e(2),e(3),e(4),a,s,c,d)
         model%mu%fitted(i)=a
         if(np>=2)model%sigma%fitted(i)=s
         if(np>=3)model%nu%fitted(i)=c
         if(np>=4)model%tau%fitted(i)=d
      end do
   end subroutine fill_result

   subroutine copy_block(b,out)
      type(rs_block_t),intent(in)::b
      type(gamlss_parameter_result_t),intent(out)::out
      out%coefficients=b%beta;out%eta=b%eta;out%covariance=b%covariance
      allocate(out%fitted(size(b%eta)));out%fitted=0.0_dp
      out%edf=real(size(b%beta),dp);out%penalty=0.0_dp;out%lambda=0.0_dp
   end subroutine copy_block
end module gamlss_correlated_rs_v05
