! Gaussian-copula joint likelihoods for continuous GAMLSS margins.
! This is a genuine joint density, distinct from correlated RS estimating equations.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_copula_v06
   use gamlss_kinds, only : dp
   use gamlss_types, only : gamlss_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf,GAMLSS_BEINF
   use gamlss_family_support, only : family_cdf,family_is_discrete
   use gamlss_special, only : normal_quantile
   use gamlss_optim, only : bfgs_minimize
   use nlme_types, only : correlation_spec,COR_NONE,COR_UNSTRUCTURED
   use nlme_gls, only : initialize_correlation
   use nlme_correlation, only : correlation_matrix,correlation_to_unconstrained, &
      correlation_from_unconstrained
   use nlme_linalg, only : solve_spd,find_group_indices,unique_integers
   use nlme_status, only : NLME_SUCCESS
   implicit none
   private
   public :: gaussian_copula_result_t,fit_gamlss_gaussian_copula

   type :: design_block_t
      real(dp),allocatable :: x(:,:),offset(:)
      integer :: p=0,start=0
   end type design_block_t

   type,public :: gaussian_copula_result_t
      type(gamlss_result_t) :: model
      type(correlation_spec) :: correlation
      real(dp),allocatable :: correlation_parameters(:)
      real(dp),allocatable :: parameter_covariance(:,:)
      real(dp),allocatable :: gaussian_scores(:)
      real(dp) :: marginal_log_likelihood=-huge(1.0_dp)
      real(dp) :: copula_log_likelihood=-huge(1.0_dp)
      real(dp) :: joint_log_likelihood=-huge(1.0_dp)
      real(dp) :: aic=huge(1.0_dp)
      integer :: optimizer_status=0
      integer :: status=0
      logical :: converged=.false.
   end type gaussian_copula_result_t

   type :: copula_context_t
      real(dp),allocatable :: y(:),time(:),coordinates(:,:)
      integer,allocatable :: group(:),levels(:)
      type(design_block_t) :: b(4)
      type(correlation_spec) :: corr_template
      integer :: family=0,np=0,n=0,nmax=0,kcorr=0,ncoef=0
      logical :: estimate_corr=.false.,has_coordinates=.false.
   end type copula_context_t

   type(copula_context_t),save :: ctx
contains

   subroutine fit_gamlss_gaussian_copula(y,x_mu,family,result,correlation,x_sigma,x_nu,x_tau, &
      offset_mu,offset_sigma,offset_nu,offset_tau,time,group,coordinates,control,max_iter,tolerance)
      real(dp),intent(in) :: y(:),x_mu(:,:)
      integer,intent(in) :: family
      type(gaussian_copula_result_t),intent(out) :: result
      type(correlation_spec),intent(in),optional :: correlation
      real(dp),intent(in),optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:)
      real(dp),intent(in),optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      real(dp),intent(in),optional :: time(:),coordinates(:,:)
      integer,intent(in),optional :: group(:)
      type(gamlss_control_t),intent(in),optional :: control
      integer,intent(in),optional :: max_iter
      real(dp),intent(in),optional :: tolerance

      type(gamlss_result_t) :: initial
      type(gamlss_control_t) :: ctl
      type(correlation_spec) :: corr
      real(dp),allocatable :: theta(:),cx(:),hinv(:,:),z(:)
      integer :: n,np,j,pos,kc,istat,niter
      real(dp) :: fval,tol,mll,cll

      call clear_context()
      n=size(y);np=family_npar(family)
      if(n<2.or.np<1.or.np>4.or.size(x_mu,1)/=n)then;result%status=1;return;end if
      if(family_is_discrete(family).or.family==GAMLSS_BEINF)then;result%status=2;return;end if
      ctx%n=n;ctx%np=np;ctx%family=family;ctx%y=y
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,ctx%b,istat)
      if(istat/=0)then;result%status=3;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,ctx%b,istat)
      if(istat/=0)then;result%status=4;return;end if
      allocate(ctx%time(n),ctx%group(n))
      ctx%time=[(real(j,dp),j=1,n)];ctx%group=1
      if(present(time))then
         if(size(time)/=n)then;result%status=5;return;end if
         ctx%time=time
      end if
      if(present(group))then
         if(size(group)/=n)then;result%status=6;return;end if
         ctx%group=group
      end if
      if(present(coordinates))then
         if(size(coordinates,1)/=n)then;result%status=7;return;end if
         ctx%coordinates=coordinates;ctx%has_coordinates=.true.
      end if
      call unique_integers(ctx%group,ctx%levels,istat)
      if(istat/=NLME_SUCCESS)then;result%status=8;return;end if
      ctx%nmax=0
      do j=1,size(ctx%levels);ctx%nmax=max(ctx%nmax,count(ctx%group==ctx%levels(j)));end do
      corr=correlation_spec();if(present(correlation))corr=correlation
      if(ctx%has_coordinates)then
         call initialize_correlation(corr,ctx%time,ctx%group,ctx%coordinates)
      else
         call initialize_correlation(corr,ctx%time,ctx%group)
      end if
      if(corr%kind==COR_UNSTRUCTURED.and.size(ctx%levels)>1)then
         do j=1,size(ctx%levels)
            if(count(ctx%group==ctx%levels(j))/=ctx%nmax)then;result%status=9;return;end if
         end do
      end if
      ctx%corr_template=corr;ctx%estimate_corr=.not.corr%fixed.and.corr%kind/=COR_NONE
      call correlation_to_unconstrained(corr,ctx%nmax,cx,istat)
      if(istat/=NLME_SUCCESS)then;result%status=10;return;end if
      kc=size(cx);ctx%kcorr=merge(kc,0,ctx%estimate_corr)

      ctl=gamlss_control_t();if(present(control))ctl=control
      call fit_initial_margin(y,x_mu,family,initial,x_sigma,x_nu,x_tau,offset_mu,offset_sigma, &
         offset_nu,offset_tau,ctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      ctx%ncoef=0
      do j=1,np
         ctx%b(j)%p=size(parameter_coefficients(initial,j))
         ctx%b(j)%start=ctx%ncoef+1;ctx%ncoef=ctx%ncoef+ctx%b(j)%p
      end do
      allocate(theta(ctx%ncoef+ctx%kcorr));pos=0
      do j=1,np
         theta(pos+1:pos+ctx%b(j)%p)=parameter_coefficients(initial,j);pos=pos+ctx%b(j)%p
      end do
      if(ctx%estimate_corr.and.kc>0)theta(ctx%ncoef+1:)=cx
      allocate(hinv(size(theta),size(theta)));hinv=0.0_dp
      niter=160;if(present(max_iter))niter=max(1,max_iter)
      tol=2.0e-6_dp;if(present(tolerance))tol=max(1.0e-10_dp,tolerance)
      call bfgs_minimize(copula_objective,theta,fval,istat,max_iter=niter,tol=tol,inverse_hessian=hinv)
      result%optimizer_status=istat
      if(istat/=0.and.istat/=1)then;result%status=30+istat;return;end if
      call evaluate_theta(theta,mll,cll,z,corr,istat)
      if(istat/=0)then;result%status=40+istat;return;end if
      call fill_model_from_theta(theta,hinv,mll+cll,result%model)
      result%correlation=corr
      if(allocated(corr%par))then;result%correlation_parameters=corr%par
      else;allocate(result%correlation_parameters(0));end if
      result%parameter_covariance=hinv;result%gaussian_scores=z
      result%marginal_log_likelihood=mll;result%copula_log_likelihood=cll
      result%joint_log_likelihood=mll+cll
      result%aic=-2.0_dp*result%joint_log_likelihood+2.0_dp*real(size(theta),dp)
      result%status=0;result%converged=(istat==0);result%model%converged=result%converged
      call clear_context()
   end subroutine fit_gamlss_gaussian_copula

   real(dp) function copula_objective(theta) result(f)
      real(dp),intent(in) :: theta(:)
      real(dp) :: mll,cll
      real(dp),allocatable :: z(:)
      type(correlation_spec) :: corr
      integer :: status
      call evaluate_theta(theta,mll,cll,z,corr,status)
      if(status/=0.or..not.(mll+cll>-huge(1.0_dp)))then
         f=huge(1.0_dp)/100.0_dp
      else
         f=-(mll+cll)
      end if
   end function copula_objective

   subroutine evaluate_theta(theta,mll,cll,z,corr,status)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: mll,cll
      real(dp),allocatable,intent(out) :: z(:)
      type(correlation_spec),intent(out) :: corr
      integer,intent(out) :: status
      real(dp),allocatable :: eta(:,:),a(:),b(:),c(:),d(:),u(:),tl(:),zl(:),r(:,:),rinvz(:),cl(:,:)
      integer,allocatable :: idx(:)
      integer :: i,j,g,pos,istat
      real(dp) :: lp,ld,pclip,e2,e3,e4
      status=0;mll=0.0_dp;cll=0.0_dp
      if(size(theta)/=ctx%ncoef+ctx%kcorr)then;status=1;return;end if
      allocate(eta(ctx%n,ctx%np),a(ctx%n),b(ctx%n),c(ctx%n),d(ctx%n),u(ctx%n),z(ctx%n))
      pos=0
      do j=1,ctx%np
         eta(:,j)=matmul(ctx%b(j)%x,theta(pos+1:pos+ctx%b(j)%p))+ctx%b(j)%offset
         pos=pos+ctx%b(j)%p
      end do
      if(ctx%estimate_corr)then
         call correlation_from_unconstrained(ctx%corr_template,ctx%nmax,theta(ctx%ncoef+1:),corr,istat)
         if(istat/=NLME_SUCCESS)then;status=2;return;end if
      else
         corr=ctx%corr_template
      end if
      do i=1,ctx%n
         e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(ctx%np>=2)e2=eta(i,2)
         if(ctx%np>=3)e3=eta(i,3)
         if(ctx%np>=4)e4=eta(i,4)
         call map_parameters(ctx%family,eta(i,1),e2,e3,e4,a(i),b(i),c(i),d(i))
         lp=family_logpdf(ctx%family,ctx%y(i),a(i),b(i),c(i),d(i))
         if(.not.(lp>-huge(1.0_dp)/10.0_dp))then;status=3;return;end if
         u(i)=family_cdf(ctx%family,ctx%y(i),a(i),b(i),c(i),d(i))
         if(u(i)<0.0_dp)then;status=4;return;end if
         pclip=min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,u(i)))
         z(i)=normal_quantile(pclip);mll=mll+lp
      end do
      do g=1,size(ctx%levels)
         call find_group_indices(ctx%group,ctx%levels(g),idx)
         if(size(idx)<=1)cycle
         allocate(tl(size(idx)),zl(size(idx)));tl=ctx%time(idx);zl=z(idx)
         if(ctx%has_coordinates)then
            allocate(cl(size(idx),size(ctx%coordinates,2)));cl=ctx%coordinates(idx,:)
            call correlation_matrix(corr,tl,r,istat,cl);deallocate(cl)
         else
            call correlation_matrix(corr,tl,r,istat)
         end if
         if(istat/=NLME_SUCCESS)then;status=5;return;end if
         call solve_spd(r,zl,rinvz,istat,ld)
         if(istat/=NLME_SUCCESS)then;status=6;return;end if
         cll=cll-0.5_dp*(ld+dot_product(zl,rinvz)-dot_product(zl,zl))
         deallocate(tl,zl,r,rinvz,idx)
      end do
   end subroutine evaluate_theta

   subroutine fill_model_from_theta(theta,hinv,jll,model)
      real(dp),intent(in) :: theta(:),hinv(:,:),jll
      type(gamlss_result_t),intent(out) :: model
      real(dp),allocatable :: eta(:,:),a(:),b(:),c(:),d(:),lp(:)
      integer :: j,pos,i,p,lo,hi
      real(dp) :: e1,e2,e3,e4
      allocate(eta(ctx%n,ctx%np),a(ctx%n),b(ctx%n),c(ctx%n),d(ctx%n),lp(ctx%n))
      pos=0
      do j=1,ctx%np
         p=ctx%b(j)%p;lo=pos+1;hi=pos+p
         eta(:,j)=matmul(ctx%b(j)%x,theta(lo:hi))+ctx%b(j)%offset
         select case(j)
         case(1)
            model%mu%coefficients=theta(lo:hi);model%mu%eta=eta(:,j);model%mu%covariance=hinv(lo:hi,lo:hi)
         case(2)
            model%sigma%coefficients=theta(lo:hi);model%sigma%eta=eta(:,j);model%sigma%covariance=hinv(lo:hi,lo:hi)
         case(3)
            model%nu%coefficients=theta(lo:hi);model%nu%eta=eta(:,j);model%nu%covariance=hinv(lo:hi,lo:hi)
         case(4)
            model%tau%coefficients=theta(lo:hi);model%tau%eta=eta(:,j);model%tau%covariance=hinv(lo:hi,lo:hi)
         end select
         pos=hi
      end do
      do i=1,ctx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(ctx%np>=2)e2=eta(i,2)
         if(ctx%np>=3)e3=eta(i,3)
         if(ctx%np>=4)e4=eta(i,4)
         call map_parameters(ctx%family,e1,e2,e3,e4,a(i),b(i),c(i),d(i))
         lp(i)=family_logpdf(ctx%family,ctx%y(i),a(i),b(i),c(i),d(i))
      end do
      model%mu%fitted=a
      if(ctx%np>=2)model%sigma%fitted=b
      if(ctx%np>=3)model%nu%fitted=c
      if(ctx%np>=4)model%tau%fitted=d
      allocate(model%residuals(ctx%n),model%case_deviance(ctx%n))
      model%residuals=ctx%y-a;model%case_deviance=-2.0_dp*lp
      model%global_deviance=-2.0_dp*jll;model%penalized_deviance=model%global_deviance
      model%df_fit=real(size(theta),dp);model%df_residual=real(ctx%n,dp)-model%df_fit
      model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(real(max(1,ctx%n),dp))*model%df_fit
      model%family=ctx%family;model%status=0;model%converged=.true.
   end subroutine fill_model_from_theta

   subroutine fit_initial_margin(y,xmu,family,fit,xs,xn,xt,om,os,on,ot,ctl)
      real(dp),intent(in) :: y(:),xmu(:,:)
      integer,intent(in) :: family
      type(gamlss_result_t),intent(out) :: fit
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:),om(:),os(:),on(:),ot(:)
      type(gamlss_control_t),intent(in) :: ctl
      call fit_gamlss_model(y,xmu,family,fit,method=GAMLSS_METHOD_RS,x_sigma=xs,x_nu=xn,x_tau=xt, &
         offset_mu=om,offset_sigma=os,offset_nu=on,offset_tau=ot,control=ctl)
   end subroutine fit_initial_margin

   function parameter_coefficients(model,j) result(v)
      type(gamlss_result_t),intent(in) :: model
      integer,intent(in) :: j
      real(dp),allocatable :: v(:)
      select case(j)
      case(1);v=model%mu%coefficients
      case(2);v=model%sigma%coefficients
      case(3);v=model%nu%coefficients
      case(4);v=model%tau%coefficients
      end select
   end function parameter_coefficients

   subroutine make_designs(n,np,xmu,xs,xn,xt,b,status)
      integer,intent(in) :: n,np
      real(dp),intent(in) :: xmu(:,:)
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:)
      type(design_block_t),intent(out) :: b(4)
      integer,intent(out) :: status
      status=0;b(1)%x=xmu
      if(np>=2)then
         if(present(xs))then;if(size(xs,1)/=n)then;status=1;return;end if;b(2)%x=xs
         else;allocate(b(2)%x(n,1));b(2)%x=1.0_dp;end if
      end if
      if(np>=3)then
         if(present(xn))then;if(size(xn,1)/=n)then;status=2;return;end if;b(3)%x=xn
         else;allocate(b(3)%x(n,1));b(3)%x=1.0_dp;end if
      end if
      if(np>=4)then
         if(present(xt))then;if(size(xt,1)/=n)then;status=3;return;end if;b(4)%x=xt
         else;allocate(b(4)%x(n,1));b(4)%x=1.0_dp;end if
      end if
   end subroutine make_designs

   subroutine make_offsets(n,np,om,os,on,ot,b,status)
      integer,intent(in) :: n,np
      real(dp),intent(in),optional :: om(:),os(:),on(:),ot(:)
      type(design_block_t),intent(inout) :: b(4)
      integer,intent(out) :: status
      integer :: j
      status=0
      do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
      if(present(om))then;if(size(om)/=n)then;status=1;return;end if;b(1)%offset=om;end if
      if(np>=2.and.present(os))then;if(size(os)/=n)then;status=2;return;end if;b(2)%offset=os;end if
      if(np>=3.and.present(on))then;if(size(on)/=n)then;status=3;return;end if;b(3)%offset=on;end if
      if(np>=4.and.present(ot))then;if(size(ot)/=n)then;status=4;return;end if;b(4)%offset=ot;end if
   end subroutine make_offsets

   subroutine clear_context()
      integer :: j
      if(allocated(ctx%y))deallocate(ctx%y)
      if(allocated(ctx%time))deallocate(ctx%time)
      if(allocated(ctx%coordinates))deallocate(ctx%coordinates)
      if(allocated(ctx%group))deallocate(ctx%group)
      if(allocated(ctx%levels))deallocate(ctx%levels)
      do j=1,4
         if(allocated(ctx%b(j)%x))deallocate(ctx%b(j)%x)
         if(allocated(ctx%b(j)%offset))deallocate(ctx%b(j)%offset)
         ctx%b(j)%p=0;ctx%b(j)%start=0
      end do
      ctx%family=0;ctx%np=0;ctx%n=0;ctx%nmax=0;ctx%kcorr=0;ctx%ncoef=0
      ctx%estimate_corr=.false.;ctx%has_coordinates=.false.
   end subroutine clear_context
end module gamlss_copula_v06
