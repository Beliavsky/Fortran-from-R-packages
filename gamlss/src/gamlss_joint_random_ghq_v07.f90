! Marginal GAMLSS random-effects likelihood by tensor Gauss-Hermite quadrature.
! Intended for low-dimensional joint random-effect blocks; the scalable
! Laplace/posterior-moment fitter remains available in gamlss_joint_random_v06.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_joint_random_ghq_v07
   use gamlss_kinds, only : dp,pi,sqrt2
   use gamlss_types, only : gamlss_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_optim, only : bfgs_minimize
   use gamlss_linalg, only : cholesky_factor
   use gamlss_joint_random_v06, only : joint_random_effects_result_t,fit_gamlss_joint_random_effects
   implicit none
   private
   public :: joint_random_ghq_result_t,fit_gamlss_joint_random_effects_ghq

   type :: gh_block_t
      real(dp),allocatable :: x(:,:),offset(:)
      integer :: pfix=0,start_fix=0,active_slot=0
   end type gh_block_t

   type,public :: joint_random_ghq_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: joint_covariance(:,:),parameter_covariance(:,:)
      real(dp),allocatable :: posterior_effects(:,:,:),group_log_likelihood(:)
      integer,allocatable :: levels(:),active_parameter_ids(:)
      logical :: active(4)=.false.
      real(dp) :: marginal_log_likelihood=-huge(1.0_dp),aic=huge(1.0_dp)
      integer :: q=0,latent_dimension=0,quadrature_order=0,optimizer_status=0,status=0
      logical :: converged=.false.
   end type joint_random_ghq_result_t

   type :: gh_context_t
      real(dp),allocatable :: y(:),weights(:),z(:,:,:)
      integer,allocatable :: group(:),levels(:),group_index(:),active_ids(:)
      type(gh_block_t) :: b(4)
      integer :: family=0,np=0,n=0,q=0,ng=0,nact=0,r=0,nfix=0,ncov=0,order=5
      logical :: active(4)=.false.
   end type gh_context_t
   type(gh_context_t),save :: gctx
contains

   subroutine fit_gamlss_joint_random_effects_ghq(y,x_mu,z_random,group,family,result,active_parameters, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control, &
      initial_covariance,quadrature_order,max_iter,tolerance)
      real(dp),intent(in) :: y(:),x_mu(:,:),z_random(:,:,:)
      integer,intent(in) :: group(:),family
      type(joint_random_ghq_result_t),intent(out) :: result
      logical,intent(in),optional :: active_parameters(4)
      real(dp),intent(in),optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional :: control
      real(dp),intent(in),optional :: initial_covariance(:,:)
      integer,intent(in),optional :: quadrature_order,max_iter
      real(dp),intent(in),optional :: tolerance
      type(gamlss_result_t) :: initial
      type(joint_random_effects_result_t) :: lap
      type(gamlss_control_t) :: ctl
      real(dp),allocatable :: sigma(:,:),theta(:),hinv(:,:),gll(:),post(:,:,:)
      integer :: n,np,q,j,a,pos,istat,niter
      real(dp) :: fval,tol,ll
      logical :: act(4)

      call clear_context();n=size(y);np=family_npar(family);q=size(z_random,2)
      if(n<2.or.np<1.or.np>4.or.q<1.or.size(x_mu,1)/=n.or.size(group)/=n)then;result%status=1;return;end if
      if(size(z_random,1)/=n.or.size(z_random,3)<np)then;result%status=2;return;end if
      act=.false.;act(1)=.true.;if(present(active_parameters))act=active_parameters
      if(any(act(np+1:4)).or..not.any(act(1:np)))then;result%status=3;return;end if
      gctx%n=n;gctx%np=np;gctx%q=q;gctx%family=family;gctx%active=act;gctx%y=y;gctx%group=group
      gctx%nact=count(act(1:np));gctx%r=gctx%nact*q
      if(gctx%r>4)then;result%status=4;return;end if
      gctx%order=5;if(present(quadrature_order))gctx%order=quadrature_order
      if(gctx%order/=5.and.gctx%order/=7)then;result%status=5;return;end if
      allocate(gctx%weights(n));gctx%weights=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;result%status=6;return;end if
         gctx%weights=weights
      end if
      gctx%z=z_random(:,:,1:np)
      call make_levels(group,gctx%levels,gctx%group_index,istat)
      if(istat/=0)then;result%status=7;return;end if
      gctx%ng=size(gctx%levels);allocate(gctx%active_ids(gctx%nact));a=0
      do j=1,np
         if(act(j))then;a=a+1;gctx%active_ids(a)=j;gctx%b(j)%active_slot=a;end if
      end do
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,gctx%b,istat)
      if(istat/=0)then;result%status=8;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,gctx%b,istat)
      if(istat/=0)then;result%status=9;return;end if
      ctl=gamlss_control_t();if(present(control))ctl=control
      call fit_initial(y,x_mu,family,initial,x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma, &
         offset_nu,offset_tau,ctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      call make_fixed_layout(initial)
      allocate(sigma(gctx%r,gctx%r));sigma=0.0_dp
      do a=1,gctx%r;sigma(a,a)=0.15_dp;end do
      if(present(initial_covariance))then
         if(any(shape(initial_covariance)/=[gctx%r,gctx%r]))then;result%status=10;return;end if
         sigma=initial_covariance
      else
         call fit_gamlss_joint_random_effects(y,x_mu,z_random,group,family,lap,active_parameters=act, &
            x_sigma=x_sigma,x_nu=x_nu,x_tau=x_tau,weights=weights,offset_mu=offset_mu,offset_sigma=offset_sigma, &
            offset_nu=offset_nu,offset_tau=offset_tau,control=ctl,max_outer=3,max_inner=45,tol_cov=1.0e-3_dp)
         if(lap%status==0.and.all(shape(lap%joint_covariance)==[gctx%r,gctx%r]))then
            sigma=lap%joint_covariance;initial=lap%model
         end if
      end if
      call covariance_to_theta(initial,sigma,theta,istat)
      if(istat/=0)then;result%status=11;return;end if
      allocate(hinv(size(theta),size(theta)));hinv=0.0_dp
      niter=100;if(present(max_iter))niter=max(1,max_iter)
      tol=5.0e-6_dp;if(present(tolerance))tol=max(1.0e-10_dp,tolerance)
      call bfgs_minimize(ghq_objective,theta,fval,istat,max_iter=niter,tol=tol,inverse_hessian=hinv)
      result%optimizer_status=istat
      if(istat/=0.and.istat/=1)then;result%status=30+istat;call clear_context();return;end if
      call evaluate_theta(theta,ll,gll,sigma,post,istat)
      if(istat/=0)then;result%status=40+istat;call clear_context();return;end if
      call fill_model(theta,hinv,ll,result%model)
      result%joint_covariance=sigma;result%parameter_covariance=hinv;result%posterior_effects=post
      result%group_log_likelihood=gll;result%levels=gctx%levels;result%active_parameter_ids=gctx%active_ids
      result%active=act;result%q=q;result%latent_dimension=gctx%r;result%quadrature_order=gctx%order
      result%marginal_log_likelihood=ll;result%aic=-2.0_dp*ll+2.0_dp*real(size(theta),dp)
      result%status=0;result%converged=(istat==0);result%model%converged=result%converged
      call clear_context()
   end subroutine fit_gamlss_joint_random_effects_ghq

   real(dp) function ghq_objective(theta) result(f)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable :: gll(:),sigma(:,:),post(:,:,:)
      real(dp) :: ll
      integer :: status
      call evaluate_theta(theta,ll,gll,sigma,post,status)
      if(status/=0.or..not.(ll>-huge(1.0_dp)/10.0_dp))then;f=huge(1.0_dp)/100.0_dp
      else;f=-ll;end if
   end function ghq_objective

   subroutine evaluate_theta(theta,ll,gll,sigma,post,status)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: ll
      real(dp),allocatable,intent(out) :: gll(:),sigma(:,:),post(:,:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: l(:,:),eta0(:,:),nodes(:),wgh(:),bvec(:),terms(:),bstore(:,:)
      integer :: g,ncomb,k,istat
      status=0;ll=0.0_dp
      call theta_to_covariance(theta,sigma,l,istat)
      if(istat/=0)then;status=1;return;end if
      call fixed_eta(theta,eta0)
      call gh_nodes(gctx%order,nodes,wgh)
      ncomb=gctx%order**gctx%r
      allocate(gll(gctx%ng),post(gctx%ng,gctx%q,4),bvec(gctx%r),terms(ncomb),bstore(gctx%r,ncomb))
      gll=0.0_dp;post=0.0_dp
      do g=1,gctx%ng
         call group_quadrature(g,eta0,l,nodes,wgh,terms,bstore,gll(g),bvec,istat)
         if(istat/=0)then;status=2;return;end if
         call store_posterior_mean(g,bvec,post);ll=ll+gll(g)
      end do
   end subroutine evaluate_theta

   subroutine group_quadrature(g,eta0,l,nodes,wgh,terms,bstore,gll,bmean,status)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),l(:,:),nodes(:),wgh(:)
      real(dp),intent(out) :: terms(:),bstore(:,:),gll,bmean(:)
      integer,intent(out) :: status
      real(dp),allocatable :: zstd(:),b(:)
      real(dp) :: logw,lc,maxterm,sw,w
      integer :: code,tmp,digit,a,ncomb
      status=0;ncomb=size(terms);allocate(zstd(gctx%r),b(gctx%r));maxterm=-huge(1.0_dp)
      do code=0,ncomb-1
         tmp=code;logw=0.0_dp
         do a=1,gctx%r
            digit=mod(tmp,gctx%order)+1;tmp=tmp/gctx%order
            zstd(a)=sqrt2*nodes(digit);logw=logw+log(wgh(digit))
         end do
         b=matmul(l,zstd);call conditional_group_loglik(g,eta0,b,lc,status)
         if(status/=0)return
         terms(code+1)=logw+lc;bstore(:,code+1)=b;maxterm=max(maxterm,terms(code+1))
      end do
      sw=0.0_dp;bmean=0.0_dp
      do code=1,ncomb
         w=exp(terms(code)-maxterm);sw=sw+w;bmean=bmean+w*bstore(:,code)
      end do
      if(sw<=0.0_dp)then;status=1;return;end if
      bmean=bmean/sw;gll=maxterm+log(sw)-0.5_dp*real(gctx%r,dp)*log(pi)
   end subroutine group_quadrature

   subroutine conditional_group_loglik(g,eta0,b,lc,status)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),b(:)
      real(dp),intent(out) :: lc
      integer,intent(out) :: status
      real(dp) :: eta(4),aa,bb,cc,dd,lp
      integer :: i,j,s,a,pos
      status=0;lc=0.0_dp
      do i=1,gctx%n
         if(gctx%group_index(i)/=g)cycle
         eta=0.0_dp;eta(1:gctx%np)=eta0(i,1:gctx%np)
         pos=0
         do s=1,gctx%nact
            j=gctx%active_ids(s)
            do a=1,gctx%q;pos=pos+1;eta(j)=eta(j)+gctx%z(i,a,j)*b(pos);end do
         end do
         call map_parameters(gctx%family,eta(1),eta(2),eta(3),eta(4),aa,bb,cc,dd)
         lp=family_logpdf(gctx%family,gctx%y(i),aa,bb,cc,dd)
         if(.not.(lp>-huge(1.0_dp)/10.0_dp))then;status=1;return;end if
         lc=lc+gctx%weights(i)*lp
      end do
   end subroutine conditional_group_loglik

   subroutine store_posterior_mean(g,b,post)
      integer,intent(in) :: g
      real(dp),intent(in) :: b(:)
      real(dp),intent(inout) :: post(:,:,:)
      integer :: s,j,a,pos
      pos=0
      do s=1,gctx%nact
         j=gctx%active_ids(s)
         do a=1,gctx%q;pos=pos+1;post(g,a,j)=b(pos);end do
      end do
   end subroutine store_posterior_mean

   subroutine fixed_eta(theta,eta)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable,intent(out) :: eta(:,:)
      integer :: j,lo,hi
      allocate(eta(gctx%n,gctx%np))
      do j=1,gctx%np
         lo=gctx%b(j)%start_fix;hi=lo+gctx%b(j)%pfix-1
         eta(:,j)=matmul(gctx%b(j)%x,theta(lo:hi))+gctx%b(j)%offset
      end do
   end subroutine fixed_eta

   subroutine make_fixed_layout(initial)
      type(gamlss_result_t),intent(in) :: initial
      real(dp),allocatable :: coef(:)
      integer :: j,pos
      pos=1
      do j=1,gctx%np
         coef=parameter_coefficients(initial,j);gctx%b(j)%pfix=size(coef);gctx%b(j)%start_fix=pos
         pos=pos+size(coef)
      end do
      gctx%nfix=pos-1;gctx%ncov=gctx%r*(gctx%r+1)/2
   end subroutine make_fixed_layout

   subroutine covariance_to_theta(initial,sigma,theta,status)
      type(gamlss_result_t),intent(in) :: initial
      real(dp),intent(in) :: sigma(:,:)
      real(dp),allocatable,intent(out) :: theta(:)
      integer,intent(out) :: status
      real(dp),allocatable :: l(:,:),coef(:)
      integer :: j,pos,i,k
      call cholesky_factor(0.5_dp*(sigma+transpose(sigma)),l,status)
      if(status/=0)return
      allocate(theta(gctx%nfix+gctx%ncov));pos=1
      do j=1,gctx%np
         coef=parameter_coefficients(initial,j);theta(pos:pos+size(coef)-1)=coef;pos=pos+size(coef)
      end do
      do i=1,gctx%r
         do k=1,i
            if(i==k)then;theta(pos)=log(max(l(i,i),1.0e-6_dp));else;theta(pos)=l(i,k);end if
            pos=pos+1
         end do
      end do
   end subroutine covariance_to_theta

   subroutine theta_to_covariance(theta,sigma,l,status)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable,intent(out) :: sigma(:,:),l(:,:)
      integer,intent(out) :: status
      integer :: i,k,pos
      status=0
      if(size(theta)/=gctx%nfix+gctx%ncov)then;status=1;return;end if
      allocate(l(gctx%r,gctx%r));l=0.0_dp;pos=gctx%nfix+1
      do i=1,gctx%r
         do k=1,i
            if(i==k)then;l(i,k)=exp(min(6.0_dp,max(-12.0_dp,theta(pos))))
            else;l(i,k)=theta(pos);end if
            pos=pos+1
         end do
      end do
      sigma=matmul(l,transpose(l))
   end subroutine theta_to_covariance

   subroutine fill_model(theta,hinv,ll,model)
      real(dp),intent(in) :: theta(:),hinv(:,:),ll
      type(gamlss_result_t),intent(out) :: model
      real(dp),allocatable :: eta(:,:),aa(:),bb(:),cc(:),dd(:),lp(:)
      integer :: j,lo,hi,i
      real(dp) :: e1,e2,e3,e4
      call fixed_eta(theta,eta);allocate(aa(gctx%n),bb(gctx%n),cc(gctx%n),dd(gctx%n),lp(gctx%n))
      do j=1,gctx%np
         lo=gctx%b(j)%start_fix;hi=lo+gctx%b(j)%pfix-1
         select case(j)
         case(1);model%mu%coefficients=theta(lo:hi);model%mu%eta=eta(:,j);model%mu%covariance=hinv(lo:hi,lo:hi)
         case(2);model%sigma%coefficients=theta(lo:hi);model%sigma%eta=eta(:,j);model%sigma%covariance=hinv(lo:hi,lo:hi)
         case(3);model%nu%coefficients=theta(lo:hi);model%nu%eta=eta(:,j);model%nu%covariance=hinv(lo:hi,lo:hi)
         case(4);model%tau%coefficients=theta(lo:hi);model%tau%eta=eta(:,j);model%tau%covariance=hinv(lo:hi,lo:hi)
         end select
      end do
      do i=1,gctx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(gctx%np>=2)e2=eta(i,2)
         if(gctx%np>=3)e3=eta(i,3)
         if(gctx%np>=4)e4=eta(i,4)
         call map_parameters(gctx%family,e1,e2,e3,e4,aa(i),bb(i),cc(i),dd(i))
         lp(i)=family_logpdf(gctx%family,gctx%y(i),aa(i),bb(i),cc(i),dd(i))
      end do
      model%mu%fitted=aa;if(gctx%np>=2)model%sigma%fitted=bb
      if(gctx%np>=3)model%nu%fitted=cc;if(gctx%np>=4)model%tau%fitted=dd
      allocate(model%residuals(gctx%n),model%case_deviance(gctx%n));model%residuals=gctx%y-aa
      model%case_deviance=-2.0_dp*lp;model%global_deviance=-2.0_dp*ll;model%penalized_deviance=model%global_deviance
      model%df_fit=real(size(theta),dp);model%df_residual=sum(gctx%weights)-model%df_fit
      model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(max(1.0_dp,sum(gctx%weights)))*model%df_fit
      model%family=gctx%family;model%status=0;model%converged=.true.
   end subroutine fill_model

   subroutine gh_nodes(order,x,w)
      integer,intent(in) :: order
      real(dp),allocatable,intent(out) :: x(:),w(:)
      if(order==5)then
         x=[-2.020182870456086_dp,-0.958572464613819_dp,0.0_dp,0.958572464613819_dp,2.020182870456086_dp]
         w=[0.019953242059046_dp,0.393619323152241_dp,0.945308720482942_dp, &
            0.393619323152241_dp,0.019953242059046_dp]
      else
         x=[-2.651961356835234_dp,-1.673551628767471_dp,-0.816287882858965_dp,0.0_dp, &
            0.816287882858965_dp,1.673551628767471_dp,2.651961356835234_dp]
         w=[0.000971781245100_dp,0.054515582819127_dp,0.425607252610128_dp,0.810264617556807_dp, &
            0.425607252610128_dp,0.054515582819127_dp,0.000971781245100_dp]
      end if
   end subroutine gh_nodes

   subroutine fit_initial(y,xmu,family,fit,xs,xn,xt,w,om,os,on,ot,ctl)
      real(dp),intent(in) :: y(:),xmu(:,:)
      integer,intent(in) :: family
      type(gamlss_result_t),intent(out) :: fit
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:),w(:),om(:),os(:),on(:),ot(:)
      type(gamlss_control_t),intent(in) :: ctl
      call fit_gamlss_model(y,xmu,family,fit,method=GAMLSS_METHOD_RS,x_sigma=xs,x_nu=xn,x_tau=xt, &
         weights=w,offset_mu=om,offset_sigma=os,offset_nu=on,offset_tau=ot,control=ctl)
   end subroutine fit_initial

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

   subroutine make_levels(group,levels,index,status)
      integer,intent(in) :: group(:)
      integer,allocatable,intent(out) :: levels(:),index(:)
      integer,intent(out) :: status
      integer,allocatable :: tmp(:)
      integer :: i,j,nlev
      status=0;allocate(tmp(size(group)));nlev=0
      do i=1,size(group)
         if(nlev==0.or..not.any(tmp(1:nlev)==group(i)))then;nlev=nlev+1;tmp(nlev)=group(i);end if
      end do
      allocate(levels(nlev),index(size(group)));levels=tmp(1:nlev)
      do i=1,size(group)
         index(i)=0
         do j=1,nlev;if(group(i)==levels(j))then;index(i)=j;exit;end if;end do
         if(index(i)==0)then;status=1;return;end if
      end do
   end subroutine make_levels

   subroutine make_designs(n,np,xmu,xs,xn,xt,b,status)
      integer,intent(in) :: n,np
      real(dp),intent(in) :: xmu(:,:)
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:)
      type(gh_block_t),intent(out) :: b(4)
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
      type(gh_block_t),intent(inout) :: b(4)
      integer,intent(out) :: status
      integer :: j
      status=0;do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
      if(present(om))then;if(size(om)/=n)then;status=1;return;end if;b(1)%offset=om;end if
      if(np>=2.and.present(os))then;if(size(os)/=n)then;status=2;return;end if;b(2)%offset=os;end if
      if(np>=3.and.present(on))then;if(size(on)/=n)then;status=3;return;end if;b(3)%offset=on;end if
      if(np>=4.and.present(ot))then;if(size(ot)/=n)then;status=4;return;end if;b(4)%offset=ot;end if
   end subroutine make_offsets

   subroutine clear_context()
      integer :: j
      if(allocated(gctx%y))deallocate(gctx%y);if(allocated(gctx%weights))deallocate(gctx%weights)
      if(allocated(gctx%z))deallocate(gctx%z);if(allocated(gctx%group))deallocate(gctx%group)
      if(allocated(gctx%levels))deallocate(gctx%levels);if(allocated(gctx%group_index))deallocate(gctx%group_index)
      if(allocated(gctx%active_ids))deallocate(gctx%active_ids)
      do j=1,4
         if(allocated(gctx%b(j)%x))deallocate(gctx%b(j)%x)
         if(allocated(gctx%b(j)%offset))deallocate(gctx%b(j)%offset)
         gctx%b(j)%pfix=0;gctx%b(j)%start_fix=0;gctx%b(j)%active_slot=0
      end do
      gctx%family=0;gctx%np=0;gctx%n=0;gctx%q=0;gctx%ng=0;gctx%nact=0;gctx%r=0
      gctx%nfix=0;gctx%ncov=0;gctx%order=5;gctx%active=.false.
   end subroutine clear_context
end module gamlss_joint_random_ghq_v07
