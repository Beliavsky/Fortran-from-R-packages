! Adaptive quasi-Monte-Carlo importance sampling for joint GAMLSS random effects.
! The proposal for each group is centered at its posterior mode and uses the
! inverse local posterior Hessian.  This avoids tensor GHQ's order**dimension
! growth and is intended for moderate-dimensional random-effect blocks.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_joint_random_ais_v08
   use gamlss_kinds, only : dp,log2pi
   use gamlss_types, only : gamlss_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_optim, only : bfgs_minimize
   use gamlss_linalg, only : cholesky_factor,invert_matrix
   use gamlss_special, only : normal_quantile
   use gamlss_joint_random_v06, only : joint_random_effects_result_t,fit_gamlss_joint_random_effects
   implicit none
   private
   public :: joint_random_ais_result_t,fit_gamlss_joint_random_effects_ais

   type :: ais_block_t
      real(dp),allocatable :: x(:,:),offset(:)
      integer :: pfix=0,start_fix=0,active_slot=0
   end type ais_block_t

   type,public :: joint_random_ais_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: joint_covariance(:,:),parameter_covariance(:,:)
      real(dp),allocatable :: posterior_effects(:,:,:),posterior_covariance(:,:,:)
      real(dp),allocatable :: group_log_likelihood(:),effective_sample_size(:)
      integer,allocatable :: levels(:),active_parameter_ids(:)
      logical :: active(4)=.false.
      real(dp) :: marginal_log_likelihood=-huge(1.0_dp),aic=huge(1.0_dp)
      real(dp) :: minimum_ess=0.0_dp,proposal_scale=1.0_dp
      integer :: q=0,latent_dimension=0,qmc_points=0,optimizer_status=0,status=0
      logical :: converged=.false.,parameters_refined=.false.
   end type joint_random_ais_result_t

   type :: ais_context_t
      real(dp),allocatable :: y(:),weights(:),z(:,:,:)
      integer,allocatable :: group(:),levels(:),group_index(:),active_ids(:)
      type(ais_block_t) :: b(4)
      integer :: family=0,np=0,n=0,q=0,ng=0,nact=0,r=0,nfix=0,ncov=0,nqmc=512
      real(dp) :: proposal_scale=1.25_dp
      logical :: active(4)=.false.
   end type ais_context_t
   type(ais_context_t),save :: actx

   integer,save :: mode_group=0
   real(dp),allocatable,save :: mode_eta0(:,:),mode_precision(:,:),mode_lprior(:,:)

contains

   subroutine fit_gamlss_joint_random_effects_ais(y,x_mu,z_random,group,family,result,active_parameters, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control,initial_covariance, &
      qmc_points,proposal_scale,refine_parameters,max_iter,tolerance)
      real(dp),intent(in) :: y(:),x_mu(:,:),z_random(:,:,:)
      integer,intent(in) :: group(:),family
      type(joint_random_ais_result_t),intent(out) :: result
      logical,intent(in),optional :: active_parameters(4),refine_parameters
      real(dp),intent(in),optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional :: control
      real(dp),intent(in),optional :: initial_covariance(:,:),proposal_scale,tolerance
      integer,intent(in),optional :: qmc_points,max_iter
      type(gamlss_result_t) :: initial
      type(joint_random_effects_result_t) :: lap
      type(gamlss_control_t) :: ctl
      real(dp),allocatable :: sigma(:,:),theta(:),hinv(:,:),gll(:),post(:,:,:),pcov(:,:,:),ess(:)
      integer :: n,np,q,j,a,istat,niter,optstat
      real(dp) :: fval,tol,ll
      logical :: act(4),refine

      call clear_context()
      n=size(y);np=family_npar(family);q=size(z_random,2)
      if(n<2.or.np<1.or.np>4.or.q<1.or.size(x_mu,1)/=n.or.size(group)/=n)then
         result%status=1;return
      end if
      if(size(z_random,1)/=n.or.size(z_random,3)<np)then;result%status=2;return;end if
      act=.false.;act(1)=.true.;if(present(active_parameters))act=active_parameters
      if(any(act(np+1:4)).or..not.any(act(1:np)))then;result%status=3;return;end if
      actx%n=n;actx%np=np;actx%q=q;actx%family=family;actx%active=act;actx%y=y;actx%group=group
      actx%nact=count(act(1:np));actx%r=actx%nact*q
      if(actx%r>16)then;result%status=4;return;end if
      actx%nqmc=512;if(present(qmc_points))actx%nqmc=max(64,qmc_points)
      if(mod(actx%nqmc,2)/=0)actx%nqmc=actx%nqmc+1
      actx%proposal_scale=1.25_dp
      if(present(proposal_scale))actx%proposal_scale=max(0.5_dp,min(3.0_dp,proposal_scale))
      allocate(actx%weights(n));actx%weights=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;result%status=5;return;end if
         actx%weights=weights
      end if
      actx%z=z_random(:,:,1:np)
      call make_levels(group,actx%levels,actx%group_index,istat)
      if(istat/=0)then;result%status=6;return;end if
      actx%ng=size(actx%levels);allocate(actx%active_ids(actx%nact));a=0
      do j=1,np
         if(act(j))then;a=a+1;actx%active_ids(a)=j;actx%b(j)%active_slot=a;end if
      end do
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,actx%b,istat)
      if(istat/=0)then;result%status=7;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,actx%b,istat)
      if(istat/=0)then;result%status=8;return;end if
      ctl=gamlss_control_t();if(present(control))ctl=control
      call fit_initial(y,x_mu,family,initial,x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma, &
         offset_nu,offset_tau,ctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      call make_fixed_layout(initial)
      allocate(sigma(actx%r,actx%r));sigma=0.0_dp
      do a=1,actx%r;sigma(a,a)=0.15_dp;end do
      if(present(initial_covariance))then
         if(any(shape(initial_covariance)/=[actx%r,actx%r]))then;result%status=9;return;end if
         sigma=initial_covariance
      else
         call fit_gamlss_joint_random_effects(y,x_mu,z_random,group,family,lap,active_parameters=act, &
            x_sigma=x_sigma,x_nu=x_nu,x_tau=x_tau,weights=weights,offset_mu=offset_mu, &
            offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau,control=ctl, &
            max_outer=4,max_inner=60,tol_cov=1.0e-3_dp)
         if(lap%status==0.and.all(shape(lap%joint_covariance)==[actx%r,actx%r]))then
            sigma=lap%joint_covariance;initial=lap%model
         end if
      end if
      call stabilize_spd(sigma,istat)
      if(istat/=0)then;result%status=10;return;end if
      call covariance_to_theta(initial,sigma,theta,istat)
      if(istat/=0)then;result%status=11;return;end if
      allocate(hinv(size(theta),size(theta)));hinv=0.0_dp
      refine=.false.;if(present(refine_parameters))refine=refine_parameters
      niter=30;if(present(max_iter))niter=max(1,max_iter)
      tol=2.0e-5_dp;if(present(tolerance))tol=max(1.0e-9_dp,tolerance)
      optstat=0
      if(refine)then
         call bfgs_minimize(ais_objective,theta,fval,istat,max_iter=niter,tol=tol,inverse_hessian=hinv)
         result%optimizer_status=istat;optstat=istat
         if(istat/=0.and.istat/=1)then;result%status=30+istat;call clear_context();return;end if
      else
         hinv=0.0_dp;istat=0;result%optimizer_status=0
      end if
      call evaluate_theta(theta,ll,gll,sigma,post,pcov,ess,istat)
      if(istat/=0)then;result%status=40+istat;call clear_context();return;end if
      call fill_model(theta,hinv,ll,result%model)
      result%joint_covariance=sigma;result%parameter_covariance=hinv
      result%posterior_effects=post;result%posterior_covariance=pcov
      result%group_log_likelihood=gll;result%effective_sample_size=ess
      result%levels=actx%levels;result%active_parameter_ids=actx%active_ids;result%active=act
      result%q=q;result%latent_dimension=actx%r;result%qmc_points=actx%nqmc
      result%proposal_scale=actx%proposal_scale;result%minimum_ess=minval(ess)
      result%marginal_log_likelihood=ll;result%aic=-2.0_dp*ll+2.0_dp*real(size(theta),dp)
      result%status=0;result%converged=(.not.refine.or.optstat==0);result%parameters_refined=refine
      result%model%converged=result%converged
      call clear_context()
   end subroutine fit_gamlss_joint_random_effects_ais

   real(dp) function ais_objective(theta) result(f)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable :: gll(:),sigma(:,:),post(:,:,:),pcov(:,:,:),ess(:)
      real(dp) :: ll
      integer :: status
      call evaluate_theta(theta,ll,gll,sigma,post,pcov,ess,status)
      if(status/=0.or..not.(ll>-huge(1.0_dp)/10.0_dp))then
         f=huge(1.0_dp)/100.0_dp
      else
         f=-ll
      end if
   end function ais_objective

   subroutine evaluate_theta(theta,ll,gll,sigma,post,pcov,ess,status)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: ll
      real(dp),allocatable,intent(out) :: gll(:),sigma(:,:),post(:,:,:),pcov(:,:,:),ess(:)
      integer,intent(out) :: status
      real(dp),allocatable :: lprior(:,:),precision(:,:),eta0(:,:),bmean(:),bcov(:,:)
      integer :: g,istat
      status=0;ll=0.0_dp
      call theta_to_covariance(theta,sigma,lprior,istat)
      if(istat/=0)then;status=1;return;end if
      call invert_matrix(sigma,precision,istat)
      if(istat/=0)then;status=2;return;end if
      call fixed_eta(theta,eta0)
      allocate(gll(actx%ng),post(actx%ng,actx%q,4),pcov(actx%ng,actx%r,actx%r),ess(actx%ng))
      allocate(bmean(actx%r),bcov(actx%r,actx%r));gll=0.0_dp;post=0.0_dp;pcov=0.0_dp;ess=0.0_dp
      do g=1,actx%ng
         call group_importance(g,eta0,lprior,precision,gll(g),bmean,bcov,ess(g),istat)
         if(istat/=0)then;status=3;return;end if
         call store_posterior_mean(g,bmean,post);pcov(g,:,:)=bcov;ll=ll+gll(g)
      end do
   end subroutine evaluate_theta

   subroutine group_importance(g,eta0,lprior,precision,gll,bmean,bcov,ess,status)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),lprior(:,:),precision(:,:)
      real(dp),intent(out) :: gll,bmean(:),bcov(:,:),ess
      integer,intent(out) :: status
      real(dp),allocatable :: mode(:),hinv(:,:),hess(:,:),pc(:,:),lprop(:,:),z(:),b(:)
      real(dp),allocatable :: lw(:),bs(:,:),wraw(:)
      real(dp) :: fmode,maxlw,sw,sw2,w,ltarget,lq,fval
      integer :: istat,k,d,m
      status=0;m=actx%nqmc
      call set_mode_context(g,eta0,precision,lprior)
      allocate(mode(actx%r),hinv(actx%r,actx%r),hess(actx%r,actx%r));mode=0.0_dp;hinv=0.0_dp
      call find_posterior_mode(mode,fmode,hess,istat)
      if(istat/=0)then;status=1;return;end if
      call invert_matrix(hess,pc,istat)
      if(istat/=0)then;status=3;return;end if
      pc=(actx%proposal_scale**2)*pc;call stabilize_spd(pc,istat)
      if(istat/=0)then;status=4;return;end if
      call cholesky_factor(pc,lprop,istat)
      if(istat/=0)then;status=5;return;end if
      allocate(z(actx%r),b(actx%r),lw(m),bs(actx%r,m),wraw(m));maxlw=-huge(1.0_dp)
      do k=1,m
         do d=1,actx%r
            z(d)=normal_quantile(qmc_normal_probability(k,d,m,g))
         end do
         b=mode+matmul(lprop,z);bs(:,k)=b
         call log_target(g,eta0,b,lprior,precision,ltarget,istat)
         if(istat/=0)then;status=6;return;end if
         lq=normal_logpdf_chol(b,mode,lprop,istat)
         if(istat/=0)then;status=7;return;end if
         lw(k)=ltarget-lq;maxlw=max(maxlw,lw(k))
      end do
      sw=0.0_dp;sw2=0.0_dp;bmean=0.0_dp
      do k=1,m
         wraw(k)=exp(lw(k)-maxlw);sw=sw+wraw(k);sw2=sw2+wraw(k)*wraw(k)
         bmean=bmean+wraw(k)*bs(:,k)
      end do
      if(sw<=0.0_dp.or.sw2<=0.0_dp)then;status=8;return;end if
      bmean=bmean/sw;bcov=0.0_dp
      do k=1,m
         b=bs(:,k)-bmean;w=wraw(k)/sw;bcov=bcov+w*outer_product(b,b)
      end do
      gll=maxlw+log(sw)-log(real(m,dp));ess=sw*sw/sw2
      fval=fmode
   end subroutine group_importance


   subroutine local_gradient(x,grad)
      real(dp),intent(in) :: x(:)
      real(dp),intent(out) :: grad(:)
      real(dp),allocatable :: xp(:),xm(:)
      real(dp) :: h
      integer :: j
      allocate(xp(size(x)),xm(size(x)))
      do j=1,size(x)
         h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(j)))
         xp=x;xm=x;xp(j)=x(j)+h;xm(j)=x(j)-h
         grad(j)=(posterior_objective(xp)-posterior_objective(xm))/(2.0_dp*h)
      end do
   end subroutine local_gradient

   subroutine local_hessian(x,hess)
      real(dp),intent(in) :: x(:)
      real(dp),intent(out) :: hess(:,:)
      real(dp),allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp) :: hi,hj,f0
      integer :: i,j,n
      n=size(x);allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n));f0=posterior_objective(x);hess=0.0_dp
      do i=1,n
         hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
         xp=x;xm=x;xp(i)=x(i)+hi;xm(i)=x(i)-hi
         hess(i,i)=(posterior_objective(xp)-2.0_dp*f0+posterior_objective(xm))/(hi*hi)
         do j=i+1,n
            hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
            xpp=x;xpm=x;xmp=x;xmm=x
            xpp(i)=x(i)+hi;xpp(j)=x(j)+hj;xpm(i)=x(i)+hi;xpm(j)=x(j)-hj
            xmp(i)=x(i)-hi;xmp(j)=x(j)+hj;xmm(i)=x(i)-hi;xmm(j)=x(j)-hj
            hess(i,j)=(posterior_objective(xpp)-posterior_objective(xpm)-posterior_objective(xmp)+ &
               posterior_objective(xmm))/(4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
   end subroutine local_hessian

   subroutine find_posterior_mode(mode,fval,hess,status)
      real(dp),intent(inout) :: mode(:)
      real(dp),intent(out) :: fval,hess(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: grad(:),hinv(:,:),step(:),trial(:)
      real(dp) :: f0,fnew,alpha,scale
      integer :: iter,istat
      allocate(grad(size(mode)),hinv(size(mode),size(mode)),step(size(mode)),trial(size(mode)))
      status=0;f0=posterior_objective(mode)
      if(.not.(f0<huge(1.0_dp)/1000.0_dp))then;status=1;return;end if
      do iter=1,18
         call local_gradient(mode,grad)
         call local_hessian(mode,hess)
         call stabilize_spd(hess,istat);if(istat/=0)then;status=2;return;end if
         call invert_matrix(hess,hinv,istat);if(istat/=0)then;status=3;return;end if
         step=-matmul(hinv,grad);scale=maxval(abs(step))
         if(scale>2.0_dp)step=step*(2.0_dp/scale)
         alpha=1.0_dp
         do
            trial=mode+alpha*step;fnew=posterior_objective(trial)
            if(fnew<=f0.or.alpha<=1.0e-4_dp)exit
            alpha=0.5_dp*alpha
         end do
         if(.not.(fnew<huge(1.0_dp)/1000.0_dp))then;status=4;return;end if
         mode=trial;f0=fnew
         if(maxval(abs(alpha*step))<2.0e-6_dp)exit
      end do
      call local_hessian(mode,hess)
      call stabilize_spd(hess,istat);if(istat/=0)then;status=5;return;end if
      fval=f0
   end subroutine find_posterior_mode

   subroutine set_mode_context(g,eta0,precision,lprior)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),precision(:,:),lprior(:,:)
      mode_group=g
      if(allocated(mode_eta0))deallocate(mode_eta0)
      if(allocated(mode_precision))deallocate(mode_precision)
      if(allocated(mode_lprior))deallocate(mode_lprior)
      mode_eta0=eta0;mode_precision=precision;mode_lprior=lprior
   end subroutine set_mode_context

   real(dp) function posterior_objective(b) result(f)
      real(dp),intent(in) :: b(:)
      real(dp) :: lc,lp
      integer :: status
      call conditional_group_loglik(mode_group,mode_eta0,b,lc,status)
      if(status/=0)then;f=huge(1.0_dp)/100.0_dp;return;end if
      lp=normal_logpdf_precision(b,mode_precision,mode_lprior,status)
      if(status/=0)then;f=huge(1.0_dp)/100.0_dp;return;end if
      f=-(lc+lp)
   end function posterior_objective

   subroutine log_target(g,eta0,b,lprior,precision,ltarget,status)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),b(:),lprior(:,:),precision(:,:)
      real(dp),intent(out) :: ltarget
      integer,intent(out) :: status
      real(dp) :: lc,lp
      call conditional_group_loglik(g,eta0,b,lc,status);if(status/=0)return
      lp=normal_logpdf_precision(b,precision,lprior,status);if(status/=0)return
      ltarget=lc+lp
   end subroutine log_target

   real(dp) function normal_logpdf_precision(x,precision,l,status) result(lp)
      real(dp),intent(in) :: x(:),precision(:,:),l(:,:)
      integer,intent(out) :: status
      real(dp) :: logdet
      integer :: i,n
      n=size(x);status=0;lp=-huge(1.0_dp)
      if(any(shape(precision)/=[n,n]).or.any(shape(l)/=[n,n]))then;status=1;return;end if
      logdet=0.0_dp
      do i=1,n;logdet=logdet+2.0_dp*log(l(i,i));end do
      lp=-0.5_dp*(real(n,dp)*log2pi+logdet+dot_product(x,matmul(precision,x)))
   end function normal_logpdf_precision

   real(dp) function normal_logpdf_chol(x,mean,l,status) result(lp)
      real(dp),intent(in) :: x(:),mean(:),l(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: z(:)
      real(dp) :: logdet
      integer :: i,j,n
      n=size(x);status=0;lp=-huge(1.0_dp)
      if(size(mean)/=n.or.any(shape(l)/=[n,n]))then;status=1;return;end if
      allocate(z(n));z=x-mean;logdet=0.0_dp
      do i=1,n
         do j=1,i-1;z(i)=z(i)-l(i,j)*z(j);end do
         if(l(i,i)<=0.0_dp)then;status=2;return;end if
         z(i)=z(i)/l(i,i);logdet=logdet+2.0_dp*log(l(i,i))
      end do
      lp=-0.5_dp*(real(n,dp)*log2pi+logdet+dot_product(z,z))
   end function normal_logpdf_chol

   real(dp) function qmc_normal_probability(k,dim,m,g) result(u)
      integer,intent(in) :: k,dim,m,g
      integer,parameter :: primes(32)=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71, &
         73,79,83,89,97,101,103,107,109,113,127,131]
      integer :: kk,base,digit,half
      real(dp) :: f,shift
      half=(m+1)/2;kk=1+mod(k-1,half);base=primes(1+mod(dim-1,size(primes)))
      u=0.0_dp;f=1.0_dp/real(base,dp)
      do while(kk>0)
         digit=mod(kk,base);u=u+f*real(digit,dp);kk=kk/base;f=f/real(base,dp)
      end do
      shift=modulo(0.6180339887498949_dp*real(dim+7*g+m,dp),1.0_dp)
      u=modulo(u+shift,1.0_dp)
      if(k>half)u=1.0_dp-u
      u=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,u))
   end function qmc_normal_probability

   pure function outer_product(a,b) result(c)
      real(dp),intent(in) :: a(:),b(:)
      real(dp) :: c(size(a),size(b))
      integer :: i,j
      do i=1,size(a);do j=1,size(b);c(i,j)=a(i)*b(j);end do;end do
   end function outer_product

   subroutine conditional_group_loglik(g,eta0,b,lc,status)
      integer,intent(in) :: g
      real(dp),intent(in) :: eta0(:,:),b(:)
      real(dp),intent(out) :: lc
      integer,intent(out) :: status
      real(dp) :: eta(4),aa,bb,cc,dd,lp
      integer :: i,j,s,a,pos
      status=0;lc=0.0_dp
      do i=1,actx%n
         if(actx%group_index(i)/=g)cycle
         eta=0.0_dp;eta(1:actx%np)=eta0(i,1:actx%np);pos=0
         do s=1,actx%nact
            j=actx%active_ids(s)
            do a=1,actx%q;pos=pos+1;eta(j)=eta(j)+actx%z(i,a,j)*b(pos);end do
         end do
         call map_parameters(actx%family,eta(1),eta(2),eta(3),eta(4),aa,bb,cc,dd)
         lp=family_logpdf(actx%family,actx%y(i),aa,bb,cc,dd)
         if(.not.(lp>-huge(1.0_dp)/10.0_dp))then;status=1;return;end if
         lc=lc+actx%weights(i)*lp
      end do
   end subroutine conditional_group_loglik

   subroutine store_posterior_mean(g,b,post)
      integer,intent(in) :: g
      real(dp),intent(in) :: b(:)
      real(dp),intent(inout) :: post(:,:,:)
      integer :: s,j,a,pos
      pos=0
      do s=1,actx%nact
         j=actx%active_ids(s)
         do a=1,actx%q;pos=pos+1;post(g,a,j)=b(pos);end do
      end do
   end subroutine store_posterior_mean

   subroutine fixed_eta(theta,eta)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable,intent(out) :: eta(:,:)
      integer :: j,lo,hi
      allocate(eta(actx%n,actx%np))
      do j=1,actx%np
         lo=actx%b(j)%start_fix;hi=lo+actx%b(j)%pfix-1
         eta(:,j)=matmul(actx%b(j)%x,theta(lo:hi))+actx%b(j)%offset
      end do
   end subroutine fixed_eta

   subroutine make_fixed_layout(initial)
      type(gamlss_result_t),intent(in) :: initial
      real(dp),allocatable :: coef(:)
      integer :: j,pos
      pos=1
      do j=1,actx%np
         coef=parameter_coefficients(initial,j);actx%b(j)%pfix=size(coef);actx%b(j)%start_fix=pos
         pos=pos+size(coef)
      end do
      actx%nfix=pos-1;actx%ncov=actx%r*(actx%r+1)/2
   end subroutine make_fixed_layout

   subroutine covariance_to_theta(initial,sigma,theta,status)
      type(gamlss_result_t),intent(in) :: initial
      real(dp),intent(in) :: sigma(:,:)
      real(dp),allocatable,intent(out) :: theta(:)
      integer,intent(out) :: status
      real(dp),allocatable :: l(:,:),coef(:)
      integer :: j,pos,i,k
      call cholesky_factor(0.5_dp*(sigma+transpose(sigma)),l,status);if(status/=0)return
      allocate(theta(actx%nfix+actx%ncov));pos=1
      do j=1,actx%np
         coef=parameter_coefficients(initial,j);theta(pos:pos+size(coef)-1)=coef;pos=pos+size(coef)
      end do
      do i=1,actx%r
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
      if(size(theta)/=actx%nfix+actx%ncov)then;status=1;return;end if
      allocate(l(actx%r,actx%r));l=0.0_dp;pos=actx%nfix+1
      do i=1,actx%r
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
      call fixed_eta(theta,eta);allocate(aa(actx%n),bb(actx%n),cc(actx%n),dd(actx%n),lp(actx%n))
      do j=1,actx%np
         lo=actx%b(j)%start_fix;hi=lo+actx%b(j)%pfix-1
         select case(j)
         case(1)
            model%mu%coefficients=theta(lo:hi);model%mu%eta=eta(:,j);model%mu%covariance=hinv(lo:hi,lo:hi)
         case(2)
            model%sigma%coefficients=theta(lo:hi);model%sigma%eta=eta(:,j)
            model%sigma%covariance=hinv(lo:hi,lo:hi)
         case(3)
            model%nu%coefficients=theta(lo:hi);model%nu%eta=eta(:,j);model%nu%covariance=hinv(lo:hi,lo:hi)
         case(4)
            model%tau%coefficients=theta(lo:hi);model%tau%eta=eta(:,j);model%tau%covariance=hinv(lo:hi,lo:hi)
         end select
      end do
      do i=1,actx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(actx%np>=2)e2=eta(i,2)
         if(actx%np>=3)e3=eta(i,3)
         if(actx%np>=4)e4=eta(i,4)
         call map_parameters(actx%family,e1,e2,e3,e4,aa(i),bb(i),cc(i),dd(i))
         lp(i)=family_logpdf(actx%family,actx%y(i),aa(i),bb(i),cc(i),dd(i))
      end do
      model%mu%fitted=aa;if(actx%np>=2)model%sigma%fitted=bb
      if(actx%np>=3)model%nu%fitted=cc;if(actx%np>=4)model%tau%fitted=dd
      allocate(model%residuals(actx%n),model%case_deviance(actx%n));model%residuals=actx%y-aa
      model%case_deviance=-2.0_dp*lp;model%global_deviance=-2.0_dp*ll
      model%penalized_deviance=model%global_deviance;model%df_fit=real(size(theta),dp)
      model%df_residual=sum(actx%weights)-model%df_fit;model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(max(1.0_dp,sum(actx%weights)))*model%df_fit
      model%family=actx%family;model%status=0;model%converged=.true.
   end subroutine fill_model

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
      type(ais_block_t),intent(out) :: b(4)
      integer,intent(out) :: status
      status=0;b(1)%x=xmu
      if(np>=2)then
         if(present(xs))then
            if(size(xs,1)/=n)then;status=1;return;end if;b(2)%x=xs
         else;allocate(b(2)%x(n,1));b(2)%x=1.0_dp;end if
      end if
      if(np>=3)then
         if(present(xn))then
            if(size(xn,1)/=n)then;status=2;return;end if;b(3)%x=xn
         else;allocate(b(3)%x(n,1));b(3)%x=1.0_dp;end if
      end if
      if(np>=4)then
         if(present(xt))then
            if(size(xt,1)/=n)then;status=3;return;end if;b(4)%x=xt
         else;allocate(b(4)%x(n,1));b(4)%x=1.0_dp;end if
      end if
   end subroutine make_designs

   subroutine make_offsets(n,np,om,os,on,ot,b,status)
      integer,intent(in) :: n,np
      real(dp),intent(in),optional :: om(:),os(:),on(:),ot(:)
      type(ais_block_t),intent(inout) :: b(4)
      integer,intent(out) :: status
      integer :: j
      status=0;do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
      if(present(om))then;if(size(om)/=n)then;status=1;return;end if;b(1)%offset=om;end if
      if(np>=2.and.present(os))then;if(size(os)/=n)then;status=2;return;end if;b(2)%offset=os;end if
      if(np>=3.and.present(on))then;if(size(on)/=n)then;status=3;return;end if;b(3)%offset=on;end if
      if(np>=4.and.present(ot))then;if(size(ot)/=n)then;status=4;return;end if;b(4)%offset=ot;end if
   end subroutine make_offsets

   subroutine stabilize_spd(a,status)
      real(dp),intent(inout) :: a(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: l(:,:)
      real(dp) :: jitter
      integer :: i,k
      a=0.5_dp*(a+transpose(a));do i=1,size(a,1);a(i,i)=max(a(i,i),1.0e-8_dp);end do
      jitter=0.0_dp
      do k=1,10
         call cholesky_factor(a,l,status);if(status==0)return
         if(jitter==0.0_dp)then;jitter=1.0e-8_dp;else;jitter=10.0_dp*jitter;end if
         do i=1,size(a,1);a(i,i)=a(i,i)+jitter;end do
         if(allocated(l))deallocate(l)
      end do
   end subroutine stabilize_spd

   subroutine clear_context()
      integer :: j
      if(allocated(actx%y))deallocate(actx%y);if(allocated(actx%weights))deallocate(actx%weights)
      if(allocated(actx%z))deallocate(actx%z);if(allocated(actx%group))deallocate(actx%group)
      if(allocated(actx%levels))deallocate(actx%levels)
      if(allocated(actx%group_index))deallocate(actx%group_index)
      if(allocated(actx%active_ids))deallocate(actx%active_ids)
      do j=1,4
         if(allocated(actx%b(j)%x))deallocate(actx%b(j)%x)
         if(allocated(actx%b(j)%offset))deallocate(actx%b(j)%offset)
         actx%b(j)%pfix=0;actx%b(j)%start_fix=0;actx%b(j)%active_slot=0
      end do
      if(allocated(mode_eta0))deallocate(mode_eta0)
      if(allocated(mode_precision))deallocate(mode_precision)
      if(allocated(mode_lprior))deallocate(mode_lprior)
      mode_group=0;actx%family=0;actx%np=0;actx%n=0;actx%q=0;actx%ng=0;actx%nact=0;actx%r=0
      actx%nfix=0;actx%ncov=0;actx%nqmc=512;actx%proposal_scale=1.25_dp;actx%active=.false.
   end subroutine clear_context
end module gamlss_joint_random_ais_v08
