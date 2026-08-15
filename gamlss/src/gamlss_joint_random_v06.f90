! Joint cross-parameter grouped random effects for GAMLSS.
! A full covariance is estimated across active parameter/random-term combinations.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_joint_random_v06
   use gamlss_kinds, only : dp
   use gamlss_types, only : gamlss_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_optim, only : bfgs_minimize
   use gamlss_linalg, only : invert_matrix,cholesky_factor
   implicit none
   private
   public :: joint_random_effects_result_t,fit_gamlss_joint_random_effects

   type :: jr_block_t
      real(dp),allocatable :: x(:,:),offset(:)
      integer :: pfix=0,start_fix=0,start_re=0,active_slot=0
   end type jr_block_t

   type,public :: joint_random_effects_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: effects(:,:,:)
      real(dp),allocatable :: joint_covariance(:,:),joint_precision(:,:)
      real(dp),allocatable :: parameter_covariance(:,:)
      integer,allocatable :: levels(:),active_parameter_ids(:)
      logical :: active(4)=.false.
      integer :: q=0,iterations=0,optimizer_status=0,status=0
      logical :: converged=.false.
   end type joint_random_effects_result_t

   type :: jr_context_t
      real(dp),allocatable :: y(:),weights(:),z(:,:,:)
      integer,allocatable :: group(:),levels(:),group_index(:),active_ids(:)
      type(jr_block_t) :: b(4)
      real(dp),allocatable :: precision(:,:)
      integer :: family=0,np=0,n=0,q=0,ng=0,nact=0,r=0,ntheta=0
      logical :: active(4)=.false.
   end type jr_context_t

   type(jr_context_t),save :: jctx
contains

   subroutine fit_gamlss_joint_random_effects(y,x_mu,z_random,group,family,result,active_parameters, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control, &
      initial_covariance,max_outer,tol_cov,max_inner,tolerance)
      real(dp),intent(in) :: y(:),x_mu(:,:),z_random(:,:,:)
      integer,intent(in) :: group(:),family
      type(joint_random_effects_result_t),intent(out) :: result
      logical,intent(in),optional :: active_parameters(4)
      real(dp),intent(in),optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional :: control
      real(dp),intent(in),optional :: initial_covariance(:,:)
      integer,intent(in),optional :: max_outer,max_inner
      real(dp),intent(in),optional :: tol_cov,tolerance

      type(gamlss_result_t) :: initial
      type(gamlss_control_t) :: ctl
      real(dp),allocatable :: theta(:),hinv(:,:),sigma(:,:),snew(:,:),pinv(:,:),effects(:,:,:)
      real(dp),allocatable :: bvec(:)
      integer :: np,n,q,ng,nouter,ninner,it,istat,j,g,a,b,ia,ib
      real(dp) :: fval,crit,tol,delta,maxdelta
      logical :: act(4)

      call clear_jctx()
      n=size(y);np=family_npar(family);q=size(z_random,2)
      if(n<2.or.np<1.or.np>4.or.q<1.or.size(x_mu,1)/=n.or.size(group)/=n)then
         result%status=1;return
      end if
      if(size(z_random,1)/=n.or.size(z_random,3)<np)then;result%status=2;return;end if
      act=.false.;act(1)=.true.;if(present(active_parameters))act=active_parameters
      if(any(act(np+1:4)).or..not.any(act(1:np)))then;result%status=3;return;end if
      jctx%n=n;jctx%np=np;jctx%q=q;jctx%family=family;jctx%active=act;jctx%y=y;jctx%group=group
      allocate(jctx%weights(n));jctx%weights=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;result%status=4;return;end if
         jctx%weights=weights
      end if
      jctx%z=z_random(:,:,1:np)
      call make_levels(group,jctx%levels,jctx%group_index,istat)
      if(istat/=0)then;result%status=5;return;end if
      jctx%ng=size(jctx%levels);ng=jctx%ng
      jctx%nact=count(act(1:np));allocate(jctx%active_ids(jctx%nact));a=0
      do j=1,np
         if(act(j))then;a=a+1;jctx%active_ids(a)=j;jctx%b(j)%active_slot=a;end if
      end do
      jctx%r=jctx%nact*q
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,jctx%b,istat)
      if(istat/=0)then;result%status=6;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,jctx%b,istat)
      if(istat/=0)then;result%status=7;return;end if

      ctl=gamlss_control_t();if(present(control))ctl=control
      call fit_initial(y,x_mu,family,initial,x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma, &
         offset_nu,offset_tau,ctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      call make_layout(initial,theta)
      allocate(sigma(jctx%r,jctx%r));sigma=0.0_dp
      do a=1,jctx%r;sigma(a,a)=0.20_dp;end do
      if(present(initial_covariance))then
         if(any(shape(initial_covariance)/=[jctx%r,jctx%r]))then;result%status=8;return;end if
         sigma=initial_covariance
      end if
      call stabilize_spd(sigma,istat)
      if(istat/=0)then;result%status=9;return;end if
      call invert_matrix(sigma,pinv,istat)
      if(istat/=0)then;result%status=10;return;end if
      jctx%precision=pinv;deallocate(pinv)

      nouter=10;if(present(max_outer))nouter=max(1,max_outer)
      ninner=120;if(present(max_inner))ninner=max(1,max_inner)
      crit=2.0e-4_dp;if(present(tol_cov))crit=max(0.0_dp,tol_cov)
      tol=2.0e-6_dp;if(present(tolerance))tol=max(1.0e-10_dp,tolerance)
      allocate(hinv(size(theta),size(theta)),effects(ng,q,4),bvec(jctx%r));effects=0.0_dp

      do it=1,nouter
         hinv=0.0_dp
         call bfgs_minimize(joint_re_objective,theta,fval,istat,max_iter=ninner,tol=tol,inverse_hessian=hinv)
         result%optimizer_status=istat
         if(istat/=0.and.istat/=1)then;result%status=30+istat;return;end if
         call extract_effects(theta,effects)
         allocate(snew(jctx%r,jctx%r));snew=0.0_dp
         do g=1,ng
            call group_effect_vector(theta,g,bvec)
            do a=1,jctx%r
               ia=theta_random_index(jctx%active_ids((a-1)/q+1),g,1+mod(a-1,q))
               do b=1,jctx%r
                  ib=theta_random_index(jctx%active_ids((b-1)/q+1),g,1+mod(b-1,q))
                  snew(a,b)=snew(a,b)+bvec(a)*bvec(b)+hinv(ia,ib)
               end do
            end do
         end do
         snew=snew/real(ng,dp);call stabilize_spd(snew,istat)
         if(istat/=0)then;result%status=11;return;end if
         delta=maxval(abs(snew-sigma)/(1.0_dp+abs(sigma)));maxdelta=delta
         sigma=0.5_dp*sigma+0.5_dp*snew;deallocate(snew)
         call stabilize_spd(sigma,istat)
         if(istat/=0)then;result%status=12;return;end if
         call invert_matrix(sigma,pinv,istat)
         if(istat/=0)then;result%status=13;return;end if
         jctx%precision=pinv;deallocate(pinv)
         result%iterations=it
         if(maxdelta<crit)then;result%converged=.true.;exit;end if
      end do
      call extract_effects(theta,effects)
      call fill_joint_model(theta,hinv,result%model)
      result%effects=effects;result%joint_covariance=sigma;result%joint_precision=jctx%precision
      result%parameter_covariance=hinv;result%levels=jctx%levels;result%active_parameter_ids=jctx%active_ids
      result%active=act;result%q=q;result%status=0;result%model%converged=result%converged
      call clear_jctx()
   end subroutine fit_gamlss_joint_random_effects

   real(dp) function joint_re_objective(theta) result(f)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable :: eta(:,:),bvec(:)
      real(dp) :: e1,e2,e3,e4,a,b,c,d,lp,pen
      integer :: i,j,g
      if(size(theta)/=jctx%ntheta)then;f=huge(1.0_dp)/100.0_dp;return;end if
      allocate(eta(jctx%n,jctx%np),bvec(jctx%r));call compute_eta(theta,eta)
      f=0.0_dp
      do i=1,jctx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(jctx%np>=2)e2=eta(i,2)
         if(jctx%np>=3)e3=eta(i,3)
         if(jctx%np>=4)e4=eta(i,4)
         call map_parameters(jctx%family,e1,e2,e3,e4,a,b,c,d)
         lp=family_logpdf(jctx%family,jctx%y(i),a,b,c,d)
         if(.not.(lp>-huge(1.0_dp)/10.0_dp))then;f=huge(1.0_dp)/100.0_dp;return;end if
         f=f-jctx%weights(i)*lp
      end do
      pen=0.0_dp
      do g=1,jctx%ng
         call group_effect_vector(theta,g,bvec)
         pen=pen+dot_product(bvec,matmul(jctx%precision,bvec))
      end do
      f=f+0.5_dp*pen
   end function joint_re_objective

   subroutine compute_eta(theta,eta)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: eta(:,:)
      integer :: i,j,g,a,idx
      do j=1,jctx%np
         eta(:,j)=matmul(jctx%b(j)%x,theta(jctx%b(j)%start_fix: &
            jctx%b(j)%start_fix+jctx%b(j)%pfix-1))+jctx%b(j)%offset
         if(.not.jctx%active(j))cycle
         do i=1,jctx%n
            g=jctx%group_index(i)
            do a=1,jctx%q
               idx=theta_random_index(j,g,a)
               eta(i,j)=eta(i,j)+jctx%z(i,a,j)*theta(idx)
            end do
         end do
      end do
   end subroutine compute_eta

   subroutine group_effect_vector(theta,g,v)
      real(dp),intent(in) :: theta(:)
      integer,intent(in) :: g
      real(dp),intent(out) :: v(:)
      integer :: s,a,j,pos
      pos=0
      do s=1,jctx%nact
         j=jctx%active_ids(s)
         do a=1,jctx%q
            pos=pos+1;v(pos)=theta(theta_random_index(j,g,a))
         end do
      end do
   end subroutine group_effect_vector

   integer function theta_random_index(j,g,a) result(idx)
      integer,intent(in) :: j,g,a
      idx=jctx%b(j)%start_re+(g-1)*jctx%q+a-1
   end function theta_random_index

   subroutine extract_effects(theta,effects)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: effects(:,:,:)
      integer :: j,g,a
      effects=0.0_dp
      do j=1,jctx%np
         if(.not.jctx%active(j))cycle
         do g=1,jctx%ng;do a=1,jctx%q
            effects(g,a,j)=theta(theta_random_index(j,g,a))
         end do;end do
      end do
   end subroutine extract_effects

   subroutine fill_joint_model(theta,hinv,model)
      real(dp),intent(in) :: theta(:),hinv(:,:)
      type(gamlss_result_t),intent(out) :: model
      real(dp),allocatable :: eta(:,:),aa(:),bb(:),cc(:),dd(:),lp(:)
      integer :: j,i,lo,hi
      real(dp) :: e1,e2,e3,e4
      allocate(eta(jctx%n,jctx%np),aa(jctx%n),bb(jctx%n),cc(jctx%n),dd(jctx%n),lp(jctx%n))
      call compute_eta(theta,eta)
      do j=1,jctx%np
         lo=jctx%b(j)%start_fix;hi=lo+jctx%b(j)%pfix-1
         select case(j)
         case(1);model%mu%coefficients=theta(lo:hi);model%mu%eta=eta(:,j);model%mu%covariance=hinv(lo:hi,lo:hi)
         case(2);model%sigma%coefficients=theta(lo:hi);model%sigma%eta=eta(:,j);model%sigma%covariance=hinv(lo:hi,lo:hi)
         case(3);model%nu%coefficients=theta(lo:hi);model%nu%eta=eta(:,j);model%nu%covariance=hinv(lo:hi,lo:hi)
         case(4);model%tau%coefficients=theta(lo:hi);model%tau%eta=eta(:,j);model%tau%covariance=hinv(lo:hi,lo:hi)
         end select
      end do
      do i=1,jctx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(jctx%np>=2)e2=eta(i,2)
         if(jctx%np>=3)e3=eta(i,3)
         if(jctx%np>=4)e4=eta(i,4)
         call map_parameters(jctx%family,e1,e2,e3,e4,aa(i),bb(i),cc(i),dd(i))
         lp(i)=family_logpdf(jctx%family,jctx%y(i),aa(i),bb(i),cc(i),dd(i))
      end do
      model%mu%fitted=aa;if(jctx%np>=2)model%sigma%fitted=bb
      if(jctx%np>=3)model%nu%fitted=cc;if(jctx%np>=4)model%tau%fitted=dd
      allocate(model%residuals(jctx%n),model%case_deviance(jctx%n))
      model%residuals=jctx%y-aa;model%case_deviance=-2.0_dp*lp
      model%global_deviance=-2.0_dp*sum(jctx%weights*lp);model%penalized_deviance=model%global_deviance
      model%df_fit=real(sum([(jctx%b(j)%pfix,j=1,jctx%np)]),dp)
      model%df_residual=sum(jctx%weights)-model%df_fit
      model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(max(1.0_dp,sum(jctx%weights)))*model%df_fit
      model%family=jctx%family;model%status=0;model%converged=.true.
   end subroutine fill_joint_model

   subroutine make_layout(initial,theta)
      type(gamlss_result_t),intent(in) :: initial
      real(dp),allocatable,intent(out) :: theta(:)
      real(dp),allocatable :: coef(:)
      integer :: j,pos,nall
      nall=0
      do j=1,jctx%np
         coef=parameter_coefficients(initial,j);jctx%b(j)%pfix=size(coef)
         nall=nall+jctx%b(j)%pfix
         if(jctx%active(j))nall=nall+jctx%ng*jctx%q
      end do
      jctx%ntheta=nall;allocate(theta(nall));theta=0.0_dp;pos=1
      do j=1,jctx%np
         coef=parameter_coefficients(initial,j);jctx%b(j)%start_fix=pos
         theta(pos:pos+size(coef)-1)=coef;pos=pos+size(coef)
         if(jctx%active(j))then;jctx%b(j)%start_re=pos;pos=pos+jctx%ng*jctx%q;end if
      end do
   end subroutine make_layout

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
         index(i)=0;do j=1,nlev;if(group(i)==levels(j))then;index(i)=j;exit;end if;end do
         if(index(i)==0)then;status=1;return;end if
      end do
   end subroutine make_levels

   subroutine make_designs(n,np,xmu,xs,xn,xt,b,status)
      integer,intent(in) :: n,np
      real(dp),intent(in) :: xmu(:,:)
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:)
      type(jr_block_t),intent(out) :: b(4)
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
      type(jr_block_t),intent(inout) :: b(4)
      integer,intent(out) :: status
      integer :: j
      status=0
      do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
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
         call cholesky_factor(a,l,status)
         if(status==0)return
         jitter=merge(1.0e-8_dp,10.0_dp*jitter,jitter==0.0_dp)
         do i=1,size(a,1);a(i,i)=a(i,i)+jitter;end do
         if(allocated(l))deallocate(l)
      end do
   end subroutine stabilize_spd

   subroutine clear_jctx()
      integer :: j
      if(allocated(jctx%y))deallocate(jctx%y)
      if(allocated(jctx%weights))deallocate(jctx%weights)
      if(allocated(jctx%z))deallocate(jctx%z)
      if(allocated(jctx%group))deallocate(jctx%group)
      if(allocated(jctx%levels))deallocate(jctx%levels)
      if(allocated(jctx%group_index))deallocate(jctx%group_index)
      if(allocated(jctx%active_ids))deallocate(jctx%active_ids)
      if(allocated(jctx%precision))deallocate(jctx%precision)
      do j=1,4
         if(allocated(jctx%b(j)%x))deallocate(jctx%b(j)%x)
         if(allocated(jctx%b(j)%offset))deallocate(jctx%b(j)%offset)
         jctx%b(j)%pfix=0;jctx%b(j)%start_fix=0;jctx%b(j)%start_re=0;jctx%b(j)%active_slot=0
      end do
      jctx%family=0;jctx%np=0;jctx%n=0;jctx%q=0;jctx%ng=0;jctx%nact=0;jctx%r=0;jctx%ntheta=0
      jctx%active=.false.
   end subroutine clear_jctx
end module gamlss_joint_random_v06
