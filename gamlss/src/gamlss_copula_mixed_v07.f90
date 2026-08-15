! Gaussian-copula likelihoods for continuous, discrete, and mixed atomic margins.
! Discrete/atomic coordinates contribute latent-normal rectangle probabilities;
! continuous coordinates contribute exact conditional Gaussian densities.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_copula_mixed_v07
   use gamlss_kinds, only : dp
   use gamlss_types, only : gamlss_result_t,gamlss_control_t,GAMLSS_METHOD_RS
   use gamlss_core, only : fit_gamlss_model
   use gamlss_fit, only : family_npar,map_parameters,family_logpdf
   use gamlss_family_support, only : family_cdf,family_cdf_left,family_is_discrete, &
      family_observation_is_atom
   use gamlss_special, only : normal_quantile,normal_pdf
   use gamlss_optim, only : bfgs_minimize
   use gamlss_mvn_v07, only : mvn_logpdf,mvn_rectangle_probability,mvn_conditional
   use nlme_types, only : correlation_spec,COR_NONE,COR_UNSTRUCTURED
   use nlme_gls, only : initialize_correlation
   use nlme_correlation, only : correlation_matrix,correlation_to_unconstrained, &
      correlation_from_unconstrained
   use nlme_linalg, only : find_group_indices,unique_integers
   use nlme_status, only : NLME_SUCCESS
   implicit none
   private
   public :: gaussian_copula_mixed_result_t,fit_gamlss_gaussian_copula_mixed

   type :: design_block_t
      real(dp),allocatable :: x(:,:),offset(:)
      integer :: p=0
   end type design_block_t

   type,public :: gaussian_copula_mixed_result_t
      type(gamlss_result_t) :: model
      type(correlation_spec) :: correlation
      real(dp),allocatable :: correlation_parameters(:),parameter_covariance(:,:)
      real(dp),allocatable :: group_log_likelihood(:)
      real(dp) :: joint_log_likelihood=-huge(1.0_dp),aic=huge(1.0_dp)
      integer :: qmc_points=0,optimizer_status=0,status=0
      logical :: converged=.false.
   end type gaussian_copula_mixed_result_t

   type :: mixed_context_t
      real(dp),allocatable :: y(:),time(:),coordinates(:,:)
      integer,allocatable :: group(:),levels(:)
      type(design_block_t) :: b(4)
      type(correlation_spec) :: corr_template
      integer :: family=0,np=0,n=0,nmax=0,kcorr=0,ncoef=0,nqmc=1024
      logical :: estimate_corr=.false.,has_coordinates=.false.
   end type mixed_context_t
   type(mixed_context_t),save :: mctx
contains

   subroutine fit_gamlss_gaussian_copula_mixed(y,x_mu,family,result,correlation,x_sigma,x_nu,x_tau, &
      offset_mu,offset_sigma,offset_nu,offset_tau,time,group,coordinates,control,max_iter,tolerance,n_qmc)
      real(dp),intent(in) :: y(:),x_mu(:,:)
      integer,intent(in) :: family
      type(gaussian_copula_mixed_result_t),intent(out) :: result
      type(correlation_spec),intent(in),optional :: correlation
      real(dp),intent(in),optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:)
      real(dp),intent(in),optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      real(dp),intent(in),optional :: time(:),coordinates(:,:)
      integer,intent(in),optional :: group(:),max_iter,n_qmc
      type(gamlss_control_t),intent(in),optional :: control
      real(dp),intent(in),optional :: tolerance
      type(gamlss_result_t) :: initial
      type(gamlss_control_t) :: ctl
      type(correlation_spec) :: corr
      real(dp),allocatable :: theta(:),cx(:),hinv(:,:),gll(:)
      integer :: n,np,j,pos,kc,istat,niter
      real(dp) :: fval,tol,jll

      call clear_context()
      n=size(y);np=family_npar(family)
      if(n<2.or.np<1.or.np>4.or.size(x_mu,1)/=n)then;result%status=1;return;end if
      if(family_is_discrete(family))then
         if(any(abs(y-real(nint(y),dp))>1.0e-9_dp))then;result%status=2;return;end if
      end if
      mctx%n=n;mctx%np=np;mctx%family=family;mctx%y=y
      if(present(n_qmc))mctx%nqmc=max(64,n_qmc)
      call make_designs(n,np,x_mu,x_sigma,x_nu,x_tau,mctx%b,istat)
      if(istat/=0)then;result%status=3;return;end if
      call make_offsets(n,np,offset_mu,offset_sigma,offset_nu,offset_tau,mctx%b,istat)
      if(istat/=0)then;result%status=4;return;end if
      allocate(mctx%time(n),mctx%group(n));mctx%time=[(real(j,dp),j=1,n)];mctx%group=1
      if(present(time))then
         if(size(time)/=n)then;result%status=5;return;end if
         mctx%time=time
      end if
      if(present(group))then
         if(size(group)/=n)then;result%status=6;return;end if
         mctx%group=group
      end if
      if(present(coordinates))then
         if(size(coordinates,1)/=n)then;result%status=7;return;end if
         mctx%coordinates=coordinates;mctx%has_coordinates=.true.
      end if
      call unique_integers(mctx%group,mctx%levels,istat)
      if(istat/=NLME_SUCCESS)then;result%status=8;return;end if
      mctx%nmax=0
      do j=1,size(mctx%levels);mctx%nmax=max(mctx%nmax,count(mctx%group==mctx%levels(j)));end do
      corr=correlation_spec();if(present(correlation))corr=correlation
      if(mctx%has_coordinates)then
         call initialize_correlation(corr,mctx%time,mctx%group,mctx%coordinates)
      else
         call initialize_correlation(corr,mctx%time,mctx%group)
      end if
      if(corr%kind==COR_UNSTRUCTURED.and.size(mctx%levels)>1)then
         do j=1,size(mctx%levels)
            if(count(mctx%group==mctx%levels(j))/=mctx%nmax)then;result%status=9;return;end if
         end do
      end if
      mctx%corr_template=corr;mctx%estimate_corr=.not.corr%fixed.and.corr%kind/=COR_NONE
      call correlation_to_unconstrained(corr,mctx%nmax,cx,istat)
      if(istat/=NLME_SUCCESS)then;result%status=10;return;end if
      kc=size(cx);mctx%kcorr=merge(kc,0,mctx%estimate_corr)
      ctl=gamlss_control_t();if(present(control))ctl=control
      call fit_initial(y,x_mu,family,initial,x_sigma,x_nu,x_tau,offset_mu,offset_sigma, &
         offset_nu,offset_tau,ctl)
      if(initial%status/=0)then;result%status=20+initial%status;return;end if
      mctx%ncoef=0
      do j=1,np;mctx%b(j)%p=size(parameter_coefficients(initial,j));mctx%ncoef=mctx%ncoef+mctx%b(j)%p;end do
      allocate(theta(mctx%ncoef+mctx%kcorr));pos=0
      do j=1,np
         theta(pos+1:pos+mctx%b(j)%p)=parameter_coefficients(initial,j);pos=pos+mctx%b(j)%p
      end do
      if(mctx%estimate_corr.and.kc>0)theta(mctx%ncoef+1:)=cx
      allocate(hinv(size(theta),size(theta)));hinv=0.0_dp
      niter=120;if(present(max_iter))niter=max(1,max_iter)
      tol=3.0e-6_dp;if(present(tolerance))tol=max(1.0e-10_dp,tolerance)
      call bfgs_minimize(mixed_objective,theta,fval,istat,max_iter=niter,tol=tol,inverse_hessian=hinv)
      result%optimizer_status=istat
      if(istat/=0.and.istat/=1)then;result%status=30+istat;call clear_context();return;end if
      call evaluate_theta(theta,jll,gll,corr,istat)
      if(istat/=0)then;result%status=40+istat;call clear_context();return;end if
      call fill_model(theta,hinv,jll,result%model)
      result%correlation=corr
      if(allocated(corr%par))then;result%correlation_parameters=corr%par
      else;allocate(result%correlation_parameters(0));end if
      result%parameter_covariance=hinv;result%group_log_likelihood=gll
      result%joint_log_likelihood=jll;result%aic=-2.0_dp*jll+2.0_dp*real(size(theta),dp)
      result%qmc_points=mctx%nqmc;result%status=0;result%converged=(istat==0)
      result%model%converged=result%converged
      call clear_context()
   end subroutine fit_gamlss_gaussian_copula_mixed

   real(dp) function mixed_objective(theta) result(f)
      real(dp),intent(in) :: theta(:)
      real(dp),allocatable :: gll(:)
      type(correlation_spec) :: corr
      real(dp) :: ll
      integer :: status
      call evaluate_theta(theta,ll,gll,corr,status)
      if(status/=0.or..not.(ll>-huge(1.0_dp)/10.0_dp))then;f=huge(1.0_dp)/100.0_dp
      else;f=-ll;end if
   end function mixed_objective

   subroutine evaluate_theta(theta,ll,gll,corr,status)
      real(dp),intent(in) :: theta(:)
      real(dp),intent(out) :: ll
      real(dp),allocatable,intent(out) :: gll(:)
      type(correlation_spec),intent(out) :: corr
      integer,intent(out) :: status
      real(dp),allocatable :: eta(:,:),aa(:),bb(:),cc(:),dd(:),tl(:),r(:,:),cl(:,:)
      integer,allocatable :: idx(:)
      integer :: i,j,g,pos,istat
      real(dp) :: e2,e3,e4
      status=0;ll=0.0_dp
      if(size(theta)/=mctx%ncoef+mctx%kcorr)then;status=1;return;end if
      allocate(eta(mctx%n,mctx%np),aa(mctx%n),bb(mctx%n),cc(mctx%n),dd(mctx%n))
      pos=0
      do j=1,mctx%np
         eta(:,j)=matmul(mctx%b(j)%x,theta(pos+1:pos+mctx%b(j)%p))+mctx%b(j)%offset
         pos=pos+mctx%b(j)%p
      end do
      if(mctx%estimate_corr)then
         call correlation_from_unconstrained(mctx%corr_template,mctx%nmax,theta(mctx%ncoef+1:),corr,istat)
         if(istat/=NLME_SUCCESS)then;status=2;return;end if
      else;corr=mctx%corr_template;end if
      do i=1,mctx%n
         e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(mctx%np>=2)e2=eta(i,2);if(mctx%np>=3)e3=eta(i,3);if(mctx%np>=4)e4=eta(i,4)
         call map_parameters(mctx%family,eta(i,1),e2,e3,e4,aa(i),bb(i),cc(i),dd(i))
      end do
      allocate(gll(size(mctx%levels)));gll=0.0_dp
      do g=1,size(mctx%levels)
         call find_group_indices(mctx%group,mctx%levels(g),idx)
         allocate(tl(size(idx)));tl=mctx%time(idx)
         if(mctx%has_coordinates)then
            allocate(cl(size(idx),size(mctx%coordinates,2)));cl=mctx%coordinates(idx,:)
            call correlation_matrix(corr,tl,r,istat,cl);deallocate(cl)
         else;call correlation_matrix(corr,tl,r,istat);end if
         if(istat/=NLME_SUCCESS)then;status=3;return;end if
         call group_mixed_loglik(idx,r,aa,bb,cc,dd,gll(g),istat)
         if(istat/=0)then;status=10+istat;return;end if
         ll=ll+gll(g);deallocate(idx,tl,r)
      end do
   end subroutine evaluate_theta

   subroutine group_mixed_loglik(idx,r,aa,bb,cc,dd,ll,status)
      integer,intent(in) :: idx(:)
      real(dp),intent(in) :: r(:,:),aa(:),bb(:),cc(:),dd(:)
      real(dp),intent(out) :: ll
      integer,intent(out) :: status
      logical,allocatable :: atom(:)
      integer,allocatable :: ic(:),id(:)
      real(dp),allocatable :: zc(:),lo(:),hi(:),md(:),cd(:,:),rcc(:,:)
      real(dp) :: lp,u,ul,prob,lpc,logphi
      integer :: i,j,n,nc,nd,istat
      n=size(idx);status=0;ll=0.0_dp;allocate(atom(n))
      do i=1,n;atom(i)=family_observation_is_atom(mctx%family,mctx%y(idx(i)));end do
      nc=count(.not.atom);nd=count(atom);allocate(ic(nc),id(nd));nc=0;nd=0
      do i=1,n
         if(atom(i))then;nd=nd+1;id(nd)=i;else;nc=nc+1;ic(nc)=i;end if
      end do
      allocate(zc(nc),lo(nd),hi(nd));ll=0.0_dp
      do i=1,nc
         j=idx(ic(i));lp=family_logpdf(mctx%family,mctx%y(j),aa(j),bb(j),cc(j),dd(j))
         u=family_cdf(mctx%family,mctx%y(j),aa(j),bb(j),cc(j),dd(j))
         if(.not.(lp>-huge(1.0_dp)/10.0_dp).or.u<0.0_dp)then;status=1;return;end if
         zc(i)=normal_quantile(min(1.0_dp-1.0e-11_dp,max(1.0e-11_dp,u)))
         logphi=log(max(normal_pdf(zc(i)),tiny(1.0_dp)));ll=ll+lp-logphi
      end do
      do i=1,nd
         j=idx(id(i));ul=family_cdf_left(mctx%family,mctx%y(j),aa(j),bb(j),cc(j),dd(j))
         u=family_cdf(mctx%family,mctx%y(j),aa(j),bb(j),cc(j),dd(j))
         if(ul<0.0_dp.or.u<0.0_dp.or.u<=ul)then;status=2;return;end if
         lo(i)=cdf_to_bound(ul,-1);hi(i)=cdf_to_bound(u,1)
      end do
      if(nc>0)then
         allocate(rcc(nc,nc))
         do i=1,nc;do j=1,nc;rcc(i,j)=r(ic(i),ic(j));end do;end do
         lpc=mvn_logpdf(zc,0.0_dp*zc,rcc,istat)
         if(istat/=0)then;status=3;return;end if
         ll=ll+lpc
      end if
      if(nd>0)then
         call mvn_conditional(r,ic,id,zc,md,cd,istat)
         if(istat/=0)then;status=4;return;end if
         prob=mvn_rectangle_probability(lo,hi,md,cd,istat,n_qmc=mctx%nqmc)
         if(istat/=0.or.prob<=0.0_dp)then;status=5;return;end if
         ll=ll+log(max(prob,tiny(1.0_dp)))
      end if
   end subroutine group_mixed_loglik

   real(dp) function cdf_to_bound(p,side) result(z)
      real(dp),intent(in) :: p
      integer,intent(in) :: side
      if(p<=0.0_dp)then;z=-huge(1.0_dp)
      else if(p>=1.0_dp)then;z=huge(1.0_dp)
      else;z=normal_quantile(p);end if
      if(side==0)z=z
   end function cdf_to_bound

   subroutine fill_model(theta,hinv,jll,model)
      real(dp),intent(in) :: theta(:),hinv(:,:),jll
      type(gamlss_result_t),intent(out) :: model
      real(dp),allocatable :: eta(:,:),a(:),b(:),c(:),d(:),lp(:)
      integer :: j,pos,i,p,lo,hi
      real(dp) :: e1,e2,e3,e4
      allocate(eta(mctx%n,mctx%np),a(mctx%n),b(mctx%n),c(mctx%n),d(mctx%n),lp(mctx%n));pos=0
      do j=1,mctx%np
         p=mctx%b(j)%p;lo=pos+1;hi=pos+p
         eta(:,j)=matmul(mctx%b(j)%x,theta(lo:hi))+mctx%b(j)%offset
         select case(j)
         case(1);model%mu%coefficients=theta(lo:hi);model%mu%eta=eta(:,j);model%mu%covariance=hinv(lo:hi,lo:hi)
         case(2);model%sigma%coefficients=theta(lo:hi);model%sigma%eta=eta(:,j);model%sigma%covariance=hinv(lo:hi,lo:hi)
         case(3);model%nu%coefficients=theta(lo:hi);model%nu%eta=eta(:,j);model%nu%covariance=hinv(lo:hi,lo:hi)
         case(4);model%tau%coefficients=theta(lo:hi);model%tau%eta=eta(:,j);model%tau%covariance=hinv(lo:hi,lo:hi)
         end select
         pos=hi
      end do
      do i=1,mctx%n
         e1=eta(i,1);e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(mctx%np>=2)e2=eta(i,2);if(mctx%np>=3)e3=eta(i,3);if(mctx%np>=4)e4=eta(i,4)
         call map_parameters(mctx%family,e1,e2,e3,e4,a(i),b(i),c(i),d(i))
         lp(i)=family_logpdf(mctx%family,mctx%y(i),a(i),b(i),c(i),d(i))
      end do
      model%mu%fitted=a;if(mctx%np>=2)model%sigma%fitted=b
      if(mctx%np>=3)model%nu%fitted=c;if(mctx%np>=4)model%tau%fitted=d
      allocate(model%residuals(mctx%n),model%case_deviance(mctx%n))
      model%residuals=mctx%y-a;model%case_deviance=-2.0_dp*lp
      model%global_deviance=-2.0_dp*jll;model%penalized_deviance=model%global_deviance
      model%df_fit=real(size(theta),dp);model%df_residual=real(mctx%n,dp)-model%df_fit
      model%aic=model%global_deviance+2.0_dp*model%df_fit
      model%sbc=model%global_deviance+log(real(max(1,mctx%n),dp))*model%df_fit
      model%family=mctx%family;model%status=0;model%converged=.true.
   end subroutine fill_model

   subroutine fit_initial(y,xmu,family,fit,xs,xn,xt,om,os,on,ot,ctl)
      real(dp),intent(in) :: y(:),xmu(:,:)
      integer,intent(in) :: family
      type(gamlss_result_t),intent(out) :: fit
      real(dp),intent(in),optional :: xs(:,:),xn(:,:),xt(:,:),om(:),os(:),on(:),ot(:)
      type(gamlss_control_t),intent(in) :: ctl
      call fit_gamlss_model(y,xmu,family,fit,method=GAMLSS_METHOD_RS,x_sigma=xs,x_nu=xn,x_tau=xt, &
         offset_mu=om,offset_sigma=os,offset_nu=on,offset_tau=ot,control=ctl)
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
      status=0;do j=1,np;allocate(b(j)%offset(n));b(j)%offset=0.0_dp;end do
      if(present(om))then;if(size(om)/=n)then;status=1;return;end if;b(1)%offset=om;end if
      if(np>=2.and.present(os))then;if(size(os)/=n)then;status=2;return;end if;b(2)%offset=os;end if
      if(np>=3.and.present(on))then;if(size(on)/=n)then;status=3;return;end if;b(3)%offset=on;end if
      if(np>=4.and.present(ot))then;if(size(ot)/=n)then;status=4;return;end if;b(4)%offset=ot;end if
   end subroutine make_offsets

   subroutine clear_context()
      integer :: j
      if(allocated(mctx%y))deallocate(mctx%y);if(allocated(mctx%time))deallocate(mctx%time)
      if(allocated(mctx%coordinates))deallocate(mctx%coordinates);if(allocated(mctx%group))deallocate(mctx%group)
      if(allocated(mctx%levels))deallocate(mctx%levels)
      do j=1,4
         if(allocated(mctx%b(j)%x))deallocate(mctx%b(j)%x)
         if(allocated(mctx%b(j)%offset))deallocate(mctx%b(j)%offset)
         mctx%b(j)%p=0
      end do
      mctx%family=0;mctx%np=0;mctx%n=0;mctx%nmax=0;mctx%kcorr=0;mctx%ncoef=0;mctx%nqmc=1024
      mctx%estimate_corr=.false.;mctx%has_coordinates=.false.
   end subroutine clear_context
end module gamlss_copula_mixed_v07
