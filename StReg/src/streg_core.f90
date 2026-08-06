! SPDX-License-Identifier: GPL-2.0-only
module streg_core
   use streg_kinds, only : dp, pi
   use streg_linalg, only : expand_vech, vech, spd_inverse_logdet, general_inverse, &
      symmetric_sqrt, least_squares, covariance_matrix, all_finite
   use streg_probability, only : anderson_darling_t, student_t_two_sided_p
   use streg_optimize, only : bfgs_minimize, nelder_mead, numerical_hessian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   integer, parameter, public :: streg_success = 0
   integer, parameter, public :: streg_invalid_input = 1
   integer, parameter, public :: streg_covariance_failure = 2
   integer, parameter, public :: streg_hessian_failure = 3

   type, public :: streg_options
      integer :: max_iter = 1000
      real(dp) :: tolerance = 1.0e-8_dp
      logical :: compute_hessian = .false.
      logical :: source_compatible = .false.
      logical :: enforce_sample_rule = .true.
      character(len=16) :: optimizer = 'bfgs'
   end type streg_options


   real(dp), allocatable, save :: objective_z(:,:), objective_trend_stack(:,:)
   integer, save :: objective_l=0, objective_c=0, objective_lag=0
   integer, save :: objective_d=0, objective_neff=0, objective_ncov=0
   real(dp), save :: objective_df=1.0_dp
   logical, save :: objective_state_ready=.false.

   type, public :: streg_fit
      integer :: status = streg_success
      character(len=240) :: message = ''
      character(len=16) :: model = ''
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: nobs = 0
      integer :: nvar = 0
      integer :: nresponse = 0
      integer :: ntrend = 0
      integer :: lag = 0
      real(dp) :: df = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: joint_scale(:,:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: innovation_scale(:,:)
      real(dp), allocatable :: precision_predictors(:,:)
      real(dp), allocatable :: conditional_factor(:)
      real(dp), allocatable :: trend(:,:)
      real(dp), allocatable :: fitted(:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: r_squared(:)
      real(dp), allocatable :: f_statistic(:)
      real(dp), allocatable :: ad_test(:,:)
      real(dp), allocatable :: var_coef(:)
      real(dp), allocatable :: coef_estimate(:)
      real(dp), allocatable :: coef_se(:)
      real(dp), allocatable :: coef_p_value(:)
      real(dp), allocatable :: var_coef_se(:)
      real(dp), allocatable :: var_coef_p_value(:)
      real(dp), allocatable :: hessian(:,:)
      real(dp), allocatable :: theta_covariance(:,:)
   end type streg_fit

   interface fit_student_lm
      module procedure stlm
   end interface
   interface fit_student_ar
      module procedure star
   end interface
   interface fit_student_dlm
      module procedure stdlm
   end interface
   interface fit_student_var
      module procedure stvar
   end interface

   public :: stlm, star, stdlm, stvar
   public :: fit_student_lm, fit_student_ar, fit_student_dlm, fit_student_var
   public :: block_top_covariance, conditional_parameters

   abstract interface
      function objective_interface(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_interface
   end interface

contains

   function stlm(y,x,trend,v,options,include_intercept,init) result(fit)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: trend(:,:)
      real(dp), intent(in), optional :: v
      type(streg_options), intent(in), optional :: options
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: init(:)
      type(streg_fit) :: fit
      real(dp), allocatable :: data(:,:),tw(:,:)
      real(dp) :: df
      type(streg_options) :: opt
      logical :: intercept,ok
      character(len=240) :: msg
      df=1.0_dp; if(present(v))df=v
      opt=streg_options(); if(present(options))opt=options
      intercept=.true.; if(present(include_intercept))intercept=include_intercept
      if(size(y)/=size(x,1) .or. size(x,2)<1)then
         call set_error(fit,'stlm','y and x must have equal rows and x must have at least one column.')
         return
      end if
      allocate(data(size(y),size(x,2)+1)); data(:,1)=y; data(:,2:)=x
      call prepare_trend(size(y),trend,intercept,tw,ok,msg)
      if(.not.ok)then; call set_error(fit,'stlm',msg); return; end if
      if(present(init))then
         call fit_joint_model(data,tw,0,1,df,opt,'stlm',fit,init)
      else
         call fit_joint_model(data,tw,0,1,df,opt,'stlm',fit)
      end if
   end function stlm

   function star(y,lag,trend,v,options,include_intercept,init) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: trend(:,:)
      real(dp), intent(in), optional :: v
      type(streg_options), intent(in), optional :: options
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: init(:)
      type(streg_fit) :: fit
      real(dp), allocatable :: data(:,:),tw(:,:)
      real(dp) :: df
      integer :: p
      type(streg_options) :: opt
      logical :: intercept,ok
      character(len=240) :: msg
      p=1; if(present(lag))p=lag
      df=1.0_dp; if(present(v))df=v
      opt=streg_options(); if(present(options))opt=options
      intercept=.true.; if(present(include_intercept))intercept=include_intercept
      if(p<1)then; call set_error(fit,'star','lag must be a positive integer.'); return; end if
      allocate(data(size(y),1)); data(:,1)=y
      call prepare_trend(size(y),trend,intercept,tw,ok,msg)
      if(.not.ok)then; call set_error(fit,'star',msg); return; end if
      if(present(init))then
         call fit_joint_model(data,tw,p,1,df,opt,'star',fit,init)
      else
         call fit_joint_model(data,tw,p,1,df,opt,'star',fit)
      end if
   end function star

   function stdlm(y,x,lag,trend,v,options,include_intercept,init) result(fit)
      real(dp), intent(in) :: y(:),x(:,:)
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: trend(:,:)
      real(dp), intent(in), optional :: v
      type(streg_options), intent(in), optional :: options
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: init(:)
      type(streg_fit) :: fit
      real(dp), allocatable :: data(:,:),tw(:,:)
      real(dp) :: df
      integer :: p
      type(streg_options) :: opt
      logical :: intercept,ok
      character(len=240) :: msg
      p=1; if(present(lag))p=lag
      df=1.0_dp; if(present(v))df=v
      opt=streg_options(); if(present(options))opt=options
      intercept=.true.; if(present(include_intercept))intercept=include_intercept
      if(p<1)then; call set_error(fit,'stdlm','lag must be a positive integer.'); return; end if
      if(size(y)/=size(x,1) .or. size(x,2)<1)then
         call set_error(fit,'stdlm','y and x must have equal rows and x must have at least one column.')
         return
      end if
      allocate(data(size(y),size(x,2)+1)); data(:,1)=y; data(:,2:)=x
      call prepare_trend(size(y),trend,intercept,tw,ok,msg)
      if(.not.ok)then; call set_error(fit,'stdlm',msg); return; end if
      if(present(init))then
         call fit_joint_model(data,tw,p,1,df,opt,'stdlm',fit,init)
      else
         call fit_joint_model(data,tw,p,1,df,opt,'stdlm',fit)
      end if
   end function stdlm

   function stvar(data,lag,trend,v,options,include_intercept,init) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: trend(:,:)
      real(dp), intent(in), optional :: v
      type(streg_options), intent(in), optional :: options
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: init(:)
      type(streg_fit) :: fit
      real(dp), allocatable :: tw(:,:)
      real(dp) :: df
      integer :: p
      type(streg_options) :: opt
      logical :: intercept,ok
      character(len=240) :: msg
      p=1; if(present(lag))p=lag
      df=1.0_dp; if(present(v))df=v
      opt=streg_options(); if(present(options))opt=options
      intercept=.true.; if(present(include_intercept))intercept=include_intercept
      if(p<1)then; call set_error(fit,'stvar','lag must be a positive integer.'); return; end if
      if(size(data,2)<1)then; call set_error(fit,'stvar','data must have at least one column.'); return; end if
      call prepare_trend(size(data,1),trend,intercept,tw,ok,msg)
      if(.not.ok)then; call set_error(fit,'stvar',msg); return; end if
      if(present(init))then
         call fit_joint_model(data,tw,p,size(data,2),df,opt,'stvar',fit,init)
      else
         call fit_joint_model(data,tw,p,size(data,2),df,opt,'stvar',fit)
      end if
   end function stvar

   subroutine prepare_trend(n,trend_in,include_intercept,trend,ok,message)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: trend_in(:,:)
      logical, intent(in) :: include_intercept
      real(dp), allocatable, intent(out) :: trend(:,:)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: message
      ok=.true.; message=''
      if(present(trend_in))then
         if(size(trend_in,1)/=n .or. size(trend_in,2)<1)then
            ok=.false.; message='data and trend must have equal rows and trend must have at least one column.'
            allocate(trend(n,0)); return
         end if
         allocate(trend(n,size(trend_in,2))); trend=trend_in
      else if(include_intercept)then
         allocate(trend(n,1)); trend=1.0_dp
      else
         allocate(trend(n,0))
      end if
   end subroutine prepare_trend

   subroutine set_error(fit,model,message)
      type(streg_fit), intent(out) :: fit
      character(len=*), intent(in) :: model,message
      fit%status=streg_invalid_input; fit%model=model; fit%message=message
   end subroutine set_error

   subroutine block_top_covariance(a,l,lag,s,status)
      real(dp), intent(in) :: a(:)
      integer, intent(in) :: l,lag
      real(dp), allocatable, intent(out) :: s(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: cmat(:,:,:),seq(:,:,:),bmat(:,:)
      integer :: nmat,nv,m,nseq,i,j,k,r1,r2,c1,c2,local_status
      logical :: ok
      local_status=streg_success
      nmat=lag+3; nv=l*(l+1)/2; m=lag+1; nseq=2*nmat+1
      if(l<1 .or. lag<0 .or. size(a)/=nmat*nv)then
         allocate(s(0,0)); local_status=streg_invalid_input
         if(present(status))status=local_status
         return
      end if
      allocate(cmat(l,l,nmat),seq(l,l,nseq),bmat(m*l,(lag+4)*l))
      do k=1,nmat
         call expand_vech(a((k-1)*nv+1:k*nv),l,cmat(:,:,k),ok)
      end do
      seq(:,:,1:nmat)=cmat
      seq(:,:,nmat+1)=cmat(:,:,1)
      seq(:,:,nmat+2:2*nmat+1)=cmat
      bmat=0.0_dp
      do i=1,m
         r1=(i-1)*l+1; r2=i*l
         do j=1,lag+4
            c1=(j-1)*l+1; c2=j*l
            bmat(r1:r2,c1:c2)=seq(:,:,i+j-1)
         end do
      end do
      s=matmul(bmat,transpose(bmat))
      s=0.5_dp*(s+transpose(s))
      if(present(status))status=local_status
   end subroutine block_top_covariance

   subroutine conditional_parameters(theta,l,lag,response_dim,df,ntrend,s,b,sigma,q,delta0,var_coef,status)
      real(dp), intent(in) :: theta(:),df
      integer, intent(in) :: l,lag,response_dim,ntrend
      real(dp), allocatable, intent(out) :: s(:,:),b(:,:),sigma(:,:),q(:,:),delta0(:,:),var_coef(:)
      integer, intent(out), optional :: status
      logical :: ok
      call decode_theta(theta,l,lag,response_dim,df,ntrend,s,b,sigma,q,delta0,var_coef,ok)
      if(present(status))status=merge(streg_success,streg_covariance_failure,ok)
   end subroutine conditional_parameters

   subroutine fit_joint_model(data,trend_base,lag,r,df,opt,model,fit,init)
      real(dp), intent(in) :: data(:,:),trend_base(:,:),df
      integer, intent(in) :: lag,r
      type(streg_options), intent(in) :: opt
      character(len=*), intent(in) :: model
      type(streg_fit), intent(out) :: fit
      real(dp), intent(in), optional :: init(:)
      real(dp), allocatable :: z(:,:),trend_stack(:,:),theta0(:),theta(:),m(:,:)
      real(dp), allocatable :: s(:,:),b(:,:),sigma(:,:),q(:,:),delta0(:,:),varcoef(:)
      real(dp) :: objective_value
      integer :: n,l,c,d,neff,ntheta,ncov,t,block,row0,col0
      logical :: ok,conv
      character(len=16) :: optimizer_name

      fit=streg_fit(); fit%model=model; fit%lag=lag; fit%df=df
      n=size(data,1); l=size(data,2); c=size(trend_base,2); d=(lag+1)*l; neff=n-lag
      fit%nobs=neff; fit%nvar=l; fit%nresponse=r; fit%ntrend=c
      if(df<=0.0_dp .or. n<=lag+2 .or. r<1 .or. r>d .or. .not.all_finite(data))then
         call set_error(fit,model,'invalid sample size, dimensions, degrees of freedom, or nonfinite data.')
         return
      end if
      if(opt%max_iter<10)then
         call set_error(fit,model,'maximum iterations must be at least 10.')
         return
      end if
      if(c>0)then
         if(.not.all_finite(trend_base))then
            call set_error(fit,model,'trend contains nonfinite values.'); return
         end if
      end if
      ncov=(lag+3)*l*(l+1)/2; ntheta=ncov+l*c
      if(2*n < 3*ncov+3*l+c)then
         if(opt%enforce_sample_rule)then
            call set_error(fit,model,'too many joint-distribution parameters for the sample size.')
            return
         else
            fit%message='warning: the upstream sample-size rule flags many parameters for this sample.'
         end if
      end if
      allocate(z(d,neff),trend_stack(c*(lag+1),neff))
      do t=1,neff
         do block=0,lag
            row0=block*l
            z(row0+1:row0+l,t)=data(t+lag-block,:)
            if(c>0)then
               col0=block*c
               trend_stack(col0+1:col0+c,t)=trend_base(t+lag-block,:)
            end if
         end do
      end do
      if(present(init))then
         if(size(init)/=ntheta)then
            call set_error(fit,model,'initial parameter vector has the wrong length.'); return
         end if
         allocate(theta0(ntheta)); theta0=init
      else
         call initial_theta(data,trend_base,lag,df,theta0,ok)
         if(.not.ok)then; call set_error(fit,model,'could not construct initial parameters.'); return; end if
      end if
      optimizer_name=lower_string(opt%optimizer)
      if(trim(optimizer_name)/='bfgs' .and. trim(optimizer_name)/='nelder-mead' .and. &
         trim(optimizer_name)/='nelder_mead')then
         call set_error(fit,model,'optimizer must be bfgs or nelder-mead.')
         return
      end if
      if(trim(optimizer_name)=='nelder-mead' .or. trim(optimizer_name)=='nelder_mead')then
         call set_objective_state(z,trend_stack,l,c,lag,d,neff,ncov,df)
         call nelder_mead(joint_objective,theta0,theta,objective_value,conv,fit%iterations,opt%max_iter,opt%tolerance)
      else
         call set_objective_state(z,trend_stack,l,c,lag,d,neff,ncov,df)
         call bfgs_minimize(joint_objective,theta0,theta,objective_value,conv,fit%iterations,opt%max_iter,opt%tolerance)
      end if
      fit%converged=conv; fit%log_likelihood=-objective_value; fit%theta=theta
      if(.not.conv)then
         if(len_trim(fit%message)>0)then
            fit%message=trim(fit%message)//' Optimizer stopped before the requested tolerance.'
         else
            fit%message='optimizer stopped before the requested tolerance; returned the best finite iterate.'
         end if
      end if
      call decode_theta(theta,l,lag,r,df,c,s,b,sigma,q,delta0,varcoef,ok)
      if(.not.ok)then
         fit%status=streg_covariance_failure; fit%message='optimized covariance is singular.'
         call clear_objective_state(); return
      end if
      fit%joint_scale=s; fit%innovation_scale=sigma; fit%precision_predictors=q
      fit%var_coef=varcoef
      call mean_matrix(theta,ncov,l,c,lag,trend_stack,m)
      call populate_fit(z,m,r,df,opt%source_compatible,b,sigma,q,delta0,fit)
      call parameter_vectors(theta,l,lag,r,df,c,fit%coef_estimate,fit%var_coef,ok)
      fit%status=streg_success
      if(opt%compute_hessian)then
         call inference(theta,l,lag,r,df,c,neff,joint_objective,fit)
      end if

      call clear_objective_state()

   end subroutine fit_joint_model

   subroutine set_objective_state(z,trend_stack,l,c,lag,d,neff,ncov,df)
      real(dp), intent(in) :: z(:,:),trend_stack(:,:),df
      integer, intent(in) :: l,c,lag,d,neff,ncov
      if(allocated(objective_z))deallocate(objective_z)
      if(allocated(objective_trend_stack))deallocate(objective_trend_stack)
      objective_z=z; objective_trend_stack=trend_stack
      objective_l=l; objective_c=c; objective_lag=lag; objective_d=d
      objective_neff=neff; objective_ncov=ncov; objective_df=df
      objective_state_ready=.true.
   end subroutine set_objective_state

   subroutine clear_objective_state()
      if(allocated(objective_z))deallocate(objective_z)
      if(allocated(objective_trend_stack))deallocate(objective_trend_stack)
      objective_state_ready=.false.
   end subroutine clear_objective_state

   function joint_objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      real(dp), allocatable :: local_s(:,:),local_inv(:,:),local_m(:,:),diff(:)
      real(dp) :: local_logdet,term,sumlog,const
      integer :: it,st
      logical :: local_ok
      if(.not.objective_state_ready .or. any(.not.ieee_is_finite(x)))then
         value=huge(1.0_dp)/100.0_dp; return
      end if
      call block_top_covariance(x(1:objective_ncov),objective_l,objective_lag,local_s,st)
      if(st/=streg_success)then; value=huge(1.0_dp)/100.0_dp; return; end if
      call spd_inverse_logdet(local_s,local_inv,local_logdet,local_ok)
      if(.not.local_ok)then; value=huge(1.0_dp)/100.0_dp+sum(x*x); return; end if
      call mean_matrix(x,objective_ncov,objective_l,objective_c,objective_lag, &
         objective_trend_stack,local_m)
      allocate(diff(objective_d)); sumlog=0.0_dp
      do it=1,objective_neff
         diff=objective_z(:,it)-local_m(:,it)
         term=dot_product(diff,matmul(local_inv,diff))/objective_df
         if(term<=-1.0_dp .or. .not.ieee_is_finite(term))then
            value=huge(1.0_dp)/100.0_dp; return
         end if
         sumlog=sumlog+log(1.0_dp+term)
      end do
      const=log_gamma(0.5_dp*(objective_df+real(objective_d,dp)))- &
         log_gamma(0.5_dp*objective_df)-0.5_dp*real(objective_d,dp)*log(pi*objective_df)
      value=-(real(objective_neff,dp)*const-0.5_dp*real(objective_neff,dp)*local_logdet- &
         0.5_dp*(objective_df+real(objective_d,dp))*sumlog)
   end function joint_objective

   subroutine initial_theta(data,trend,lag,df,theta,ok)
      real(dp), intent(in) :: data(:,:),trend(:,:),df
      integer, intent(in) :: lag
      real(dp), allocatable, intent(out) :: theta(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: coef(:,:),resid(:,:),cov(:,:),root(:,:),vv(:)
      real(dp) :: scale_factor,small
      integer :: l,c,ncov,nv,k,last0
      l=size(data,2); c=size(trend,2); nv=l*(l+1)/2; ncov=(lag+3)*nv
      allocate(theta(ncov+l*c)); theta=0.0_dp
      if(c>0)then
         call least_squares(trend,data,coef,ok)
         if(.not.ok)return
         resid=data-matmul(trend,coef)
         theta(ncov+1:)=reshape(transpose(coef),[l*c])
      else
         allocate(resid(size(data,1),l)); resid=data
      end if
      call covariance_matrix(resid,cov,ok); if(.not.ok)return
      scale_factor=1.0_dp
      if(df>2.0_dp)scale_factor=(df-2.0_dp)/df
      cov=scale_factor*cov
      do k=1,l
         cov(k,k)=cov(k,k)+sqrt(epsilon(1.0_dp))*max(1.0_dp,cov(k,k))
      end do
      call symmetric_sqrt(cov,root,ok); if(.not.ok)return
      small=1.0e-4_dp*max(1.0_dp,maxval(abs(root)))
      do k=1,ncov
         theta(k)=small*sin(real(17*k+3,dp))
      end do
      vv=vech(root); last0=(lag+2)*nv
      theta(last0+1:last0+nv)=vv
      ok=.true.
   end subroutine initial_theta

   subroutine mean_matrix(theta,ncov,l,c,lag,trend_stack,m)
      real(dp), intent(in) :: theta(:),trend_stack(:,:)
      integer, intent(in) :: ncov,l,c,lag
      real(dp), allocatable, intent(out) :: m(:,:)
      real(dp), allocatable :: mu(:,:)
      integer :: block,r0,c0
      allocate(m((lag+1)*l,size(trend_stack,2))); m=0.0_dp
      if(c==0)return
      mu=reshape(theta(ncov+1:ncov+l*c),[l,c])
      do block=0,lag
         r0=block*l; c0=block*c
         m(r0+1:r0+l,:)=matmul(mu,trend_stack(c0+1:c0+c,:))
      end do
   end subroutine mean_matrix

   subroutine decode_theta(theta,l,lag,r,df,c,s,b,sigma,q,delta0,varcoef,ok)
      real(dp), intent(in) :: theta(:),df
      integer, intent(in) :: l,lag,r,c
      real(dp), allocatable, intent(out) :: s(:,:),b(:,:),sigma(:,:),q(:,:),delta0(:,:),varcoef(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: sxx(:,:),invxx(:,:),syx(:,:),mu(:,:),mux(:,:),vv(:)
      real(dp) :: logdet,denom
      integer :: d,p,ncov,g,var,st
      d=(lag+1)*l; p=d-r; ncov=(lag+3)*l*(l+1)/2
      ok=size(theta)==ncov+l*c .and. p>0
      if(.not.ok)then
         allocate(s(0,0),b(0,0),sigma(0,0),q(0,0),delta0(0,0),varcoef(0)); return
      end if
      call block_top_covariance(theta(1:ncov),l,lag,s,st)
      if(st/=streg_success)then; ok=.false.; return; end if
      sxx=s(r+1:d,r+1:d); syx=s(1:r,r+1:d)
      call spd_inverse_logdet(sxx,invxx,logdet,ok); if(.not.ok)return
      allocate(b(p,r)); b=matmul(invxx,transpose(syx))
      sigma=s(1:r,1:r)-matmul(syx,matmul(invxx,transpose(syx)))
      sigma=0.5_dp*(sigma+transpose(sigma)); q=invxx/df
      allocate(delta0(r,c)); delta0=0.0_dp
      if(c>0)then
         mu=reshape(theta(ncov+1:ncov+l*c),[l,c]); allocate(mux(p,c))
         do g=r+1,d
            var=mod(g-1,l)+1
            mux(g-r,:)=mu(var,:)
         end do
         delta0=mu(1:r,:)-matmul(transpose(b),mux)
      end if
      vv=vech(invxx); allocate(varcoef(1+size(vv)))
      denom=df+real(p,dp)-2.0_dp
      if(abs(denom)<=sqrt(epsilon(1.0_dp)))then
         varcoef=huge(1.0_dp)
      else
         varcoef(1)=df*sigma(1,1)/denom
         varcoef(2:)=sigma(1,1)*vv/denom
      end if
   end subroutine decode_theta

   subroutine populate_fit(z,m,r,df,source,b,sigma,q,delta0,fit)
      real(dp), intent(in) :: z(:,:),m(:,:),df
      integer, intent(in) :: r
      logical, intent(in) :: source
      real(dp), intent(in) :: b(:,:),sigma(:,:),q(:,:),delta0(:,:)
      type(streg_fit), intent(inout) :: fit
      real(dp), allocatable :: diff(:),standardized(:)
      real(dp) :: actual_mean,sse,sst,denom,scale,ad_df
      integer :: t,k,p,c,neff,ncoef,dfadd
      p=size(b,1); c=size(delta0,2); neff=size(z,2)
      allocate(fit%beta(r,c+p),fit%conditional_factor(neff),fit%trend(neff,r), &
         fit%fitted(neff,r),fit%residuals(neff,r),fit%r_squared(r), &
         fit%f_statistic(r),fit%ad_test(r,2))
      if(c>0)fit%beta(:,1:c)=delta0
      fit%beta(:,c+1:c+p)=transpose(b)
      allocate(diff(p))
      do t=1,neff
         diff=z(r+1:,t)-m(r+1:,t)
         fit%conditional_factor(t)=1.0_dp+dot_product(diff,matmul(q,diff))
         fit%trend(t,:)=m(1:r,t)-matmul(transpose(b),m(r+1:,t))
         fit%fitted(t,:)=fit%trend(t,:)+matmul(transpose(b),z(r+1:,t))
         fit%residuals(t,:)=z(1:r,t)-fit%fitted(t,:)
      end do
      ncoef=c+p
      do k=1,r
         actual_mean=sum(z(k,:))/real(neff,dp)
         sse=sum((z(k,:)-fit%fitted(:,k))**2)
         sst=sum((z(k,:)-actual_mean)**2)
         if(sst>0.0_dp)then
            fit%r_squared(k)=1.0_dp-sse/sst
         else
            fit%r_squared(k)=0.0_dp
         end if
         if(ncoef>1 .and. neff>ncoef .and. fit%r_squared(k)<1.0_dp)then
            fit%f_statistic(k)=(fit%r_squared(k)/(1.0_dp-fit%r_squared(k)))* &
               real(neff-ncoef,dp)/real(ncoef-1,dp)
         else
            fit%f_statistic(k)=0.0_dp
         end if
         dfadd=merge(r,p,source); denom=df+real(dfadd,dp)-2.0_dp
         allocate(standardized(neff))
         if(denom>0.0_dp .and. sigma(k,k)>0.0_dp)then
            do t=1,neff
               scale=(df/denom)*sigma(k,k)*fit%conditional_factor(t)
               standardized(t)=fit%residuals(t,k)/sqrt(scale)
            end do
            ad_df=df+real(dfadd,dp)
            call anderson_darling_t(standardized,ad_df,fit%ad_test(k,1),fit%ad_test(k,2))
         else
            fit%ad_test(k,:)=[0.0_dp,1.0_dp]
         end if
         deallocate(standardized)
      end do
   end subroutine populate_fit

   subroutine parameter_vectors(theta,l,lag,r,df,c,coefvec,varvec,ok)
      real(dp), intent(in) :: theta(:),df
      integer, intent(in) :: l,lag,r,c
      real(dp), allocatable, intent(out) :: coefvec(:),varvec(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: s(:,:),b(:,:),sigma(:,:),q(:,:),delta0(:,:),vv(:)
      integer :: n1,n2
      call decode_theta(theta,l,lag,r,df,c,s,b,sigma,q,delta0,varvec,ok)
      if(.not.ok)then; allocate(coefvec(0)); return; end if
      vv=vech(sigma); n1=size(delta0); n2=size(b)
      allocate(coefvec(n1+n2+size(vv)))
      if(n1>0)coefvec(1:n1)=reshape(delta0,[n1])
      coefvec(n1+1:n1+n2)=reshape(b,[n2])
      coefvec(n1+n2+1:)=vv
   end subroutine parameter_vectors

   subroutine inference(theta,l,lag,r,df,c,neff,objective,fit)
      real(dp), intent(in) :: theta(:),df
      integer, intent(in) :: l,lag,r,c,neff
      procedure(objective_interface) :: objective
      type(streg_fit), intent(inout) :: fit
      real(dp), allocatable :: h(:,:),covt(:,:),cp(:),cm(:),vp(:),vm(:),jc(:,:),jv(:,:)
      real(dp) :: xp(size(theta)),xm(size(theta)),step,ridge
      logical :: okp,okm,invok
      integer :: i,nc,nv,nt
      call numerical_hessian(objective,theta,h)
      ridge=1.0e-7_dp*max(1.0_dp,maxval(abs(h)))
      call general_inverse(h,covt,invok,ridge)
      fit%hessian=h
      if(.not.invok)then
         fit%status=streg_hessian_failure
         fit%message=trim(fit%message)//' numerical Hessian could not be regularized.'
         return
      end if
      fit%theta_covariance=covt; nt=size(theta)
      nc=size(fit%coef_estimate); nv=size(fit%var_coef)
      allocate(jc(nc,nt),jv(nv,nt))
      do i=1,nt
         step=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(theta(i)))
         xp=theta; xm=theta; xp(i)=xp(i)+step; xm(i)=xm(i)-step
         call parameter_vectors(xp,l,lag,r,df,c,cp,vp,okp)
         call parameter_vectors(xm,l,lag,r,df,c,cm,vm,okm)
         if(okp.and.okm)then
            jc(:,i)=(cp-cm)/(2.0_dp*step); jv(:,i)=(vp-vm)/(2.0_dp*step)
         else
            jc(:,i)=0.0_dp; jv(:,i)=0.0_dp
         end if
      end do
      allocate(fit%coef_se(nc),fit%coef_p_value(nc),fit%var_coef_se(nv),fit%var_coef_p_value(nv))
      do i=1,nc
         fit%coef_se(i)=sqrt(max(0.0_dp,dot_product(jc(i,:),matmul(covt,jc(i,:)))))
         if(fit%coef_se(i)>0.0_dp)then
            fit%coef_p_value(i)=student_t_two_sided_p(fit%coef_estimate(i)/fit%coef_se(i),real(neff,dp))
         else
            fit%coef_p_value(i)=1.0_dp
         end if
      end do
      do i=1,nv
         fit%var_coef_se(i)=sqrt(max(0.0_dp,dot_product(jv(i,:),matmul(covt,jv(i,:)))))
         if(fit%var_coef_se(i)>0.0_dp)then
            fit%var_coef_p_value(i)=student_t_two_sided_p(fit%var_coef(i)/fit%var_coef_se(i),real(neff,dp))
         else
            fit%var_coef_p_value(i)=1.0_dp
         end if
      end do
   end subroutine inference

   pure function lower_string(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i,k
      out=text
      do i=1,len(text)
         k=iachar(out(i:i))
         if(k>=iachar('A').and.k<=iachar('Z'))out(i:i)=achar(k+32)
      end do
   end function lower_string

end module streg_core
