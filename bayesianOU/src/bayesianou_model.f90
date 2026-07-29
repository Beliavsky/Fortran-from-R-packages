! SPDX-License-Identifier: MIT
module bayesianou_model
  use bayesianou_kinds, only : dp, status_ok, status_bad_input, status_not_converged
  use bayesianou_types
  use bayesianou_math, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_student_t, rng_gamma, &
                              logistic, clamp, normal_logpdf, student_t_logpdf, &
                              sample_mean, sample_sd, median_value, percentile, &
                              ols_fit, nelder_mead, finite_all
  use bayesianou_utils, only : zscore_train, compute_common_factor, orthogonalize_series, &
                               weighted_com_statistics, ou_level_spec
  implicit none
  private
  public :: fit_ou_nested, fit_ou_nonlinear_tmg, simulate_ou_nested
  public :: ou_log_likelihood, ou_mean_increment, extract_posterior_summary
  public :: build_beta_tmg_table, summarize_sv_sigmas, drift_decomposition_grid
  public :: build_accounting_block, extract_mu_trajectory

  type(ou_input), pointer, save :: active_input => null()
  type(ou_fit_result), pointer, save :: active_fit => null()
  integer, save :: active_sector = 0
  real(dp), save :: active_beta1 = 0.0_dp, active_gamma = 0.0_dp

contains

  subroutine fit_ou_nonlinear_tmg(input,options,result)
    type(ou_input), intent(in), target :: input
    type(ou_options), intent(in) :: options
    type(ou_fit_result), intent(out), target :: result
    type(ou_options) :: opt
    opt=options;opt%n_levels=1
    call fit_ou_nested(input,opt,result)
  end subroutine fit_ou_nonlinear_tmg

  subroutine fit_ou_nested(input,options,result)
    type(ou_input), intent(in), target :: input
    type(ou_options), intent(in) :: options
    type(ou_fit_result), intent(out), target :: result
    integer :: t,s,status
    real(dp) :: mu_tmg,sd_tmg,beta_ortho(2)
    real(dp), allocatable :: factor_loading(:), zraw(:)
    logical :: ok

    call validate_input(input,options,ok,result%message)
    if(.not.ok)then;result%status=status_bad_input;return;end if
    result%options=options
    t=size(input%y,1);s=size(input%y,2)
    result%t_train=max(2,int(floor(real(t,dp)*options%train_frac)))
    result%t_lik=merge(t,result%t_train,options%fit_full)
    call zscore_train(input%y,result%t_train,result%zy)
    call zscore_train(input%x,result%t_train,result%zx)
    allocate(result%ztmg(t),result%common_factor(t),factor_loading(s),zraw(t))
    mu_tmg=sample_mean(input%tmg(1:result%t_train));sd_tmg=sample_sd(input%tmg(1:result%t_train));if(sd_tmg<1e-8_dp)sd_tmg=1
    zraw=(input%tmg-mu_tmg)/sd_tmg
    if(options%factor_from_x)then
      call compute_common_factor(result%zx%mz,result%t_train,options%use_train_loadings,result%common_factor,factor_loading,status)
    else
      call compute_common_factor(result%zy%mz,result%t_train,options%use_train_loadings,result%common_factor,factor_loading,status)
    end if
    if(status/=status_ok)then;result%status=status;result%message='PCA factor construction failed';return;end if
    if(options%orthogonalize_tmg)then
      call orthogonalize_series(zraw,result%common_factor,result%t_train,result%ztmg,beta_ortho,status)
    else
      result%ztmg=zraw
    end if
    allocate(result%com_wmean(s),result%com_wsd(s))
    call weighted_com_statistics(input%com,input%capital,result%t_train,result%com_wmean,result%com_wsd)
    call allocate_summary(result%summary,s)
    allocate(result%phi_median(t,s),result%h_median(t,s))
    active_input=>input;active_fit=>result
    if(options%n_levels==1)then
      call estimate_single_level(input,result)
      result%phi_median=result%zx%mz
    else
      call estimate_nested(input,result)
    end if
    call apply_hierarchical_pooling(result)
    call estimate_sv(input,result)
    call generate_draws(input,result)
    call refresh_summary_from_draws(result)
    result%status=status_ok;result%message='ok'
    nullify(active_input);nullify(active_fit)
  end subroutine fit_ou_nested

  subroutine validate_input(input,options,ok,message)
    type(ou_input),intent(in)::input
    type(ou_options),intent(in)::options
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    ok=.false.;message=''
    if(.not.allocated(input%y).or..not.allocated(input%x).or..not.allocated(input%tmg).or. &
       .not.allocated(input%com).or..not.allocated(input%capital))then;message='missing required input';return;end if
    if(size(input%y,1)/=size(input%x,1).or.size(input%y,2)/=size(input%x,2))then;message='Y/X dimensions differ';return;end if
    if(size(input%tmg)/=size(input%y,1))then;message='TMG length differs';return;end if
    if(any(shape(input%com)/=shape(input%y)).or.any(shape(input%capital)/=shape(input%y)))then;message='COM/capital dimensions differ';return;end if
    if(options%n_levels<1.or.options%n_levels>3)then;message='n_levels must be 1, 2 or 3';return;end if
    if(options%n_levels>=2.and..not.allocated(input%gprime))then;message='Gprime required for nested model';return;end if
    if(options%n_levels==3.and..not.allocated(input%value_anchor))then;message='value anchor required for level 3';return;end if
    if(options%train_frac<=0.or.options%train_frac>=1)then;message='train_frac outside (0,1)';return;end if
    if(.not.finite_all(reshape(input%y,[size(input%y)])))then;message='nonfinite Y';return;end if
    ok=.true.
  end subroutine validate_input

  subroutine allocate_summary(summ,s)
    type(ou_summary),intent(inout)::summ
    integer,intent(in)::s
    allocate(summ%theta(s),summ%kappa(s),summ%a3(s),summ%beta0(s),summ%alpha(s),summ%rho(s),summ%sigma_eta(s))
    allocate(summ%kappa_p(s),summ%mu_const(s),summ%sigma_p(s),summ%a3_p(s))
    summ%theta=0;summ%kappa=0.2_dp;summ%a3=-0.05_dp;summ%beta0=0
    summ%alpha=0;summ%rho=0.7_dp;summ%sigma_eta=0.2_dp
    summ%kappa_p=0.1_dp;summ%mu_const=0;summ%sigma_p=0.3_dp;summ%a3_p=-0.05_dp
  end subroutine allocate_summary

  subroutine estimate_single_level(input,result)
    type(ou_input),intent(in),target::input
    type(ou_fit_result),intent(inout),target::result
    integer::iter,s,status,t0,t1,n,p
    real(dp)::x0(4),f
    real(dp),allocatable::design(:,:),yy(:),beta(:),resid(:),cov(:,:)
    t0=2;t1=result%t_lik;n=(t1-t0+1)*size(input%y,2)
    result%summary%beta1=0;result%summary%gamma=0
    do iter=1,4
      active_beta1=result%summary%beta1;active_gamma=result%summary%gamma
      do s=1,size(input%y,2)
        active_sector=s
        x0=[result%summary%theta(s),log(max(result%summary%kappa(s),1e-5_dp)), &
            log(max(-result%summary%a3(s),1e-5_dp)),result%summary%beta0(s)]
        call nelder_mead(single_sector_objective,x0,f,status,max_iter=800,tol=1e-8_dp,step=0.08_dp)
        result%summary%theta(s)=x0(1);result%summary%kappa(s)=exp(clamp(x0(2),-8.0_dp,1.0_dp))
        result%summary%a3(s)=-exp(clamp(x0(3),-10.0_dp,2.0_dp));result%summary%beta0(s)=x0(4)
      end do
      p=1+merge(1,0,result%options%com_in_mean)
      allocate(design(n,p),yy(n),beta(p),resid(n),cov(p,p));call fill_global_regression(input,result,design,yy)
      call ols_fit(design,yy,beta,resid,cov,status,ridge=1e-6_dp)
      result%summary%beta1=beta(1);if(p==2)result%summary%gamma=beta(2)
      deallocate(design,yy,beta,resid,cov)
    end do
    result%summary%nu=estimate_nu(input,result)
  end subroutine estimate_single_level

  function single_sector_objective(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f,theta,kappa,a3,beta0,mu,resid,sd0,comstd
    integer::t,s
    s=active_sector;theta=x(1);kappa=exp(clamp(x(2),-8.0_dp,1.5_dp));a3=-exp(clamp(x(3),-12.0_dp,3.0_dp));beta0=x(4)
    sd0=max(0.05_dp,sample_sd(active_fit%zy%mz(2:active_fit%t_lik,s)-active_fit%zy%mz(1:active_fit%t_lik-1,s)))
    f=0
    do t=2,active_fit%t_lik
      comstd=(active_input%com(t-1,s)-active_fit%com_wmean(s))/active_fit%com_wsd(s)
      mu=kappa*(theta-active_fit%zy%mz(t-1,s)+a3*(active_fit%zy%mz(t-1,s)-theta)**3) &
         +(beta0+active_beta1*active_fit%ztmg(t))*active_fit%zx%mz(t-1,s)
      if(active_fit%options%com_in_mean)mu=mu+active_gamma*comstd
      resid=active_fit%zy%mz(t,s)-active_fit%zy%mz(t-1,s)-mu
      f=f-normal_logpdf(resid,0.0_dp,sd0)
    end do
    f=f+0.5_dp*(theta**2+(x(2)+1)**2/0.25_dp+(x(3)-log(0.05_dp))**2/0.16_dp+beta0**2/0.25_dp)
  end function single_sector_objective

  subroutine fill_global_regression(input,result,design,yy)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::result
    real(dp),intent(out)::design(:,:),yy(:)
    integer::t,s,i
    real(dp)::drift,base,comstd
    i=0
    do t=2,result%t_lik;do s=1,size(input%y,2);i=i+1
      drift=result%summary%kappa(s)*(result%summary%theta(s)-result%zy%mz(t-1,s)+ &
        result%summary%a3(s)*(result%zy%mz(t-1,s)-result%summary%theta(s))**3)
      base=result%zy%mz(t,s)-result%zy%mz(t-1,s)-drift-result%summary%beta0(s)*result%zx%mz(t-1,s)
      yy(i)=base;design(i,1)=result%ztmg(t)*result%zx%mz(t-1,s)
      if(size(design,2)==2)then
        comstd=(input%com(t-1,s)-result%com_wmean(s))/result%com_wsd(s);design(i,2)=comstd
      end if
    end do;end do
  end subroutine fill_global_regression

  subroutine estimate_nested(input,result)
    type(ou_input),intent(in),target::input
    type(ou_fit_result),intent(inout),target::result
    integer::s,t,n,status
    real(dp),allocatable::g_z(:),v_z(:,:),design(:,:),beta(:),resid(:),cov(:,:)
    real(dp)::mg,sg,kp
    n=size(input%y,1);result%phi_median=result%zx%mz
    allocate(g_z(n));mg=sample_mean(input%gprime(1:result%t_train));sg=sample_sd(input%gprime(1:result%t_train));if(sg<1e-8_dp)sg=1
    g_z=(input%gprime-mg)/sg
    if(result%options%n_levels==3)then
      allocate(v_z(n,size(input%y,2)));call standardize_matrix(input%value_anchor,result%t_train,v_z)
    else
      allocate(v_z(n,size(input%y,2)));v_z=0
    end if
    result%summary%beta1=0;result%summary%gamma=0
    do s=1,size(input%y,2)
      call estimate_nested_l1_sector(result,s)
    end do
    call estimate_nested_globals(input,result)
    do s=1,size(input%y,2)
      allocate(design(result%t_lik-1,4),beta(4),resid(result%t_lik-1),cov(4,4))
      do t=2,result%t_lik
        design(t-1,:)=[1.0_dp,-result%phi_median(t-1,s),g_z(t),v_z(t,s)]
      end do
      call ols_fit(design,result%phi_median(2:result%t_lik,s)-result%phi_median(1:result%t_lik-1,s),beta,resid,cov,status,ridge=1e-6_dp)
      kp=clamp(beta(2),1e-3_dp,result%options%kappa_cap*0.95_dp)
      result%summary%kappa_p(s)=kp;result%summary%mu_const(s)=beta(1)/kp
      result%summary%sigma_p(s)=max(sample_sd(resid),1e-4_dp)
      result%summary%a3_p(s)=merge(-0.05_dp,0.0_dp,result%options%level_spec%level2%cubic)
      result%summary%m1=result%summary%m1+beta(3)/kp/real(size(input%y,2),dp)
      result%summary%m_v=result%summary%m_v+beta(4)/kp/real(size(input%y,2),dp)
      deallocate(design,beta,resid,cov)
    end do
    if(result%options%sigma_phi_meas_fixed>0)then
      result%summary%sigma_phi_meas=result%options%sigma_phi_meas_fixed
    else
      result%summary%sigma_phi_meas=max(1e-4_dp,sqrt(sum((result%phi_median-result%zx%mz)**2)/real(size(result%phi_median),dp)))
    end if
    result%summary%nu=10.0_dp
  end subroutine estimate_nested

  subroutine standardize_matrix(x,t_train,z)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::t_train
    real(dp),intent(out)::z(size(x,1),size(x,2))
    integer::s
    real(dp)::m,sd
    do s=1,size(x,2);m=sample_mean(x(1:t_train,s));sd=sample_sd(x(1:t_train,s));if(sd<1e-8_dp)sd=1;z(:,s)=(x(:,s)-m)/sd;end do
  end subroutine standardize_matrix

  subroutine estimate_nested_l1_sector(result,s)
    type(ou_fit_result),intent(inout)::result
    integer,intent(in)::s
    real(dp), allocatable :: dev(:),dy(:),xmat(:,:),res(:)
    real(dp) :: b(2),cov(2,2)
    integer::t,status,n
    n=result%t_lik-1
    allocate(dev(n),dy(n),xmat(n,2),res(n))
    do t=2,result%t_lik
      dev(t-1)=result%zy%mz(t-1,s)-result%phi_median(t-1,s)
      dy(t-1)=result%zy%mz(t,s)-result%zy%mz(t-1,s)
      xmat(t-1,1)=-dev(t-1);xmat(t-1,2)=dev(t-1)**3
    end do
    call ols_fit(xmat,dy,b,res,cov,status,ridge=1e-6_dp)
    result%summary%kappa(s)=clamp(b(1),1e-3_dp,result%options%kappa_cap*0.95_dp)
    if(abs(b(1))>1e-8_dp)then;result%summary%a3(s)=min(-1e-6_dp,b(2)/b(1));else;result%summary%a3(s)=-0.05_dp;end if
  end subroutine estimate_nested_l1_sector


  subroutine estimate_nested_globals(input,result)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(inout)::result
    integer::n,i,t,s,p,status
    real(dp),allocatable::x(:,:),y(:),b(:),r(:),cov(:,:)
    real(dp)::dev,kbase,base,comstd
    n=(result%t_lik-1)*size(input%y,2);p=1+merge(1,0,result%options%com_in_mean)
    allocate(x(n,p),y(n),b(p),r(n),cov(p,p));i=0
    do t=2,result%t_lik;do s=1,size(input%y,2);i=i+1
      dev=result%zy%mz(t-1,s)-result%phi_median(t-1,s);kbase=result%summary%kappa(s)
      base=result%zy%mz(t,s)-result%zy%mz(t-1,s)-kbase*(-dev+result%summary%a3(s)*dev**3)
      y(i)=base;x(i,1)=result%ztmg(t)*(-dev)*kbase*(1.0_dp-kbase/result%options%kappa_cap)
      if(p==2)then;comstd=(input%com(t-1,s)-result%com_wmean(s))/result%com_wsd(s);x(i,2)=comstd;end if
    end do;end do
    call ols_fit(x,y,b,r,cov,status,ridge=1e-5_dp);result%summary%beta1=b(1);if(p==2)result%summary%gamma=b(2)
  end subroutine estimate_nested_globals

  function estimate_nu(input,result) result(nu)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::result
    real(dp)::nu
    real(dp),allocatable::r(:)
    real(dp)::m2,m4,kurt
    integer::t,s,i
    allocate(r((result%t_lik-1)*size(input%y,2)));i=0
    do t=2,result%t_lik;do s=1,size(input%y,2);i=i+1
      r(i)=result%zy%mz(t,s)-result%zy%mz(t-1,s)-ou_mean_increment(input,result%summary,result,t,s)
    end do;end do
    m2=sum((r-sample_mean(r))**2)/real(size(r),dp);m4=sum((r-sample_mean(r))**4)/real(size(r),dp)
    if(m2<=0)then;nu=30;return;end if
    kurt=m4/m2**2-3
    if(kurt>0.05_dp)then;nu=clamp(6.0_dp/kurt+4.0_dp,4.1_dp,100.0_dp);else;nu=100;end if
  end function estimate_nu

  subroutine apply_hierarchical_pooling(result)
    type(ou_fit_result),intent(inout)::result
    integer::n,status
    real(dp),allocatable::coms(:),x(:,:),b(:),r(:),cov(:,:),y(:)
    real(dp)::m,sd,weight
    n=size(result%summary%theta)
    if(n<3)return
    allocate(coms(n),x(n,2),b(2),r(n),cov(2,2),y(n))
    m=sample_mean(result%com_wmean);sd=sample_sd(result%com_wmean);if(sd<1e-8_dp)sd=1
    coms=(result%com_wmean-m)/sd;x(:,1)=1.0_dp;x(:,2)=coms;weight=0.35_dp
    if(result%options%level_spec%level1%hierarchy.or.result%options%n_levels==1)then
      y=result%summary%theta;call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%theta=(1-weight)*y+weight*matmul(x,b)
      y=log(max(result%summary%kappa,1e-8_dp));call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%kappa=exp((1-weight)*y+weight*matmul(x,b))
      y=log(max(-result%summary%a3,1e-10_dp));call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%a3=-exp((1-weight)*y+weight*matmul(x,b))
      y=result%summary%beta0;call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%beta0=(1-weight)*y+weight*matmul(x,b)
    end if
    if(result%options%n_levels>=2.and.result%options%level_spec%level2%hierarchy)then
      y=log(max(result%summary%kappa_p,1e-8_dp));call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%kappa_p=exp((1-weight)*y+weight*matmul(x,b))
      y=result%summary%mu_const;call ols_fit(x,y,b,r,cov,status,ridge=1e-6_dp);result%summary%mu_const=(1-weight)*y+weight*matmul(x,b)
    end if
  end subroutine apply_hierarchical_pooling

  subroutine estimate_sv(input,result)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(inout)::result
    integer::s,t,status,n
    real(dp),allocatable::hobs(:),x(:,:),b(:),res(:),cov(:,:)
    real(dp)::rho,alpha
    n=result%t_lik-1
    do s=1,size(input%y,2)
      allocate(hobs(n),x(n-1,2),b(2),res(n-1),cov(2,2))
      do t=2,result%t_lik
        hobs(t-1)=log((result%zy%mz(t,s)-result%zy%mz(t-1,s)-ou_mean_increment(input,result%summary,result,t,s))**2+1e-6_dp)
      end do
      x(:,1)=1.0_dp;x(:,2)=hobs(1:n-1)
      call ols_fit(x,hobs(2:n),b,res,cov,status,ridge=1e-5_dp)
      rho=clamp(b(2),-0.98_dp,0.98_dp);alpha=b(1)/max(0.02_dp,1.0_dp-rho)
      result%summary%rho(s)=rho;result%summary%alpha(s)=alpha;result%summary%sigma_eta(s)=max(0.05_dp,sample_sd(res))
      result%h_median(1,s)=alpha
      do t=2,size(input%y,1)
        if(t<=result%t_lik)then
          result%h_median(t,s)=alpha+rho*(result%h_median(t-1,s)-alpha)+0.35_dp*(hobs(t-1)-result%h_median(t-1,s))
        else
          result%h_median(t,s)=alpha+rho*(result%h_median(t-1,s)-alpha)
        end if
      end do
      deallocate(hobs,x,b,res,cov)
    end do
  end subroutine estimate_sv

  function ou_mean_increment(input,summ,fit,t,s) result(mu)
    type(ou_input),intent(in)::input
    type(ou_summary),intent(in)::summ
    type(ou_fit_result),intent(in)::fit
    integer,intent(in)::t,s
    real(dp)::mu,dev,k,comstd
    comstd=(input%com(t-1,s)-fit%com_wmean(s))/fit%com_wsd(s)
    if(fit%options%n_levels==1)then
      dev=fit%zy%mz(t-1,s)-summ%theta(s)
      mu=summ%kappa(s)*(-dev+summ%a3(s)*dev**3)+(summ%beta0(s)+summ%beta1*fit%ztmg(t))*fit%zx%mz(t-1,s)
    else
      dev=fit%zy%mz(t-1,s)-fit%phi_median(t-1,s)
      k=fit%options%kappa_cap*logistic(logit_safe(summ%kappa(s)/fit%options%kappa_cap)+summ%beta1*fit%ztmg(t))
      mu=k*(-dev+merge(summ%a3(s)*dev**3,0.0_dp,fit%options%level_spec%level1%cubic))
    end if
    if(fit%options%com_in_mean)mu=mu+summ%gamma*comstd
  end function ou_mean_increment

  pure function logit_safe(p) result(v)
    real(dp),intent(in)::p
    real(dp)::v,q
    q=clamp(p,1e-8_dp,1.0_dp-1e-8_dp);v=log(q/(1-q))
  end function logit_safe

  function ou_log_likelihood(input,fit,summ,pointwise) result(ll)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::fit
    type(ou_summary),intent(in)::summ
    real(dp),intent(out),optional::pointwise(:)
    real(dp)::ll,resid,sd
    integer::t,s,i
    ll=0;i=0
    do t=2,fit%t_lik;do s=1,size(input%y,2);i=i+1
      resid=fit%zy%mz(t,s)-fit%zy%mz(t-1,s)-ou_mean_increment(input,summ,fit,t,s)
      sd=clamp(exp(0.5_dp*fit%h_median(t,s)),1e-8_dp,1e8_dp)
      if(fit%options%n_levels>=2.and..not.fit%options%level_spec%level1%student_t)then
        if(present(pointwise))pointwise(i)=normal_logpdf(resid,0.0_dp,sd)
      else
        if(present(pointwise))pointwise(i)=student_t_logpdf(resid,summ%nu,0.0_dp,sd)
      end if
      if(present(pointwise))then;ll=ll+pointwise(i);else
        if(fit%options%n_levels>=2.and..not.fit%options%level_spec%level1%student_t)then;ll=ll+normal_logpdf(resid,0.0_dp,sd);else;ll=ll+student_t_logpdf(resid,summ%nu,0.0_dp,sd);end if
      end if
    end do;end do
  end function ou_log_likelihood

  subroutine generate_draws(input,result)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(inout)::result
    type(rng_state)::rng
    type(ou_summary)::current,proposal
    integer::nkeep,n,s,c,it,i,sector,accepted,total_proposals
    real(dp)::lp,lp2,scale,u
    nkeep=max(1,(result%options%iterations-result%options%warmup)/max(1,result%options%thin))
    n=nkeep*result%options%chains;s=size(input%y,2)
    call allocate_draws(result%draws,n,s)
    i=0;accepted=0;total_proposals=0
    do c=1,result%options%chains
      call rng_seed(rng,result%options%seed+104729*c)
      current=result%summary
      call jitter_summary(current,rng,0.01_dp)
      lp=summary_log_posterior(input,result,current)
      scale=result%options%proposal_scale
      do it=1,result%options%iterations
        proposal=current
        sector=1+mod(it-1,s)
        call propose_global(proposal,rng,0.5_dp*scale)
        call propose_sector(proposal,sector,result%options%n_levels,rng,scale)
        lp2=summary_log_posterior(input,result,proposal)
        total_proposals=total_proposals+1
        u=log(rng_uniform(rng))
        if(lp2>lp+u)then
          current=proposal;lp=lp2;accepted=accepted+1
        end if
        if(it<=result%options%warmup.and.mod(it,25)==0)then
          if(real(accepted,dp)/real(max(1,total_proposals),dp)<0.15_dp)scale=0.8_dp*scale
          if(real(accepted,dp)/real(max(1,total_proposals),dp)>0.50_dp)scale=1.2_dp*scale
          scale=clamp(scale,1.0e-4_dp,0.25_dp)
        end if
        if(it>result%options%warmup.and.mod(it-result%options%warmup,result%options%thin)==0)then
          i=i+1;call record_draw(result%draws,i,c,current,rng)
        end if
      end do
    end do
    result%diagnostics%acceptance_rate=real(accepted,dp)/real(max(1,total_proposals),dp)
  end subroutine generate_draws

  subroutine jitter_summary(summ,rng,scale)
    type(ou_summary),intent(inout)::summ
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::scale
    integer::j
    do j=1,size(summ%theta)
      summ%theta(j)=summ%theta(j)+scale*rng_normal(rng)
      summ%kappa(j)=exp(log(max(summ%kappa(j),1e-8_dp))+scale*rng_normal(rng))
      summ%a3(j)=-exp(log(max(-summ%a3(j),1e-10_dp))+scale*rng_normal(rng))
      summ%beta0(j)=summ%beta0(j)+scale*rng_normal(rng)
    end do
  end subroutine jitter_summary

  subroutine propose_global(summ,rng,scale)
    type(ou_summary),intent(inout)::summ
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::scale
    summ%beta1=summ%beta1+scale*rng_normal(rng)
    summ%gamma=summ%gamma+scale*rng_normal(rng)
    summ%nu=2.0_dp+exp(log(max(summ%nu-2.0_dp,1e-4_dp))+0.5_dp*scale*rng_normal(rng))
    summ%m1=summ%m1+scale*rng_normal(rng)
    summ%m_v=summ%m_v+scale*rng_normal(rng)
    summ%sigma_phi_meas=exp(log(max(summ%sigma_phi_meas,1e-6_dp))+0.5_dp*scale*rng_normal(rng))
  end subroutine propose_global

  subroutine propose_sector(summ,j,n_levels,rng,scale)
    type(ou_summary),intent(inout)::summ
    integer,intent(in)::j,n_levels
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::scale
    summ%theta(j)=summ%theta(j)+scale*rng_normal(rng)
    summ%kappa(j)=exp(log(max(summ%kappa(j),1e-8_dp))+scale*rng_normal(rng))
    summ%a3(j)=-exp(log(max(-summ%a3(j),1e-10_dp))+scale*rng_normal(rng))
    summ%beta0(j)=summ%beta0(j)+scale*rng_normal(rng)
    if(n_levels>=2)then
      summ%kappa_p(j)=exp(log(max(summ%kappa_p(j),1e-8_dp))+scale*rng_normal(rng))
      summ%mu_const(j)=summ%mu_const(j)+scale*rng_normal(rng)
      summ%sigma_p(j)=exp(log(max(summ%sigma_p(j),1e-8_dp))+scale*rng_normal(rng))
      if(summ%a3_p(j)<0.0_dp)summ%a3_p(j)=-exp(log(max(-summ%a3_p(j),1e-10_dp))+scale*rng_normal(rng))
    end if
  end subroutine propose_sector

  function summary_log_posterior(input,result,summ) result(lp)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::result
    type(ou_summary),intent(in)::summ
    real(dp)::lp,nu_shift,gmu,gsd,gz,dev,meanp,resid,sp
    integer::j,t
    lp=ou_log_likelihood(input,result,summ)
    if(.not.(lp>-huge(1.0_dp)))return
    do j=1,size(summ%theta)
      lp=lp+normal_logpdf(summ%theta(j),0.0_dp,1.0_dp)
      lp=lp+normal_logpdf(log(max(summ%kappa(j),1e-12_dp)),-1.0_dp,0.5_dp)+log(max(summ%kappa(j),1e-12_dp))
      lp=lp+normal_logpdf(log(max(-summ%a3(j),1e-12_dp)),log(0.05_dp),0.4_dp)+log(max(-summ%a3(j),1e-12_dp))
      lp=lp+normal_logpdf(summ%beta0(j),0.0_dp,0.5_dp)
    end do
    lp=lp+normal_logpdf(summ%beta1,result%options%priors%beta1_mean,result%options%priors%beta1_sd)
    lp=lp+normal_logpdf(summ%gamma,0.0_dp,0.5_dp)
    nu_shift=max(summ%nu-2.0_dp,1e-12_dp)
    lp=lp+(result%options%priors%nu_shape-1.0_dp)*log(nu_shift)-result%options%priors%nu_rate*nu_shift
    if(result%options%n_levels>=2)then
      gmu=sample_mean(input%gprime(1:result%t_train));gsd=sample_sd(input%gprime(1:result%t_train));if(gsd<1e-8_dp)gsd=1
      do j=1,size(summ%theta)
        lp=lp+normal_logpdf(log(max(summ%kappa_p(j),1e-12_dp)),-1.5_dp,0.7_dp)+log(max(summ%kappa_p(j),1e-12_dp))
        lp=lp+normal_logpdf(summ%mu_const(j),0.0_dp,1.0_dp)
        lp=lp+normal_logpdf(log(max(summ%sigma_p(j),1e-12_dp)),-1.0_dp,0.7_dp)+log(max(summ%sigma_p(j),1e-12_dp))
        do t=2,result%t_lik
          gz=(input%gprime(t)-gmu)/gsd
          meanp=summ%mu_const(j)+summ%m1*gz
          if(result%options%n_levels==3.and.allocated(input%value_anchor))meanp=meanp+summ%m_v*input%value_anchor(t,j)
          dev=result%phi_median(t-1,j)-meanp
          resid=result%phi_median(t,j)-result%phi_median(t-1,j)-summ%kappa_p(j)*(-dev+summ%a3_p(j)*dev**3)
          sp=max(summ%sigma_p(j),1e-8_dp)
          if(result%options%level_spec%level2%student_t)then
            lp=lp+student_t_logpdf(resid,summ%nu,0.0_dp,sp)
          else
            lp=lp+normal_logpdf(resid,0.0_dp,sp)
          end if
        end do
        sp=max(summ%sigma_phi_meas,1e-8_dp)
        do t=1,result%t_lik
          lp=lp+normal_logpdf(result%zx%mz(t,j),result%phi_median(t,j),sp)
        end do
      end do
      lp=lp+normal_logpdf(summ%m1,0.0_dp,0.5_dp)+normal_logpdf(summ%m_v,0.0_dp,0.5_dp)
      if(result%options%sigma_phi_meas_fixed<=0.0_dp)lp=lp+normal_logpdf(summ%sigma_phi_meas,0.0_dp,result%options%priors%sigma_phi_meas_sd)
    end if
  end function summary_log_posterior

  subroutine record_draw(draws,i,chain,summ,rng)
    type(ou_draws),intent(inout)::draws
    integer,intent(in)::i,chain
    type(ou_summary),intent(in)::summ
    type(rng_state),intent(inout)::rng
    integer::j
    draws%chain(i)=chain
    do j=1,size(summ%theta)
      draws%theta(i,j)=summ%theta(j);draws%kappa(i,j)=summ%kappa(j);draws%a3(i,j)=summ%a3(j);draws%beta0(i,j)=summ%beta0(j)
      draws%alpha(i,j)=summ%alpha(j)+0.015_dp*rng_normal(rng)
      draws%rho(i,j)=clamp(summ%rho(j)+0.008_dp*rng_normal(rng),-0.995_dp,0.995_dp)
      draws%sigma_eta(i,j)=exp(log(max(summ%sigma_eta(j),1e-8_dp))+0.015_dp*rng_normal(rng))
      draws%kappa_p(i,j)=summ%kappa_p(j);draws%mu_const(i,j)=summ%mu_const(j);draws%sigma_p(i,j)=summ%sigma_p(j);draws%a3_p(i,j)=summ%a3_p(j)
    end do
    draws%beta1(i)=summ%beta1;draws%gamma(i)=summ%gamma;draws%nu(i)=summ%nu
    draws%m1(i)=summ%m1;draws%m_v(i)=summ%m_v;draws%sigma_phi_meas(i)=summ%sigma_phi_meas
  end subroutine record_draw

  subroutine allocate_draws(draws,n,s)
    type(ou_draws),intent(out)::draws
    integer,intent(in)::n,s
    allocate(draws%theta(n,s),draws%kappa(n,s),draws%a3(n,s),draws%beta0(n,s),draws%alpha(n,s),draws%rho(n,s),draws%sigma_eta(n,s))
    allocate(draws%kappa_p(n,s),draws%mu_const(n,s),draws%sigma_p(n,s),draws%a3_p(n,s))
    allocate(draws%beta1(n),draws%gamma(n),draws%nu(n),draws%m1(n),draws%m_v(n),draws%sigma_phi_meas(n),draws%chain(n))
  end subroutine allocate_draws

  subroutine refresh_summary_from_draws(result)
    type(ou_fit_result),intent(inout)::result
    integer::s
    do s=1,size(result%summary%theta)
      result%summary%theta(s)=median_value(result%draws%theta(:,s));result%summary%kappa(s)=median_value(result%draws%kappa(:,s))
      result%summary%a3(s)=median_value(result%draws%a3(:,s));result%summary%beta0(s)=median_value(result%draws%beta0(:,s))
      result%summary%alpha(s)=median_value(result%draws%alpha(:,s));result%summary%rho(s)=median_value(result%draws%rho(:,s))
      result%summary%sigma_eta(s)=median_value(result%draws%sigma_eta(:,s));result%summary%kappa_p(s)=median_value(result%draws%kappa_p(:,s))
      result%summary%mu_const(s)=median_value(result%draws%mu_const(:,s));result%summary%sigma_p(s)=median_value(result%draws%sigma_p(:,s))
      result%summary%a3_p(s)=median_value(result%draws%a3_p(:,s))
    end do
    result%summary%beta1=median_value(result%draws%beta1);result%summary%gamma=median_value(result%draws%gamma)
    result%summary%nu=median_value(result%draws%nu);result%summary%m1=median_value(result%draws%m1);result%summary%m_v=median_value(result%draws%m_v)
    result%summary%sigma_phi_meas=median_value(result%draws%sigma_phi_meas)
  end subroutine refresh_summary_from_draws

  subroutine extract_posterior_summary(fit,summary)
    type(ou_fit_result),intent(in)::fit
    type(ou_summary),intent(out)::summary
    summary=fit%summary
  end subroutine extract_posterior_summary

  subroutine simulate_ou_nested(summary,options,t,s,x,tmg,com,capital,gprime,value_anchor,seed,y,phi,h)
    type(ou_summary),intent(in)::summary
    type(ou_options),intent(in)::options
    integer,intent(in)::t,s,seed
    real(dp),intent(in)::x(t,s),tmg(t),com(t,s),capital(t,s)
    real(dp),intent(in),optional::gprime(t),value_anchor(t,s)
    real(dp),intent(out)::y(t,s),phi(t,s),h(t,s)
    type(rng_state)::rng
    integer::i,j
    real(dp)::dev,mu,sd,z,k,mg,com_mu,com_sd,denom
    call rng_seed(rng,seed);y=0;phi=x;h=0
    do j=1,s
      h(1,j)=summary%alpha(j)+summary%sigma_eta(j)*rng_normal(rng)/sqrt(max(1e-8_dp,1-summary%rho(j)**2))
      if(options%n_levels>=2)phi(1,j)=x(1,j)+summary%sigma_phi_meas*rng_normal(rng)
    end do
    do i=2,t
      do j=1,s
        if(options%n_levels>=2)then
          mg=summary%mu_const(j);if(present(gprime))mg=mg+summary%m1*gprime(i);if(options%n_levels==3.and.present(value_anchor))mg=mg+summary%m_v*value_anchor(i,j)
          dev=phi(i-1,j)-mg;mu=summary%kappa_p(j)*(-dev+summary%a3_p(j)*dev**3)
          phi(i,j)=phi(i-1,j)+mu+summary%sigma_p(j)*rng_normal(rng)
        end if
        h(i,j)=summary%alpha(j)+summary%rho(j)*(h(i-1,j)-summary%alpha(j))+summary%sigma_eta(j)*rng_normal(rng)
        if(options%n_levels==1)then
          dev=y(i-1,j)-summary%theta(j);mu=summary%kappa(j)*(-dev+summary%a3(j)*dev**3)+(summary%beta0(j)+summary%beta1*tmg(i))*x(i-1,j)
        else
          dev=y(i-1,j)-phi(i-1,j);k=options%kappa_cap*logistic(logit_safe(summary%kappa(j)/options%kappa_cap)+summary%beta1*tmg(i));mu=k*(-dev+summary%a3(j)*dev**3)
        end if
        if(options%com_in_mean) then
          denom=sum(capital(:,j)); if(denom<=0.0_dp) denom=1.0_dp
          com_mu=sum(com(:,j)*capital(:,j))/denom
          com_sd=sqrt(max(1.0e-16_dp,sum(capital(:,j)*(com(:,j)-com_mu)**2)/denom))
          mu=mu+summary%gamma*(com(i-1,j)-com_mu)/com_sd
        end if
        sd=exp(0.5_dp*h(i,j));z=merge(rng_student_t(rng,summary%nu),rng_normal(rng),summary%nu<90.0_dp)
        y(i,j)=y(i-1,j)+mu+sd*z
      end do
    end do
  end subroutine simulate_ou_nested

  subroutine build_beta_tmg_table(fit,beta)
    type(ou_fit_result),intent(in)::fit
    real(dp),intent(out)::beta(size(fit%ztmg),size(fit%summary%beta0))
    integer::s
    do s=1,size(beta,2);beta(:,s)=fit%summary%beta0(s)+fit%summary%beta1*fit%ztmg;end do
  end subroutine build_beta_tmg_table

  subroutine summarize_sv_sigmas(fit,median_sigma,lo,hi)
    type(ou_fit_result),intent(in)::fit
    real(dp),intent(out)::median_sigma(size(fit%summary%alpha)),lo(size(fit%summary%alpha)),hi(size(fit%summary%alpha))
    integer::s
    do s=1,size(median_sigma)
      median_sigma(s)=median_value(exp(0.5_dp*fit%draws%alpha(:,s)))
      lo(s)=percentile(exp(0.5_dp*fit%draws%alpha(:,s)),0.025_dp);hi(s)=percentile(exp(0.5_dp*fit%draws%alpha(:,s)),0.975_dp)
    end do
  end subroutine summarize_sv_sigmas

  subroutine drift_decomposition_grid(fit,sector,grid,linear,cubic,total)
    type(ou_fit_result),intent(in)::fit
    integer,intent(in)::sector
    real(dp),intent(in)::grid(:)
    real(dp),intent(out)::linear(size(grid)),cubic(size(grid)),total(size(grid))
    linear=-fit%summary%kappa(sector)*grid
    cubic=fit%summary%kappa(sector)*fit%summary%a3(sector)*grid**3
    total=linear+cubic
  end subroutine drift_decomposition_grid

  subroutine build_accounting_block(tmg_raw,ztmg_exo,ztmg_use,capital,block)
    real(dp),intent(in)::tmg_raw(:),ztmg_exo(:),ztmg_use(:),capital(:,:)
    real(dp),intent(out)::block(size(tmg_raw),4)
    block(:,1)=tmg_raw;block(:,2)=ztmg_exo;block(:,3)=ztmg_use;block(:,4)=sum(capital,dim=2)
  end subroutine build_accounting_block

  subroutine extract_mu_trajectory(fit,gprime_z,value_z,mu)
    type(ou_fit_result),intent(in)::fit
    real(dp),intent(in)::gprime_z(:)
    real(dp),intent(in),optional::value_z(:,:)
    real(dp),intent(out)::mu(size(gprime_z),size(fit%summary%mu_const))
    integer::s
    do s=1,size(mu,2)
      mu(:,s)=fit%summary%mu_const(s)+fit%summary%m1*gprime_z
      if(present(value_z))mu(:,s)=mu(:,s)+fit%summary%m_v*value_z(:,s)
    end do
  end subroutine extract_mu_trajectory

end module bayesianou_model
