! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_gls
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp, pi_dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR, NLME_MAX_ITER
  use nlme_types, only : correlation_spec, variance_spec, gls_result, nlme_control, &
       COR_NONE, COR_AR1, COR_CAR1, COR_ARMA, COR_COMPOUND_SYMM, COR_EXPONENTIAL, &
       COR_GAUSSIAN, COR_LINEAR, COR_RATIO, COR_SPHERICAL, COR_UNSTRUCTURED, &
       VAR_CONSTANT, VAR_IDENT
  use nlme_correlation, only : correlation_parameter_count, correlation_to_unconstrained, &
       correlation_from_unconstrained
  use nlme_variance, only : variance_parameter_count, variance_to_unconstrained, variance_from_unconstrained
  use nlme_covariance, only : build_residual_covariance, default_group_vector
  use nlme_linalg, only : solve_spd, inverse_spd
  use nlme_optimize, only : nelder_mead
  implicit none
  private
  public :: fit_gls, gls_log_likelihood, initialize_correlation, initialize_variance

  type :: gls_context
    real(dp), allocatable :: y(:),x(:,:),time(:),var_covariate(:),coordinates(:,:)
    integer, allocatable :: group(:),var_group(:)
    type(correlation_spec) :: corr_template
    type(variance_spec) :: var_template
    integer :: kc=0,kv=0,nlevels=1,max_group_size=1
    logical :: reml=.true.
  end type gls_context
contains
  subroutine initialize_correlation(spec,time,group,coordinates)
    type(correlation_spec), intent(inout) :: spec
    real(dp), intent(in) :: time(:)
    integer, intent(in) :: group(:)
    real(dp), intent(in), optional :: coordinates(:,:)
    integer :: k,nmax,g
    real(dp) :: span
    nmax=0
    do g=minval(group),maxval(group)
      nmax=max(nmax,count(group==g))
    end do
    k=correlation_parameter_count(spec,nmax)
    if (allocated(spec%par)) then
      if (size(spec%par)==k) return
      deallocate(spec%par)
    end if
    allocate(spec%par(k))
    spec%par=0.0_dp
    select case(spec%kind)
    case(COR_AR1)
      spec%par=0.1_dp
    case(COR_CAR1)
      spec%par=0.5_dp
    case(COR_ARMA)
      spec%par=0.0_dp
    case(COR_COMPOUND_SYMM)
      spec%par=0.05_dp
    case(COR_EXPONENTIAL,COR_GAUSSIAN,COR_LINEAR,COR_RATIO,COR_SPHERICAL)
      if (present(coordinates)) then
        span=maxval(coordinates)-minval(coordinates)
      else
        span=maxval(time)-minval(time)
      end if
      spec%par(1)=max(1.0_dp,0.25_dp*span)
      if (spec%nugget) spec%par(2)=0.9_dp
    case(COR_UNSTRUCTURED)
      spec%par=0.0_dp
    end select
  end subroutine initialize_correlation

  subroutine initialize_variance(spec,var_group)
    type(variance_spec), intent(inout) :: spec
    integer, intent(in) :: var_group(:)
    integer :: k,nlevels
    nlevels=max(1,maxval(var_group))
    k=variance_parameter_count(spec,nlevels)
    if (allocated(spec%par)) then
      if (size(spec%par)==k) return
      deallocate(spec%par)
    end if
    allocate(spec%par(k))
    spec%par=0.0_dp
  end subroutine initialize_variance

  subroutine fit_gls(y,x,result,correlation,variance,time,group,var_covariate,var_group,coordinates,control)
    real(dp), intent(in) :: y(:),x(:,:)
    type(gls_result), intent(out) :: result
    type(correlation_spec), intent(in), optional :: correlation
    type(variance_spec), intent(in), optional :: variance
    real(dp), intent(in), optional :: time(:),var_covariate(:),coordinates(:,:)
    integer, intent(in), optional :: group(:),var_group(:)
    type(nlme_control), intent(in), optional :: control
    type(gls_context) :: ctx
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    type(nlme_control) :: ctl
    real(dp), allocatable :: xc(:),xv(:),par(:)
    real(dp) :: objective
    integer :: n,nlevels,status,it

    n=size(y)
    if (size(x,1)/=n .or. size(x,2)<1 .or. n<=size(x,2)) then
      result%status=NLME_DIMENSION_ERROR
      return
    end if
    ctl=nlme_control()
    if (present(control)) ctl=control
    corr=correlation_spec()
    if (present(correlation)) corr=correlation
    var=variance_spec()
    if (present(variance)) var=variance
    allocate(ctx%y(n),ctx%x(n,size(x,2)),ctx%time(n),ctx%var_covariate(n))
    ctx%y=y
    ctx%x=x
    if (present(time)) then
      if (size(time)/=n) then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      ctx%time=time
    else
      ctx%time=[(real(it,dp),it=1,n)]
    end if
    if (present(group)) then
      if (size(group)/=n) then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      allocate(ctx%group(n))
      ctx%group=group
    else
      call default_group_vector(n,ctx%group)
    end if
    if (present(var_covariate)) then
      if (size(var_covariate)/=n) then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      ctx%var_covariate=var_covariate
    else
      ctx%var_covariate=abs(y)
    end if
    if (present(var_group)) then
      if (size(var_group)/=n) then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      allocate(ctx%var_group(n))
      ctx%var_group=var_group
    else
      allocate(ctx%var_group(n))
      ctx%var_group=1
    end if
    if (present(coordinates)) then
      if (size(coordinates,1)/=n) then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      allocate(ctx%coordinates(n,size(coordinates,2)))
      ctx%coordinates=coordinates
    else
      allocate(ctx%coordinates(n,0))
    end if
    call initialize_correlation(corr,ctx%time,ctx%group,ctx%coordinates)
    call initialize_variance(var,ctx%var_group)
    ctx%corr_template=corr
    ctx%var_template=var
    ctx%reml=ctl%reml
    nlevels=max(1,maxval(ctx%var_group))
    ctx%nlevels=nlevels
    ctx%max_group_size=1
    do it=1,size(ctx%group)
      ctx%max_group_size=max(ctx%max_group_size,count(ctx%group==ctx%group(it)))
    end do
    ctx%kc=merge(0,correlation_parameter_count(corr,ctx%max_group_size),corr%fixed)
    ctx%kv=merge(0,variance_parameter_count(var,nlevels),var%fixed)
    allocate(par(ctx%kc+ctx%kv))
    par=0.0_dp
    if (ctx%kc>0) then
      call correlation_to_unconstrained(corr,ctx%max_group_size,xc,status)
      if (status/=NLME_SUCCESS) then
      result%status=status
      return
      end if
      par(:ctx%kc)=xc
    end if
    if (ctx%kv>0) then
      call variance_to_unconstrained(var,nlevels,xv,status)
      if (status/=NLME_SUCCESS) then
      result%status=status
      return
      end if
      par(ctx%kc+1:)=xv
    end if
    if (size(par)>0 .and. ctl%optimize_covariance) then
      call nelder_mead(gls_objective,ctx,par,objective,status,it,ctl%max_iter,ctl%tolerance,ctl%step,ctl%verbose)
    else
      objective=gls_objective(par,ctx)
      status=NLME_SUCCESS
      it=0
    end if
    call fill_gls_result(par,ctx,result,status,it)
  end subroutine fit_gls

  function gls_objective(par,context) result(value)
    real(dp), intent(in) :: par(:)
    class(*), intent(inout) :: context
    real(dp) :: value
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    integer :: status
    select type(ctx=>context)
    type is(gls_context)
      call unpack_specs(par,ctx,corr,var,status)
      if (status/=NLME_SUCCESS) then
      value=huge(1.0_dp)/100.0_dp
      return
      end if
      call evaluate_gls(ctx,corr,var,value,status)
      if (status/=NLME_SUCCESS .or. .not.ieee_is_finite(value)) value=huge(1.0_dp)/100.0_dp
    class default
      value=huge(1.0_dp)/100.0_dp
    end select
  end function gls_objective

  subroutine unpack_specs(par,ctx,corr,var,status)
    real(dp), intent(in) :: par(:)
    type(gls_context), intent(in) :: ctx
    type(correlation_spec), intent(out) :: corr
    type(variance_spec), intent(out) :: var
    integer, intent(out) :: status
    if (ctx%kc>0) then
      call correlation_from_unconstrained(ctx%corr_template,ctx%max_group_size,par(:ctx%kc),corr,status)
      if (status/=NLME_SUCCESS) return
    else
      corr=ctx%corr_template
    end if
    if (ctx%kv>0) then
      call variance_from_unconstrained(ctx%var_template,ctx%nlevels,par(ctx%kc+1:),var,status)
    else
      var=ctx%var_template
      status=NLME_SUCCESS
    end if
  end subroutine unpack_specs

  subroutine evaluate_gls(ctx,corr,var,nll,status,beta,beta_cov,sigma,covariance,residuals,loglik)
    type(gls_context), intent(in) :: ctx
    type(correlation_spec), intent(in) :: corr
    type(variance_spec), intent(in) :: var
    real(dp), intent(out) :: nll
    integer, intent(out) :: status
    real(dp), allocatable, intent(out), optional :: beta(:),beta_cov(:,:),covariance(:,:),residuals(:)
    real(dp), intent(out), optional :: sigma,loglik
    real(dp), allocatable :: v(:,:),vinvx(:,:),vinvy(:),xtvinvx(:,:),xtvinvy(:),b(:),binv(:,:),r(:),vinvr(:)
    real(dp) :: logdetv,logdetx,rss,s2,ll
    integer :: n,p,df
    n=size(ctx%y)
    p=size(ctx%x,2)
    if (size(ctx%coordinates,2)>0) then
      call build_residual_covariance(corr,var,ctx%time,ctx%group,ctx%var_covariate,ctx%var_group,v,status,ctx%coordinates)
    else
      call build_residual_covariance(corr,var,ctx%time,ctx%group,ctx%var_covariate,ctx%var_group,v,status)
    end if
    if (status/=NLME_SUCCESS) then
    nll=huge(1.0_dp)
    return
    end if
    call solve_spd(v,ctx%x,vinvx,status,logdetv)
    if (status/=NLME_SUCCESS) then
    nll=huge(1.0_dp)
    return
    end if
    call solve_spd(v,ctx%y,vinvy,status)
    if (status/=NLME_SUCCESS) then
    nll=huge(1.0_dp)
    return
    end if
    xtvinvx=matmul(transpose(ctx%x),vinvx)
    xtvinvy=matmul(transpose(ctx%x),vinvy)
    call solve_spd(xtvinvx,xtvinvy,b,status,logdetx)
    if (status/=NLME_SUCCESS) then
    nll=huge(1.0_dp)
    return
    end if
    r=ctx%y-matmul(ctx%x,b)
    call solve_spd(v,r,vinvr,status)
    if (status/=NLME_SUCCESS) then
    nll=huge(1.0_dp)
    return
    end if
    rss=dot_product(r,vinvr)
    df=merge(n-p,n,ctx%reml)
    if (rss<=0.0_dp .or. df<=0) then
    status=NLME_INVALID_ARGUMENT
    nll=huge(1.0_dp)
    return
    end if
    s2=rss/real(df,dp)
    if (ctx%reml) then
      ll=-0.5_dp*(real(df,dp)*(log(2.0_dp*pi_dp)+1.0_dp+log(s2))+logdetv+logdetx)
    else
      ll=-0.5_dp*(real(n,dp)*(log(2.0_dp*pi_dp)+1.0_dp+log(s2))+logdetv)
    end if
    nll=-ll
    status=NLME_SUCCESS
    if (present(beta)) then
    allocate(beta(size(b)))
    beta=b
    end if
    if (present(beta_cov)) then
      call inverse_spd(xtvinvx,binv,status)
      if (status/=NLME_SUCCESS) return
      allocate(beta_cov(p,p))
      beta_cov=s2*binv
    end if
    if (present(sigma)) sigma=sqrt(s2)
    if (present(covariance)) then
    allocate(covariance(n,n))
    covariance=s2*v
    end if
    if (present(residuals)) then
    allocate(residuals(n))
    residuals=r
    end if
    if (present(loglik)) loglik=ll
  end subroutine evaluate_gls

  subroutine fill_gls_result(par,ctx,result,opt_status,iterations)
    real(dp), intent(in) :: par(:)
    type(gls_context), intent(in) :: ctx
    type(gls_result), intent(out) :: result
    integer, intent(in) :: opt_status,iterations
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    real(dp) :: nll
    integer :: status,k
    call unpack_specs(par,ctx,corr,var,status)
    if (status/=NLME_SUCCESS) then
    result%status=status
    return
    end if
    call evaluate_gls(ctx,corr,var,nll,status,result%beta,result%beta_cov,result%sigma,result%covariance,&
         result%residuals,result%log_likelihood)
    if (status/=NLME_SUCCESS) then
    result%status=status
    return
    end if
    allocate(result%fitted(size(ctx%y)))
    result%fitted=ctx%y-result%residuals
    if (allocated(corr%par)) then
      allocate(result%correlation_parameters(size(corr%par)))
      result%correlation_parameters=corr%par
    else
      allocate(result%correlation_parameters(0))
    end if
    if (allocated(var%par)) then
      allocate(result%variance_parameters(size(var%par)))
      result%variance_parameters=var%par
    else
      allocate(result%variance_parameters(0))
    end if
    k=size(result%beta)+1+ctx%kc+ctx%kv
    result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(k,dp)
    result%bic=-2.0_dp*result%log_likelihood+log(real(size(ctx%y),dp))*real(k,dp)
    result%iterations=iterations
    result%converged=(opt_status==NLME_SUCCESS .or. size(par)==0)
    result%status=merge(NLME_SUCCESS,opt_status,result%converged)
    if (opt_status==NLME_MAX_ITER) result%status=NLME_MAX_ITER
  end subroutine fill_gls_result

  subroutine gls_log_likelihood(y,x,covariance,reml,loglik,beta,sigma,status)
    real(dp), intent(in) :: y(:),x(:,:),covariance(:,:)
    logical, intent(in) :: reml
    real(dp), intent(out) :: loglik,sigma
    real(dp), allocatable, intent(out) :: beta(:)
    integer, intent(out) :: status
    real(dp), allocatable :: vinvx(:,:),vinvy(:),xtx(:,:),xty(:),r(:),vr(:)
    real(dp) :: ldv,ldx,rss,s2
    integer :: n,p,df
    n=size(y)
    p=size(x,2)
    if (size(x,1)/=n .or. any(shape(covariance)/=[n,n])) then
      allocate(beta(0))
      status=NLME_DIMENSION_ERROR
      return
    end if
    call solve_spd(covariance,x,vinvx,status,ldv)
    if(status/=NLME_SUCCESS)return
    call solve_spd(covariance,y,vinvy,status)
    if(status/=NLME_SUCCESS)return
    xtx=matmul(transpose(x),vinvx)
    xty=matmul(transpose(x),vinvy)
    call solve_spd(xtx,xty,beta,status,ldx)
    if(status/=NLME_SUCCESS)return
    r=y-matmul(x,beta)
    call solve_spd(covariance,r,vr,status)
    if(status/=NLME_SUCCESS)return
    rss=dot_product(r,vr)
    df=merge(n-p,n,reml)
    s2=rss/real(df,dp)
    sigma=sqrt(s2)
    if(reml)then
      loglik=-0.5_dp*(real(df,dp)*(log(2.0_dp*pi_dp)+1.0_dp+log(s2))+ldv+ldx)
    else
      loglik=-0.5_dp*(real(n,dp)*(log(2.0_dp*pi_dp)+1.0_dp+log(s2))+ldv)
    end if
  end subroutine gls_log_likelihood
end module nlme_gls
