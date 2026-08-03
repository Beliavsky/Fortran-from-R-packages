! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_lme
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp, pi_dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR, NLME_MAX_ITER
  use nlme_types, only : correlation_spec, variance_spec, pd_spec, lme_result, nlme_control, &
       PD_IDENT, PD_DIAG, PD_LOG_CHOL, PD_COMPOUND_SYMM
  use nlme_pdmat, only : pd_matrix, pd_parameter_count
  use nlme_correlation, only : correlation_parameter_count, correlation_to_unconstrained, &
       correlation_from_unconstrained
  use nlme_variance, only : variance_parameter_count, variance_to_unconstrained, variance_from_unconstrained
  use nlme_covariance, only : build_residual_covariance, default_group_vector, group_level_indices
  use nlme_gls, only : initialize_correlation, initialize_variance
  use nlme_linalg, only : solve_spd, inverse_spd, find_group_indices
  use nlme_optimize, only : nelder_mead
  implicit none
  private
  public :: fit_lme, initialize_pd

  type :: lme_context
    real(dp), allocatable :: y(:),x(:,:),z(:,:),time(:),var_covariate(:),coordinates(:,:)
    integer, allocatable :: group(:),var_group(:),levels(:),group_index(:)
    type(correlation_spec) :: corr_template
    type(variance_spec) :: var_template
    type(pd_spec) :: pd_template
    integer :: kp=0,ks=1,kc=0,kv=0,nlevels=1,max_group_size=1
    logical :: reml=.true.,sigma_fixed=.false.
    real(dp) :: fixed_sigma=1.0_dp
  end type lme_context
contains
  subroutine initialize_pd(spec,q,y)
    type(pd_spec), intent(inout) :: spec
    integer, intent(in) :: q
    real(dp), intent(in) :: y(:)
    integer :: k,i,j
    real(dp) :: sy
    spec%dim=q
    k=pd_parameter_count(spec)
    if (allocated(spec%par)) then
      if (size(spec%par)==k) return
      deallocate(spec%par)
    end if
    allocate(spec%par(k)); spec%par=0.0_dp
    sy=sqrt(max(1.0e-8_dp,sum((y-sum(y)/real(size(y),dp))**2)/real(max(1,size(y)-1),dp)))
    select case(spec%kind)
    case(PD_IDENT)
      spec%par(1)=log(max(1.0e-4_dp,0.35_dp*sy))
    case(PD_DIAG)
      spec%par=log(max(1.0e-4_dp,0.35_dp*sy))
    case(PD_LOG_CHOL)
      k=0
      do i=1,q
        k=k+1; spec%par(k)=log(max(1.0e-4_dp,0.35_dp*sy))
        do j=1,i-1
          k=k+1; spec%par(k)=0.0_dp
        end do
      end do
    case(PD_COMPOUND_SYMM)
      spec%par(1)=log(max(1.0e-4_dp,0.35_dp*sy)); spec%par(2)=0.0_dp
    end select
  end subroutine initialize_pd

  subroutine fit_lme(y,x,z,group,result,random,correlation,variance,time,var_covariate,var_group,coordinates,control,fixed_sigma)
    real(dp), intent(in) :: y(:),x(:,:),z(:,:)
    integer, intent(in) :: group(:)
    type(lme_result), intent(out) :: result
    type(pd_spec), intent(in), optional :: random
    type(correlation_spec), intent(in), optional :: correlation
    type(variance_spec), intent(in), optional :: variance
    real(dp), intent(in), optional :: time(:),var_covariate(:),coordinates(:,:)
    integer, intent(in), optional :: var_group(:)
    type(nlme_control), intent(in), optional :: control
    real(dp), intent(in), optional :: fixed_sigma
    type(lme_context) :: ctx
    type(pd_spec) :: pds
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    type(nlme_control) :: ctl
    real(dp), allocatable :: xp(:),xc(:),xv(:),par(:)
    real(dp) :: objective,sy
    integer :: n,q,status,it,pos

    n=size(y); q=size(z,2)
    if (size(x,1)/=n .or. size(z,1)/=n .or. size(group)/=n .or. size(x,2)<1 .or. q<1 .or. n<=size(x,2)) then
      result%status=NLME_DIMENSION_ERROR; return
    end if
    ctl=nlme_control(); if(present(control))ctl=control
    pds=pd_spec(); pds%kind=PD_LOG_CHOL; if(present(random))pds=random
    corr=correlation_spec(); if(present(correlation))corr=correlation
    var=variance_spec(); if(present(variance))var=variance
    allocate(ctx%y(n),ctx%x(n,size(x,2)),ctx%z(n,q),ctx%group(n),ctx%time(n),ctx%var_covariate(n))
    ctx%y=y; ctx%x=x; ctx%z=z; ctx%group=group
    if(present(time))then
      if(size(time)/=n)then;result%status=NLME_DIMENSION_ERROR;return;end if
      ctx%time=time
    else
      ctx%time=[(real(it,dp),it=1,n)]
    end if
    if(present(var_covariate))then
      if(size(var_covariate)/=n)then;result%status=NLME_DIMENSION_ERROR;return;end if
      ctx%var_covariate=var_covariate
    else
      ctx%var_covariate=abs(y)
    end if
    if(present(var_group))then
      if(size(var_group)/=n)then;result%status=NLME_DIMENSION_ERROR;return;end if
      allocate(ctx%var_group(n));ctx%var_group=var_group
    else
      allocate(ctx%var_group(n));ctx%var_group=1
    end if
    if(present(coordinates))then
      if(size(coordinates,1)/=n)then;result%status=NLME_DIMENSION_ERROR;return;end if
      allocate(ctx%coordinates(n,size(coordinates,2)));ctx%coordinates=coordinates
    else
      allocate(ctx%coordinates(n,0))
    end if
    call group_level_indices(ctx%group,ctx%levels,ctx%group_index,status)
    if(status/=NLME_SUCCESS)then;result%status=status;return;end if
    call initialize_pd(pds,q,y); call initialize_correlation(corr,ctx%time,ctx%group,ctx%coordinates)
    call initialize_variance(var,ctx%var_group)
    ctx%pd_template=pds;ctx%corr_template=corr;ctx%var_template=var;ctx%reml=ctl%reml
    ctx%nlevels=max(1,maxval(ctx%var_group));ctx%max_group_size=1
    do it=1,size(ctx%levels);ctx%max_group_size=max(ctx%max_group_size,count(ctx%group==ctx%levels(it)));end do
    ctx%kp=merge(0,pd_parameter_count(pds),pds%fixed)
    ctx%ks=1
    if(present(fixed_sigma))then
      if(fixed_sigma<=0.0_dp)then;result%status=NLME_INVALID_ARGUMENT;return;end if
      ctx%sigma_fixed=.true.;ctx%fixed_sigma=fixed_sigma;ctx%ks=0
    end if
    ctx%kc=merge(0,correlation_parameter_count(corr,ctx%max_group_size),corr%fixed)
    ctx%kv=merge(0,variance_parameter_count(var,ctx%nlevels),var%fixed)
    allocate(par(ctx%kp+ctx%ks+ctx%kc+ctx%kv));pos=0
    if(ctx%kp>0)then;xp=pds%par;par(1:ctx%kp)=xp;pos=ctx%kp;end if
    if(ctx%ks>0)then
      sy=sqrt(max(1.0e-8_dp,sum((y-sum(y)/real(n,dp))**2)/real(max(1,n-1),dp)))
      par(pos+1)=log(max(1.0e-4_dp,0.7_dp*sy));pos=pos+1
    end if
    if(ctx%kc>0)then
      call correlation_to_unconstrained(corr,ctx%max_group_size,xc,status)
      if(status/=NLME_SUCCESS)then;result%status=status;return;end if
      par(pos+1:pos+ctx%kc)=xc;pos=pos+ctx%kc
    end if
    if(ctx%kv>0)then
      call variance_to_unconstrained(var,ctx%nlevels,xv,status)
      if(status/=NLME_SUCCESS)then;result%status=status;return;end if
      par(pos+1:)=xv
    end if
    if(size(par)>0 .and. ctl%optimize_covariance)then
      call nelder_mead(lme_objective,ctx,par,objective,status,it,ctl%max_iter,ctl%tolerance,ctl%step,ctl%verbose)
    else
      objective=lme_objective(par,ctx);status=NLME_SUCCESS;it=0
    end if
    call fill_lme_result(par,ctx,result,status,it)
  end subroutine fit_lme

  function lme_objective(par,context) result(value)
    real(dp), intent(in) :: par(:)
    class(*), intent(inout) :: context
    real(dp) :: value
    type(pd_spec) :: pds
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    real(dp) :: sigma
    integer :: status
    select type(ctx=>context)
    type is(lme_context)
      call unpack_lme(par,ctx,pds,sigma,corr,var,status)
      if(status/=NLME_SUCCESS)then;value=huge(1.0_dp)/100.0_dp;return;end if
      call evaluate_lme(ctx,pds,sigma,corr,var,value,status)
      if(status/=NLME_SUCCESS .or. .not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
    class default
      value=huge(1.0_dp)/100.0_dp
    end select
  end function lme_objective

  subroutine unpack_lme(par,ctx,pds,sigma,corr,var,status)
    real(dp), intent(in) :: par(:)
    type(lme_context), intent(in) :: ctx
    type(pd_spec), intent(out) :: pds
    real(dp), intent(out) :: sigma
    type(correlation_spec), intent(out) :: corr
    type(variance_spec), intent(out) :: var
    integer, intent(out) :: status
    integer :: pos
    pds=ctx%pd_template;pos=0
    if(ctx%kp>0)then
      if(allocated(pds%par))deallocate(pds%par);allocate(pds%par(ctx%kp));pds%par=par(1:ctx%kp);pos=ctx%kp
    end if
    if(ctx%ks>0)then;sigma=exp(max(-30.0_dp,min(30.0_dp,par(pos+1))));pos=pos+1;else;sigma=ctx%fixed_sigma;end if
    if(ctx%kc>0)then
      call correlation_from_unconstrained(ctx%corr_template,ctx%max_group_size,par(pos+1:pos+ctx%kc),corr,status)
      if(status/=NLME_SUCCESS)return;pos=pos+ctx%kc
    else
      corr=ctx%corr_template
    end if
    if(ctx%kv>0)then
      call variance_from_unconstrained(ctx%var_template,ctx%nlevels,par(pos+1:),var,status)
    else
      var=ctx%var_template;status=NLME_SUCCESS
    end if
  end subroutine unpack_lme

  subroutine build_lme_covariance(ctx,gmat,sigma,corr,var,v,status,rbase)
    type(lme_context), intent(in) :: ctx
    real(dp), intent(in) :: gmat(:,:),sigma
    type(correlation_spec), intent(in) :: corr
    type(variance_spec), intent(in) :: var
    real(dp), allocatable, intent(out) :: v(:,:)
    integer, intent(out) :: status
    real(dp), allocatable, intent(out), optional :: rbase(:,:)
    real(dp), allocatable :: r(:,:),zi(:,:),block(:,:)
    integer, allocatable :: idx(:)
    integer :: g,n
    n=size(ctx%y)
    if(size(ctx%coordinates,2)>0)then
      call build_residual_covariance(corr,var,ctx%time,ctx%group,ctx%var_covariate,ctx%var_group,r,status,ctx%coordinates)
    else
      call build_residual_covariance(corr,var,ctx%time,ctx%group,ctx%var_covariate,ctx%var_group,r,status)
    end if
    if(status/=NLME_SUCCESS)then;allocate(v(0,0));return;end if
    allocate(v(n,n));v=sigma*sigma*r
    do g=1,size(ctx%levels)
      call find_group_indices(ctx%group,ctx%levels(g),idx)
      zi=ctx%z(idx,:);block=matmul(matmul(zi,gmat),transpose(zi))
      v(idx,idx)=v(idx,idx)+block
    end do
    if(present(rbase))then;allocate(rbase(n,n));rbase=r;end if
  end subroutine build_lme_covariance

  subroutine evaluate_lme(ctx,pds,sigma,corr,var,nll,status,beta,beta_cov,gmat_out,v_out,random_effects,&
       fitted_marginal,fitted_conditional,resid_marginal,resid_conditional,loglik)
    type(lme_context), intent(in) :: ctx
    type(pd_spec), intent(in) :: pds
    real(dp), intent(in) :: sigma
    type(correlation_spec), intent(in) :: corr
    type(variance_spec), intent(in) :: var
    real(dp), intent(out) :: nll
    integer, intent(out) :: status
    real(dp), allocatable, intent(out), optional :: beta(:),beta_cov(:,:),gmat_out(:,:),v_out(:,:),random_effects(:,:),&
         fitted_marginal(:),fitted_conditional(:),resid_marginal(:),resid_conditional(:)
    real(dp), intent(out), optional :: loglik
    real(dp), allocatable :: gmat(:,:),v(:,:),vinvx(:,:),vinvy(:),xtx(:,:),xty(:),b(:),bcov(:,:),r(:),vr(:)
    real(dp), allocatable :: bre(:,:),vi_r(:),zi(:,:),vi(:,:),ri(:),fitc(:)
    integer, allocatable :: idx(:)
    real(dp) :: ldv,ldx,quad,ll
    integer :: n,p,gg
    call pd_matrix(pds,gmat,status);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    call build_lme_covariance(ctx,gmat,sigma,corr,var,v,status);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    call solve_spd(v,ctx%x,vinvx,status,ldv);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    call solve_spd(v,ctx%y,vinvy,status);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    xtx=matmul(transpose(ctx%x),vinvx);xty=matmul(transpose(ctx%x),vinvy)
    call solve_spd(xtx,xty,b,status,ldx);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    r=ctx%y-matmul(ctx%x,b);call solve_spd(v,r,vr,status);if(status/=NLME_SUCCESS)then;nll=huge(1.0_dp);return;end if
    quad=dot_product(r,vr);n=size(ctx%y);p=size(ctx%x,2)
    if(ctx%reml)then
      ll=-0.5_dp*(real(n-p,dp)*log(2.0_dp*pi_dp)+ldv+ldx+quad)
    else
      ll=-0.5_dp*(real(n,dp)*log(2.0_dp*pi_dp)+ldv+quad)
    end if
    nll=-ll;status=NLME_SUCCESS
    if(present(beta))then;allocate(beta(size(b)));beta=b;end if
    if(present(beta_cov))then
      call inverse_spd(xtx,bcov,status);if(status/=NLME_SUCCESS)return
      allocate(beta_cov(size(b),size(b)));beta_cov=bcov
    end if
    if(present(gmat_out))then;allocate(gmat_out(size(gmat,1),size(gmat,2)));gmat_out=gmat;end if
    if(present(v_out))then;allocate(v_out(n,n));v_out=v;end if
    if(present(resid_marginal))then;allocate(resid_marginal(n));resid_marginal=r;end if
    if(present(fitted_marginal))then;allocate(fitted_marginal(n));fitted_marginal=matmul(ctx%x,b);end if
    if(present(random_effects) .or. present(fitted_conditional) .or. present(resid_conditional))then
      allocate(bre(size(ctx%levels),size(ctx%z,2)));bre=0.0_dp;allocate(fitc(n));fitc=matmul(ctx%x,b)
      do gg=1,size(ctx%levels)
        call find_group_indices(ctx%group,ctx%levels(gg),idx)
        zi=ctx%z(idx,:);vi=v(idx,idx);ri=r(idx)
        call solve_spd(vi,ri,vi_r,status);if(status/=NLME_SUCCESS)return
        bre(gg,:)=matmul(gmat,matmul(transpose(zi),vi_r))
        fitc(idx)=fitc(idx)+matmul(zi,bre(gg,:))
      end do
      if(present(random_effects))then;allocate(random_effects(size(bre,1),size(bre,2)));random_effects=bre;end if
      if(present(fitted_conditional))then;allocate(fitted_conditional(n));fitted_conditional=fitc;end if
      if(present(resid_conditional))then;allocate(resid_conditional(n));resid_conditional=ctx%y-fitc;end if
    end if
    if(present(loglik))loglik=ll
  end subroutine evaluate_lme

  subroutine fill_lme_result(par,ctx,result,opt_status,iterations)
    real(dp), intent(in) :: par(:)
    type(lme_context), intent(in) :: ctx
    type(lme_result), intent(out) :: result
    integer, intent(in) :: opt_status,iterations
    type(pd_spec) :: pds
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    real(dp) :: sigma,nll
    integer :: status,k
    call unpack_lme(par,ctx,pds,sigma,corr,var,status);if(status/=NLME_SUCCESS)then;result%status=status;return;end if
    call evaluate_lme(ctx,pds,sigma,corr,var,nll,status,result%beta,result%beta_cov,result%random_covariance,&
         result%covariance,result%random_effects,result%fitted_marginal,result%fitted_conditional,&
         result%residuals_marginal,result%residuals_conditional,result%log_likelihood)
    if(status/=NLME_SUCCESS)then;result%status=status;return;end if
    allocate(result%group_levels(size(ctx%levels)));result%group_levels=ctx%levels
    if(allocated(corr%par))then;allocate(result%correlation_parameters(size(corr%par)));result%correlation_parameters=corr%par
    else;allocate(result%correlation_parameters(0));end if
    if(allocated(var%par))then;allocate(result%variance_parameters(size(var%par)));result%variance_parameters=var%par
    else;allocate(result%variance_parameters(0));end if
    result%sigma=sigma;k=size(result%beta)+ctx%kp+ctx%ks+ctx%kc+ctx%kv
    result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(k,dp)
    result%bic=-2.0_dp*result%log_likelihood+log(real(size(ctx%y),dp))*real(k,dp)
    result%iterations=iterations;result%converged=(opt_status==NLME_SUCCESS .or. size(par)==0)
    result%status=merge(NLME_SUCCESS,opt_status,result%converged)
    if(opt_status==NLME_MAX_ITER)result%status=NLME_MAX_ITER
  end subroutine fill_lme_result
end module nlme_lme
