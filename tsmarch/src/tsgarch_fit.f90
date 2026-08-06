! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_fit_module
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ghyp_kinds, only : dp
  use tsd_types, only : parameter_specification
  use tsd_fit, only : distribution_modelspec
  use tsd_optimize, only : nelder_mead, parameters_to_unconstrained, unconstrained_to_parameters
  use tsd_math, only : invert_matrix
  use tsgarch_types
  use tsgarch_model, only : filter_garch, initialize_parameters, persistence, validate_specification
  implicit none
  private

  type :: fit_context
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: vreg(:, :)
    logical :: has_vreg = .false.
    type(garch_spec) :: spec
    type(garch_parameters) :: template
    real(dp), allocatable :: lower(:), upper(:)
    real(dp) :: penalty_scale = 1.0e8_dp
  end type fit_context

  public :: estimate_garch, pack_parameters, unpack_parameters, parameter_bounds
  public :: numerical_inference, garch_loglikelihood

contains

  integer function parameter_count(spec, nreg) result(n)
    type(garch_spec), intent(in) :: spec
    integer, intent(in) :: nreg
    character(len=:), allocatable :: model
    integer :: total_ab
    model = trim(spec%model)
    n = 0
    if (spec%constant) n = n + 1
    if (.not. spec%variance_targeting .and. model /= 'ewma') n = n + 1
    total_ab = spec%p + spec%q
    if (model == 'igarch') then
      n = n + max(0, total_ab - 1)
    else if (model == 'ewma') then
      n = n + 1
    else
      n = n + spec%p + spec%q
    end if
    select case(model)
    case('gjrgarch','egarch','aparch','fgarch')
      n = n + spec%p
    end select
    if (model == 'fgarch') n = n + spec%p
    if (model == 'aparch' .or. model == 'fgarch') n = n + 1
    if (model == 'cgarch') n = n + 2
    n = n + max(0,nreg)
    select case(trim(spec%distribution))
    case('snorm')
      n=n+1
    case('std','ged')
      n=n+1
    case('sstd','sged','nig','jsu','ghst')
      n=n+2
    case('gh')
      n=n+3
    end select
  end function parameter_count

  subroutine pack_parameters(spec, par, packed, names)
    type(garch_spec), intent(in) :: spec
    type(garch_parameters), intent(in) :: par
    real(dp), allocatable, intent(out) :: packed(:)
    character(len=20), allocatable, intent(out), optional :: names(:)
    integer :: n, k, i, total_ab, omit_index
    character(len=:), allocatable :: model
    n = parameter_count(spec,size(par%xi))
    allocate(packed(n))
    if (present(names)) allocate(names(n))
    k=0
    model=trim(spec%model)
    if(spec%constant) then
      k=k+1
      packed(k)=par%mu
      if(present(names))names(k)='mu'
    end if
    if(.not.spec%variance_targeting .and. model/='ewma') then
      k=k+1
      packed(k)=par%omega
      if(present(names))names(k)='omega'
    end if
    total_ab=spec%p+spec%q
    if(model=='igarch') then
      omit_index=total_ab
      do i=1,spec%p
        if(i==omit_index) cycle
        k=k+1
        packed(k)=par%alpha(i)
        if(present(names))write(names(k),'("alpha",i0)')i
      end do
      do i=1,spec%q
        if(spec%p+i==omit_index) cycle
        k=k+1
        packed(k)=par%beta(i)
        if(present(names))write(names(k),'("beta",i0)')i
      end do
    else if(model=='ewma') then
      k=k+1
      packed(k)=par%beta(1)
      if(present(names))names(k)='lambda'
    else
      do i=1,spec%p
        k=k+1
        packed(k)=par%alpha(i)
        if(present(names))write(names(k),'("alpha",i0)')i
      end do
      do i=1,spec%q
        k=k+1
        packed(k)=par%beta(i)
        if(present(names))write(names(k),'("beta",i0)')i
      end do
    end if
    select case(model)
    case('gjrgarch','egarch','aparch','fgarch')
      do i=1,spec%p
        k=k+1
        packed(k)=par%gamma(i)
        if(present(names))write(names(k),'("gamma",i0)')i
      end do
    end select
    if(model=='fgarch') then
      do i=1,spec%p
        k=k+1
        packed(k)=par%eta(i)
        if(present(names))write(names(k),'("eta",i0)')i
      end do
    end if
    if(model=='aparch'.or.model=='fgarch') then
      k=k+1
      packed(k)=par%delta
      if(present(names))names(k)='delta'
    end if
    if(model=='cgarch') then
      k=k+1
      packed(k)=par%rho
      if(present(names))names(k)='rho'
      k=k+1
      packed(k)=par%phi
      if(present(names))names(k)='phi'
    end if
    do i=1,size(par%xi)
      k=k+1
      packed(k)=par%xi(i)
      if(present(names))write(names(k),'("xi",i0)')i
    end do
    select case(trim(spec%distribution))
    case('snorm')
      k=k+1
      packed(k)=par%dist%skew
      if(present(names))names(k)='skew'
    case('std','ged')
      k=k+1
      packed(k)=par%dist%shape
      if(present(names))names(k)='shape'
    case('sstd','sged','nig','jsu','ghst')
      k=k+1
      packed(k)=par%dist%skew
      if(present(names))names(k)='skew'
      k=k+1
      packed(k)=par%dist%shape
      if(present(names))names(k)='shape'
    case('gh')
      k=k+1
      packed(k)=par%dist%skew
      if(present(names))names(k)='skew'
      k=k+1
      packed(k)=par%dist%shape
      if(present(names))names(k)='shape'
      k=k+1
      packed(k)=par%dist%lambda
      if(present(names))names(k)='lambda'
    end select
  end subroutine pack_parameters

  subroutine unpack_parameters(spec, packed, template, par, status)
    type(garch_spec), intent(in) :: spec
    real(dp), intent(in) :: packed(:)
    type(garch_parameters), intent(in) :: template
    type(garch_parameters), intent(out) :: par
    integer, intent(out) :: status
    integer :: k,i,total_ab,omit_index
    real(dp) :: remainder
    character(len=:),allocatable :: model
    par=template
    k=0
    model=trim(spec%model)
    status=tsg_invalid_argument
    if(size(packed)/=parameter_count(spec,size(template%xi))) return
    if(spec%constant) then
    k=k+1
    par%mu=packed(k)
    end if
    if(.not.spec%variance_targeting.and.model/='ewma') then
    k=k+1
    par%omega=packed(k)
    end if
    total_ab=spec%p+spec%q
    if(model=='igarch') then
      omit_index=total_ab
      par%alpha=0.0_dp
      par%beta=0.0_dp
      do i=1,spec%p
        if(i==omit_index)cycle
        k=k+1
        par%alpha(i)=packed(k)
      end do
      do i=1,spec%q
        if(spec%p+i==omit_index)cycle
        k=k+1
        par%beta(i)=packed(k)
      end do
      remainder=1.0_dp-sum(par%alpha)-sum(par%beta)
      if(spec%q>0)then
      par%beta(spec%q)=remainder
      else if(spec%p>0)then
      par%alpha(spec%p)=remainder
      end if
    else if(model=='ewma') then
      k=k+1
      par%beta(1)=packed(k)
      par%alpha(1)=1.0_dp-packed(k)
      par%omega=0.0_dp
    else
      do i=1,spec%p
      k=k+1
      par%alpha(i)=packed(k)
      end do
      do i=1,spec%q
      k=k+1
      par%beta(i)=packed(k)
      end do
    end if
    select case(model)
    case('gjrgarch','egarch','aparch','fgarch')
      do i=1,spec%p
      k=k+1
      par%gamma(i)=packed(k)
      end do
    end select
    if(model=='fgarch')then
      do i=1,spec%p
      k=k+1
      par%eta(i)=packed(k)
      end do
    end if
    if(model=='aparch'.or.model=='fgarch')then
    k=k+1
    par%delta=packed(k)
    end if
    if(model=='cgarch')then
    k=k+1
    par%rho=packed(k)
    k=k+1
    par%phi=packed(k)
    end if
    do i=1,size(par%xi)
    k=k+1
    par%xi(i)=packed(k)
    end do
    select case(trim(spec%distribution))
    case('snorm')
      k=k+1
      par%dist%skew=packed(k)
    case('std','ged')
      k=k+1
      par%dist%shape=packed(k)
    case('sstd','sged','nig','jsu','ghst')
      k=k+1
      par%dist%skew=packed(k)
      k=k+1
      par%dist%shape=packed(k)
    case('gh')
      k=k+1
      par%dist%skew=packed(k)
      k=k+1
      par%dist%shape=packed(k)
      k=k+1
      par%dist%lambda=packed(k)
    end select
    status=tsg_success
  end subroutine unpack_parameters

  subroutine parameter_bounds(y,spec,template,lower,upper)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::template
    real(dp),allocatable,intent(out)::lower(:),upper(:)
    type(parameter_specification)::ds
    real(dp),allocatable::p(:)
    integer::n,k,i,total_ab,omit_index
    real(dp)::m,v,sd
    character(len=:),allocatable::model
    call pack_parameters(spec,template,p)
    n=size(p)
    allocate(lower(n),upper(n))
    k=0
    model=trim(spec%model)
    m=sum(y)/real(max(1,size(y)),dp)
    v=sum((y-m)**2)/real(max(1,size(y)),dp)
    sd=sqrt(max(v,1.0e-8_dp))
    if(spec%constant)then
    k=k+1
    lower(k)=m-10.0_dp*sd
    upper(k)=m+10.0_dp*sd
    end if
    if(.not.spec%variance_targeting.and.model/='ewma')then
      k=k+1
      if(spec%multiplicative.or.model=='egarch')then
      lower(k)=-20.0_dp
      upper(k)=20.0_dp
      else
      lower(k)=1.0e-12_dp
      upper(k)=max(100.0_dp*v,1.0_dp)
      end if
    end if
    total_ab=spec%p+spec%q
    if(model=='igarch')then
      omit_index=total_ab
      do i=1,spec%p
      if(i==omit_index)cycle
      k=k+1
      lower(k)=1.0e-8_dp
      upper(k)=0.999_dp
      end do
      do i=1,spec%q
      if(spec%p+i==omit_index)cycle
      k=k+1
      lower(k)=1.0e-8_dp
      upper(k)=0.999_dp
      end do
    else if(model=='ewma')then
      k=k+1
      lower(k)=0.001_dp
      upper(k)=0.9999_dp
    else
      do i=1,spec%p
      k=k+1
      lower(k)=0.0_dp
      upper(k)=0.999_dp
      end do
      do i=1,spec%q
      k=k+1
      lower(k)=0.0_dp
      upper(k)=0.999_dp
      end do
    end if
    select case(model)
    case('gjrgarch')
      do i=1,spec%p
      k=k+1
      lower(k)=-0.999_dp
      upper(k)=0.999_dp
      end do
    case('egarch')
      do i=1,spec%p
      k=k+1
      lower(k)=-3.0_dp
      upper(k)=3.0_dp
      end do
    case('aparch','fgarch')
      do i=1,spec%p
      k=k+1
      lower(k)=-0.999_dp
      upper(k)=0.999_dp
      end do
    end select
    if(model=='fgarch')then
      do i=1,spec%p
      k=k+1
      lower(k)=-3.0_dp
      upper(k)=3.0_dp
      end do
    end if
    if(model=='aparch'.or.model=='fgarch')then
    k=k+1
    lower(k)=0.2_dp
    upper(k)=4.0_dp
    end if
    if(model=='cgarch')then
      k=k+1
      lower(k)=0.01_dp
      upper(k)=0.999_dp
      k=k+1
      lower(k)=-0.999_dp
      upper(k)=0.999_dp
    end if
    do i=1,size(template%xi)
    k=k+1
    lower(k)=-20.0_dp
    upper(k)=20.0_dp
    end do
    ds=distribution_modelspec([ -1.0_dp,1.0_dp ],spec%distribution)
    select case(trim(spec%distribution))
    case('snorm')
      k=k+1
      lower(k)=ds%lower%skew
      upper(k)=ds%upper%skew
    case('std','ged')
      k=k+1
      lower(k)=ds%lower%shape
      upper(k)=ds%upper%shape
    case('sstd','sged','nig','jsu','ghst')
      k=k+1
      lower(k)=ds%lower%skew
      upper(k)=ds%upper%skew
      k=k+1
      lower(k)=ds%lower%shape
      upper(k)=ds%upper%shape
    case('gh')
      k=k+1
      lower(k)=ds%lower%skew
      upper(k)=ds%upper%skew
      k=k+1
      lower(k)=ds%lower%shape
      upper(k)=ds%upper%shape
      k=k+1
      lower(k)=ds%lower%lambda
      upper(k)=ds%upper%lambda
    end select
  end subroutine parameter_bounds

  function estimate_garch(y,specification,start,vreg,options) result(fit)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::specification
    type(garch_parameters),intent(in),optional::start
    real(dp),intent(in),optional::vreg(:,:)
    type(fit_options),intent(in),optional::options
    type(garch_fit)::fit
    type(garch_spec)::spec
    type(garch_parameters)::initial,estimated
    type(fit_options)::opt
    type(fit_context)::ctx
    real(dp),allocatable::natural(:),z(:),lower(:),upper(:)
    real(dp)::fval
    integer::status,opt_status,iters,nreg
    character(len=240)::message
    spec=specification
    call validate_specification(spec,status,message)
    if(status/=tsg_success)then
    fit%message=message
    return
    end if
    if(size(y)<max(20,2*max(spec%p,spec%q)+5))then
    fit%message='too few observations'
    return
    end if
    nreg=0
    if(present(vreg))then
      if(size(vreg,1)/=size(y))then
      fit%message='variance regressor row count mismatch'
      return
      end if
      nreg=size(vreg,2)
    end if
    if(present(start))then
    initial=start
    else
    initial=initialize_parameters(y,spec,nreg)
    end if
    if(.not.allocated(initial%xi))then
    allocate(initial%xi(nreg))
    initial%xi=0.0_dp
    end if
    if(size(initial%xi)/=nreg)then
    fit%message='starting variance coefficients do not conform'
    return
    end if
    call pack_parameters(spec,initial,natural,fit%parameter_names)
    call parameter_bounds(y,spec,initial,lower,upper)
    call parameters_to_unconstrained(natural,lower,upper,z,status)
    if(status/=0)then
    fit%message='could not transform starting parameters'
    return
    end if
    ctx%y=y
    ctx%spec=spec
    ctx%template=initial
    ctx%lower=lower
    ctx%upper=upper
    if(present(vreg))then
    ctx%vreg=vreg
    ctx%has_vreg=.true.
    end if
    opt=fit_options()
    if(present(options))opt=options
    call nelder_mead(garch_objective,ctx,z,fval,opt_status,iters,opt%max_iterations,opt%tolerance,opt%simplex_scale)
    call unconstrained_to_parameters(z,lower,upper,natural)
    call unpack_parameters(spec,natural,initial,estimated,status)
    fit%spec=spec
    fit%parameters=estimated
    fit%packed_parameters=natural
    fit%npars=size(natural)
    fit%iterations=iters
    if(ctx%has_vreg)then
    fit%filtered=filter_garch(y,spec,estimated,ctx%vreg)
    else
    fit%filtered=filter_garch(y,spec,estimated)
    end if
    if(fit%filtered%status/=tsg_success)then
    fit%message=fit%filtered%message
    fit%status=tsg_numerical_failure
    return
    end if
    fit%log_likelihood=fit%filtered%log_likelihood
    fit%aic=-2.0_dp*fit%log_likelihood+2.0_dp*real(fit%npars,dp)
    fit%bic=-2.0_dp*fit%log_likelihood+log(real(size(y),dp))*real(fit%npars,dp)
    fit%status=merge(tsg_success,tsg_no_convergence,opt_status==0)
    if(fit%status==tsg_success)then
    fit%message='ok'
    else
    fit%message='optimizer reached its iteration limit'
    end if
    if(opt%compute_inference)then
      if(ctx%has_vreg)then
      call numerical_inference(y,spec,estimated,fit,ctx%vreg)
      else
      call numerical_inference(y,spec,estimated,fit)
      end if
    end if
  end function estimate_garch

  real(dp) function garch_objective(z,context) result(value)
    real(dp),intent(in)::z(:)
    class(*),intent(inout)::context
    real(dp),allocatable::natural(:)
    type(garch_parameters)::par
    type(garch_filter_result)::filtered
    integer::status
    real(dp)::pen
    select type(ctx=>context)
    type is(fit_context)
      allocate(natural(size(z)))
      call unconstrained_to_parameters(z,ctx%lower,ctx%upper,natural)
      call unpack_parameters(ctx%spec,natural,ctx%template,par,status)
      if(status/=tsg_success)then
      value=huge(1.0_dp)/100.0_dp
      return
      end if
      pen=constraint_penalty(ctx%spec,par,ctx%penalty_scale)
      if(pen>1.0e50_dp)then
      value=pen
      return
      end if
      if(ctx%has_vreg)then
      filtered=filter_garch(ctx%y,ctx%spec,par,ctx%vreg)
      else
      filtered=filter_garch(ctx%y,ctx%spec,par)
      end if
      if(filtered%status/=tsg_success)then
      value=huge(1.0_dp)/100.0_dp
      else
      value=-filtered%log_likelihood+pen
      end if
    class default
      value=huge(1.0_dp)/100.0_dp
    end select
  end function garch_objective

  real(dp) function constraint_penalty(spec,par,scale) result(value)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in)::scale
    real(dp)::p,viol,beta_sum
    integer::i
    value=0.0_dp
    p=persistence(spec,par)
    if(trim(spec%model)/='igarch'.and.trim(spec%model)/='ewma')then
      viol=max(0.0_dp,p-spec%stationarity_limit)
      value=value+scale*viol*viol
    end if
    if(trim(spec%model)=='igarch')then
      if(any(par%alpha<0.0_dp).or.any(par%beta<0.0_dp))value=huge(1.0_dp)/100.0_dp
    end if
    if(trim(spec%model)=='gjrgarch')then
      do i=1,size(par%alpha)
      viol=max(0.0_dp,-par%alpha(i)-par%gamma(i))
      value=value+scale*viol*viol
      end do
    end if
    if(trim(spec%model)=='cgarch')then
      beta_sum=sum(par%beta)
      viol=max(0.0_dp,sum(par%alpha)+beta_sum-par%rho)
      value=value+scale*viol*viol
      viol=max(0.0_dp,par%rho-spec%stationarity_limit)
      value=value+scale*viol*viol
      viol=max(0.0_dp,par%phi-beta_sum)
      value=value+scale*viol*viol
    end if
  end function constraint_penalty

  real(dp) function garch_loglikelihood(y,spec,par,vreg) result(value)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in),optional::vreg(:,:)
    type(garch_filter_result)::filtered
    if(present(vreg))then
    filtered=filter_garch(y,spec,par,vreg)
    else
    filtered=filter_garch(y,spec,par)
    end if
    if(filtered%status==tsg_success)then
    value=filtered%log_likelihood
    else
    value=-huge(1.0_dp)
    end if
  end function garch_loglikelihood

  subroutine numerical_inference(y,spec,par,fit,vreg)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    type(garch_fit),intent(inout)::fit
    real(dp),intent(in),optional::vreg(:,:)
    real(dp),allocatable::theta(:),lower(:),upper(:),hessian(:,:),inv(:,:),scores(:,:),lp(:),lm(:),lpp(:),lpm(:),lmp(:),lmm(:)
    type(garch_parameters)::pp,pm,p2
    type(garch_filter_result)::fp,fm,fpp,fpm,fmp,fmm
    integer::npar,n,i,j,status
    real(dp)::hi,hj,f0
    logical::ok
    call pack_parameters(spec,par,theta)
    call parameter_bounds(y,spec,par,lower,upper)
    npar=size(theta)
    n=size(y)
    allocate(hessian(npar,npar),scores(n,npar))
    hessian=0.0_dp
    scores=0.0_dp
    f0=-garch_loglikelihood(y,spec,par,vreg)
    do i=1,npar
      hi=max(1.0e-5_dp,1.0e-4_dp*max(1.0_dp,abs(theta(i))))
      hi=min(hi,0.24_dp*max(upper(i)-lower(i),1.0e-8_dp))
      allocate(lp(npar),lm(npar))
      lp=theta
      lm=theta
      lp(i)=min(upper(i)-1e-10_dp,theta(i)+hi)
      lm(i)=max(lower(i)+1e-10_dp,theta(i)-hi)
      hi=0.5_dp*(lp(i)-lm(i))
      call unpack_parameters(spec,lp,par,pp,status)
      call unpack_parameters(spec,lm,par,pm,status)
      if(present(vreg))then
      fp=filter_garch(y,spec,pp,vreg)
      fm=filter_garch(y,spec,pm,vreg)
      else
      fp=filter_garch(y,spec,pp)
      fm=filter_garch(y,spec,pm)
      end if
      if(fp%status==tsg_success.and.fm%status==tsg_success.and.hi>0.0_dp)then
        scores(:,i)=(fp%loglik_vector-fm%loglik_vector)/(2.0_dp*hi)
        hessian(i,i)=(-fp%log_likelihood-2.0_dp*f0-fm%log_likelihood)/(hi*hi)
      else
        hessian(i,i)=1.0e8_dp
      end if
      deallocate(lp,lm)
    end do
    do i=1,npar
      hi=max(1.0e-5_dp,1.0e-4_dp*max(1.0_dp,abs(theta(i))))
      do j=i+1,npar
        hj=max(1.0e-5_dp,1.0e-4_dp*max(1.0_dp,abs(theta(j))))
        allocate(lpp(npar),lpm(npar),lmp(npar),lmm(npar))
        lpp=theta
        lpm=theta
        lmp=theta
        lmm=theta
        lpp(i)=min(upper(i)-1e-10_dp,theta(i)+hi)
        lpp(j)=min(upper(j)-1e-10_dp,theta(j)+hj)
        lpm(i)=lpp(i)
        lpm(j)=max(lower(j)+1e-10_dp,theta(j)-hj)
        lmp(i)=max(lower(i)+1e-10_dp,theta(i)-hi)
        lmp(j)=lpp(j)
        lmm(i)=lmp(i)
        lmm(j)=lpm(j)
        hi=0.5_dp*(lpp(i)-lmp(i))
        hj=0.5_dp*(lpp(j)-lpm(j))
        call unpack_parameters(spec,lpp,par,p2,status)
        fpp=filter_optional(y,spec,p2,vreg)
        call unpack_parameters(spec,lpm,par,p2,status)
        fpm=filter_optional(y,spec,p2,vreg)
        call unpack_parameters(spec,lmp,par,p2,status)
        fmp=filter_optional(y,spec,p2,vreg)
        call unpack_parameters(spec,lmm,par,p2,status)
        fmm=filter_optional(y,spec,p2,vreg)
        if(min(fpp%status,fpm%status,fmp%status,fmm%status)==tsg_success.and.hi>0.0_dp.and.hj>0.0_dp)then
          hessian(i,j)=(-fpp%log_likelihood+fpm%log_likelihood+fmp%log_likelihood-fmm%log_likelihood)/(4.0_dp*hi*hj)
          hessian(j,i)=hessian(i,j)
        end if
        deallocate(lpp,lpm,lmp,lmm)
      end do
    end do
    call invert_matrix(hessian,inv,ok)
    fit%hessian=hessian
    fit%scores=scores
    if(ok)then
      fit%covariance=inv
      allocate(fit%standard_errors(npar))
      do i=1,npar
      fit%standard_errors(i)=sqrt(max(inv(i,i),0.0_dp))
      end do
    else
      allocate(fit%covariance(npar,npar),fit%standard_errors(npar))
      fit%covariance=0.0_dp
      fit%standard_errors=huge(1.0_dp)
    end if
  end subroutine numerical_inference

  function filter_optional(y,spec,par,vreg) result(out)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in),optional::vreg(:,:)
    type(garch_filter_result)::out
    if(present(vreg))then
    out=filter_garch(y,spec,par,vreg)
    else
    out=filter_garch(y,spec,par)
    end if
  end function filter_optional

end module tsgarch_fit_module
