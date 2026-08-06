! SPDX-License-Identifier: GPL-2.0-only
module tsd_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ghyp_kinds, only : dp
  use tsd_types, only : distribution_parameters, parameter_specification, distribution_fit, &
                        valid_distribution, distribution_id, dist_norm, dist_std, dist_snorm, &
                        dist_sstd, dist_ged, dist_sged, dist_nig, dist_gh, dist_jsu, dist_ghst, &
                        tsd_success, tsd_invalid_argument, tsd_singular
  use tsd_math, only : mean_value, sample_sd, invert_matrix
  use tsd_distributions, only : ddist
  use tsd_optimize, only : nelder_mead, parameters_to_unconstrained, unconstrained_to_parameters
  implicit none
  private

  type :: likelihood_context
    real(dp), allocatable :: y(:)
    type(parameter_specification) :: spec
    real(dp), allocatable :: lower(:), upper(:)
  end type likelihood_context

  public :: distribution_modelspec, distribution_bounds, estimate_distribution
  public :: log_likelihood, information_covariance

contains

  function distribution_modelspec(y, distribution) result(spec)
    real(dp), intent(in) :: y(:)
    character(len=*), intent(in), optional :: distribution
    type(parameter_specification) :: spec
    character(len=8) :: name
    real(dp) :: s
    name='norm'
    if(present(distribution)) name=adjustl(distribution)
    spec%distribution=name
    if(size(y)>0)then
      spec%parameters%mu=mean_value(y)
      s=sample_sd(y)
      if(.not.ieee_is_finite(s).or.s<=0.0_dp)s=max(1.0_dp,abs(spec%parameters%mu)*0.1_dp)
    else
      s=1.0_dp
    end if
    spec%parameters%sigma=s
    spec%lower%mu=-huge(1.0_dp)
    spec%upper%mu=huge(1.0_dp)
    spec%lower%sigma=1.0e-14_dp
    spec%upper%sigma=max(10.0_dp*s,1.0e-6_dp)
    spec%lower%skew=0.0_dp
    spec%upper%skew=0.0_dp
    spec%lower%shape=0.0_dp
    spec%upper%shape=0.0_dp
    spec%lower%lambda=-6.0_dp
    spec%upper%lambda=6.0_dp
    spec%estimate_mu=.true.
    spec%estimate_sigma=.true.
    spec%estimate_skew=.false.
    spec%estimate_shape=.false.
    spec%estimate_lambda=.false.
    select case(distribution_id(name))
    case(dist_norm)
      spec%parameters%skew=0.0_dp
      spec%parameters%shape=0.0_dp
      spec%parameters%lambda=-0.5_dp
    case(dist_ged)
      spec%parameters%shape=2.0_dp
      spec%lower%shape=0.1_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_std)
      spec%parameters%shape=4.0_dp
      spec%lower%shape=2.01_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_snorm)
      spec%parameters%skew=0.5_dp
      spec%lower%skew=0.1_dp
      spec%upper%skew=10.0_dp
      spec%estimate_skew=.true.
    case(dist_sged)
      spec%parameters%skew=1.0_dp
      spec%lower%skew=0.01_dp
      spec%upper%skew=30.0_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=2.0_dp
      spec%lower%shape=0.1_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_sstd)
      spec%parameters%skew=1.0_dp
      spec%lower%skew=0.01_dp
      spec%upper%skew=30.0_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=4.0_dp
      spec%lower%shape=2.01_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_nig)
      spec%parameters%skew=0.2_dp
      spec%lower%skew=-0.99_dp
      spec%upper%skew=0.99_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=0.4_dp
      spec%lower%shape=0.01_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_gh)
      spec%parameters%skew=0.2_dp
      spec%lower%skew=-0.99_dp
      spec%upper%skew=0.99_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=2.0_dp
      spec%lower%shape=0.25_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
      spec%parameters%lambda=-0.5_dp
      spec%lower%lambda=-30.0_dp
      spec%upper%lambda=30.0_dp
      spec%estimate_lambda=.true.
    case(dist_jsu)
      spec%parameters%skew=0.0_dp
      spec%lower%skew=-20.0_dp
      spec%upper%skew=20.0_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=1.0_dp
      spec%lower%shape=0.1_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case(dist_ghst)
      spec%parameters%skew=0.2_dp
      spec%lower%skew=-80.0_dp
      spec%upper%skew=80.0_dp
      spec%estimate_skew=.true.
      spec%parameters%shape=8.2_dp
      spec%lower%shape=4.01_dp
      spec%upper%shape=100.0_dp
      spec%estimate_shape=.true.
    case default
      spec%distribution='invalid'
    end select
  end function distribution_modelspec

  subroutine distribution_bounds(distribution, lower, upper, status)
    character(len=*), intent(in) :: distribution
    type(distribution_parameters), intent(out) :: lower, upper
    integer, intent(out), optional :: status
    type(parameter_specification) :: spec
    real(dp) :: dummy(2)
    dummy=[-1.0_dp,1.0_dp]
    spec=distribution_modelspec(dummy,distribution)
    lower=spec%lower
    upper=spec%upper
    if(present(status))then
      if(valid_distribution(distribution))then
      status=tsd_success
      else
      status=tsd_invalid_argument
      end if
    end if
  end subroutine distribution_bounds

  function estimate_distribution(y, specification, max_iterations, tolerance, use_hessian) result(fit)
    real(dp), intent(in) :: y(:)
    type(parameter_specification), intent(in), optional :: specification
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: use_hessian
    type(distribution_fit) :: fit
    type(likelihood_context) :: ctx
    type(parameter_specification) :: spec
    real(dp), allocatable :: p0(:),lo(:),hi(:),z(:),popt(:)
    real(dp) :: fval,tol
    integer :: status,iters,maxit,k
    logical :: hess

    if(size(y)<2 .or. any(.not.ieee_is_finite(y)))then
      fit%status=tsd_invalid_argument
      fit%message='at least two finite observations are required'
      return
    end if
    if(present(specification))then
    spec=specification
    else
    spec=distribution_modelspec(y,'norm')
    end if
    if(.not.valid_distribution(spec%distribution))then
      fit%status=tsd_invalid_argument
      fit%message='invalid distribution'
      return
    end if
    call pack_parameters(spec,spec%parameters,p0,lo,hi)
    k=size(p0)
    if(k==0)then
      fit%parameters=spec%parameters
      fit%distribution=spec%distribution
      fit%nobs=size(y)
      fit%degrees_of_freedom=0
      fit%negative_log_likelihood=-log_likelihood(y,spec%distribution,spec%parameters)
      fit%log_likelihood=-fit%negative_log_likelihood
      fit%aic=-2.0_dp*fit%log_likelihood
      fit%bic=-2.0_dp*fit%log_likelihood
      fit%status=tsd_success
      fit%message='no free parameters'
      return
    end if
    call parameters_to_unconstrained(p0,lo,hi,z,status)
    if(status/=tsd_success)then
    fit%status=status
    fit%message='invalid parameter bounds'
    return
    end if
    ctx%y=y
    ctx%spec=spec
    ctx%lower=lo
    ctx%upper=hi
    maxit=1000
    if(present(max_iterations))maxit=max_iterations
    tol=1.0e-7_dp
    if(present(tolerance))tol=tolerance
    call nelder_mead(likelihood_objective,ctx,z,fval,status,iters,maxit,tol,0.15_dp)
    allocate(popt(k))
    call unconstrained_to_parameters(z,lo,hi,popt)
    call unpack_parameters(spec,popt,fit%parameters)
    fit%distribution=spec%distribution
    fit%estimated=[spec%estimate_mu,spec%estimate_sigma,spec%estimate_skew,spec%estimate_shape,spec%estimate_lambda]
    fit%negative_log_likelihood=fval
    fit%log_likelihood=-fval
    fit%nobs=size(y)
    fit%degrees_of_freedom=k
    fit%iterations=iters
    fit%aic=2.0_dp*real(k,dp)-2.0_dp*fit%log_likelihood
    fit%bic=log(real(size(y),dp))*real(k,dp)-2.0_dp*fit%log_likelihood
    fit%status=status
    if(status==tsd_success)then
    fit%message='converged'
    else
    fit%message='optimizer reached its iteration limit'
    end if
    hess=.true.
    if(present(use_hessian))hess=use_hessian
    call numerical_inference(y,spec,fit%parameters,fit,hess)
  end function estimate_distribution

  real(dp) function log_likelihood(y,distribution,parameters) result(value)
    real(dp),intent(in)::y(:)
    character(len=*),intent(in)::distribution
    type(distribution_parameters),intent(in)::parameters
    integer::i
    real(dp)::ld
    value=0.0_dp
    do i=1,size(y)
      ld=ddist(distribution,y(i),parameters,.true.)
      if(.not.ieee_is_finite(ld))then
      value=-huge(1.0_dp)/100.0_dp
      return
      end if
      value=value+ld
    end do
  end function log_likelihood

  real(dp) function likelihood_objective(z,context) result(value)
    real(dp),intent(in)::z(:)
    class(*),intent(inout)::context
    real(dp),allocatable::p(:)
    type(distribution_parameters)::pars
    select type(ctx=>context)
    type is(likelihood_context)
      allocate(p(size(z)))
      call unconstrained_to_parameters(z,ctx%lower,ctx%upper,p)
      call unpack_parameters(ctx%spec,p,pars)
      value=-log_likelihood(ctx%y,ctx%spec%distribution,pars)
    class default
      value=huge(1.0_dp)/100.0_dp
    end select
  end function likelihood_objective

  subroutine numerical_inference(y,spec,pars,fit,compute_hessian)
    real(dp),intent(in)::y(:)
    type(parameter_specification),intent(in)::spec
    type(distribution_parameters),intent(in)::pars
    type(distribution_fit),intent(inout)::fit
    logical,intent(in)::compute_hessian
    real(dp),allocatable::p(:),lo(:),hi(:),pp(:),pm(:),fplus(:),fminus(:),opg(:,:),invh(:,:),invo(:,:)
    real(dp)::h1,h2,f0,fpp,fpm,fmp,fmm
    integer::i,j,k,n
    logical::okh,oko
    call pack_parameters(spec,pars,p,lo,hi)
    k=size(p)
    n=size(y)
    allocate(fit%gradient(k),fit%scores(n,k))
    fit%gradient=0.0_dp
    fit%scores=0.0_dp
    allocate(pp(k),pm(k),fplus(n),fminus(n))
    do j=1,k
      h1=step_size(p(j),lo(j),hi(j))
      pp=p
      pm=p
      pp(j)=min(hi(j)-epsilon(1.0_dp),p(j)+h1)
      pm(j)=max(lo(j)+epsilon(1.0_dp),p(j)-h1)
      h1=0.5_dp*(pp(j)-pm(j))
      call observation_logdensities(y,spec,pp,fplus)
      call observation_logdensities(y,spec,pm,fminus)
      fit%scores(:,j)=(fplus-fminus)/(2.0_dp*h1)
      fit%gradient(j)=-sum(fit%scores(:,j))
    end do
    allocate(opg(k,k))
    opg=matmul(transpose(fit%scores),fit%scores)
    allocate(fit%covariance_opg(k,k))
    call invert_matrix(opg,invo,oko)
    if(oko)then
    fit%covariance_opg=invo
    else
    fit%covariance_opg=0.0_dp
    end if
    if(.not.compute_hessian)then
      allocate(fit%hessian(0,0),fit%covariance_hessian(0,0),fit%covariance_qmle(0,0))
      return
    end if
    allocate(fit%hessian(k,k))
    fit%hessian=0.0_dp
    f0=-log_likelihood(y,spec%distribution,pars)
    do i=1,k
      h1=step_size(p(i),lo(i),hi(i))
      pp=p
      pm=p
      pp(i)=min(hi(i)-epsilon(1.0_dp),p(i)+h1)
      pm(i)=max(lo(i)+epsilon(1.0_dp),p(i)-h1)
      h1=0.5_dp*(pp(i)-pm(i))
      fpp=objective_natural(y,spec,pp)
      fmm=objective_natural(y,spec,pm)
      fit%hessian(i,i)=(fpp-2.0_dp*f0+fmm)/(h1*h1)
      do j=i+1,k
        h2=step_size(p(j),lo(j),hi(j))
        pp=p
        pp(i)=min(hi(i)-epsilon(1.0_dp),p(i)+h1)
        pp(j)=min(hi(j)-epsilon(1.0_dp),p(j)+h2)
        pm=p
        pm(i)=min(hi(i)-epsilon(1.0_dp),p(i)+h1)
        pm(j)=max(lo(j)+epsilon(1.0_dp),p(j)-h2)
        fpm=objective_natural(y,spec,pm)
        fpp=objective_natural(y,spec,pp)
        pp=p
        pp(i)=max(lo(i)+epsilon(1.0_dp),p(i)-h1)
        pp(j)=min(hi(j)-epsilon(1.0_dp),p(j)+h2)
        fmp=objective_natural(y,spec,pp)
        pm=p
        pm(i)=max(lo(i)+epsilon(1.0_dp),p(i)-h1)
        pm(j)=max(lo(j)+epsilon(1.0_dp),p(j)-h2)
        fmm=objective_natural(y,spec,pm)
        fit%hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*h1*h2)
        fit%hessian(j,i)=fit%hessian(i,j)
      end do
    end do
    allocate(fit%covariance_hessian(k,k),fit%covariance_qmle(k,k))
    call invert_matrix(fit%hessian,invh,okh)
    if(okh)then
      fit%covariance_hessian=invh
      fit%covariance_qmle=matmul(invh,matmul(opg,invh))
    else
      fit%covariance_hessian=0.0_dp
      fit%covariance_qmle=0.0_dp
      if(fit%status==tsd_success)then
      fit%status=tsd_singular
      fit%message='fit converged but numerical Hessian is singular'
      end if
    end if
  end subroutine numerical_inference

  subroutine information_covariance(fit,method,covariance,status)
    type(distribution_fit),intent(in)::fit
    character(len=*),intent(in)::method
    real(dp),allocatable,intent(out)::covariance(:,:)
    integer,intent(out),optional::status
    select case(trim(adjustl(method)))
    case('H','h','hessian')
      covariance=fit%covariance_hessian
    case('OP','op','opg')
      covariance=fit%covariance_opg
    case('QMLE','qmle','sandwich')
      covariance=fit%covariance_qmle
    case default
      allocate(covariance(0,0))
      if(present(status))status=tsd_invalid_argument
      return
    end select
    if(present(status))status=tsd_success
  end subroutine information_covariance

  subroutine observation_logdensities(y,spec,p,values)
    real(dp),intent(in)::y(:),p(:)
    type(parameter_specification),intent(in)::spec
    real(dp),intent(out)::values(:)
    type(distribution_parameters)::pars
    integer::i
    call unpack_parameters(spec,p,pars)
    do i=1,size(y)
    values(i)=ddist(spec%distribution,y(i),pars,.true.)
    end do
  end subroutine observation_logdensities

  real(dp) function objective_natural(y,spec,p) result(value)
    real(dp),intent(in)::y(:),p(:)
    type(parameter_specification),intent(in)::spec
    type(distribution_parameters)::pars
    call unpack_parameters(spec,p,pars)
    value=-log_likelihood(y,spec%distribution,pars)
  end function objective_natural

  pure real(dp) function step_size(x,lo,hi) result(h)
    real(dp),intent(in)::x,lo,hi
    h=max(1.0e-5_dp*max(1.0_dp,abs(x)),1.0e-7_dp)
    h=min(h,0.2_dp*max(hi-lo,1.0e-6_dp))
  end function step_size

  subroutine pack_parameters(spec,pars,p,lo,hi)
    type(parameter_specification),intent(in)::spec
    type(distribution_parameters),intent(in)::pars
    real(dp),allocatable,intent(out)::p(:),lo(:),hi(:)
    real(dp)::allp(5),alllo(5),allhi(5)
    logical::mask(5)
    allp=[pars%mu,pars%sigma,pars%skew,pars%shape,pars%lambda]
    alllo=[spec%lower%mu,spec%lower%sigma,spec%lower%skew,spec%lower%shape,spec%lower%lambda]
    allhi=[spec%upper%mu,spec%upper%sigma,spec%upper%skew,spec%upper%shape,spec%upper%lambda]
    mask=[spec%estimate_mu,spec%estimate_sigma,spec%estimate_skew,spec%estimate_shape,spec%estimate_lambda]
    p=pack(allp,mask)
    lo=pack(alllo,mask)
    hi=pack(allhi,mask)
  end subroutine pack_parameters

  subroutine unpack_parameters(spec,p,pars)
    type(parameter_specification),intent(in)::spec
    real(dp),intent(in)::p(:)
    type(distribution_parameters),intent(out)::pars
    integer::j
    pars=spec%parameters
    j=0
    if(spec%estimate_mu)then
    j=j+1
    pars%mu=p(j)
    end if
    if(spec%estimate_sigma)then
    j=j+1
    pars%sigma=p(j)
    end if
    if(spec%estimate_skew)then
    j=j+1
    pars%skew=p(j)
    end if
    if(spec%estimate_shape)then
    j=j+1
    pars%shape=p(j)
    end if
    if(spec%estimate_lambda)then
    j=j+1
    pars%lambda=p(j)
    end if
  end subroutine unpack_parameters

end module tsd_fit
