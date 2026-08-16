! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types
  use fitdistrplus_math, only : sample_mean, sample_variance, type7_quantile, &
    weighted_quantile, numerical_hessian, covariance_from_hessian
  use fitdistrplus_optimize, only : nelder_mead, parameters_to_unconstrained, &
    unconstrained_to_parameters, transformation_jacobian
  implicit none
  private

  integer, parameter :: objective_mle = 1
  integer, parameter :: objective_censored_mle = 2
  integer, parameter :: objective_mme = 3
  integer, parameter :: objective_qme = 4
  integer, parameter :: objective_mge = 5
  integer, parameter :: objective_mse = 6

  type :: fit_objective_context
    integer :: objective_kind = 0
    integer :: metric = 0
    integer :: divergence = phi_kl
    real(dp) :: power = 2.0_dp
    logical :: total_likelihood = .false.
    type(distribution_model) :: dist
    real(dp), allocatable :: data(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp), allocatable :: sorted(:)
    real(dp), allocatable :: sorted_weights(:)
    real(dp), allocatable :: empirical(:)
    real(dp), allocatable :: probabilities(:)
    integer, allocatable :: orders(:)
    real(dp), allocatable :: censored_left(:)
    real(dp), allocatable :: censored_right(:)
  end type fit_objective_context

  public :: fitdist, fitdist_auto, fitdistcens, prefit
  public :: mledist, mledist_censored, mmedist, qmedist, mgedist, msedist
  public :: default_start, log_likelihood

contains

  subroutine fitdist(data, dist, method, start, result, control, lower, upper, &
      weights, orders, probs, gof, phidiv, power)
    real(dp), intent(in) :: data(:), start(:)
    type(distribution_model), intent(in) :: dist
    integer, intent(in) :: method
    type(fit_result), intent(out) :: result
    type(fit_control), intent(in), optional :: control
    real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
    integer, intent(in), optional :: orders(:), gof, phidiv
    real(dp), intent(in), optional :: probs(:), power

    select case(method)
    case(method_mle)
      call mledist(data,dist,start,result,control,lower,upper,weights)
    case(method_mme)
      if (.not.present(orders)) then
        call invalid_result(result,"MME requires moment orders")
      else
        call mmedist(data,dist,start,orders,result,control,lower,upper,weights)
      end if
    case(method_qme)
      if (.not.present(probs)) then
        call invalid_result(result,"QME requires probabilities")
      else
        call qmedist(data,dist,start,probs,result,control,lower,upper,weights)
      end if
    case(method_mge)
      if (present(gof)) then
        call mgedist(data,dist,start,result,gof,control,lower,upper)
      else
        call mgedist(data,dist,start,result,gof_cvm,control,lower,upper)
      end if
    case(method_mse)
      if (present(phidiv) .and. present(power)) then
        call msedist(data,dist,start,result,phidiv,power,control,lower,upper,weights)
      else if (present(phidiv)) then
        call msedist(data,dist,start,result,phidiv,control=control,lower=lower, &
          upper=upper,weights=weights)
      else
        call msedist(data,dist,start,result,phi_kl,control=control,lower=lower, &
          upper=upper,weights=weights)
      end if
    case default
      call invalid_result(result,"unknown fitting method")
    end select
  end subroutine fitdist

  subroutine fitdist_auto(data,dist,method,result,control,lower,upper,weights,orders,probs,gof,phidiv,power)
    real(dp),intent(in)::data(:)
    type(distribution_model),intent(in)::dist
    integer,intent(in)::method
    type(fit_result),intent(out)::result
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:),probs(:),power
    integer,intent(in),optional::orders(:),gof,phidiv
    real(dp),allocatable::start(:)
    integer::status
    call default_start(data,dist,start,status)
    if(status/=fit_success)then
      call invalid_result(result,"automatic starting values are unavailable");return
    end if
    call fitdist(data,dist,method,start,result,control,lower,upper,weights,orders,probs,gof,phidiv,power)
  end subroutine fitdist_auto

  subroutine prefit(data,dist,method,feasible_parameters,result,control,lower,upper,weights,orders,probs,gof)
    real(dp),intent(in)::data(:),feasible_parameters(:)
    type(distribution_model),intent(in)::dist
    integer,intent(in)::method
    type(fit_result),intent(out)::result
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:),probs(:)
    integer,intent(in),optional::orders(:),gof
    call fitdist(data,dist,method,feasible_parameters,result,control,lower,upper,weights,orders,probs,gof)
  end subroutine prefit

  subroutine fitdistcens(sample,dist,start,result,control,lower,upper,weights)
    type(censored_sample),intent(in)::sample
    type(distribution_model),intent(in)::dist
    real(dp),intent(in)::start(:)
    type(fit_result),intent(out)::result
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:)
    call mledist_censored(sample,dist,start,result,control,lower,upper,weights)
  end subroutine fitdistcens

  subroutine mledist(data,dist,start,result,control,lower,upper,weights)
    real(dp), intent(in) :: data(:), start(:)
    type(distribution_model), intent(in) :: dist
    type(fit_result), intent(out) :: result
    type(fit_control), intent(in), optional :: control
    real(dp), intent(in), optional :: lower(:),upper(:),weights(:)
    real(dp), allocatable :: lo(:),hi(:),z(:),par(:),w(:),hess(:,:),covz(:,:),jac(:,:)
    real(dp) :: fval, loglik
    type(fit_control) :: ctl
    type(fit_objective_context) :: context
    integer :: status,it,stcov

    call initialize_fit(result,dist,"mle",size(data))
    if (.not.associated(dist%logpdf) .or. size(data)<2 .or. size(start)/=dist%npar) then
      call invalid_result(result,"invalid data, start, or density callback"); return
    end if
    if (any(.not.ieee_is_finite(data))) then
      call invalid_result(result,"data must be finite"); return
    end if
    call prepare_control(control,ctl)
    call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success) then
      call invalid_result(result,"invalid parameter bounds"); return
    end if
    call prepare_weights(data,weights,w,status)
    if(status/=fit_success) then
      call invalid_result(result,"invalid weights"); return
    end if

    if (trim(dist%name)=="uniform") then
      allocate(par(2)); par=[minval(data),maxval(data)]
      if(par(2)<=par(1)) then
        call invalid_result(result,"uniform sample has zero range"); return
      end if
      loglik=log_likelihood(data,dist,par,w)
      call finish_fit(result,par,loglik,-loglik/real(size(data),dp),fit_success,0)
      result%message="closed-form uniform MLE"
      allocate(result%covariance(0,0),result%standard_error(0))
      return
    end if

    call parameters_to_unconstrained(start,lo,hi,z,status)
    if(status/=fit_success) then
      call invalid_result(result,"start is incompatible with bounds"); return
    end if
    call set_base_context(context,objective_mle,dist,data,w,lo,hi)
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=log_likelihood(data,dist,par,w)
    call finish_fit(result,par,loglik,fval,status,it)
    if(ctl%calculate_vcov .and. status==fit_success) then
      context%total_likelihood=.true.
      call numerical_hessian(fit_objective,context,z,hess)
      call covariance_from_hessian(hess,covz,stcov)
      if(stcov==fit_success) then
        call transformation_jacobian(z,lo,hi,jac)
        result%covariance=matmul(jac,matmul(covz,transpose(jac)))
        allocate(result%standard_error(dist%npar))
        result%standard_error=sqrt(max(diagonal(result%covariance),0.0_dp))
      else
        allocate(result%covariance(0,0),result%standard_error(0))
      end if
    else
      allocate(result%covariance(0,0),result%standard_error(0))
    end if
  end subroutine mledist

  subroutine mledist_censored(sample,dist,start,result,control,lower,upper,weights)
    type(censored_sample), intent(in) :: sample
    type(distribution_model), intent(in) :: dist
    real(dp), intent(in) :: start(:)
    type(fit_result), intent(out) :: result
    type(fit_control), intent(in), optional :: control
    real(dp), intent(in), optional :: lower(:),upper(:),weights(:)
    real(dp), allocatable :: lo(:),hi(:),z(:),par(:),w(:),hess(:,:),covz(:,:),jac(:,:)
    real(dp) :: fval,loglik
    type(fit_control) :: ctl
    type(fit_objective_context) :: context
    integer :: status,it,stcov,n

    n=size(sample%left)
    call initialize_fit(result,dist,"mle-censored",n)
    if(n<2 .or. size(sample%right)/=n .or. size(start)/=dist%npar .or. &
       .not.associated(dist%logpdf) .or. .not.associated(dist%cdf)) then
      call invalid_result(result,"invalid censored sample or distribution callbacks"); return
    end if
    if(any(sample%left>sample%right)) then
      call invalid_result(result,"left censoring endpoint exceeds right endpoint"); return
    end if
    call prepare_control(control,ctl)
    call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success) then
      call invalid_result(result,"invalid parameter bounds"); return
    end if
    call prepare_weights(sample%left,weights,w,status,allow_infinite=.true.)
    if(status/=fit_success) then
      call invalid_result(result,"invalid weights"); return
    end if
    call parameters_to_unconstrained(start,lo,hi,z,status)
    if(status/=fit_success) then
      call invalid_result(result,"start is incompatible with bounds"); return
    end if
    call set_base_context(context,objective_censored_mle,dist,sample%left,w,lo,hi)
    context%censored_left=sample%left;context%censored_right=sample%right
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=censored_log_likelihood_arrays(sample%left,sample%right,dist,par,w)
    call finish_fit(result,par,loglik,fval,status,it)
    if(ctl%calculate_vcov .and. status==fit_success) then
      context%total_likelihood=.true.
      call numerical_hessian(fit_objective,context,z,hess)
      call covariance_from_hessian(hess,covz,stcov)
      if(stcov==fit_success) then
        call transformation_jacobian(z,lo,hi,jac)
        result%covariance=matmul(jac,matmul(covz,transpose(jac)))
        allocate(result%standard_error(dist%npar))
        result%standard_error=sqrt(max(diagonal(result%covariance),0.0_dp))
      else
        allocate(result%covariance(0,0),result%standard_error(0))
      end if
    else
      allocate(result%covariance(0,0),result%standard_error(0))
    end if
  end subroutine mledist_censored

  subroutine mmedist(data,dist,start,orders,result,control,lower,upper,weights)
    real(dp),intent(in)::data(:),start(:)
    type(distribution_model),intent(in)::dist
    integer,intent(in)::orders(:)
    type(fit_result),intent(out)::result
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:)
    real(dp),allocatable::lo(:),hi(:),z(:),par(:),w(:),emp(:)
    real(dp)::fval,loglik
    type(fit_control)::ctl
    type(fit_objective_context)::context
    integer::status,it,j

    call initialize_fit(result,dist,"mme",size(data))
    if(.not.associated(dist%raw_moment) .or. size(start)/=dist%npar .or. &
       size(orders)/=dist%npar .or. any(orders<1)) then
      call invalid_result(result,"MME requires one valid moment order per free parameter");return
    end if
    call prepare_control(control,ctl);call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid bounds");return;end if
    call prepare_weights(data,weights,w,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid weights");return;end if
    allocate(emp(size(orders)))
    do j=1,size(orders);emp(j)=sum(w*data**orders(j))/sum(w);end do
    call parameters_to_unconstrained(start,lo,hi,z,status)
    call set_base_context(context,objective_mme,dist,data,w,lo,hi)
    context%orders=orders;context%empirical=emp
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=log_likelihood(data,dist,par,w)
    call finish_fit(result,par,loglik,fval,status,it)
    allocate(result%covariance(0,0),result%standard_error(0))
  end subroutine mmedist

  subroutine qmedist(data,dist,start,probs,result,control,lower,upper,weights)
    real(dp),intent(in)::data(:),start(:),probs(:)
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(out)::result
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:)
    real(dp),allocatable::lo(:),hi(:),z(:),par(:),w(:),emp(:)
    real(dp)::fval,loglik
    type(fit_control)::ctl
    type(fit_objective_context)::context
    integer::status,it,j

    call initialize_fit(result,dist,"qme",size(data))
    if(.not.associated(dist%quantile) .or. size(start)/=dist%npar .or. &
       size(probs)/=dist%npar .or. any(probs<=0.0_dp) .or. any(probs>=1.0_dp)) then
      call invalid_result(result,"QME requires one interior probability per free parameter");return
    end if
    call prepare_control(control,ctl);call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid bounds");return;end if
    call prepare_weights(data,weights,w,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid weights");return;end if
    allocate(emp(size(probs)))
    do j=1,size(probs)
      if(present(weights))then
        emp(j)=weighted_quantile(data,w,probs(j))
      else
        emp(j)=type7_quantile(data,probs(j))
      end if
    end do
    call parameters_to_unconstrained(start,lo,hi,z,status)
    call set_base_context(context,objective_qme,dist,data,w,lo,hi)
    context%probabilities=probs;context%empirical=emp
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=log_likelihood(data,dist,par,w)
    call finish_fit(result,par,loglik,fval,status,it)
    allocate(result%covariance(0,0),result%standard_error(0))
  end subroutine qmedist

  subroutine mgedist(data,dist,start,result,gof,control,lower,upper)
    real(dp),intent(in)::data(:),start(:)
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(out)::result
    integer,intent(in)::gof
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:)
    real(dp),allocatable::lo(:),hi(:),z(:),par(:),sorted(:),unitw(:)
    real(dp)::fval,loglik
    type(fit_control)::ctl
    type(fit_objective_context)::context
    integer::status,it

    call initialize_fit(result,dist,"mge",size(data))
    if(.not.associated(dist%cdf) .or. size(start)/=dist%npar .or. &
       gof<gof_cvm .or. gof>gof_ad2)then
      call invalid_result(result,"invalid MGE arguments");return
    end if
    call prepare_control(control,ctl);call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid bounds");return;end if
    sorted=data;call sort_local(sorted);allocate(unitw(size(data)));unitw=1.0_dp
    call parameters_to_unconstrained(start,lo,hi,z,status)
    call set_base_context(context,objective_mge,dist,data,unitw,lo,hi)
    context%sorted=sorted;context%metric=gof
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=log_likelihood(data,dist,par,unitw)
    call finish_fit(result,par,loglik,fval,status,it)
    allocate(result%covariance(0,0),result%standard_error(0))
  end subroutine mgedist

  subroutine msedist(data,dist,start,result,phidiv,power,control,lower,upper,weights)
    real(dp),intent(in)::data(:),start(:)
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(out)::result
    integer,intent(in)::phidiv
    real(dp),intent(in),optional::power
    type(fit_control),intent(in),optional::control
    real(dp),intent(in),optional::lower(:),upper(:),weights(:)
    real(dp),allocatable::lo(:),hi(:),z(:),par(:),sorted(:),w(:),sortedw(:)
    real(dp)::fval,loglik,pow
    type(fit_control)::ctl
    type(fit_objective_context)::context
    integer::status,it,i,j

    call initialize_fit(result,dist,"mse",size(data))
    if(.not.associated(dist%cdf) .or. size(start)/=dist%npar .or. &
       phidiv<phi_kl .or. phidiv>phi_v)then
      call invalid_result(result,"invalid MSE arguments");return
    end if
    pow=2.0_dp;if(present(power))pow=power
    call prepare_control(control,ctl);call prepare_bounds(dist,lower,upper,lo,hi,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid bounds");return;end if
    call prepare_weights(data,weights,w,status)
    if(status/=fit_success)then;call invalid_result(result,"invalid weights");return;end if
    sorted=data;sortedw=w
    do i=2,size(sorted)
      j=i
      do while(j>1 .and. sorted(j)<sorted(j-1))
        call swap_real(sorted(j),sorted(j-1));call swap_real(sortedw(j),sortedw(j-1));j=j-1
      end do
    end do
    call parameters_to_unconstrained(start,lo,hi,z,status)
    call set_base_context(context,objective_mse,dist,data,w,lo,hi)
    context%sorted=sorted;context%sorted_weights=sortedw
    context%divergence=phidiv;context%power=pow
    call nelder_mead(fit_objective,context,z,fval,status,it,ctl%max_iterations, &
      ctl%tolerance,ctl%simplex_scale)
    allocate(par(dist%npar));call unconstrained_to_parameters(z,lo,hi,par)
    loglik=log_likelihood(data,dist,par,w)
    call finish_fit(result,par,loglik,fval,status,it)
    allocate(result%covariance(0,0),result%standard_error(0))
  end subroutine msedist

  function fit_objective(z,context_any) result(value)
    real(dp),intent(in)::z(:)
    class(*),intent(inout)::context_any
    real(dp)::value
    real(dp)::par(size(z)),scale
    integer::i,n,k,npos
    real(dp),allocatable::p(:),spacings(:),ww(:)
    logical,allocatable::positive(:)

    value=huge(1.0_dp)/100.0_dp
    select type(context=>context_any)
    type is(fit_objective_context)
      call unconstrained_to_parameters(z,context%lower,context%upper,par)
      select case(context%objective_kind)
      case(objective_mle)
        value=-log_likelihood(context%data,context%dist,par,context%weights)
        if(.not.context%total_likelihood)value=value/sum(context%weights)
      case(objective_censored_mle)
        value=-censored_log_likelihood_arrays(context%censored_left,context%censored_right, &
          context%dist,par,context%weights)
        if(.not.context%total_likelihood)value=value/sum(context%weights)
      case(objective_mme)
        value=0.0_dp
        do k=1,size(context%orders)
          scale=context%dist%raw_moment(context%orders(k),par)
          if(.not.ieee_is_finite(scale))then;value=huge(1.0_dp)/100.0_dp;return;end if
          value=value+(context%empirical(k)-scale)**2
        end do
        value=value/real(size(context%orders),dp)
      case(objective_qme)
        value=0.0_dp
        do k=1,size(context%probabilities)
          scale=context%dist%quantile(context%probabilities(k),par)
          if(.not.ieee_is_finite(scale))then;value=huge(1.0_dp)/100.0_dp;return;end if
          value=value+(context%empirical(k)-scale)**2
        end do
        value=value/real(size(context%probabilities),dp)
      case(objective_mge)
        n=size(context%sorted);allocate(p(n))
        do i=1,n;p(i)=context%dist%cdf(context%sorted(i),par);end do
        if(any(p<=0.0_dp) .or. any(p>=1.0_dp))then;value=huge(1.0_dp)/100.0_dp;return;end if
        select case(context%metric)
        case(gof_cvm)
          value=1.0_dp/(12.0_dp*real(n*n,dp))+ &
            sum((p-[(real(2*i-1,dp)/(2.0_dp*n),i=1,n)])**2)/real(n,dp)
        case(gof_ks)
          value=maxval(max(abs(p-[(real(i,dp)/n,i=1,n)]), &
            abs(p-[(real(i-1,dp)/n,i=1,n)])))
        case(gof_ad)
          value=-1.0_dp-sum([(real(2*i-1,dp)*(log(p(i))+ &
            log(1.0_dp-p(n+1-i))),i=1,n)])/real(n*n,dp)
        case(gof_adr)
          value=0.5_dp-2.0_dp*sum(p)/real(n,dp)- &
            sum([(real(2*i-1,dp)*log(1.0_dp-p(n+1-i)),i=1,n)])/real(n*n,dp)
        case(gof_adl)
          value=-1.5_dp+2.0_dp*sum(p)/real(n,dp)- &
            sum([(real(2*i-1,dp)*log(p(i)),i=1,n)])/real(n*n,dp)
        case(gof_ad2r)
          value=2.0_dp*sum(log(1.0_dp-p))/real(n,dp)+ &
            sum([(real(2*i-1,dp)/(1.0_dp-p(n+1-i)),i=1,n)])/real(n*n,dp)
        case(gof_ad2l)
          value=2.0_dp*sum(log(p))/real(n,dp)+ &
            sum([(real(2*i-1,dp)/p(i),i=1,n)])/real(n*n,dp)
        case(gof_ad2)
          value=2.0_dp*sum(log(p*(1.0_dp-p)))/real(n,dp)+ &
            sum([(real(2*i-1,dp)*(1.0_dp/p(i)+ &
            1.0_dp/(1.0_dp-p(n+1-i))),i=1,n)])/real(n*n,dp)
        case default
          value=huge(1.0_dp)/100.0_dp
        end select
      case(objective_mse)
        n=size(context%sorted);allocate(spacings(n+1),ww(n+1),positive(n+1))
        scale=0.0_dp
        do k=1,n
          value=context%dist%cdf(context%sorted(k),par)
          spacings(k)=value-scale;scale=value
        end do
        spacings(n+1)=1.0_dp-scale;ww(1)=context%sorted_weights(1)
        ww(2:)=context%sorted_weights;positive=spacings>0.0_dp;npos=count(positive)
        if(npos==0)then;value=huge(1.0_dp)/100.0_dp;return;end if
        select case(context%divergence)
        case(phi_kl)
          value=-sum(pack(ww*log(max(spacings,tiny(1.0_dp))),positive))/real(npos,dp)
        case(phi_j)
          value=-sum(pack(ww*log(max(spacings,tiny(1.0_dp)))*(1.0_dp-spacings),positive))/real(npos,dp)
        case(phi_r)
          value=-sign(1.0_dp,1.0_dp-context%power)* &
            sum(pack(ww*spacings**context%power,positive))/real(npos,dp)
        case(phi_h)
          value=sum(pack(ww*abs(1.0_dp-spacings**(1.0_dp/context%power))** &
            context%power,positive))/real(npos,dp)
        case(phi_v)
          value=sum(pack(ww*abs(1.0_dp-spacings)**context%power,positive))/real(npos,dp)
        end select
      case default
        value=huge(1.0_dp)/100.0_dp
      end select
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
    class default
      value=huge(1.0_dp)/100.0_dp
    end select
  end function fit_objective

  subroutine set_base_context(context,kind,dist,data,weights,lower,upper)
    type(fit_objective_context),intent(out)::context
    integer,intent(in)::kind
    type(distribution_model),intent(in)::dist
    real(dp),intent(in)::data(:),weights(:),lower(:),upper(:)
    context%objective_kind=kind;context%dist=dist;context%data=data
    context%weights=weights;context%lower=lower;context%upper=upper
  end subroutine set_base_context

  subroutine default_start(data,dist,start,status)
    real(dp),intent(in)::data(:)
    type(distribution_model),intent(in)::dist
    real(dp),allocatable,intent(out)::start(:)
    integer,intent(out)::status
    real(dp)::m,v,s,mlog,vlog,q1,q3,common
    if(size(data)<2 .or. any(.not.ieee_is_finite(data)))then
      allocate(start(0));status=fit_invalid_argument;return
    end if
    m=sample_mean(data);v=sample_variance(data,.false.);s=sqrt(max(v,tiny(1.0_dp)))
    q1=type7_quantile(data,0.25_dp);q3=type7_quantile(data,0.75_dp)
    select case(trim(dist%name))
    case("normal");start=[m,max(s,1.0e-6_dp)]
    case("lognormal")
      if(any(data<=0.0_dp))then;allocate(start(0));status=fit_invalid_argument;return;end if
      mlog=sum(log(data))/real(size(data),dp);vlog=sum((log(data)-mlog)**2)/real(size(data),dp)
      start=[mlog,sqrt(max(vlog,1.0e-12_dp))]
    case("exponential")
      if(m<=0.0_dp)then;allocate(start(0));status=fit_invalid_argument;return;end if
      start=[1.0_dp/m]
    case("gamma")
      if(m<=0.0_dp .or. v<=0.0_dp)then;allocate(start(0));status=fit_invalid_argument;return;end if
      start=[m*m/v,m/v]
    case("weibull")
      if(any(data<=0.0_dp))then;allocate(start(0));status=fit_invalid_argument;return;end if
      mlog=sum(log(data))/real(size(data),dp);vlog=sum((log(data)-mlog)**2)/real(size(data),dp)
      common=max(sqrt(vlog),0.05_dp);start=[1.2_dp/common,exp(mlog+0.572_dp*common/1.2_dp)]
    case("uniform");start=[minval(data),maxval(data)]
    case("logistic");start=[type7_quantile(data,0.5_dp),max(s*sqrt(3.0_dp)/acos(-1.0_dp),1.0e-6_dp)]
    case("cauchy");start=[type7_quantile(data,0.5_dp),max(0.5_dp*(q3-q1),1.0e-6_dp)]
    case("beta")
      if(any(data<=0.0_dp) .or. any(data>=1.0_dp) .or. v<=0.0_dp)then
        allocate(start(0));status=fit_invalid_argument;return
      end if
      common=max(m*(1.0_dp-m)/v-1.0_dp,0.1_dp)
      start=[max(m*common,0.1_dp),max((1.0_dp-m)*common,0.1_dp)]
    case("poisson");start=[max(m,1.0e-6_dp)]
    case("geometric");start=[min(max(1.0_dp/(1.0_dp+m),1.0e-6_dp),1.0_dp-1.0e-6_dp)]
    case("negative-binomial")
      common=merge(m*m/(v-m),100.0_dp,v>m+1.0e-10_dp)
      start=[max(common,1.0e-3_dp),max(m,1.0e-6_dp)]
    case default
      allocate(start(0));status=fit_not_supported;return
    end select
    status=fit_success
  end subroutine default_start

  function log_likelihood(data,dist,par,weights) result(value)
    real(dp),intent(in)::data(:),par(:),weights(:)
    type(distribution_model),intent(in)::dist
    real(dp)::value,term
    integer::i
    value=0.0_dp
    if(size(weights)/=size(data) .or. .not.associated(dist%logpdf))then
      value=-huge(1.0_dp);return
    end if
    do i=1,size(data)
      term=dist%logpdf(data(i),par)
      if(.not.ieee_is_finite(term) .or. term<=-0.5_dp*huge(1.0_dp))then
        value=-huge(1.0_dp);return
      end if
      value=value+weights(i)*term
    end do
  end function log_likelihood

  function censored_log_likelihood_arrays(left,right,dist,par,weights) result(value)
    real(dp),intent(in)::left(:),right(:),par(:),weights(:)
    type(distribution_model),intent(in)::dist
    real(dp)::value,term,pl,pr,tol
    integer::i
    value=0.0_dp
    do i=1,size(left)
      tol=epsilon(1.0_dp)*max(1.0_dp,abs(left(i)),abs(right(i)))
      if(ieee_is_finite(left(i)) .and. ieee_is_finite(right(i)) .and. &
         abs(left(i)-right(i))<=tol)then
        term=dist%logpdf(left(i),par)
      else if(.not.ieee_is_finite(left(i)))then
        pr=dist%cdf(right(i),par);term=log(max(pr,tiny(1.0_dp)))
      else if(.not.ieee_is_finite(right(i)))then
        pl=dist%cdf(left(i),par);term=log(max(1.0_dp-pl,tiny(1.0_dp)))
      else
        pl=dist%cdf(left(i),par);pr=dist%cdf(right(i),par)
        term=log(max(pr-pl,tiny(1.0_dp)))
      end if
      if(.not.ieee_is_finite(term))then;value=-huge(1.0_dp);return;end if
      value=value+weights(i)*term
    end do
  end function censored_log_likelihood_arrays

  subroutine prepare_control(control,ctl)
    type(fit_control),intent(in),optional::control
    type(fit_control),intent(out)::ctl
    ctl=fit_control();if(present(control))ctl=control
  end subroutine prepare_control

  subroutine prepare_bounds(dist,lower,upper,lo,hi,status)
    type(distribution_model),intent(in)::dist
    real(dp),intent(in),optional::lower(:),upper(:)
    real(dp),allocatable,intent(out)::lo(:),hi(:)
    integer,intent(out)::status
    if(allocated(dist%default_lower))then
      lo=dist%default_lower;hi=dist%default_upper
    else
      allocate(lo(dist%npar),hi(dist%npar));lo=-huge(1.0_dp);hi=huge(1.0_dp)
    end if
    if(present(lower))then
      if(size(lower)/=dist%npar)then;status=fit_invalid_argument;return;end if
      lo=lower
    end if
    if(present(upper))then
      if(size(upper)/=dist%npar)then;status=fit_invalid_argument;return;end if
      hi=upper
    end if
    if(any(lo>=hi))then;status=fit_invalid_argument;else;status=fit_success;end if
  end subroutine prepare_bounds

  subroutine prepare_weights(data,weights,w,status,allow_infinite)
    real(dp),intent(in)::data(:)
    real(dp),intent(in),optional::weights(:)
    real(dp),allocatable,intent(out)::w(:)
    integer,intent(out)::status
    logical,intent(in),optional::allow_infinite
    logical::allowinf
    allowinf=.false.;if(present(allow_infinite))allowinf=allow_infinite
    allocate(w(size(data)))
    if(present(weights))then
      if(size(weights)/=size(data) .or. any(weights<0.0_dp) .or. sum(weights)<=0.0_dp)then
        w=0.0_dp;status=fit_invalid_argument;return
      end if
      w=weights
    else
      w=1.0_dp
    end if
    if(.not.allowinf .and. any(.not.ieee_is_finite(data)))then;status=fit_invalid_argument;return;end if
    status=fit_success
  end subroutine prepare_weights

  subroutine initialize_fit(result,dist,method,n)
    type(fit_result),intent(out)::result
    type(distribution_model),intent(in)::dist
    character(len=*),intent(in)::method
    integer,intent(in)::n
    result%distribution=dist%name;result%method=method;result%nobs=n
    result%message="not fitted";result%convergence=fit_invalid_argument
  end subroutine initialize_fit

  subroutine finish_fit(result,par,loglik,objective,status,iterations)
    type(fit_result),intent(inout)::result
    real(dp),intent(in)::par(:),loglik,objective
    integer,intent(in)::status,iterations
    result%estimate=par;result%log_likelihood=loglik;result%objective=objective
    result%iterations=iterations;result%convergence=status
    result%aic=-2.0_dp*loglik+2.0_dp*real(size(par),dp)
    result%bic=-2.0_dp*loglik+log(real(max(result%nobs,1),dp))*real(size(par),dp)
    if(status==fit_success)then;result%message="converged";else;result%message="optimizer did not converge";end if
  end subroutine finish_fit

  subroutine invalid_result(result,message)
    type(fit_result),intent(inout)::result
    character(len=*),intent(in)::message
    result%convergence=fit_invalid_argument;result%message=message
    if(.not.allocated(result%estimate))allocate(result%estimate(0))
    if(.not.allocated(result%covariance))allocate(result%covariance(0,0))
    if(.not.allocated(result%standard_error))allocate(result%standard_error(0))
  end subroutine invalid_result

  pure function diagonal(a) result(value)
    real(dp),intent(in)::a(:,:)
    real(dp)::value(min(size(a,1),size(a,2)))
    integer::i
    do i=1,size(value);value(i)=a(i,i);end do
  end function diagonal

  subroutine sort_local(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x);key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_local

  subroutine swap_real(a,b)
    real(dp),intent(inout)::a,b
    real(dp)::tmp
    tmp=a;a=b;b=tmp
  end subroutine swap_real

end module fitdistrplus_fit
