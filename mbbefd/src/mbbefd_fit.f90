! SPDX-License-Identifier: GPL-2.0-only
module mbbefd_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mbbefd_kinds, only : dp
  use mbbefd_math, only : nan_dp, invert_matrix, type7_quantile
  use mbbefd_distributions
  use fitdistrplus_types, only : distribution_model, fit_result, fit_control, fit_success
  use fitdistrplus_distributions, only : make_beta
  use fitdistrplus_fit, only : mledist
  use fitdistrplus_optimize, only : nelder_mead
  use alabama, only : alabama_result_t, alabama_outer_control_t, alabama_inner_control_t, &
    constr_optim_nl, al_success
  implicit none
  private

  integer, parameter :: dist_oiunif=1, dist_oistpareto=2, dist_oibeta=3, &
    dist_oigbeta=4, dist_mbbefd_ab=5, dist_mbbefd_gb=6

  type, public :: dr_fit_result
    character(len=16) :: distribution = ""
    character(len=8) :: method = ""
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: correlation(:,:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: nobs = 0
    integer :: convergence = 1
    character(len=160) :: message = "not fitted"
  end type dr_fit_result

  type, public :: dr_boot_result
    real(dp), allocatable :: estimates(:,:)
    integer, allocatable :: convergence(:)
    real(dp), allocatable :: confidence_interval(:,:) ! median, 2.5%, 97.5%
    integer :: nsim = 0
    integer :: successful = 0
    integer :: status = 1
  end type dr_boot_result

  type :: tlmme_context_t
    integer :: family = 0
    real(dp) :: p1 = 0.0_dp
    real(dp) :: emp(3) = 0.0_dp
    integer :: nmom = 0
  end type tlmme_context_t

  real(dp), allocatable, save :: constrained_obs(:)
  integer, save :: constrained_family = 0
  logical, save :: constrained_tlmme = .false.
  real(dp), save :: constrained_emp_mean = 0.0_dp
  real(dp), save :: constrained_emp_tl = 0.0_dp

  public :: fit_dr, boot_dr, simulate_dr

contains

  subroutine fit_dr(data, dist, result, method, start)
    real(dp), intent(in) :: data(:)
    character(len=*), intent(in) :: dist
    type(dr_fit_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: start(:)
    character(len=16) :: key
    character(len=8) :: meth

    call initialize_dr_result(result, dist, size(data))
    if (size(data) < 2 .or. any(.not.ieee_is_finite(data)) .or. &
        any(data < 0.0_dp) .or. any(data > 1.0_dp)) then
      result%message = "data must contain at least two finite observations in [0,1]"
      return
    end if
    if (trim(dist) == "MBBEFD") then
      key = "mbbefd_gb"
      result%distribution = "mbbefd_gb"
    else
      key = lower_ascii(trim(dist))
    end if
    meth = "mle"
    if (present(method)) meth = lower_ascii(trim(method))
    if (trim(meth) /= "mle" .and. trim(meth) /= "tlmme") then
      result%message = "method must be 'mle' or 'tlmme'"
      return
    end if
    result%method = trim(meth)

    select case (trim(key))
    case ("oiunif")
      call fit_oiunif(data, result)
    case ("oistpareto")
      call fit_one_inflated(data, dist_oistpareto, meth, result, start)
    case ("oibeta")
      call fit_one_inflated(data, dist_oibeta, meth, result, start)
    case ("oigbeta")
      call fit_one_inflated(data, dist_oigbeta, meth, result, start)
    case ("mbbefd")
      call fit_mbbefd_constrained(data, dist_mbbefd_ab, meth, result, start)
    case ("mbbefd_gb", "mbbefd-gb", "mbbefd2", "gb")
      call fit_mbbefd_constrained(data, dist_mbbefd_gb, meth, result, start)
    case default
      result%message = "unknown destruction-rate distribution"
    end select
  end subroutine fit_dr

  subroutine fit_oiunif(data, result)
    real(dp), intent(in) :: data(:)
    type(dr_fit_result), intent(inout) :: result
    real(dp) :: p1
    integer :: n
    n=size(data); p1=etl(data)
    allocate(result%estimate(1), result%covariance(1,1), result%standard_error(1), &
      result%correlation(1,1))
    result%estimate=[p1]
    result%covariance(1,1)=p1*(1.0_dp-p1)/real(n,dp)
    result%standard_error(1)=sqrt(max(result%covariance(1,1),0.0_dp))
    result%correlation(1,1)=1.0_dp
    result%log_likelihood=loglik_dispatch(data,dist_oiunif,result%estimate)
    call finish_information(result)
    result%objective=-result%log_likelihood/real(n,dp)
    result%convergence=0; result%message="closed-form one-inflated uniform fit"
  end subroutine fit_oiunif

  subroutine fit_one_inflated(data, family, method, result, start)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: family
    character(len=*), intent(in) :: method
    type(dr_fit_result), intent(inout) :: result
    real(dp), intent(in), optional :: start(:)
    real(dp), allocatable :: base(:), start_base(:)
    real(dp) :: p1
    integer :: i,j,nbase

    p1=etl(data); nbase=count(data /= 1.0_dp)
    if (nbase < 2) then
      result%message="at least two observations below 1 are required"
      return
    end if
    allocate(base(nbase)); j=0
    do i=1,size(data)
      if(data(i)/=1.0_dp)then;j=j+1;base(j)=data(i);end if
    end do

    select case(family)
    case(dist_oistpareto); allocate(start_base(1)); start_base=[1.0_dp]
    case(dist_oibeta)
      allocate(start_base(2)); call beta_start(base,start_base)
    case(dist_oigbeta)
      allocate(start_base(3)); call gbeta_start(base,start_base)
    end select
    if(present(start)) then
      if(size(start)>=size(start_base)) start_base=start(1:size(start_base))
    end if

    if(trim(method)=="mle") then
      call fit_base_mle(base,family,start_base,result)
    else
      call fit_base_tlmme(data,family,p1,start_base,result)
    end if
    if(result%convergence/=0) return
    call append_p1(result,p1,data,family)
  end subroutine fit_one_inflated

  subroutine fit_base_mle(base, family, start, result)
    real(dp),intent(in)::base(:),start(:)
    integer,intent(in)::family
    type(dr_fit_result),intent(inout)::result
    type(distribution_model)::model
    type(fit_result)::fr
    type(fit_control)::ctl
    real(dp),allocatable::lower(:),upper(:)

    ctl%max_iterations=5000;ctl%tolerance=1.0e-10_dp;ctl%calculate_vcov=.true.
    select case(family)
    case(dist_oibeta)
      call make_beta(model)
    case(dist_oistpareto)
      call make_stpareto_model(model)
    case(dist_oigbeta)
      call make_gbeta_model(model)
    end select
    allocate(lower(model%npar),upper(model%npar))
    lower=1.0e-8_dp;upper=huge(1.0_dp)
    call mledist(base,model,start,fr,ctl,lower,upper)
    result%convergence=fr%convergence;result%message=fr%message
    if(fr%convergence/=fit_success)return
    result%estimate=fr%estimate
    if(allocated(fr%covariance))result%covariance=fr%covariance
    if(allocated(fr%standard_error))result%standard_error=fr%standard_error
  end subroutine fit_base_mle

  subroutine fit_base_tlmme(data,family,p1,start,result)
    real(dp),intent(in)::data(:),p1,start(:)
    integer,intent(in)::family
    type(dr_fit_result),intent(inout)::result
    type(tlmme_context_t)::ctx
    real(dp),allocatable::z(:)
    real(dp)::fval
    integer::st,it,k
    ctx%family=family;ctx%p1=p1;ctx%nmom=size(start)
    do k=1,min(3,size(start));ctx%emp(k)=sum(data**real(k,dp))/real(size(data),dp);end do
    allocate(z(size(start)));z=log(max(start,1.0e-8_dp))
    call nelder_mead(tlmme_base_objective,ctx,z,fval,st,it,5000,1.0e-10_dp,0.2_dp)
    result%convergence=st;result%objective=fval
    if(st/=fit_success)then;result%message="TLMMe optimization did not converge";return;end if
    result%estimate=exp(z);result%message="TLMMe fit converged"
    allocate(result%covariance(0,0),result%standard_error(0),result%correlation(0,0))
  end subroutine fit_base_tlmme

  function tlmme_base_objective(z,context) result(v)
    real(dp),intent(in)::z(:)
    class(*),intent(inout)::context
    real(dp)::v,par(3),th(3)
    integer::k
    v=huge(1.0_dp);par=0.0_dp;par(1:size(z))=exp(min(z,700.0_dp))
    select type(ctx=>context)
    type is(tlmme_context_t)
      th=0.0_dp
      do k=1,ctx%nmom
        select case(ctx%family)
        case(dist_oistpareto);th(k)=moistpareto(real(k,dp),par(1),ctx%p1)
        case(dist_oibeta);th(k)=moibeta(real(k,dp),par(1),par(2),ctx%p1)
        case(dist_oigbeta);th(k)=moigbeta(real(k,dp),par(1),par(2),par(3),ctx%p1)
        end select
      end do
      v=sum((th(1:ctx%nmom)-ctx%emp(1:ctx%nmom))**2)
      v=v+(ctx%p1-etl_value(ctx%p1))**2 ! zero, kept to mirror TLM objective
    end select
  end function tlmme_base_objective

  pure function etl_value(p) result(v)
    real(dp),intent(in)::p;real(dp)::v;v=p
  end function etl_value

  subroutine append_p1(result,p1,data,family)
    type(dr_fit_result),intent(inout)::result
    real(dp),intent(in)::p1,data(:)
    integer,intent(in)::family
    real(dp),allocatable::par(:),oldcov(:,:),oldse(:)
    integer::k,n
    n=size(data); k=size(result%estimate)
    allocate(par(k+1));par(1:k)=result%estimate;par(k+1)=p1
    call move_alloc(par,result%estimate)
    if(allocated(result%covariance))then
      oldcov=result%covariance;deallocate(result%covariance)
      allocate(result%covariance(k+1,k+1));result%covariance=0.0_dp
      if(size(oldcov,1)==k)result%covariance(1:k,1:k)=oldcov
      result%covariance(k+1,k+1)=p1*(1.0_dp-p1)/real(n,dp)
    else
      allocate(result%covariance(k+1,k+1));result%covariance=0.0_dp
      result%covariance(k+1,k+1)=p1*(1.0_dp-p1)/real(n,dp)
    end if
    if(allocated(result%standard_error))then
      oldse=result%standard_error;deallocate(result%standard_error)
      allocate(result%standard_error(k+1));result%standard_error=0.0_dp
      if(size(oldse)==k)result%standard_error(1:k)=oldse
      result%standard_error(k+1)=sqrt(max(result%covariance(k+1,k+1),0.0_dp))
    else
      allocate(result%standard_error(k+1));result%standard_error=sqrt(max(diagonal(result%covariance),0.0_dp))
    end if
    call make_correlation(result%covariance,result%correlation)
    result%log_likelihood=loglik_dispatch(data,family,result%estimate)
    call finish_information(result)
    result%objective=-result%log_likelihood/real(n,dp)
  end subroutine append_p1

  subroutine fit_mbbefd_constrained(data,family,method,result,start)
    real(dp),intent(in)::data(:)
    integer,intent(in)::family
    character(len=*),intent(in)::method
    type(dr_fit_result),intent(inout)::result
    real(dp),intent(in),optional::start(:)
    type(alabama_result_t)::r1,r2
    type(alabama_outer_control_t)::oc
    type(alabama_inner_control_t)::ic
    real(dp)::p1(2),p2(2),bestll
    real(dp),allocatable::h(:,:),info(:,:),cov(:,:)
    integer::st

    constrained_obs=data;constrained_family=family;constrained_tlmme=(trim(method)=="tlmme")
    constrained_emp_mean=sum(data)/real(size(data),dp);constrained_emp_tl=etl(data)
    oc%itmax=60;oc%eps=1.0e-8_dp;oc%trace=.false.;oc%kkt2_check=.false.
    ic%max_iterations=500;ic%reltol=1.0e-10_dp
    if(family==dist_mbbefd_ab)then
      p1=[-0.5_dp,2.0_dp];p2=[0.5_dp,0.2_dp]
      if(present(start).and.size(start)>=2)then
        if(start(1)<0.0_dp.and.start(1)>-1.0_dp.and.start(2)>1.0_dp)p1=start(1:2)
        if(start(1)>0.0_dp.and.start(2)>0.0_dp.and.start(2)<1.0_dp)p2=start(1:2)
      end if
      if(constrained_tlmme)then
        call constr_optim_nl(p1,constrained_objective,r1,hin=constraint_ab1,hin_jac=jac_ab1, &
          control_outer=oc,control_inner=ic)
        call constr_optim_nl(p2,constrained_objective,r2,hin=constraint_ab2,hin_jac=jac_ab2, &
          control_outer=oc,control_inner=ic)
      else
        call constr_optim_nl(p1,constrained_objective,r1,gr=constrained_gradient, &
          hin=constraint_ab1,hin_jac=jac_ab1,control_outer=oc,control_inner=ic)
        call constr_optim_nl(p2,constrained_objective,r2,gr=constrained_gradient, &
          hin=constraint_ab2,hin_jac=jac_ab2,control_outer=oc,control_inner=ic)
      end if
    else
      p1=[max(2.0_dp,1.0_dp/max(constrained_emp_tl,1.0e-3_dp)),2.0_dp]
      p2=[p1(1),min(0.5_dp/p1(1),0.25_dp)]
      if(present(start).and.size(start)>=2)then
        if(start(1)>1.0_dp.and.start(2)>1.0_dp)p1=start(1:2)
        if(start(1)>1.0_dp.and.start(2)>0.0_dp.and.start(2)<1.0_dp)p2=start(1:2)
      end if
      if(constrained_tlmme)then
        call constr_optim_nl(p1,constrained_objective,r1,hin=constraint_gb1,hin_jac=jac_gb1, &
          control_outer=oc,control_inner=ic)
        call constr_optim_nl(p2,constrained_objective,r2,hin=constraint_gb2,hin_jac=jac_gb2, &
          control_outer=oc,control_inner=ic)
      else
        call constr_optim_nl(p1,constrained_objective,r1,gr=constrained_gradient, &
          hin=constraint_gb1,hin_jac=jac_gb1,control_outer=oc,control_inner=ic)
        call constr_optim_nl(p2,constrained_objective,r2,gr=constrained_gradient, &
          hin=constraint_gb2,hin_jac=jac_gb2,control_outer=oc,control_inner=ic)
      end if
    end if

    if(r1%convergence/=al_success.and.r2%convergence/=al_success)then
      result%message="both constrained regions failed to converge";return
    end if
    if(r2%convergence/=al_success.or.(r1%convergence==al_success.and.r1%value<=r2%value))then
      result%estimate=r1%par;result%objective=r1%value
    else
      result%estimate=r2%par;result%objective=r2%value
    end if
    result%convergence=0;result%message="constrained fit converged"
    result%log_likelihood=loglik_dispatch(data,family,result%estimate)
    call finish_information(result)

    if(trim(method)=="mle")then
      if(family==dist_mbbefd_ab)then
        h=hessloglik_mbbefd(data,result%estimate(1),result%estimate(2));info=-h
      else
        call numerical_info_gb(data,result%estimate,info)
      end if
      call invert_matrix(info,cov,st)
      if(st==0)then
        result%covariance=cov
        allocate(result%standard_error(2));result%standard_error=sqrt(max(diagonal(cov),0.0_dp))
        call make_correlation(cov,result%correlation)
      else
        allocate(result%covariance(0,0),result%standard_error(0),result%correlation(0,0))
      end if
    else
      allocate(result%covariance(0,0),result%standard_error(0),result%correlation(0,0))
    end if
    bestll=result%log_likelihood
    if(.not.ieee_is_finite(bestll))result%message="fit converged but log-likelihood is non-finite"
  end subroutine fit_mbbefd_constrained

  function constrained_objective(par) result(v)
    real(dp),intent(in)::par(:);real(dp)::v,m,t
    if(.not.allocated(constrained_obs))then;v=huge(1.0_dp);return;end if
    if(constrained_tlmme)then
      if(constrained_family==dist_mbbefd_ab)then
        m=mmbbefd(1.0_dp,par(1),par(2));t=tlmbbefd(par(1),par(2))
      else
        m=mmbbefd_gb(1.0_dp,par(1),par(2));t=tlmbbefd_gb(par(1),par(2))
      end if
      if(.not.ieee_is_finite(m).or..not.ieee_is_finite(t))then;v=huge(1.0_dp);return;end if
      v=(m-constrained_emp_mean)**2+(t-constrained_emp_tl)**2
    else
      if(constrained_family==dist_mbbefd_ab)then
        v=-loglik_mbbefd(constrained_obs,par(1),par(2))/real(size(constrained_obs),dp)
      else
        v=-loglik_mbbefd_gb(constrained_obs,par(1),par(2))/real(size(constrained_obs),dp)
      end if
      if(.not.ieee_is_finite(v))v=huge(1.0_dp)
    end if
  end function constrained_objective

  subroutine constrained_gradient(par,g)
    real(dp),intent(in)::par(:);real(dp),intent(out)::g(:);real(dp)::tmp(2)
    if(constrained_tlmme)then;g=0.0_dp;return;end if
    if(constrained_family==dist_mbbefd_ab)then
      tmp=gradloglik_mbbefd(constrained_obs,par(1),par(2))
    else
      tmp=gradloglik_mbbefd_gb(constrained_obs,par(1),par(2))
    end if
    g=-tmp/real(size(constrained_obs),dp)
  end subroutine constrained_gradient

  function constraint_ab1(x) result(v)
    real(dp),intent(in)::x(:);real(dp),allocatable::v(:)
    allocate(v(4));v=[x(1)+1.0_dp,-x(1),x(2)-1.0_dp,x(1)*(1.0_dp-x(2))]
  end function constraint_ab1
  function jac_ab1(x) result(j)
    real(dp),intent(in)::x(:);real(dp),allocatable::j(:,:)
    allocate(j(4,2));j(1,:)=[1.0_dp,0.0_dp];j(2,:)=[-1.0_dp,0.0_dp]; &
      j(3,:)=[0.0_dp,1.0_dp];j(4,:)=[1.0_dp-x(2),-x(1)]
  end function jac_ab1
  function constraint_ab2(x) result(v)
    real(dp),intent(in)::x(:);real(dp),allocatable::v(:)
    allocate(v(4));v=[x(1),x(2),1.0_dp-x(2),x(1)*(1.0_dp-x(2))]
  end function constraint_ab2
  function jac_ab2(x) result(j)
    real(dp),intent(in)::x(:);real(dp),allocatable::j(:,:)
    allocate(j(4,2));j(1,:)=[1.0_dp,0.0_dp];j(2,:)=[0.0_dp,1.0_dp]; &
      j(3,:)=[0.0_dp,-1.0_dp];j(4,:)=[1.0_dp-x(2),-x(1)]
  end function jac_ab2
  function constraint_gb1(x) result(v)
    real(dp),intent(in)::x(:);real(dp),allocatable::v(:)
    allocate(v(3));v=[x(1)-1.0_dp,x(2)-1.0_dp,x(1)*x(2)-1.0_dp]
  end function constraint_gb1
  function jac_gb1(x) result(j)
    real(dp),intent(in)::x(:);real(dp),allocatable::j(:,:)
    allocate(j(3,2));j(1,:)=[1.0_dp,0.0_dp];j(2,:)=[0.0_dp,1.0_dp];j(3,:)=[x(2),x(1)]
  end function jac_gb1
  function constraint_gb2(x) result(v)
    real(dp),intent(in)::x(:);real(dp),allocatable::v(:)
    allocate(v(4));v=[x(1)-1.0_dp,1.0_dp-x(2),x(2),1.0_dp-x(1)*x(2)]
  end function constraint_gb2
  function jac_gb2(x) result(j)
    real(dp),intent(in)::x(:);real(dp),allocatable::j(:,:)
    allocate(j(4,2));j(1,:)=[1.0_dp,0.0_dp];j(2,:)=[0.0_dp,-1.0_dp]; &
      j(3,:)=[0.0_dp,1.0_dp];j(4,:)=[-x(2),-x(1)]
  end function jac_gb2

  subroutine numerical_info_gb(data,par,info)
    real(dp),intent(in)::data(:),par(:)
    real(dp),allocatable,intent(out)::info(:,:)
    real(dp)::h1,h2,f0,fpp,fpm,fmp,fmm,p(2)
    h1=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(par(1)))
    h2=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(par(2)))
    allocate(info(2,2));p=par;f0=loglik_mbbefd_gb(data,p(1),p(2))
    p=par;p(1)=par(1)+h1;fpp=loglik_mbbefd_gb(data,p(1),p(2))
    p=par;p(1)=par(1)-h1;fpm=loglik_mbbefd_gb(data,p(1),p(2))
    info(1,1)=-(fpp-2.0_dp*f0+fpm)/(h1*h1)
    p=par;p(2)=par(2)+h2;fpp=loglik_mbbefd_gb(data,p(1),p(2))
    p=par;p(2)=par(2)-h2;fpm=loglik_mbbefd_gb(data,p(1),p(2))
    info(2,2)=-(fpp-2.0_dp*f0+fpm)/(h2*h2)
    p=par+[h1,h2];fpp=loglik_mbbefd_gb(data,p(1),p(2))
    p=par+[h1,-h2];fpm=loglik_mbbefd_gb(data,p(1),p(2))
    p=par+[-h1,h2];fmp=loglik_mbbefd_gb(data,p(1),p(2))
    p=par+[-h1,-h2];fmm=loglik_mbbefd_gb(data,p(1),p(2))
    info(1,2)=-(fpp-fpm-fmp+fmm)/(4.0_dp*h1*h2);info(2,1)=info(1,2)
  end subroutine numerical_info_gb

  subroutine boot_dr(fit,data,nsim,result,parametric)
    type(dr_fit_result),intent(in)::fit
    real(dp),intent(in)::data(:)
    integer,intent(in)::nsim
    type(dr_boot_result),intent(out)::result
    logical,intent(in),optional::parametric
    logical::param
    real(dp),allocatable::sample(:),vals(:)
    type(dr_fit_result)::fr
    integer::j,k,np,st
    real(dp)::u
    param=.true.;if(present(parametric))param=parametric
    np=size(fit%estimate);result%nsim=nsim
    if(nsim<1.or.np<1.or.size(data)<2.or.fit%convergence/=0)then;result%status=1;return;end if
    allocate(result%estimates(np,nsim),result%convergence(nsim),sample(size(data)))
    result%estimates=nan_dp();result%convergence=1
    do j=1,nsim
      if(param)then
        call simulate_dr(sample,fit%distribution,fit%estimate,st)
        if(st/=0)cycle
      else
        do k=1,size(sample);call random_number(u);sample(k)=data(1+min(int(u*size(data)),size(data)-1));end do
      end if
      call fit_dr(sample,fit%distribution,fr,fit%method,fit%estimate)
      result%convergence(j)=fr%convergence
      if(fr%convergence==0.and.size(fr%estimate)==np)result%estimates(:,j)=fr%estimate
    end do
    result%successful=count(result%convergence==0)
    allocate(result%confidence_interval(np,3));result%confidence_interval=nan_dp()
    do k=1,np
      call successful_values(result%estimates(k,:),result%convergence,vals)
      if(size(vals)>0)then
        result%confidence_interval(k,1)=type7_quantile(vals,0.5_dp)
        result%confidence_interval(k,2)=type7_quantile(vals,0.025_dp)
        result%confidence_interval(k,3)=type7_quantile(vals,0.975_dp)
      end if
    end do
    result%status=merge(0,1,result%successful>0)
  end subroutine boot_dr

  subroutine simulate_dr(x,dist,par,status)
    real(dp),intent(out)::x(:);character(len=*),intent(in)::dist;real(dp),intent(in)::par(:)
    integer,intent(out)::status
    character(len=16)::key
    key=lower_ascii(trim(dist));status=0
    select case(trim(key))
    case("oiunif");if(size(par)>=1)then;call roiunif(x,par(1));else;status=1;end if
    case("oistpareto");if(size(par)>=2)then;call roistpareto(x,par(1),par(2));else;status=1;end if
    case("oibeta");if(size(par)>=3)then;call roibeta(x,par(1),par(2),par(3));else;status=1;end if
    case("oigbeta");if(size(par)>=4)then;call roigbeta(x,par(1),par(2),par(3),par(4));else;status=1;end if
    case("mbbefd");if(size(par)>=2)then;call rmbbefd(x,par(1),par(2));else;status=1;end if
    case("mbbefd_gb","mbbefd-gb","mbbefd2","gb")
      if(size(par)>=2)then;call rmbbefd_gb(x,par(1),par(2));else;status=1;end if
    case default;status=1
    end select
  end subroutine simulate_dr

  subroutine make_stpareto_model(model)
    type(distribution_model),intent(out)::model
    model%name="stpareto";model%npar=1;model%discrete=.false.
    allocate(model%parameter_names(1),model%default_lower(1),model%default_upper(1))
    model%parameter_names=[character(len=24)::"a"]
    model%default_lower=[1.0e-8_dp];model%default_upper=[huge(1.0_dp)]
    model%logpdf=>stpareto_logpdf_cb;model%raw_moment=>stpareto_moment_cb
  end subroutine make_stpareto_model
  function stpareto_logpdf_cb(x,par) result(v)
    real(dp),intent(in)::x,par(:);real(dp)::v;v=dstpareto(x,par(1),.true.)
  end function stpareto_logpdf_cb
  function stpareto_moment_cb(order,par) result(v)
    integer,intent(in)::order;real(dp),intent(in)::par(:);real(dp)::v;v=mstpareto(real(order,dp),par(1))
  end function stpareto_moment_cb

  subroutine make_gbeta_model(model)
    type(distribution_model),intent(out)::model
    model%name="gbeta";model%npar=3;model%discrete=.false.
    allocate(model%parameter_names(3),model%default_lower(3),model%default_upper(3))
    model%parameter_names=[character(len=24)::"shape0","shape1","shape2"]
    model%default_lower=1.0e-8_dp;model%default_upper=huge(1.0_dp)
    model%logpdf=>gbeta_logpdf_cb;model%raw_moment=>gbeta_moment_cb
  end subroutine make_gbeta_model
  function gbeta_logpdf_cb(x,par) result(v)
    real(dp),intent(in)::x,par(:);real(dp)::v;v=dgbeta(x,par(1),par(2),par(3),.true.)
  end function gbeta_logpdf_cb
  function gbeta_moment_cb(order,par) result(v)
    integer,intent(in)::order;real(dp),intent(in)::par(:);real(dp)::v;v=mgbeta(real(order,dp),par(1),par(2),par(3))
  end function gbeta_moment_cb

  subroutine beta_start(x,par)
    real(dp),intent(in)::x(:);real(dp),intent(out)::par(2);real(dp)::m,v,aux
    m=sum(x)/real(size(x),dp);v=sum((x-m)**2)/real(size(x),dp)
    if(v>0.0_dp.and.m>0.0_dp.and.m<1.0_dp)then
      aux=m*(1.0_dp-m)/v-1.0_dp
      if(aux>0.0_dp)then;par=[m*aux,(1.0_dp-m)*aux];return;end if
    end if
    par=[1.0_dp,1.0_dp]
  end subroutine beta_start
  subroutine gbeta_start(x,par)
    real(dp),intent(in)::x(:);real(dp),intent(out)::par(3);real(dp)::bp(2)
    call beta_start(x,bp);par=[1.0_dp,bp]
  end subroutine gbeta_start

  function loglik_dispatch(data,family,par) result(v)
    real(dp),intent(in)::data(:),par(:);integer,intent(in)::family;real(dp)::v,d
    integer::i
    v=0.0_dp
    do i=1,size(data)
      select case(family)
      case(dist_oiunif);d=doiunif(data(i),par(1))
      case(dist_oistpareto);d=doistpareto(data(i),par(1),par(2))
      case(dist_oibeta);d=doibeta(data(i),par(1),par(2),par(3))
      case(dist_oigbeta);d=doigbeta(data(i),par(1),par(2),par(3),par(4))
      case(dist_mbbefd_ab);d=dmbbefd(data(i),par(1),par(2))
      case(dist_mbbefd_gb);d=dmbbefd_gb(data(i),par(1),par(2))
      end select
      if(d<=0.0_dp.or..not.ieee_is_finite(d))then;v=-huge(1.0_dp);return;end if
      v=v+log(d)
    end do
  end function loglik_dispatch

  subroutine finish_information(result)
    type(dr_fit_result),intent(inout)::result
    integer::k
    k=size(result%estimate)
    result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(k,dp)
    result%bic=-2.0_dp*result%log_likelihood+real(k,dp)*log(real(result%nobs,dp))
  end subroutine finish_information

  subroutine initialize_dr_result(result,dist,n)
    type(dr_fit_result),intent(out)::result;character(len=*),intent(in)::dist;integer,intent(in)::n
    result%distribution=lower_ascii(trim(dist));result%nobs=n;result%message="not fitted"
  end subroutine initialize_dr_result

  subroutine make_correlation(cov,cor)
    real(dp),intent(in)::cov(:,:);real(dp),allocatable,intent(out)::cor(:,:)
    integer::i,j,n;n=size(cov,1);allocate(cor(n,n));cor=0.0_dp
    do i=1,n;do j=1,n
      if(cov(i,i)>0.0_dp.and.cov(j,j)>0.0_dp)cor(i,j)=cov(i,j)/sqrt(cov(i,i)*cov(j,j))
    end do;end do
  end subroutine make_correlation

  pure function diagonal(a) result(d)
    real(dp),intent(in)::a(:,:);real(dp)::d(min(size(a,1),size(a,2)));integer::i
    do i=1,size(d);d(i)=a(i,i);end do
  end function diagonal

  subroutine successful_values(x,conv,y)
    real(dp),intent(in)::x(:);integer,intent(in)::conv(:);real(dp),allocatable,intent(out)::y(:)
    integer::i,j,n;n=count(conv==0);allocate(y(n));j=0
    do i=1,size(x);if(conv(i)==0.and.ieee_is_finite(x(i)))then;j=j+1;y(j)=x(i);end if;end do
    if(j<n)y=y(:j)
  end subroutine successful_values

  pure function lower_ascii(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,c
    out=s;do i=1,len(s);c=iachar(out(i:i));if(c>=65.and.c<=90)out(i:i)=achar(c+32);end do
  end function lower_ascii

end module mbbefd_fit
