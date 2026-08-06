! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_model
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ghyp_kinds, only : dp, pi
  use tsd_types, only : distribution_parameters, valid_distribution, canonical_distribution_name
  use tsd_distributions, only : ddist, pdist, qdist
  use tsgarch_types
  implicit none
  private

  public :: initialize_parameters, validate_specification, filter_garch
  public :: persistence, unconditional_variance, distribution_power_moment
  public :: news_impact, backcast_variance, model_equation, effective_omega

contains

  pure function lower_string(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, k
    value = text
    do i = 1, len(text)
      k = iachar(value(i:i))
      if (k >= iachar('A') .and. k <= iachar('Z')) value(i:i) = achar(k + 32)
    end do
  end function lower_string

  subroutine validate_specification(spec, status, message)
    type(garch_spec), intent(inout) :: spec
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=12) :: model
    model = trim(lower_string(spec%model))
    spec%model = model
    spec%distribution = canonical_distribution_name(spec%distribution)
    status = tsg_invalid_argument
    message = ''
    select case (trim(model))
    case ('garch','gjrgarch','aparch','egarch','fgarch','cgarch','igarch','ewma')
    case default
      message = 'unknown GARCH model'
      return
    end select
    if (.not. valid_distribution(spec%distribution)) then
      message = 'unknown innovation distribution'
      return
    end if
    if (spec%p < 0 .or. spec%q < 0) then
      message = 'orders must be nonnegative'
      return
    end if
    if (trim(model) == 'cgarch' .and. (spec%p < 1 .or. spec%q < 1)) then
      message = 'component GARCH requires p and q of at least one'
      return
    end if
    if (trim(model) == 'ewma') then
      spec%p = 1
      spec%q = 1
      spec%variance_targeting = .false.
      spec%multiplicative = .false.
    end if
    if (trim(model) == 'igarch') spec%variance_targeting = .false.
    if (trim(model) == 'egarch') spec%multiplicative = .false.
    if (spec%variance_targeting .and. spec%multiplicative) spec%multiplicative = .false.
    if (spec%backcast_lambda < 0.0_dp .or. spec%backcast_lambda > 1.0_dp) then
      message = 'backcast lambda must be in [0,1]'
      return
    end if
    if (spec%sample_n < 1) then
      message = 'sample_n must be positive'
      return
    end if
    if (spec%stationarity_limit <= 0.0_dp .or. spec%stationarity_limit > 1.0_dp) then
      message = 'stationarity limit must be in (0,1]'
      return
    end if
    status = tsg_success
    message = 'ok'
  end subroutine validate_specification

  function initialize_parameters(y, spec, nreg) result(par)
    real(dp), intent(in) :: y(:)
    type(garch_spec), intent(in) :: spec
    integer, intent(in), optional :: nreg
    type(garch_parameters) :: par
    integer :: nr
    real(dp) :: v, m
    nr = 0
    if (present(nreg)) nr = max(0,nreg)
    allocate(par%alpha(max(spec%p,0)), par%beta(max(spec%q,0)))
    allocate(par%gamma(max(spec%p,0)), par%eta(max(spec%p,0)), par%xi(nr))
    par%alpha = 0.0_dp
    par%beta = 0.0_dp
    par%gamma = 0.0_dp
    par%eta = 0.0_dp
    par%xi = 0.0_dp
    if (size(y) > 0) then
      m = sum(y)/real(size(y),dp)
      v = sum((y-m)**2)/real(max(1,size(y)),dp)
    else
      m = 0.0_dp
      v = 1.0_dp
    end if
    v = max(v,1.0e-8_dp)
    if (spec%constant) par%mu = m
    par%omega = 0.05_dp*v
    if (spec%p > 0) par%alpha = 0.08_dp/real(spec%p,dp)
    if (spec%q > 0) par%beta = 0.88_dp/real(spec%q,dp)
    select case(trim(spec%model))
    case('egarch')
      par%omega = log(v)*max(0.05_dp,1.0_dp-sum(par%beta))
      par%gamma = 0.05_dp
    case('gjrgarch')
      par%gamma = 0.04_dp/real(max(1,spec%p),dp)
    case('aparch')
      par%delta = 2.0_dp
      par%gamma = 0.05_dp
    case('fgarch')
      par%delta = 2.0_dp
      par%gamma = 0.05_dp
      par%eta = 0.0_dp
    case('cgarch')
      par%rho = 0.95_dp
      par%phi = 0.03_dp
      if (spec%p > 0) par%alpha = 0.04_dp/real(spec%p,dp)
      if (spec%q > 0) par%beta = 0.80_dp/real(spec%q,dp)
      par%omega = (1.0_dp-par%rho)*v
    case('igarch')
      if (spec%p > 0) par%alpha = 0.08_dp/real(spec%p,dp)
      if (spec%q > 0) par%beta = 0.92_dp/real(spec%q,dp)
      call normalize_igarch(par%alpha,par%beta)
    case('ewma')
      par%alpha = 0.06_dp
      par%beta = 0.94_dp
      par%omega = 0.0_dp
    end select
    select case(trim(spec%distribution))
    case('norm')
      par%dist%skew=0.0_dp
      par%dist%shape=0.0_dp
      par%dist%lambda=-0.5_dp
    case('std')
      par%dist%shape=6.0_dp
      par%dist%skew=0.0_dp
    case('snorm')
      par%dist%skew=1.2_dp
      par%dist%shape=0.0_dp
    case('sstd')
      par%dist%skew=1.1_dp
      par%dist%shape=6.0_dp
    case('ged')
      par%dist%shape=2.0_dp
      par%dist%skew=0.0_dp
    case('sged')
      par%dist%skew=1.1_dp
      par%dist%shape=2.0_dp
    case('nig')
      par%dist%skew=0.0_dp
      par%dist%shape=1.0_dp
    case('gh')
      par%dist%skew=0.0_dp
      par%dist%shape=1.0_dp
      par%dist%lambda=-0.5_dp
    case('jsu')
      par%dist%skew=0.0_dp
      par%dist%shape=1.5_dp
    case('ghst')
      par%dist%skew=0.0_dp
      par%dist%shape=6.0_dp
    end select
    par%dist%mu = 0.0_dp
    par%dist%sigma = 1.0_dp
  end function initialize_parameters

  subroutine normalize_igarch(alpha,beta)
    real(dp), intent(inout) :: alpha(:),beta(:)
    real(dp) :: s
    s = sum(max(alpha,0.0_dp))+sum(max(beta,0.0_dp))
    if (s <= tiny(1.0_dp)) then
      if (size(alpha)>0) alpha=0.05_dp/real(size(alpha),dp)
      if (size(beta)>0) beta=0.95_dp/real(size(beta),dp)
    else
      alpha=max(alpha,0.0_dp)/s
      beta=max(beta,0.0_dp)/s
    end if
  end subroutine normalize_igarch

  real(dp) function backcast_variance(y, lambda, delta) result(value)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in), optional :: lambda, delta
    real(dp) :: lam, d, avg
    integer :: i, n
    lam=0.7_dp
    if(present(lambda))lam=lambda
    d=2.0_dp
    if(present(delta))d=delta
    n=size(y)
    if(n==0)then
    value=1.0_dp
    return
    end if
    avg=sum(abs(y)**d)/real(n,dp)
    value=(lam**n)*avg
    do i=1,n
      value=value+(1.0_dp-lam)*lam**(i-1)*abs(y(i))**d
    end do
  end function backcast_variance

  real(dp) function initial_variance(residuals,spec,delta) result(value)
    real(dp),intent(in)::residuals(:)
    type(garch_spec),intent(in)::spec
    real(dp),intent(in),optional::delta
    real(dp)::d
    integer::m
    d=2.0_dp
    if(present(delta))d=delta
    select case(spec%initialization)
    case(init_sample)
      m=min(size(residuals),max(1,spec%sample_n))
      value=sum(abs(residuals(1:m))**2)/real(m,dp)
    case(init_backcast)
      value=backcast_variance(residuals,spec%backcast_lambda,2.0_dp)
    case default
      value=sum(residuals**2)/real(max(1,size(residuals)),dp)
    end select
    value=max(value,1.0e-10_dp)
    if(abs(d-2.0_dp)>1.0e-12_dp)value=value**(d/2.0_dp)
  end function initial_variance

  real(dp) function distribution_power_moment(distribution, pars, gamma, delta, eta, mode) result(value)
    character(len=*),intent(in)::distribution
    type(distribution_parameters),intent(in)::pars
    real(dp),intent(in),optional::gamma,delta,eta
    integer,intent(in),optional::mode
    integer,parameter::nq=64
    real(dp)::g,d,e,p,z,s
    integer::i,md
    g=0.0_dp
    if(present(gamma))g=gamma
    d=1.0_dp
    if(present(delta))d=delta
    e=0.0_dp
    if(present(eta))e=eta
    md=1
    if(present(mode))md=mode
    if(trim(distribution)=='norm')then
      if(md==1 .and. abs(d-1.0_dp)<1.0e-12_dp .and. abs(g)<1.0e-12_dp .and. abs(e)<1.0e-12_dp)then
        value=sqrt(2.0_dp/pi)
        return
      else if(md==2)then
        value=0.5_dp
        return
      else if(md==1 .and. abs(d-2.0_dp)<1.0e-12_dp .and. abs(e)<1.0e-12_dp)then
        value=1.0_dp+g*g
        return
      end if
    end if
    s=0.0_dp
    do i=1,nq
      p=(real(i,dp)-0.5_dp)/real(nq,dp)
      z=qdist(distribution,p,pars)
      select case(md)
      case(2)
        if(z<=0.0_dp)s=s+1.0_dp
      case default
        s=s+(abs(z-e)-g*(z-e))**d
      end select
    end do
    value=s/real(nq,dp)
  end function distribution_power_moment

  real(dp) function persistence(spec,par) result(value)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp)::k,pneg
    integer :: i
    select case(trim(spec%model))
    case('egarch')
      value=sum(par%beta)
    case('gjrgarch')
      pneg=distribution_power_moment(spec%distribution,par%dist,mode=2)
      value=sum(par%alpha)+sum(par%beta)+pneg*sum(par%gamma)
    case('aparch')
      k=0.0_dp
      if(size(par%alpha)>0)then
        k=sum([(distribution_power_moment(spec%distribution,par%dist, &
          par%gamma(i),par%delta,0.0_dp,1),i=1,size(par%alpha))]*par%alpha)
      end if
      value=sum(par%beta)+k
    case('fgarch')
      k=0.0_dp
      if(size(par%alpha)>0)then
        k=sum([(distribution_power_moment(spec%distribution,par%dist, &
          par%gamma(i),par%delta,par%eta(i),1),i=1,size(par%alpha))]*par%alpha)
      end if
      value=sum(par%beta)+k
    case('cgarch')
      value=max(par%rho,sum(par%alpha)+sum(par%beta))
    case('ewma','igarch')
      value=1.0_dp
    case default
      value=sum(par%alpha)+sum(par%beta)
    end select
  end function persistence

  real(dp) function effective_omega(y,spec,par,vreg) result(value)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in),optional::vreg(:,:)
    real(dp)::m,v,p,avgv
    integer :: i
    m=par%mu
    v=sum((y-m)**2)/real(max(1,size(y)),dp)
    p=persistence(spec,par)
    avgv=0.0_dp
    if(present(vreg).and.size(par%xi)>0) avgv=sum(matmul(transpose(vreg),[(1.0_dp/real(size(y),dp),i=1,size(y))])*par%xi)
    if(.not.spec%variance_targeting)then
      value=par%omega
      return
    end if
    select case(trim(spec%model))
    case('egarch')
      value=(1.0_dp-sum(par%beta))*log(max(v,1.0e-12_dp))-avgv
    case('aparch','fgarch')
      value=max(v,1.0e-12_dp)**(par%delta/2.0_dp)*max(1.0_dp-p,1.0e-8_dp)-avgv
    case('cgarch')
      value=max(v,1.0e-12_dp)*max(1.0_dp-par%rho,1.0e-8_dp)-avgv
    case default
      value=max(v,1.0e-12_dp)*max(1.0_dp-p,1.0e-8_dp)-avgv
    end select
  end function effective_omega

  real(dp) function unconditional_variance(spec,par,mean_vreg) result(value)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in),optional::mean_vreg
    real(dp)::p,num,v
    p=persistence(spec,par)
    num=par%omega
    if(present(mean_vreg))num=num+mean_vreg
    if(spec%multiplicative)num=exp(num)
    select case(trim(spec%model))
    case('egarch')
      if(abs(1.0_dp-sum(par%beta))<1.0e-10_dp)then
      value=huge(1.0_dp)
      else
      value=exp(num/(1.0_dp-sum(par%beta)))
      end if
    case('aparch','fgarch')
      if(p>=1.0_dp)then
      value=huge(1.0_dp)
      else
      v=num/(1.0_dp-p)
      value=max(v,0.0_dp)**(2.0_dp/par%delta)
      end if
    case('cgarch')
      if(par%rho>=1.0_dp)then
      value=huge(1.0_dp)
      else
      value=num/(1.0_dp-par%rho)
      end if
    case('igarch','ewma')
      value=huge(1.0_dp)
    case default
      if(p>=1.0_dp)then
      value=huge(1.0_dp)
      else
      value=num/(1.0_dp-p)
      end if
    end select
  end function unconditional_variance

  function filter_garch(y,spec,par,vreg,initial_state) result(out)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in),optional::vreg(:,:)
    real(dp),intent(in),optional::initial_state(:)
    type(garch_filter_result)::out
    integer::n,m,i,j
    real(dp)::iv,base,vterm,kappa,powv,ll
    real(dp),allocatable::a(:),b(:),power_sigma(:),log_variance(:)
    type(distribution_parameters)::dpar
    logical::bad
    n=size(y)
    m=max(spec%p,spec%q)
    allocate(out%sigma(n),out%variance(n),out%residuals(n),out%standardized_residuals(n),out%loglik_vector(n))
    allocate(out%permanent_component(n),out%transitory_component(n))
    out%sigma=0.0_dp
    out%variance=0.0_dp
    out%residuals=y-par%mu
    out%standardized_residuals=0.0_dp
    out%loglik_vector=0.0_dp
    out%permanent_component=0.0_dp
    out%transitory_component=0.0_dp
    if(n<max(2,m+1))then
    out%message='series is too short'
    return
    end if
    if(present(vreg))then
      if(size(vreg,1)/=n.or.size(vreg,2)/=size(par%xi))then
      out%message='variance regressor dimensions do not conform'
      return
      end if
    end if
    allocate(a(size(par%alpha)),b(size(par%beta)))
    a=par%alpha
    b=par%beta
    if(trim(spec%model)=='igarch')call normalize_igarch(a,b)
    if(trim(spec%model)=='ewma')then
      if(size(b)==0)then
      out%message='EWMA requires beta'
      return
      end if
      b(1)=min(max(b(1),1.0e-6_dp),0.999999_dp)
      a=0.0_dp
      if(size(a)>0)a(1)=1.0_dp-b(1)
    end if
    iv=initial_variance(out%residuals,spec,merge(par%delta,2.0_dp,trim(spec%model)=='aparch'.or.trim(spec%model)=='fgarch'))
    if(present(initial_state))then
      if(size(initial_state)>0)iv=max(initial_state(1),1.0e-12_dp)
    end if
    out%effective_omega=effective_omega(y,spec,par,vreg)
    allocate(power_sigma(n),log_variance(n))
    power_sigma=0.0_dp
    log_variance=0.0_dp
    do i=1,m
      out%variance(i)=merge(iv**(2.0_dp/par%delta),iv,trim(spec%model)=='aparch'.or.trim(spec%model)=='fgarch')
      if(trim(spec%model)=='aparch'.or.trim(spec%model)=='fgarch')power_sigma(i)=iv
      if(trim(spec%model)=='egarch')log_variance(i)=log(max(out%variance(i),1.0e-12_dp))
      if(trim(spec%model)=='cgarch')then
      out%permanent_component(i)=out%variance(i)
      out%transitory_component(i)=0.0_dp
      end if
      out%sigma(i)=sqrt(max(out%variance(i),1.0e-12_dp))
      out%standardized_residuals(i)=out%residuals(i)/out%sigma(i)
    end do
    if(m==0)then
      m=0
    end if
    kappa=distribution_power_moment(spec%distribution,par%dist,gamma=0.0_dp,delta=1.0_dp,eta=0.0_dp,mode=1)
    bad=.false.
    do i=m+1,n
      vterm=0.0_dp
      if(present(vreg).and.size(par%xi)>0)vterm=dot_product(vreg(i,:),par%xi)
      select case(trim(spec%model))
      case('garch','igarch','ewma')
        base=out%effective_omega+vterm
        if(spec%multiplicative)base=exp(base)
        out%variance(i)=base
        do j=1,min(spec%p,i-1)
        out%variance(i)=out%variance(i)+a(j)*out%residuals(i-j)**2
        end do
        do j=1,min(spec%q,i-1)
        out%variance(i)=out%variance(i)+b(j)*out%variance(i-j)
        end do
      case('gjrgarch')
        base=out%effective_omega+vterm
        if(spec%multiplicative)base=exp(base)
        out%variance(i)=base
        do j=1,min(spec%p,i-1)
          out%variance(i)=out%variance(i)+a(j)*out%residuals(i-j)**2
          if(out%residuals(i-j)<=0.0_dp)out%variance(i)=out%variance(i)+par%gamma(j)*out%residuals(i-j)**2
        end do
        do j=1,min(spec%q,i-1)
        out%variance(i)=out%variance(i)+b(j)*out%variance(i-j)
        end do
      case('egarch')
        log_variance(i)=out%effective_omega+vterm
        do j=1,min(spec%p,i-1)
          log_variance(i)=log_variance(i)+a(j)*out%standardized_residuals(i-j) &
            +par%gamma(j)*(abs(out%standardized_residuals(i-j))-kappa)
        end do
        do j=1,min(spec%q,i-1)
        log_variance(i)=log_variance(i)+b(j)*log_variance(i-j)
        end do
        out%variance(i)=exp(min(max(log_variance(i),-700.0_dp),700.0_dp))
      case('aparch')
        base=out%effective_omega+vterm
        if(spec%multiplicative)base=exp(base)
        power_sigma(i)=base
        do j=1,min(spec%p,i-1)
          power_sigma(i)=power_sigma(i)+a(j)*(abs(out%residuals(i-j))-par%gamma(j)*out%residuals(i-j))**par%delta
        end do
        do j=1,min(spec%q,i-1)
        power_sigma(i)=power_sigma(i)+b(j)*power_sigma(i-j)
        end do
        out%variance(i)=max(power_sigma(i),1.0e-14_dp)**(2.0_dp/par%delta)
      case('fgarch')
        base=out%effective_omega+vterm
        if(spec%multiplicative)base=exp(base)
        power_sigma(i)=base
        do j=1,min(spec%p,i-1)
          powv=max(out%variance(i-j),1.0e-14_dp)**(par%delta/2.0_dp)
          power_sigma(i)=power_sigma(i)+a(j)*powv* &
            (abs(out%standardized_residuals(i-j)-par%eta(j))-par%gamma(j)* &
            (out%standardized_residuals(i-j)-par%eta(j)))**par%delta
        end do
        do j=1,min(spec%q,i-1)
        power_sigma(i)=power_sigma(i)+b(j)*power_sigma(i-j)
        end do
        out%variance(i)=max(power_sigma(i),1.0e-14_dp)**(2.0_dp/par%delta)
      case('cgarch')
        base=out%effective_omega+vterm
        if(spec%multiplicative)base=exp(base)
        out%permanent_component(i)=base+par%rho*out%permanent_component(i-1)+par%phi*(out%residuals(i-1)**2-out%variance(i-1))
        out%transitory_component(i)=0.0_dp
        do j=1,min(spec%p,i-1)
          out%transitory_component(i)=out%transitory_component(i) &
            +a(j)*out%transitory_component(i-j) &
            +a(j)*(out%residuals(i-j)**2-out%variance(i-j))
        end do
        do j=1,min(spec%q,i-1)
        out%transitory_component(i)=out%transitory_component(i)+b(j)*out%transitory_component(i-j)
        end do
        out%variance(i)=out%permanent_component(i)+out%transitory_component(i)
      end select
      if(.not.ieee_is_finite(out%variance(i)).or.out%variance(i)<=1.0e-14_dp.or.out%variance(i)>1.0e100_dp)then
      bad=.true.
      exit
      end if
      out%sigma(i)=sqrt(out%variance(i))
      out%standardized_residuals(i)=out%residuals(i)/out%sigma(i)
    end do
    if(bad)then
    out%message='variance recursion became nonpositive or nonfinite'
    return
    end if
    dpar=par%dist
    dpar%mu=0.0_dp
    dpar%sigma=1.0_dp
    ll=0.0_dp
    do i=1,n
      out%loglik_vector(i)=ddist(spec%distribution,out%standardized_residuals(i),dpar,.true.)-log(out%sigma(i))
      if(.not.ieee_is_finite(out%loglik_vector(i)))then
      out%message='nonfinite innovation likelihood'
      return
      end if
      ll=ll+out%loglik_vector(i)
    end do
    out%log_likelihood=ll
    out%nobs=n
    out%persistence=persistence(spec,par)
    out%unconditional_variance=unconditional_variance(spec,par)
    out%status=tsg_success
    out%message='ok'
  end function filter_garch

  function news_impact(epsilon,spec,par) result(variance)
    real(dp),intent(in)::epsilon(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),allocatable::variance(:)
    real(dp)::uv,kappa,base,z
    integer::i
    allocate(variance(size(epsilon)))
    uv=unconditional_variance(spec,par)
    if(.not.ieee_is_finite(uv))uv=1.0_dp
    base=par%omega
    if(spec%multiplicative)base=exp(base)
    kappa=distribution_power_moment(spec%distribution,par%dist,gamma=0.0_dp,delta=1.0_dp,eta=0.0_dp,mode=1)
    do i=1,size(epsilon)
      select case(trim(spec%model))
      case('garch')
        variance(i)=base+merge(par%alpha(1)*epsilon(i)**2,0.0_dp,size(par%alpha)>0)+merge(par%beta(1)*uv,0.0_dp,size(par%beta)>0)
      case('gjrgarch')
        variance(i)=base+par%alpha(1)*epsilon(i)**2+par%beta(1)*uv
        if(epsilon(i)<=0.0_dp)variance(i)=variance(i)+par%gamma(1)*epsilon(i)**2
      case('egarch')
        z=epsilon(i)/sqrt(uv)
        variance(i)=exp(par%omega+par%alpha(1)*z+par%gamma(1)*(abs(z)-kappa)+par%beta(1)*log(uv))
      case('aparch')
        variance(i)=(base+par%alpha(1)* &
          (abs(epsilon(i))-par%gamma(1)*epsilon(i))**par%delta &
          +par%beta(1)*uv**(par%delta/2.0_dp))**(2.0_dp/par%delta)
      case('fgarch')
        z=epsilon(i)/sqrt(uv)
        variance(i)=(base+par%alpha(1)*uv**(par%delta/2.0_dp)* &
          (abs(z-par%eta(1))-par%gamma(1)*(z-par%eta(1)))**par%delta &
          +par%beta(1)*uv**(par%delta/2.0_dp))**(2.0_dp/par%delta)
      case('cgarch')
        variance(i)=base+par%rho*uv+par%phi*(epsilon(i)**2-uv)+par%alpha(1)*(epsilon(i)**2-uv)
      case default
        variance(i)=huge(1.0_dp)
      end select
    end do
  end function news_impact

  function model_equation(spec) result(text)
    type(garch_spec),intent(in)::spec
    character(len=:),allocatable::text
    select case(trim(spec%model))
    case('garch')
    text='sigma_t^2 = omega + sum alpha_i epsilon_{t-i}^2 + sum beta_j sigma_{t-j}^2'
    case('gjrgarch')
    text='sigma_t^2 = omega + sum (alpha_i + gamma_i I[e<0]) epsilon_{t-i}^2 + sum beta_j sigma_{t-j}^2'
    case('egarch')
    text='log sigma_t^2 = omega + sum alpha_i z_{t-i} + gamma_i(|z_{t-i}|-E|z|) + sum beta_j log sigma_{t-j}^2'
    case('aparch')
    text='sigma_t^delta = omega + sum alpha_i(|epsilon|-gamma_i epsilon)^delta + sum beta_j sigma_{t-j}^delta'
    case('fgarch')
    text='sigma_t^delta = omega + sum alpha_i sigma_{t-i}^delta(|z-eta_i|-gamma_i(z-eta_i))^delta + sum beta_j sigma_{t-j}^delta'
    case('cgarch')
    text='sigma_t^2 = q_t + h_t with permanent and transitory component recursions'
    case('igarch')
    text='GARCH recursion constrained to unit persistence'
    case('ewma')
    text='sigma_t^2 = (1-lambda) epsilon_{t-1}^2 + lambda sigma_{t-1}^2'
    case default
    text='unknown model'
    end select
  end function model_equation

end module tsgarch_model
