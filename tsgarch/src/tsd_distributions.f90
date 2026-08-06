! SPDX-License-Identifier: GPL-2.0-only
module tsd_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
  use ghyp_kinds, only : dp, i8, pi, sqrt_two_pi
  use ghyp_special, only : log_bessel_k, bessel_k
  use ghyp_rng, only : rng_state, seed_rng, uniform_rng, normal_rng, gamma_rng
  use ghyp_model, only : ghyp_model_type, ghyp_ad
  use ghyp_distribution, only : ghyp_density => dghyp, ghyp_cdf => pghyp, &
                                ghyp_quantile => qghyp, ghyp_random_one => rghyp_one
  use tsd_types, only : distribution_parameters, distribution_id, &
                        dist_norm, dist_std, dist_snorm, dist_sstd, dist_ged, &
                        dist_sged, dist_nig, dist_gh, dist_jsu, dist_ghst
  use tsd_math, only : normal_pdf, normal_cdf, normal_quantile, student_pdf, student_cdf, &
                       student_quantile, regularized_gamma_p, gamma_quantile, signum, &
                       heaviside, clamp_probability
  implicit none
  private

  public :: rng_state, seed_rng
  public :: ddist, pdist, qdist, rdist
  public :: dstd, pstd, qstd, rstd
  public :: dsnorm, psnorm, qsnorm, rsnorm
  public :: dsstd, psstd, qsstd, rsstd
  public :: dged, pged, qged, rged
  public :: dsged, psged, qsged, rsged
  public :: djsu, pjsu, qjsu, rjsu
  public :: dghyp_raw, pghyp_raw, qghyp_raw, rghyp_raw
  public :: dghyp, pghyp, qghyp, rghyp
  public :: dnig, pnig, qnig, rnig
  public :: dgh, pgh, qgh, rgh
  public :: dghst, pghst, qghst, rghst
  public :: paramgh, paramghst, nigtransform, ghyptransform

contains

  real(dp) function ddist(distribution, x, parameters, log_density) result(value)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in) :: x
    type(distribution_parameters), intent(in), optional :: parameters
    logical, intent(in), optional :: log_density
    type(distribution_parameters) :: p
    logical :: lg
    p = distribution_parameters()
    if (present(parameters)) p = parameters
    lg = .false.
    if (present(log_density)) lg = log_density
    select case(distribution_id(distribution))
    case(dist_norm)
      value = normal_location_density(x,p%mu,p%sigma,lg)
    case(dist_std)
      value = dstd(x,p%mu,p%sigma,p%shape,lg)
    case(dist_snorm)
      value = dsnorm(x,p%mu,p%sigma,p%skew,lg)
    case(dist_sstd)
      value = dsstd(x,p%mu,p%sigma,p%skew,p%shape,lg)
    case(dist_ged)
      value = dged(x,p%mu,p%sigma,p%shape,lg)
    case(dist_sged)
      value = dsged(x,p%mu,p%sigma,p%skew,p%shape,lg)
    case(dist_nig)
      value = dnig(x,p%mu,p%sigma,p%skew,p%shape,lg)
    case(dist_gh)
      value = dgh(x,p%mu,p%sigma,p%skew,p%shape,p%lambda,lg)
    case(dist_jsu)
      value = djsu(x,p%mu,p%sigma,p%skew,p%shape,lg)
    case(dist_ghst)
      value = dghst(x,p%mu,p%sigma,p%skew,p%shape,lg)
    case default
      value = ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function ddist

  real(dp) function pdist(distribution, q, parameters, lower_tail, log_probability) result(value)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in) :: q
    type(distribution_parameters), intent(in), optional :: parameters
    logical, intent(in), optional :: lower_tail, log_probability
    type(distribution_parameters) :: p
    logical :: lower, logp
    p=distribution_parameters()
    if(present(parameters))p=parameters
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    logp=.false.
    if(present(log_probability))logp=log_probability
    select case(distribution_id(distribution))
    case(dist_norm)
      value=normal_location_cdf(q,p%mu,p%sigma)
    case(dist_std)
      value=pstd(q,p%mu,p%sigma,p%shape)
    case(dist_snorm)
      value=psnorm(q,p%mu,p%sigma,p%skew)
    case(dist_sstd)
      value=psstd(q,p%mu,p%sigma,p%skew,p%shape)
    case(dist_ged)
      value=pged(q,p%mu,p%sigma,p%shape)
    case(dist_sged)
      value=psged(q,p%mu,p%sigma,p%skew,p%shape)
    case(dist_nig)
      value=pnig(q,p%mu,p%sigma,p%skew,p%shape)
    case(dist_gh)
      value=pgh(q,p%mu,p%sigma,p%skew,p%shape,p%lambda)
    case(dist_jsu)
      value=pjsu(q,p%mu,p%sigma,p%skew,p%shape)
    case(dist_ghst)
      value=pghst(q,p%mu,p%sigma,p%skew,p%shape)
    case default
      value=ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end select
    if(.not.lower)value=1.0_dp-value
    value=clamp_probability(value)
    if(logp)value=log(max(value,tiny(1.0_dp)))
  end function pdist

  real(dp) function qdist(distribution, probability, parameters, lower_tail, log_probability) result(value)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in) :: probability
    type(distribution_parameters), intent(in), optional :: parameters
    logical, intent(in), optional :: lower_tail, log_probability
    type(distribution_parameters) :: p
    real(dp) :: pr
    logical :: lower, logp
    p=distribution_parameters()
    if(present(parameters))p=parameters
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    logp=.false.
    if(present(log_probability))logp=log_probability
    pr=probability
    if(logp)pr=exp(pr)
    if(.not.lower)pr=1.0_dp-pr
    select case(distribution_id(distribution))
    case(dist_norm)
      value=p%mu+p%sigma*normal_quantile(pr)
    case(dist_std)
      value=qstd(pr,p%mu,p%sigma,p%shape)
    case(dist_snorm)
      value=qsnorm(pr,p%mu,p%sigma,p%skew)
    case(dist_sstd)
      value=qsstd(pr,p%mu,p%sigma,p%skew,p%shape)
    case(dist_ged)
      value=qged(pr,p%mu,p%sigma,p%shape)
    case(dist_sged)
      value=qsged(pr,p%mu,p%sigma,p%skew,p%shape)
    case(dist_nig)
      value=qnig(pr,p%mu,p%sigma,p%skew,p%shape)
    case(dist_gh)
      value=qgh(pr,p%mu,p%sigma,p%skew,p%shape,p%lambda)
    case(dist_jsu)
      value=qjsu(pr,p%mu,p%sigma,p%skew,p%shape)
    case(dist_ghst)
      value=qghst(pr,p%mu,p%sigma,p%skew,p%shape)
    case default
      value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function qdist

  function rdist(distribution, n, rng, parameters) result(x)
    character(len=*), intent(in) :: distribution
    integer, intent(in) :: n
    type(rng_state), intent(inout) :: rng
    type(distribution_parameters), intent(in), optional :: parameters
    real(dp), allocatable :: x(:)
    type(distribution_parameters) :: p
    integer :: i
    p=distribution_parameters()
    if(present(parameters))p=parameters
    allocate(x(max(n,0)))
    do i=1,n
      select case(distribution_id(distribution))
      case(dist_norm)
      x(i)=p%mu+p%sigma*normal_rng(rng)
      case(dist_std)
      x(i)=rstd_one(rng,p%mu,p%sigma,p%shape)
      case(dist_snorm)
      x(i)=rsnorm_one(rng,p%mu,p%sigma,p%skew)
      case(dist_sstd)
      x(i)=rsstd_one(rng,p%mu,p%sigma,p%skew,p%shape)
      case(dist_ged)
      x(i)=rged_one(rng,p%mu,p%sigma,p%shape)
      case(dist_sged)
      x(i)=rsged_one(rng,p%mu,p%sigma,p%skew,p%shape)
      case(dist_nig)
      x(i)=rnig_one(rng,p%mu,p%sigma,p%skew,p%shape)
      case(dist_gh)
      x(i)=rgh_one(rng,p%mu,p%sigma,p%skew,p%shape,p%lambda)
      case(dist_jsu)
      x(i)=rjsu_one(rng,p%mu,p%sigma,p%skew,p%shape)
      case(dist_ghst)
      x(i)=rghst_one(rng,p%mu,p%sigma,p%skew,p%shape)
      case default
      x(i)=ieee_value(0.0_dp,ieee_quiet_nan)
      end select
    end do
  end function rdist

  pure real(dp) function normal_location_density(x,mu,sigma,log_density) result(value)
    real(dp),intent(in)::x,mu,sigma
    logical,intent(in)::log_density
    if(sigma<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    if(log_density)then
      value=-log(sigma*sqrt_two_pi)-0.5_dp*((x-mu)/sigma)**2
    else
      value=normal_pdf((x-mu)/sigma)/sigma
    end if
  end function normal_location_density

  pure real(dp) function normal_location_cdf(x,mu,sigma) result(value)
    real(dp),intent(in)::x,mu,sigma
    if(sigma<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    value=normal_cdf((x-mu)/sigma)
    end if
  end function normal_location_cdf

  real(dp) function dstd(x,mu,sigma,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,nu,z,scale
    logical::lg
    m=0.0_dp
    s=1.0_dp
    nu=5.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    if(present(log_density))lg=log_density
    if(s<=0.0_dp.or.nu<=2.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    scale=sqrt(nu/(nu-2.0_dp))
    z=(x-m)/s
    value=scale*student_pdf(z*scale,nu)/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dstd

  real(dp) function pstd(q,mu,sigma,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp)::m,s,nu,scale
    m=0.0_dp
    s=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    if(s<=0.0_dp.or.nu<=2.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    scale=sqrt(nu/(nu-2.0_dp))
    value=student_cdf((q-m)/s*scale,nu)
  end function pstd

  real(dp) function qstd(p,mu,sigma,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp)::m,s,nu,scale
    m=0.0_dp
    s=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    if(s<=0.0_dp.or.nu<=2.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    scale=sqrt(nu/(nu-2.0_dp))
    value=m+s*student_quantile(p,nu)/scale
  end function qstd

  function rstd(n,rng,mu,sigma,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,nu
    integer::i
    m=0.0_dp
    s=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rstd_one(rng,m,s,nu)
    end do
  end function rstd

  real(dp) function rstd_one(rng,mu,sigma,nu) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,nu
    real(dp)::chi
    if(nu<=2.0_dp.or.sigma<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    chi=gamma_rng(0.5_dp*nu,2.0_dp,rng)
    value=mu+sigma*(normal_rng(rng)/sqrt(chi/nu))/sqrt(nu/(nu-2.0_dp))
  end function rstd_one

  subroutine fs_standardization(base_m1,xi,center,scale)
    real(dp),intent(in)::base_m1,xi
    real(dp),intent(out)::center,scale
    center=base_m1*(xi-1.0_dp/xi)
    scale=sqrt((1.0_dp-base_m1*base_m1)*(xi*xi+1.0_dp/(xi*xi))+2.0_dp*base_m1*base_m1-1.0_dp)
  end subroutine fs_standardization

  real(dp) function dsnorm(x,mu,sigma,skew,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew
    logical,intent(in),optional::log_density
    real(dp)::m,s,xi,z,xxi,c,sc,g
    logical::lg
    m=0.0_dp
    s=1.0_dp
    xi=1.5_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(log_density))lg=log_density
    if(s<=0.0_dp.or.xi<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    call fs_standardization(sqrt(2.0_dp/pi),xi,c,sc)
    z=((x-m)/s)*sc+c
    if(abs(z)<=tiny(1.0_dp))then
    xxi=1.0_dp
    else if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=g*normal_pdf(z/xxi)*sc/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dsnorm

  real(dp) function psnorm(q,mu,sigma,skew) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew
    real(dp)::m,s,xi,z,xxi,c,sc,g
    m=0.0_dp
    s=1.0_dp
    xi=1.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(s<=0.0_dp.or.xi<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    call fs_standardization(sqrt(2.0_dp/pi),xi,c,sc)
    z=((q-m)/s)*sc+c
    if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=heaviside(z,0.0_dp)-signum(z)*g*xxi*normal_cdf(-abs(z)/xxi)
    value=clamp_probability(value)
  end function psnorm

  real(dp) function qsnorm(p,mu,sigma,skew) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew
    real(dp)::m,s,xi,z,xxi,c,sc,g,tmp
    m=0.0_dp
    s=1.0_dp
    xi=1.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(p<=0.0_dp)then
    value=-ieee_value(0.0_dp,ieee_positive_inf)
    return
    else if(p>=1.0_dp)then
    value=ieee_value(0.0_dp,ieee_positive_inf)
    return
    end if
    call fs_standardization(sqrt(2.0_dp/pi),xi,c,sc)
    g=2.0_dp/(xi+1.0_dp/xi)
    z=p-1.0_dp/(1.0_dp+xi*xi)
    xxi=xi**signum(z)
    tmp=(heaviside(z,0.0_dp)-signum(z)*p)/(g*xxi)
    value=m+s*((-signum(z)*normal_quantile(tmp)*xxi-c)/sc)
  end function qsnorm

  function rsnorm(n,rng,mu,sigma,skew) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew
    real(dp),allocatable::x(:)
    real(dp)::m,s,xi
    integer::i
    m=0.0_dp
    s=1.0_dp
    xi=1.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rsnorm_one(rng,m,s,xi)
    end do
  end function rsnorm

  real(dp) function rsnorm_one(rng,mu,sigma,xi) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,xi
    real(dp)::weight,z,xx,rr,c,sc
    if(xi<=0.0_dp.or.sigma<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    weight=xi/(xi+1.0_dp/xi)
    z=uniform_rng(rng)-weight
    if(z<0.0_dp)then
    xx=1.0_dp/xi
    else
    xx=xi
    end if
    rr=-abs(normal_rng(rng))/xx*merge(-1.0_dp,1.0_dp,z<0.0_dp)
    call fs_standardization(sqrt(2.0_dp/pi),xi,c,sc)
    value=mu+sigma*(rr-c)/sc
  end function rsnorm_one

  real(dp) function dsstd(x,mu,sigma,skew,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,xi,nu,m1,c,sc,z,xxi,g
    logical::lg
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=5.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    if(present(log_density))lg=log_density
    if(nu<=2.0_dp.or.xi<=0.0_dp.or.s<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    m1=2.0_dp*sqrt(nu-2.0_dp)/(nu-1.0_dp)/exp(log_gamma(0.5_dp)+log_gamma(nu/2.0_dp)-log_gamma(0.5_dp+nu/2.0_dp))
    call fs_standardization(m1,xi,c,sc)
    z=((x-m)/s)*sc+c
    if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=g*dstd(z/xxi,0.0_dp,1.0_dp,nu)*sc/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dsstd

  real(dp) function psstd(q,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,xi,nu,m1,c,sc,z,xxi,g
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    m1=2.0_dp*sqrt(nu-2.0_dp)/(nu-1.0_dp)/exp(log_gamma(0.5_dp)+log_gamma(nu/2.0_dp)-log_gamma(0.5_dp+nu/2.0_dp))
    call fs_standardization(m1,xi,c,sc)
    z=((q-m)/s)*sc+c
    if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=heaviside(z,0.0_dp)-signum(z)*g*xxi*pstd(-abs(z)/xxi,0.0_dp,1.0_dp,nu)
    value=clamp_probability(value)
  end function psstd

  real(dp) function qsstd(p,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,xi,nu,m1,c,sc,g,z,xxi,tmp
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    m1=2.0_dp*sqrt(nu-2.0_dp)/(nu-1.0_dp)/exp(log_gamma(0.5_dp)+log_gamma(nu/2.0_dp)-log_gamma(0.5_dp+nu/2.0_dp))
    call fs_standardization(m1,xi,c,sc)
    g=2.0_dp/(xi+1.0_dp/xi)
    z=p-1.0_dp/(1.0_dp+xi*xi)
    xxi=xi**signum(z)
    tmp=(heaviside(z,0.0_dp)-signum(z)*p)/(g*xxi)
    value=m+s*((-signum(z)*qstd(tmp,0.0_dp,1.0_dp,nu)*xxi-c)/sc)
  end function qsstd

  function rsstd(n,rng,mu,sigma,skew,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,xi,nu
    integer::i
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=5.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rsstd_one(rng,m,s,xi,nu)
    end do
  end function rsstd

  real(dp) function rsstd_one(rng,mu,sigma,xi,nu) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,xi,nu
    real(dp)::weight,z,xx,rr,m1,c,sc
    weight=xi/(xi+1.0_dp/xi)
    z=uniform_rng(rng)-weight
    if(z<0.0_dp)then
    xx=1.0_dp/xi
    else
    xx=xi
    end if
    rr=-abs(rstd_one(rng,0.0_dp,1.0_dp,nu))/xx*merge(-1.0_dp,1.0_dp,z<0.0_dp)
    m1=2.0_dp*sqrt(nu-2.0_dp)/(nu-1.0_dp)/exp(log_gamma(0.5_dp)+log_gamma(nu/2.0_dp)-log_gamma(0.5_dp+nu/2.0_dp))
    call fs_standardization(m1,xi,c,sc)
    value=mu+sigma*(rr-c)/sc
  end function rsstd_one

  pure real(dp) function ged_lambda(nu) result(lambda)
    real(dp),intent(in)::nu
    lambda=sqrt((0.5_dp**(2.0_dp/nu))*gamma(1.0_dp/nu)/gamma(3.0_dp/nu))
  end function ged_lambda

  real(dp) function dged(x,mu,sigma,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,nu,lam,z,g
    logical::lg
    m=0.0_dp
    s=1.0_dp
    nu=2.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    if(present(log_density))lg=log_density
    if(s<=0.0_dp.or.nu<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    lam=ged_lambda(nu)
    z=(x-m)/s
    g=nu/(lam*2.0_dp**(1.0_dp+1.0_dp/nu)*gamma(1.0_dp/nu))
    value=g*exp(-0.5_dp*abs(z/lam)**nu)/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dged

  real(dp) function pged(q,mu,sigma,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp)::m,s,nu,lam,z,pg
    m=0.0_dp
    s=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    lam=ged_lambda(nu)
    z=(q-m)/s
    pg=regularized_gamma_p(1.0_dp/nu,0.5_dp*abs(z/lam)**nu)
    value=0.5_dp+0.5_dp*signum(z)*pg
    value=clamp_probability(value)
  end function pged

  real(dp) function qged(p,mu,sigma,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp)::m,s,nu,lam,y,q
    m=0.0_dp
    s=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    if(p<=0.0_dp)then
    value=-ieee_value(0.0_dp,ieee_positive_inf)
    return
    else if(p>=1.0_dp)then
    value=ieee_value(0.0_dp,ieee_positive_inf)
    return
    end if
    y=2.0_dp*p-1.0_dp
    lam=ged_lambda(nu)
    q=lam*(2.0_dp*gamma_quantile(abs(y),1.0_dp/nu))**(1.0_dp/nu)*signum(y)
    value=m+s*q
  end function qged

  function rged(n,rng,mu,sigma,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,nu
    integer::i
    m=0.0_dp
    s=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(shape))nu=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rged_one(rng,m,s,nu)
    end do
  end function rged

  real(dp) function rged_one(rng,mu,sigma,nu) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,nu
    real(dp)::sgn
    sgn=merge(-1.0_dp,1.0_dp,uniform_rng(rng)<0.5_dp)
    value=mu+sigma*ged_lambda(nu)*(2.0_dp*gamma_rng(1.0_dp/nu,1.0_dp,rng))**(1.0_dp/nu)*sgn
  end function rged_one

  pure real(dp) function ged_m1(nu) result(m1)
    real(dp),intent(in)::nu
    m1=2.0_dp**(1.0_dp/nu)*ged_lambda(nu)*gamma(2.0_dp/nu)/gamma(1.0_dp/nu)
  end function ged_m1

  real(dp) function dsged(x,mu,sigma,skew,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,xi,nu,c,sc,z,xxi,g
    logical::lg
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=2.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    if(present(log_density))lg=log_density
    call fs_standardization(ged_m1(nu),xi,c,sc)
    z=((x-m)/s)*sc+c
    if(abs(z)<=tiny(1.0_dp))then
    xxi=1.0_dp
    else if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=g*dged(z/xxi,0.0_dp,1.0_dp,nu)*sc/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dsged

  real(dp) function psged(q,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,xi,nu,c,sc,z,xxi,g
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    call fs_standardization(ged_m1(nu),xi,c,sc)
    z=((q-m)/s)*sc+c
    if(z<0.0_dp)then
    xxi=1.0_dp/xi
    else
    xxi=xi
    end if
    g=2.0_dp/(xi+1.0_dp/xi)
    value=heaviside(z,0.0_dp)-signum(z)*g*xxi*pged(-abs(z)/xxi,0.0_dp,1.0_dp,nu)
    value=clamp_probability(value)
  end function psged

  real(dp) function qsged(p,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,xi,nu,c,sc,g,z,xxi,tmp
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    call fs_standardization(ged_m1(nu),xi,c,sc)
    g=2.0_dp/(xi+1.0_dp/xi)
    z=p-1.0_dp/(1.0_dp+xi*xi)
    xxi=xi**signum(z)
    tmp=(heaviside(z,0.0_dp)-signum(z)*p)/(g*xxi)
    value=m+s*((-signum(z)*qged(tmp,0.0_dp,1.0_dp,nu)*xxi-c)/sc)
  end function qsged

  function rsged(n,rng,mu,sigma,skew,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,xi,nu
    integer::i
    m=0.0_dp
    s=1.0_dp
    xi=1.0_dp
    nu=2.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))xi=skew
    if(present(shape))nu=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rsged_one(rng,m,s,xi,nu)
    end do
  end function rsged

  real(dp) function rsged_one(rng,mu,sigma,xi,nu) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,xi,nu
    real(dp)::weight,z,xx,rr,c,sc
    weight=xi/(xi+1.0_dp/xi)
    z=uniform_rng(rng)-weight
    if(z<0.0_dp)then
    xx=1.0_dp/xi
    else
    xx=xi
    end if
    rr=-abs(rged_one(rng,0.0_dp,1.0_dp,nu))/xx*merge(-1.0_dp,1.0_dp,z<0.0_dp)
    call fs_standardization(ged_m1(nu),xi,c,sc)
    value=mu+sigma*(rr-c)/sc
  end function rsged_one

  subroutine jsu_constants(skew,shape,rtau,w,omega,c)
    real(dp),intent(in)::skew,shape
    real(dp),intent(out)::rtau,w,omega,c
    rtau=1.0_dp/shape
    if(rtau<1.0e-7_dp)then
    w=1.0_dp+rtau*rtau
    else
    w=exp(rtau*rtau)
    end if
    omega=-skew*rtau
    c=sqrt(1.0_dp/(0.5_dp*(w-1.0_dp)*(w*cosh(2.0_dp*omega)+1.0_dp)))
  end subroutine jsu_constants

  real(dp) function djsu(x,mu,sigma,skew,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,nu,tau,rtau,w,omega,c,z,r,logpdf
    logical::lg
    m=0.0_dp
    s=1.0_dp
    nu=0.0_dp
    tau=1.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))nu=skew
    if(present(shape))tau=shape
    if(present(log_density))lg=log_density
    call jsu_constants(nu,tau,rtau,w,omega,c)
    z=((x-m)/s-c*sqrt(w)*sinh(omega))/c
    r=-nu+asinh(z)/rtau
    logpdf=-log(c)-log(rtau)-0.5_dp*log(z*z+1.0_dp)-log(sqrt_two_pi)-0.5_dp*r*r-log(s)
    if(lg)then
    value=logpdf
    else
    value=exp(logpdf)
    end if
  end function djsu

  real(dp) function pjsu(q,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,nu,tau,rtau,w,omega,c,z,r
    m=0.0_dp
    s=1.0_dp
    nu=0.0_dp
    tau=1.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))nu=skew
    if(present(shape))tau=shape
    call jsu_constants(nu,tau,rtau,w,omega,c)
    z=(q-(m+c*s*sqrt(w)*sinh(omega)))/(c*s)
    r=-nu+asinh(z)/rtau
    value=normal_cdf(r)
  end function pjsu

  real(dp) function qjsu(p,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,nu,tau,rtau,w,omega,c,z
    m=0.0_dp
    s=1.0_dp
    nu=0.0_dp
    tau=1.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))nu=skew
    if(present(shape))tau=shape
    call jsu_constants(nu,tau,rtau,w,omega,c)
    z=sinh(rtau*(normal_quantile(p)+nu))
    value=m+s*(c*sqrt(w)*sinh(omega)+c*z)
  end function qjsu

  function rjsu(n,rng,mu,sigma,skew,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,nu,tau
    integer::i
    m=0.0_dp
    s=1.0_dp
    nu=0.0_dp
    tau=1.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))nu=skew
    if(present(shape))tau=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rjsu_one(rng,m,s,nu,tau)
    end do
  end function rjsu

  real(dp) function rjsu_one(rng,mu,sigma,skew,shape) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,skew,shape
    value=qjsu(uniform_rng(rng),mu,sigma,skew,shape)
  end function rjsu_one

  subroutine paramgh(rho,zeta,lambda,alpha,beta_value,delta,mu,status)
    real(dp),intent(in)::rho,zeta,lambda
    real(dp),intent(out)::alpha,beta_value,delta,mu
    integer,intent(out),optional::status
    real(dp)::rho2,zeta2,kappa,dkappa
    if(abs(rho)>=1.0_dp.or.zeta<=0.0_dp)then
      alpha=0.0_dp
      beta_value=0.0_dp
      delta=0.0_dp
      mu=0.0_dp
      if(present(status))status=1
      return
    end if
    rho2=1.0_dp-rho*rho
    zeta2=zeta*zeta
    kappa=exp(log_bessel_k(lambda+1.0_dp,zeta)-log_bessel_k(lambda,zeta))/zeta
    dkappa=exp(log_bessel_k(lambda+2.0_dp,zeta)-log_bessel_k(lambda+1.0_dp,zeta))/zeta-kappa
    alpha=sqrt(zeta2*kappa/rho2*(1.0_dp+rho*rho*zeta2*dkappa/rho2))
    beta_value=alpha*rho
    delta=zeta/(alpha*sqrt(rho2))
    mu=-beta_value*delta*delta*kappa
    if(present(status))status=0
  end subroutine paramgh

  subroutine paramghst(betabar,nu,beta_value,delta,mu,status)
    real(dp),intent(in)::betabar,nu
    real(dp),intent(out)::beta_value,delta,mu
    integer,intent(out),optional::status
    if(nu<=4.0_dp)then
    beta_value=0.0_dp
    delta=0.0_dp
    mu=0.0_dp
    if(present(status))status=1
    return
    end if
    delta=sqrt(1.0_dp/((2.0_dp*betabar*betabar)/((nu-2.0_dp)**2*(nu-4.0_dp))+1.0_dp/(nu-2.0_dp)))
    beta_value=betabar/delta
    mu=-beta_value*delta*delta/(nu-2.0_dp)
    if(present(status))status=0
  end subroutine paramghst

  function gh_model_raw(alpha,beta_value,delta,mu,lambda) result(model)
    real(dp),intent(in)::alpha,beta_value,delta,mu,lambda
    type(ghyp_model_type)::model
    real(dp)::bv(1),mv(1),dm(1,1)
    bv(1)=beta_value
    mv(1)=mu
    dm(1,1)=1.0_dp
    model=ghyp_ad(lambda,alpha,delta,bv,mv,dm)
  end function gh_model_raw

  real(dp) function dghyp_raw(x,alpha,beta_value,delta,mu,lambda,log_density) result(value)
    real(dp),intent(in)::x,alpha,beta_value,delta,mu,lambda
    logical,intent(in),optional::log_density
    type(ghyp_model_type)::model
    logical::lg
    lg=.false.
    if(present(log_density))lg=log_density
    model=gh_model_raw(alpha,beta_value,delta,mu,lambda)
    value=ghyp_density(x,model)
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dghyp_raw

  real(dp) function pghyp_raw(q,alpha,beta_value,delta,mu,lambda) result(value)
    real(dp),intent(in)::q,alpha,beta_value,delta,mu,lambda
    type(ghyp_model_type)::model
    model=gh_model_raw(alpha,beta_value,delta,mu,lambda)
    value=ghyp_cdf(q,model)
  end function pghyp_raw

  real(dp) function qghyp_raw(p,alpha,beta_value,delta,mu,lambda) result(value)
    real(dp),intent(in)::p,alpha,beta_value,delta,mu,lambda
    type(ghyp_model_type)::model
    model=gh_model_raw(alpha,beta_value,delta,mu,lambda)
    value=ghyp_quantile(p,model)
  end function qghyp_raw

  function rghyp_raw(n,rng,alpha,beta_value,delta,mu,lambda) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::alpha,beta_value,delta,mu,lambda
    real(dp),allocatable::x(:),z(:)
    type(ghyp_model_type)::model
    integer::i
    allocate(x(max(n,0)))
    model=gh_model_raw(alpha,beta_value,delta,mu,lambda)
    do i=1,n
    z=ghyp_random_one(model,rng)
    x(i)=z(1)
    end do
  end function rghyp_raw

  subroutine standardized_gh_model(skew,shape,lambda,model,status)
    real(dp),intent(in)::skew,shape,lambda
    type(ghyp_model_type),intent(out)::model
    integer,intent(out)::status
    real(dp)::a,b,d,m
    call paramgh(skew,shape,lambda,a,b,d,m,status)
    if(status==0)model=gh_model_raw(a,b,d,m,lambda)
  end subroutine standardized_gh_model

  real(dp) function dgh(x,mu,sigma,skew,shape,lambda,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape,lambda
    logical,intent(in),optional::log_density
    real(dp)::m,s,rho,zeta,lam
    type(ghyp_model_type)::model
    integer::status
    logical::lg
    m=0.0_dp
    s=1.0_dp
    rho=0.0_dp
    zeta=3.0_dp
    lam=-0.5_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))rho=skew
    if(present(shape))zeta=shape
    if(present(lambda))lam=lambda
    if(present(log_density))lg=log_density
    call standardized_gh_model(rho,zeta,lam,model,status)
    if(status/=0.or.s<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    value=ghyp_density((x-m)/s,model)/s
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dgh

  real(dp) function pgh(q,mu,sigma,skew,shape,lambda) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape,lambda
    real(dp)::m,s,rho,zeta,lam
    type(ghyp_model_type)::model
    integer::status
    m=0.0_dp
    s=1.0_dp
    rho=0.0_dp
    zeta=3.0_dp
    lam=-0.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))rho=skew
    if(present(shape))zeta=shape
    if(present(lambda))lam=lambda
    call standardized_gh_model(rho,zeta,lam,model,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    value=ghyp_cdf((q-m)/s,model)
    end if
  end function pgh

  real(dp) function qgh(p,mu,sigma,skew,shape,lambda) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape,lambda
    real(dp)::m,s,rho,zeta,lam
    type(ghyp_model_type)::model
    integer::status
    m=0.0_dp
    s=1.0_dp
    rho=0.0_dp
    zeta=3.0_dp
    lam=-0.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))rho=skew
    if(present(shape))zeta=shape
    if(present(lambda))lam=lambda
    call standardized_gh_model(rho,zeta,lam,model,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    value=m+s*ghyp_quantile(p,model)
    end if
  end function qgh

  function rgh(n,rng,mu,sigma,skew,shape,lambda) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape,lambda
    real(dp),allocatable::x(:),z(:)
    real(dp)::m,s,rho,zeta,lam
    type(ghyp_model_type)::model
    integer::status,i
    m=0.0_dp
    s=1.0_dp
    rho=0.0_dp
    zeta=3.0_dp
    lam=-0.5_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))rho=skew
    if(present(shape))zeta=shape
    if(present(lambda))lam=lambda
    allocate(x(max(n,0)))
    call standardized_gh_model(rho,zeta,lam,model,status)
    do i=1,n
    z=ghyp_random_one(model,rng)
    x(i)=m+s*z(1)
    end do
  end function rgh

  real(dp) function rgh_one(rng,mu,sigma,skew,shape,lambda) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,skew,shape,lambda
    real(dp),allocatable::z(:)
    type(ghyp_model_type)::model
    integer::status
    call standardized_gh_model(skew,shape,lambda,model,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    z=ghyp_random_one(model,rng)
    value=mu+sigma*z(1)
    end if
  end function rgh_one

  real(dp) function dnig(x,mu,sigma,skew,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,r,z
    logical::lg
    m=0.0_dp
    s=1.0_dp
    r=0.0_dp
    z=3.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))r=skew
    if(present(shape))z=shape
    if(present(log_density))lg=log_density
    value=dgh(x,m,s,r,z,-0.5_dp,lg)
  end function dnig
  real(dp) function pnig(q,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,r,z
    m=0.0_dp
    s=1.0_dp
    r=0.0_dp
    z=3.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))r=skew
    if(present(shape))z=shape
    value=pgh(q,m,s,r,z,-0.5_dp)
  end function pnig
  real(dp) function qnig(p,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,r,z
    m=0.0_dp
    s=1.0_dp
    r=0.0_dp
    z=3.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))r=skew
    if(present(shape))z=shape
    value=qgh(p,m,s,r,z,-0.5_dp)
  end function qnig
  function rnig(n,rng,mu,sigma,skew,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,r,z
    m=0.0_dp
    s=1.0_dp
    r=0.0_dp
    z=3.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))r=skew
    if(present(shape))z=shape
    x=rgh(n,rng,m,s,r,z,-0.5_dp)
  end function rnig
  real(dp) function rnig_one(rng,mu,sigma,skew,shape) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,skew,shape
    value=rgh_one(rng,mu,sigma,skew,shape,-0.5_dp)
  end function rnig_one

  subroutine ghst_model(skew,shape,model,status)
    real(dp),intent(in)::skew,shape
    type(ghyp_model_type),intent(out)::model
    integer,intent(out)::status
    real(dp)::b,d,m,a
    call paramghst(skew,shape,b,d,m,status)
    if(status/=0)return
    a=abs(b)+1.0e-10_dp
    model=gh_model_raw(a,b,d,m,-0.5_dp*shape)
  end subroutine ghst_model

  real(dp) function dghst(x,mu,sigma,skew,shape,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma,skew,shape
    logical,intent(in),optional::log_density
    real(dp)::m,s,bbar,nu,b,d,m0,res,arg,logpdf
    integer::status
    logical::lg
    m=0.0_dp
    s=1.0_dp
    bbar=1.0_dp
    nu=8.0_dp
    lg=.false.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))bbar=skew
    if(present(shape))nu=shape
    if(present(log_density))lg=log_density
    if(abs(bbar)<1.0e-12_dp)then
      if(bbar<0.0_dp)then
        bbar=-1.0e-12_dp
      else
        bbar=1.0e-12_dp
      end if
    end if
    call paramghst(bbar,nu,b,d,m0,status)
    if(status/=0.or.s<=0.0_dp)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    res=(x-m)/s-m0
    arg=sqrt(b*b*(d*d+res*res))
    logpdf=0.5_dp*(1.0_dp-nu)*log(2.0_dp)+nu*log(d)+0.5_dp*(nu+1.0_dp)*log(abs(b))+ &
      log_bessel_k(0.5_dp*(nu+1.0_dp),arg)+b*res-log_gamma(0.5_dp*nu)-0.5_dp*log(pi)- &
      0.25_dp*(nu+1.0_dp)*log(d*d+res*res)-log(s)
    if(lg)then
    value=logpdf
    else
    value=exp(logpdf)
    end if
  end function dghst

  real(dp) function pghst(q,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,bbar,nu
    type(ghyp_model_type)::model
    integer::status
    m=0.0_dp
    s=1.0_dp
    bbar=1.0_dp
    nu=8.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))bbar=skew
    if(present(shape))nu=shape
    call ghst_model(bbar,nu,model,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    value=ghyp_cdf((q-m)/s,model)
    end if
  end function pghst

  real(dp) function qghst(p,mu,sigma,skew,shape) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp)::m,s,bbar,nu
    type(ghyp_model_type)::model
    integer::status
    m=0.0_dp
    s=1.0_dp
    bbar=1.0_dp
    nu=8.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))bbar=skew
    if(present(shape))nu=shape
    call ghst_model(bbar,nu,model,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    else
    value=m+s*ghyp_quantile(p,model)
    end if
  end function qghst

  function rghst(n,rng,mu,sigma,skew,shape) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::mu,sigma,skew,shape
    real(dp),allocatable::x(:)
    real(dp)::m,s,bbar,nu
    integer::i
    m=0.0_dp
    s=1.0_dp
    bbar=1.0_dp
    nu=8.0_dp
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(skew))bbar=skew
    if(present(shape))nu=shape
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=rghst_one(rng,m,s,bbar,nu)
    end do
  end function rghst

  real(dp) function rghst_one(rng,mu,sigma,skew,shape) result(value)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::mu,sigma,skew,shape
    real(dp)::b,d,m0,w
    integer::status
    call paramghst(skew,shape,b,d,m0,status)
    if(status/=0)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    w=1.0_dp/gamma_rng(0.5_dp*shape,2.0_dp/(d*d),rng)
    value=mu+sigma*(m0+b*w+sqrt(w)*normal_rng(rng))
  end function rghst_one

  real(dp) function dghyp(x,alpha,beta_value,delta,mu,lambda,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::alpha,beta_value,delta,mu,lambda
    logical,intent(in),optional::log_density
    real(dp)::a,b,d,m,l
    a=1.0_dp
    b=0.0_dp
    d=1.0_dp
    m=0.0_dp
    l=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta_value))b=beta_value
    if(present(delta))d=delta
    if(present(mu))m=mu
    if(present(lambda))l=lambda
    value=dghyp_raw(x,a,b,d,m,l,log_density)
  end function dghyp

  real(dp) function pghyp(q,alpha,beta_value,delta,mu,lambda,lower_tail,log_probability) result(value)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::alpha,beta_value,delta,mu,lambda
    logical,intent(in),optional::lower_tail,log_probability
    real(dp)::a,b,d,m,l
    logical::lower,logp
    a=1.0_dp
    b=0.0_dp
    d=1.0_dp
    m=0.0_dp
    l=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta_value))b=beta_value
    if(present(delta))d=delta
    if(present(mu))m=mu
    if(present(lambda))l=lambda
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    logp=.false.
    if(present(log_probability))logp=log_probability
    value=pghyp_raw(q,a,b,d,m,l)
    if(.not.lower)value=1.0_dp-value
    if(logp)value=log(max(value,tiny(1.0_dp)))
  end function pghyp

  real(dp) function qghyp(p,alpha,beta_value,delta,mu,lambda,lower_tail,log_probability) result(value)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::alpha,beta_value,delta,mu,lambda
    logical,intent(in),optional::lower_tail,log_probability
    real(dp)::a,b,d,m,l,pr
    logical::lower,logp
    a=1.0_dp
    b=0.0_dp
    d=1.0_dp
    m=0.0_dp
    l=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta_value))b=beta_value
    if(present(delta))d=delta
    if(present(mu))m=mu
    if(present(lambda))l=lambda
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    logp=.false.
    if(present(log_probability))logp=log_probability
    pr=p
    if(logp)pr=exp(pr)
    if(.not.lower)pr=1.0_dp-pr
    value=qghyp_raw(pr,a,b,d,m,l)
  end function qghyp

  function rghyp(n,rng,alpha,beta_value,delta,mu,lambda) result(x)
    integer,intent(in)::n
    type(rng_state),intent(inout)::rng
    real(dp),intent(in),optional::alpha,beta_value,delta,mu,lambda
    real(dp),allocatable::x(:)
    real(dp)::a,b,d,m,l
    a=1.0_dp
    b=0.0_dp
    d=1.0_dp
    m=0.0_dp
    l=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta_value))b=beta_value
    if(present(delta))d=delta
    if(present(mu))m=mu
    if(present(lambda))l=lambda
    x=rghyp_raw(n,rng,a,b,d,m,l)
  end function rghyp

  subroutine nigtransform(mu,sigma,skew,shape,out,status)
    real(dp),intent(in)::mu,sigma,skew,shape
    real(dp),intent(out)::out(4)
    integer,intent(out),optional::status
    real(dp)::a,b,d,m0
    integer::st
    call paramgh(skew,shape,-0.5_dp,a,b,d,m0,st)
    out=[m0*sigma+mu,d*sigma,b/sigma,a/sigma]
    if(present(status))status=st
  end subroutine nigtransform

  subroutine ghyptransform(mu,sigma,skew,shape,lambda,out,status)
    real(dp),intent(in)::mu,sigma,skew,shape,lambda
    real(dp),intent(out)::out(4)
    integer,intent(out),optional::status
    real(dp)::a,b,d,m0
    integer::st
    call paramgh(skew,shape,lambda,a,b,d,m0,st)
    out=[m0*sigma+mu,d*sigma,b/sigma,a/sigma]
    if(present(status))status=st
  end subroutine ghyptransform

end module tsd_distributions
