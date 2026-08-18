! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_distributions
  use runuran_kinds, only : dp, pi, sqrt2pi, huge_dp, clamp
  use runuran_math, only : normal_pdf, normal_cdf, normal_quantile, reg_beta, reg_gamma_p, &
    gamma_cdf, gamma_quantile, beta_cdf, beta_quantile, chisq_cdf, chisq_quantile, &
    student_t_cdf, student_t_quantile, f_cdf, f_quantile, adaptive_integral, &
    log_bessel_k_nu, log_gamma_abs_complex, log1p_safe, expm1_safe
  use runuran_rng, only : rng_state, rng_uniform, rng_normal, rng_exponential, rng_gamma, rng_beta, rng_poisson, &
    rng_binomial, rng_geometric, rng_negative_binomial, rng_cauchy, rng_chisq, rng_student_t, &
    rng_f, rng_inverse_gaussian
  implicit none
  private

  integer, parameter :: D_CUSTOM=0,D_NORMAL=1,D_BETA=2,D_CAUCHY=3,D_CHI=4,D_CHISQ=5,D_EXP=6,D_F=7
  integer, parameter :: D_FRECHET=8,D_GAMMA=9,D_GUMBEL=10,D_IG=11,D_LAPLACE=12,D_LNORM=13,D_LOGIS=14
  integer, parameter :: D_LOMAX=15,D_PARETO=16,D_POWEREXP=17,D_RAYLEIGH=18,D_SLASH=19,D_T=20,D_TRIANG=21
  integer, parameter :: D_WEIBULL=22,D_BURR=23,D_GIG=24,D_GIGA=25,D_HYPERB=26,D_GHYP=27,D_VG=28,D_MEIX=29
  integer, parameter :: D_PLANCK=30
  integer, parameter :: Q_CUSTOM=0,Q_BINOM=1,Q_GEOM=2,Q_HYPER=3,Q_LOGARITHMIC=4,Q_NBINOM=5,Q_POIS=6,Q_ZIPF=7,Q_TABLE=8

  abstract interface
    function cont_callback(x,params) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: params(:)
      real(dp) :: y
    end function cont_callback
    function discr_callback(k,params) result(y)
      import dp
      integer, intent(in) :: k
      real(dp), intent(in) :: params(:)
      real(dp) :: y
    end function discr_callback
  end interface

  type, public :: continuous_distribution
    integer :: id=D_CUSTOM
    character(len=32) :: name='custom'
    real(dp) :: params(8)=0.0_dp
    integer :: nparams=0
    real(dp) :: support_lb=-huge_dp, support_ub=huge_dp
    real(dp) :: lb=-huge_dp, ub=huge_dp
    real(dp) :: center=0.0_dp
    procedure(cont_callback), pointer, nopass :: pdf_cb=>null(), cdf_cb=>null(), logpdf_cb=>null(), dpdf_cb=>null()
  contains
    procedure :: pdf=>cont_pdf
    procedure :: logpdf=>cont_logpdf
    procedure :: dpdf=>cont_dpdf
    procedure :: dlogpdf=>cont_dlogpdf
    procedure :: cdf=>cont_cdf
    procedure :: quantile=>cont_quantile
    procedure :: sample=>cont_sample
    procedure :: sample_n=>cont_sample_n
  end type continuous_distribution

  type, public :: discrete_distribution
    integer :: id=Q_CUSTOM
    character(len=32) :: name='custom'
    real(dp) :: params(8)=0.0_dp
    integer :: nparams=0
    integer :: support_lb=0, support_ub=huge(1)
    integer :: lb=0, ub=huge(1)
    real(dp), allocatable :: probvec(:)
    procedure(discr_callback), pointer, nopass :: pmf_cb=>null(), cdf_cb=>null()
  contains
    procedure :: pmf=>discr_pmf
    procedure :: cdf=>discr_cdf
    procedure :: quantile=>discr_quantile
    procedure :: sample=>discr_sample
    procedure :: sample_n=>discr_sample_n
  end type discrete_distribution

  public :: ud_continuous, ud_continuous_cdf, ud_continuous_logpdf
  public :: ud_discrete, ud_probability_vector
  public :: udnorm, udbeta, udcauchy, udchi, udchisq, udexp, udf, udfrechet, udgamma
  public :: udgumbel, udig, udlaplace, udlnorm, udlogis, udlomax, udpareto, udpowerexp
  public :: udrayleigh, udslash, udt, udweibull, udburr, udgig, udgiga, udhyperbolic
  public :: udghyp, udvg, udmeixner, udplanck, udtriang
  public :: udbinom, udgeom, udhyper, udlogarithmic, udnbinom, udpois, udzipf

contains
  function ud_continuous(pdf,cdf,lb,ub,params,logpdf,dpdf,name) result(d)
    procedure(cont_callback) :: pdf
    procedure(cont_callback), optional :: cdf,logpdf,dpdf
    real(dp), intent(in) :: lb,ub
    real(dp), intent(in), optional :: params(:)
    character(len=*), intent(in), optional :: name
    type(continuous_distribution) :: d
    d%id=D_CUSTOM
    d%pdf_cb=>pdf
    d%lb=lb
    d%ub=ub
    d%support_lb=lb
    d%support_ub=ub
    if(present(cdf)) d%cdf_cb=>cdf
    if(present(logpdf)) d%logpdf_cb=>logpdf
    if(present(dpdf)) d%dpdf_cb=>dpdf
    if(present(params))then
    d%nparams=min(size(params),8)
    d%params(1:d%nparams)=params(1:d%nparams)
    end if
    if(present(name)) d%name=name
  end function ud_continuous

  function ud_continuous_cdf(cdf,lb,ub,params,name) result(d)
    procedure(cont_callback) :: cdf
    real(dp), intent(in) :: lb,ub
    real(dp), intent(in), optional :: params(:)
    character(len=*), intent(in), optional :: name
    type(continuous_distribution) :: d
    d%id=D_CUSTOM
    d%cdf_cb=>cdf
    d%lb=lb
    d%ub=ub
    d%support_lb=lb
    d%support_ub=ub
    if(present(params))then
      d%nparams=min(size(params),8)
      d%params(1:d%nparams)=params(1:d%nparams)
    end if
    if(present(name))d%name=name
  end function ud_continuous_cdf

  function ud_continuous_logpdf(logpdf,lb,ub,params,name) result(d)
    procedure(cont_callback) :: logpdf
    real(dp), intent(in) :: lb,ub
    real(dp), intent(in), optional :: params(:)
    character(len=*), intent(in), optional :: name
    type(continuous_distribution) :: d
    d%id=D_CUSTOM
    d%logpdf_cb=>logpdf
    d%lb=lb
    d%ub=ub
    d%support_lb=lb
    d%support_ub=ub
    if(present(params))then
      d%nparams=min(size(params),8)
      d%params(1:d%nparams)=params(1:d%nparams)
    end if
    if(present(name))d%name=name
  end function ud_continuous_logpdf

  function ud_discrete(pmf,cdf,lb,ub,params,name) result(d)
    procedure(discr_callback) :: pmf
    procedure(discr_callback), optional :: cdf
    integer,intent(in)::lb,ub
    real(dp),intent(in),optional::params(:)
    character(len=*),intent(in),optional::name
    type(discrete_distribution)::d
    d%id=Q_CUSTOM
    d%pmf_cb=>pmf
    d%lb=lb
    d%ub=ub
    d%support_lb=lb
    d%support_ub=ub
    if(present(cdf))d%cdf_cb=>cdf
    if(present(params))then
    d%nparams=min(size(params),8)
    d%params(1:d%nparams)=params(1:d%nparams)
    end if
    if(present(name))d%name=name
  end function ud_discrete

  function ud_probability_vector(pv,lb) result(d)
    real(dp), intent(in) :: pv(:)
    integer, intent(in), optional :: lb
    type(discrete_distribution) :: d
    integer :: l
    real(dp) :: z
    l=0
    if(present(lb))l=lb
    if(any(pv<0.0_dp))error stop 'ud_probability_vector: negative probability'
    z=sum(pv)
    if(z<=0.0_dp)error stop 'ud_probability_vector: zero probability sum'
    d%id=Q_TABLE
    d%name='probability_vector'
    d%support_lb=l
    d%support_ub=l+size(pv)-1
    d%lb=d%support_lb
    d%ub=d%support_ub
    allocate(d%probvec(size(pv)))
    d%probvec=pv/z
  end function ud_probability_vector

  function make_cont(id,name,p,sl,su,lb,ub,center) result(d)
    integer,intent(in)::id
    character(len=*),intent(in)::name
    real(dp),intent(in)::p(:),sl,su
    real(dp),intent(in),optional::lb,ub,center
    type(continuous_distribution)::d
    d%id=id
    d%name=name
    d%nparams=min(size(p),8)
    d%params(1:d%nparams)=p(1:d%nparams)
    d%support_lb=sl
    d%support_ub=su
    d%lb=sl
    d%ub=su
    if(present(lb))d%lb=max(sl,lb)
    if(present(ub))d%ub=min(su,ub)
    if(present(center))d%center=center
  end function make_cont

  function udnorm(mean,sd,lb,ub) result(d)
    real(dp),intent(in),optional::mean,sd,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=-huge_dp
    u=huge_dp
    if(present(mean))m=mean
    if(present(sd))s=sd
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_NORMAL,'normal',[m,s],-huge_dp,huge_dp,l,u,m)
  end function
  function udbeta(shape1,shape2,lb,ub) result(d)
    real(dp),intent(in)::shape1,shape2
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=1.0_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_BETA,'beta',[shape1,shape2],0.0_dp,1.0_dp,l,u,shape1/(shape1+shape2))
  end function
  function udcauchy(location,scale,lb,ub) result(d)
    real(dp),intent(in),optional::location,scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=-huge_dp
    u=huge_dp
    if(present(location))m=location
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_CAUCHY,'cauchy',[m,s],-huge_dp,huge_dp,l,u,m)
  end function
  function udchi(df,lb,ub) result(d)
    real(dp),intent(in)::df
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_CHI,'chi',[df],0.0_dp,huge_dp,l,u,sqrt(max(df-1.0_dp,0.0_dp)))
  end function
  function udchisq(df,lb,ub) result(d)
    real(dp),intent(in)::df
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_CHISQ,'chisquare',[df],0.0_dp,huge_dp,l,u,max(df-2.0_dp,0.0_dp))
  end function
  function udexp(rate,lb,ub) result(d)
    real(dp),intent(in),optional::rate,lb,ub
    type(continuous_distribution)::d
    real(dp)::r,l,u
    r=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(rate))r=rate
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_EXP,'exponential',[r],0.0_dp,huge_dp,l,u,0.0_dp)
  end function
  function udf(df1,df2,lb,ub) result(d)
    real(dp),intent(in)::df1,df2
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_F,'F',[df1,df2],0.0_dp,huge_dp,l,u,1.0_dp)
  end function
  function udfrechet(shape,location,scale,lb,ub) result(d)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::location,scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    if(present(location))m=location
    if(present(scale))s=scale
    l=m
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_FRECHET,'frechet',[shape,m,s],m,huge_dp,l,u,m+s*(shape/(shape+1.0_dp))**(1.0_dp/shape))
  end function
  function udgamma(shape,scale,lb,ub) result(d)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::s,l,u
    s=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_GAMMA,'gamma',[shape,s],0.0_dp,huge_dp,l,u,max(shape-1.0_dp,0.0_dp)*s)
  end function
  function udgumbel(location,scale,lb,ub) result(d)
    real(dp),intent(in),optional::location,scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=-huge_dp
    u=huge_dp
    if(present(location))m=location
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_GUMBEL,'gumbel',[m,s],-huge_dp,huge_dp,l,u,m)
  end function
  function udig(mu,lambda,lb,ub) result(d)
    real(dp),intent(in)::mu,lambda
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_IG,'inverse_gaussian',[mu,lambda],0.0_dp,huge_dp,l,u, &
      mu*(sqrt(1.0_dp+9.0_dp*mu*mu/(4.0_dp*lambda*lambda))-3.0_dp*mu/(2.0_dp*lambda)))
  end function
  function udlaplace(location,scale,lb,ub) result(d)
    real(dp),intent(in),optional::location,scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=-huge_dp
    u=huge_dp
    if(present(location))m=location
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_LAPLACE,'laplace',[m,s],-huge_dp,huge_dp,l,u,m)
  end function
  function udlnorm(meanlog,sdlog,lb,ub) result(d)
    real(dp),intent(in),optional::meanlog,sdlog,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(meanlog))m=meanlog
    if(present(sdlog))s=sdlog
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_LNORM,'lognormal',[m,s],0.0_dp,huge_dp,l,u,exp(m-s*s))
  end function
  function udlogis(location,scale,lb,ub) result(d)
    real(dp),intent(in),optional::location,scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::m,s,l,u
    m=0.0_dp
    s=1.0_dp
    l=-huge_dp
    u=huge_dp
    if(present(location))m=location
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_LOGIS,'logistic',[m,s],-huge_dp,huge_dp,l,u,m)
  end function
  function udlomax(shape,scale,lb,ub) result(d)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::s,l,u
    s=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_LOMAX,'lomax',[shape,s],0.0_dp,huge_dp,l,u,0.0_dp)
  end function
  function udpareto(k,a,lb,ub) result(d)
    real(dp),intent(in)::k,a
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=k
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_PARETO,'pareto',[k,a],k,huge_dp,l,u,k)
  end function
  function udpowerexp(shape,lb,ub) result(d)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_POWEREXP,'powerexponential',[shape],-huge_dp,huge_dp,l,u,0.0_dp)
  end function
  function udrayleigh(scale,lb,ub) result(d)
    real(dp),intent(in),optional::scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::s,l,u
    s=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_RAYLEIGH,'rayleigh',[s],0.0_dp,huge_dp,l,u,s)
  end function
  function udslash(lb,ub) result(d)
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_SLASH,'slash',[0.0_dp],-huge_dp,huge_dp,l,u,0.0_dp)
  end function
  function udt(df,lb,ub) result(d)
    real(dp),intent(in)::df
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_T,'student_t',[df],-huge_dp,huge_dp,l,u,0.0_dp)
  end function
  function udtriang(a,m,b,lb,ub) result(d)
    real(dp),intent(in)::a,m,b
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=a
    u=b
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_TRIANG,'triangular',[a,m,b],a,b,l,u,m)
  end function
  function udweibull(shape,scale,lb,ub) result(d)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::scale,lb,ub
    type(continuous_distribution)::d
    real(dp)::s,l,u,m
    s=1.0_dp
    l=0.0_dp
    u=huge_dp
    if(present(scale))s=scale
    if(present(lb))l=lb
    if(present(ub))u=ub
    m=0.0_dp
    if(shape>1.0_dp)m=s*((shape-1.0_dp)/shape)**(1.0_dp/shape)
    d=make_cont(D_WEIBULL,'weibull',[shape,s],0.0_dp,huge_dp,l,u,m)
  end function
  function udburr(a,b,lb,ub) result(d)
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_BURR,'burr',[a,b],0.0_dp,huge_dp,l,u,0.0_dp)
  end function
  function udgig(theta,psi,chi,lb,ub) result(d)
    real(dp),intent(in)::theta,psi,chi
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u,m
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    m=((theta-1.0_dp)+sqrt((theta-1.0_dp)**2+psi*chi))/psi
    d=make_cont(D_GIG,'gig2',[theta,psi,chi],0.0_dp,huge_dp,l,u,m)
  end function
  function udgiga(theta,omega,eta,lb,ub) result(d)
    real(dp),intent(in)::theta,omega
    real(dp),intent(in),optional::eta,lb,ub
    type(continuous_distribution)::d
    real(dp)::e,l,u,m
    e=1.0_dp
    if(present(eta))e=eta
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    m=e*(sqrt(omega*omega+(theta-1.0_dp)**2)+(theta-1.0_dp))/omega
    d=make_cont(D_GIGA,'gig',[theta,omega,e],0.0_dp,huge_dp,l,u,m)
  end function
  function udhyperbolic(alpha,beta,delta,mu,lb,ub) result(d)
    real(dp),intent(in)::alpha,beta,delta,mu
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u,m
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    m=mu+delta*beta/sqrt(alpha*alpha-beta*beta)
    d=make_cont(D_HYPERB,'hyperbolic',[alpha,beta,delta,mu],-huge_dp,huge_dp,l,u,m)
  end function
  function udghyp(lambda,alpha,beta,delta,mu,lb,ub) result(d)
    real(dp),intent(in)::lambda,alpha,beta,delta,mu
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_GHYP,'ghyp',[lambda,alpha,beta,delta,mu],-huge_dp,huge_dp,l,u,mu)
  end function
  function udvg(lambda,alpha,beta,mu,lb,ub) result(d)
    real(dp),intent(in)::lambda,alpha,beta,mu
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u,c,g
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    g=sqrt(alpha*alpha-beta*beta)
    c=mu+2.0_dp*beta*lambda/(g*g)
    d=make_cont(D_VG,'variance_gamma',[lambda,alpha,beta,mu],-huge_dp,huge_dp,l,u,c)
  end function
  function udmeixner(alpha,beta,delta,mu,lb,ub) result(d)
    real(dp),intent(in)::alpha,beta,delta,mu
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=-huge_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_MEIX,'meixner',[alpha,beta,delta,mu],-huge_dp,huge_dp,l,u,mu)
  end function
  function udplanck(a,lb,ub) result(d)
    real(dp),intent(in)::a
    real(dp),intent(in),optional::lb,ub
    type(continuous_distribution)::d
    real(dp)::l,u
    l=0.0_dp
    u=huge_dp
    if(present(lb))l=lb
    if(present(ub))u=ub
    d=make_cont(D_PLANCK,'planck',[a],0.0_dp,huge_dp,l,u,max(a,1.0_dp))
  end function

  pure real(dp) function zeta_s(s) result(z)
    real(dp),intent(in)::s
    integer::k
    real(dp)::term
    z=0.0_dp
    do k=1,200000
      term=real(k,dp)**(-s)
      z=z+term
      if(term<1.0e-15_dp*max(z,1.0_dp))exit
    end do
    ! Euler tail correction.
    z=z+real(k,dp)**(1.0_dp-s)/(s-1.0_dp)-0.5_dp*real(k,dp)**(-s)
  end function zeta_s

  real(dp) function base_cont_logpdf(d,x) result(lp)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    real(dp)::p1,p2,p3,p4,p5,z,y,nu,gam,logc
    if(x<d%support_lb .or. x>d%support_ub)then
    lp=-huge_dp
    return
    end if
    p1=d%params(1)
    p2=d%params(2)
    p3=d%params(3)
    p4=d%params(4)
    p5=d%params(5)
    select case(d%id)
    case(D_CUSTOM)
      if(associated(d%logpdf_cb))then
      lp=d%logpdf_cb(x,d%params(1:d%nparams))
      return
      end if
      if(associated(d%pdf_cb))then
      lp=log(max(d%pdf_cb(x,d%params(1:d%nparams)),tiny(1.0_dp)))
      return
      end if
      if(associated(d%cdf_cb))then
      y=sqrt(epsilon(1.0_dp))*(1.0_dp+abs(x))
      z=max(d%support_lb,x-y)
      nu=min(d%support_ub,x+y)
      if(nu>z)then
      lp=log(max((d%cdf_cb(nu,d%params(1:d%nparams))- &
        d%cdf_cb(z,d%params(1:d%nparams)))/(nu-z),tiny(1.0_dp)))
      else
      lp=-huge_dp
      end if
      return
      end if
      lp=-huge_dp
    case(D_NORMAL)
    lp=-log(p2*sqrt2pi)-0.5_dp*((x-p1)/p2)**2
    case(D_BETA)
      if(x<=0.0_dp.or.x>=1.0_dp)then
      lp=-huge_dp
      else
      lp=(p1-1.0_dp)*log(x)+(p2-1.0_dp)*log1p_safe(-x)-log_gamma(p1)-log_gamma(p2)+log_gamma(p1+p2)
      end if
    case(D_CAUCHY)
    lp=-log(pi*p2)-log(1.0_dp+((x-p1)/p2)**2)
    case(D_CHI)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=(p1-1.0_dp)*log(x)-0.5_dp*x*x-(0.5_dp*p1-1.0_dp)*log(2.0_dp)-log_gamma(0.5_dp*p1)
      end if
    case(D_CHISQ)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=(0.5_dp*p1-1.0_dp)*log(x)-0.5_dp*x-0.5_dp*p1*log(2.0_dp)-log_gamma(0.5_dp*p1)
      end if
    case(D_EXP)
    if(x<0.0_dp)then
    lp=-huge_dp
    else
    lp=log(p1)-p1*x
    end if
    case(D_F)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=0.5_dp*p1*log(p1/p2)+(0.5_dp*p1-1.0_dp)*log(x) &
        -0.5_dp*(p1+p2)*log1p_safe(p1*x/p2)-log_gamma(0.5_dp*p1) &
        -log_gamma(0.5_dp*p2)+log_gamma(0.5_dp*(p1+p2))
      end if
    case(D_FRECHET)
      z=(x-p2)/p3
      if(z<=0.0_dp)then
      lp=-huge_dp
      else
      lp=log(p1/p3)+(-1.0_dp-p1)*log(z)-z**(-p1)
      end if
    case(D_GAMMA)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=(p1-1.0_dp)*log(x)-x/p2-log_gamma(p1)-p1*log(p2)
      end if
    case(D_GUMBEL)
    z=(x-p1)/p2
    lp=-log(p2)-z-exp(-z)
    case(D_IG)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=0.5_dp*log(p2/(2.0_dp*pi*x**3))-p2*(x-p1)**2/(2.0_dp*p1*p1*x)
      end if
    case(D_LAPLACE)
    lp=-log(2.0_dp*p2)-abs(x-p1)/p2
    case(D_LNORM)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=-log(x*p2*sqrt2pi)-0.5_dp*((log(x)-p1)/p2)**2
      end if
    case(D_LOGIS)
    z=(x-p1)/p2
    lp=-log(p2)-z-2.0_dp*log1p_safe(exp(-z))
    case(D_LOMAX)
      if(x<0.0_dp)then
      lp=-huge_dp
      else
      lp=log(p1/p2)-(p1+1.0_dp)*log1p_safe(x/p2)
      end if
    case(D_PARETO)
      if(x<p1)then
      lp=-huge_dp
      else
      lp=log(p2)+p2*log(p1)-(p2+1.0_dp)*log(x)
      end if
    case(D_POWEREXP)
    lp=-abs(x)**p1-log(2.0_dp)-log_gamma(1.0_dp+1.0_dp/p1)
    case(D_RAYLEIGH)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=log(x)-2.0_dp*log(p1)-0.5_dp*(x/p1)**2
      end if
    case(D_SLASH)
      if(abs(x)<1.0e-8_dp)then
      lp=log(0.5_dp/sqrt2pi)
      else
      lp=log((1.0_dp-exp(-0.5_dp*x*x))/(x*x*sqrt2pi))
      end if
    case(D_T)
    lp=log_gamma(0.5_dp*(p1+1.0_dp))-log_gamma(0.5_dp*p1)-0.5_dp*log(p1*pi)-0.5_dp*(p1+1.0_dp)*log1p_safe(x*x/p1)
    case(D_TRIANG)
      if(x<p1.or.x>p3)then
      lp=-huge_dp
      else if(x<=p2)then
      lp=log(2.0_dp*(x-p1)/((p3-p1)*(p2-p1)))
      else
      lp=log(2.0_dp*(p3-x)/((p3-p1)*(p3-p2)))
      end if
    case(D_WEIBULL)
      if(x<0.0_dp)then
      lp=-huge_dp
      else if(x<=tiny(1.0_dp).and.p1<1.0_dp)then
      lp=log(huge_dp)
      else
      lp=log(p1/p2)+(p1-1.0_dp)*log(max(x/p2,tiny(1.0_dp)))-(x/p2)**p1
      end if
    case(D_BURR)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=log(p1*(p2-1.0_dp))+(p1-1.0_dp)*log(x)-p2*log1p_safe(x**p1)
      end if
    case(D_GIG)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      logc=0.5_dp*p1*log(p2/p3)-log(2.0_dp)-log_bessel_k_nu(sqrt(p2*p3),p1)
      lp=logc+(p1-1.0_dp)*log(x)-0.5_dp*(p3/x+p2*x)
      end if
    case(D_GIGA)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      logc=-log(2.0_dp)-p1*log(p3)-log_bessel_k_nu(p2,p1)
      lp=logc+(p1-1.0_dp)*log(x)-0.5_dp*p2*(x/p3+p3/x)
      end if
    case(D_HYPERB)
      gam=sqrt(p1*p1-p2*p2)
      logc=log(gam)-log(2.0_dp*p1*p3)-log_bessel_k_nu(p3*gam,1.0_dp)
      lp=logc-p1*sqrt(p3*p3+(x-p4)**2)+p2*(x-p4)
    case(D_GHYP)
      gam=sqrt(p2*p2-p3*p3)
      nu=p1-0.5_dp
      y=sqrt(p4*p4+(x-p5)**2)
      logc=-0.5_dp*log(2.0_dp*pi)+p1*log(gam/p4)-nu*log(p2)-log_bessel_k_nu(p4*gam,p1)
      lp=logc+log_bessel_k_nu(p2*y,nu)+nu*log(y)+p3*(x-p5)
    case(D_VG)
      nu=p1-0.5_dp
      y=abs(x-p4)
      logc=p1*log(p2*p2-p3*p3)-0.5_dp*log(pi)-nu*log(2.0_dp*p2)-log_gamma(p1)
      if(y<1.0e-12_dp.and.nu>0.0_dp)then
      lp=logc-log(2.0_dp)+log_gamma(nu)+nu*log(2.0_dp/p2)+p3*(x-p4)
      else
      lp=logc+log_bessel_k_nu(p2*y,nu)+nu*log(max(y,tiny(1.0_dp)))+p3*(x-p4)
      end if
    case(D_MEIX)
      y=(x-p4)/p1
      logc=2.0_dp*p3*log(2.0_dp*cos(p2/2.0_dp))-log(2.0_dp*p1*pi)-log_gamma(2.0_dp*p3)
      lp=logc+p2*y+2.0_dp*log_gamma_abs_complex(p3,y)
    case(D_PLANCK)
      if(x<=0.0_dp)then
      lp=-huge_dp
      else
      lp=p1*log(x)-log(max(expm1_safe(x),tiny(1.0_dp)))-log_gamma(p1+1.0_dp)-log(zeta_s(p1+1.0_dp))
      end if
    case default
    lp=-huge_dp
    end select
  end function base_cont_logpdf

  real(dp) function base_cont_pdf(d,x) result(y)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    real(dp)::lp
    lp=base_cont_logpdf(d,x)
    if(lp < log(tiny(1.0_dp)))then
    y=0.0_dp
    else if(lp>log(huge_dp))then
    y=huge_dp
    else
    y=exp(lp)
    end if
  end function

  real(dp) function cont_integrand(d,t,mode) result(y)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::t
    integer,intent(in)::mode
    real(dp)::z,c
    select case(mode)
    case(0)
      y=base_cont_pdf(d,t)
    case(1)
      if(t>=1.0_dp)then
        y=0.0_dp
      else
        z=t/(1.0_dp-t)
        y=base_cont_pdf(d,z)/(1.0_dp-t)**2
      end if
    case(2)
      if(t<=0.0_dp .or. t>=1.0_dp)then
        y=0.0_dp
      else
        z=tan(pi*(t-0.5_dp))
        c=cos(pi*(t-0.5_dp))
        y=base_cont_pdf(d,z)*pi/(c*c)
      end if
    case default
      y=0.0_dp
    end select
  end function cont_integrand

  recursive real(dp) function cont_simpson(d,mode,a,b,fa,fb,fm,s,eps,depth) result(v)
    class(continuous_distribution),intent(in)::d
    integer,intent(in)::mode,depth
    real(dp),intent(in)::a,b,fa,fb,fm,s,eps
    real(dp)::m,lm,rm,fl,fr,sl,sr
    m=0.5_dp*(a+b)
    lm=0.5_dp*(a+m)
    rm=0.5_dp*(m+b)
    fl=cont_integrand(d,lm,mode)
    fr=cont_integrand(d,rm,mode)
    sl=(m-a)*(fa+4.0_dp*fl+fm)/6.0_dp
    sr=(b-m)*(fm+4.0_dp*fr+fb)/6.0_dp
    if(depth<=0 .or. abs(sl+sr-s)<=15.0_dp*eps)then
      v=sl+sr+(sl+sr-s)/15.0_dp
    else
      v=cont_simpson(d,mode,a,m,fa,fm,fl,sl,0.5_dp*eps,depth-1)+ &
        cont_simpson(d,mode,m,b,fm,fb,fr,sr,0.5_dp*eps,depth-1)
    end if
  end function cont_simpson

  real(dp) function cont_adaptive_integral(d,a,b,mode,tol) result(v)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::a,b,tol
    integer,intent(in)::mode
    real(dp)::fa,fb,fm,s
    fa=cont_integrand(d,a,mode)
    fb=cont_integrand(d,b,mode)
    fm=cont_integrand(d,0.5_dp*(a+b),mode)
    s=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    v=cont_simpson(d,mode,a,b,fa,fb,fm,s,tol,25)
  end function cont_adaptive_integral

  real(dp) function numeric_base_cdf(d,x) result(p)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    real(dp)::tmax
    if(x<=d%support_lb)then
      p=0.0_dp
      return
    end if
    if(x>=d%support_ub)then
      p=1.0_dp
      return
    end if
    if(d%support_lb > -0.5_dp*huge_dp .and. d%support_ub < 0.5_dp*huge_dp)then
      p=cont_adaptive_integral(d,d%support_lb,x,0,1.0e-9_dp)
    else if(d%support_lb>=0.0_dp .and. d%support_ub>0.5_dp*huge_dp)then
      tmax=x/(1.0_dp+x)
      p=cont_adaptive_integral(d,0.0_dp,tmax,1,1.0e-9_dp)
    else
      tmax=0.5_dp+atan(x)/pi
      p=cont_adaptive_integral(d,1.0e-10_dp,tmax,2,1.0e-9_dp)
    end if
    p=clamp(p,0.0_dp,1.0_dp)
  end function numeric_base_cdf

  real(dp) function base_cont_cdf(d,x) result(p)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    real(dp)::p1,p2,p3,z,a
    if(x<=d%support_lb)then
    p=0.0_dp
    return
    end if
    if(x>=d%support_ub)then
    p=1.0_dp
    return
    end if
    p1=d%params(1)
    p2=d%params(2)
    p3=d%params(3)
    select case(d%id)
    case(D_CUSTOM)
    if(associated(d%cdf_cb))then
    p=d%cdf_cb(x,d%params(1:d%nparams))
    else
    p=numeric_base_cdf(d,x)
    end if
    case(D_NORMAL)
    p=normal_cdf((x-p1)/p2)
    case(D_BETA)
    p=beta_cdf(x,p1,p2)
    case(D_CAUCHY)
    p=0.5_dp+atan((x-p1)/p2)/pi
    case(D_CHI)
    p=chisq_cdf(x*x,p1)
    case(D_CHISQ)
    p=chisq_cdf(x,p1)
    case(D_EXP)
    p=1.0_dp-exp(-p1*x)
    case(D_F)
    p=f_cdf(x,p1,p2)
    case(D_FRECHET)
    z=(x-p2)/p3
    if(z<=0.0_dp)then
    p=0.0_dp
    else
    p=exp(-z**(-p1))
    end if
    case(D_GAMMA)
    p=gamma_cdf(x,p1,p2)
    case(D_GUMBEL)
    p=exp(-exp(-(x-p1)/p2))
    case(D_IG)
    if(x<=0.0_dp)then
    p=0.0_dp
    else
    a=sqrt(p2/x)
    p=normal_cdf(a*(x/p1-1.0_dp))+exp(2.0_dp*p2/p1)*normal_cdf(-a*(x/p1+1.0_dp))
    end if
    case(D_LAPLACE)
    if(x<p1)then
    p=0.5_dp*exp((x-p1)/p2)
    else
    p=1.0_dp-0.5_dp*exp(-(x-p1)/p2)
    end if
    case(D_LNORM)
    if(x<=0.0_dp)then
    p=0.0_dp
    else
    p=normal_cdf((log(x)-p1)/p2)
    end if
    case(D_LOGIS)
    p=1.0_dp/(1.0_dp+exp(-(x-p1)/p2))
    case(D_LOMAX)
    if(x<0.0_dp)then
    p=0.0_dp
    else
    p=1.0_dp-(1.0_dp+x/p2)**(-p1)
    end if
    case(D_PARETO)
    if(x<p1)then
    p=0.0_dp
    else
    p=1.0_dp-(p1/x)**p2
    end if
    case(D_POWEREXP)
    z=0.5_dp*reg_gamma_p(1.0_dp/p1,abs(x)**p1)
    if(x<0.0_dp)then
    p=0.5_dp-z
    else
    p=0.5_dp+z
    end if
    case(D_RAYLEIGH)
    if(x<0.0_dp)then
    p=0.0_dp
    else
    p=1.0_dp-exp(-0.5_dp*(x/p1)**2)
    end if
    case(D_T)
    p=student_t_cdf(x,p1)
    case(D_TRIANG)
    if(x<=p1)then
    p=0.0_dp
    else if(x<=p2)then
    p=(x-p1)**2/((p3-p1)*(p2-p1))
    else if(x<p3)then
    p=1.0_dp-(p3-x)**2/((p3-p1)*(p3-p2))
    else
    p=1.0_dp
    end if
    case(D_WEIBULL)
    if(x<0.0_dp)then
    p=0.0_dp
    else
    p=1.0_dp-exp(-(x/p2)**p1)
    end if
    case(D_BURR)
    if(x<0.0_dp)then
    p=0.0_dp
    else
    p=1.0_dp-(1.0_dp+x**p1)**(1.0_dp-p2)
    end if
    case default
    p=numeric_base_cdf(d,x)
    end select
    p=clamp(p,0.0_dp,1.0_dp)
  end function base_cont_cdf

  real(dp) function trunc_norm_cont(d) result(z)
    class(continuous_distribution),intent(in)::d
    z=base_cont_cdf(d,d%ub)-base_cont_cdf(d,d%lb)
    if(z<=0.0_dp)z=1.0_dp
  end function trunc_norm_cont

  real(dp) function cont_pdf(self,x) result(y)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::x
    if(x<self%lb.or.x>self%ub)then
    y=0.0_dp
    else
    y=base_cont_pdf(self,x)/trunc_norm_cont(self)
    end if
  end function
  real(dp) function cont_logpdf(self,x) result(y)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::x
    if(x<self%lb.or.x>self%ub)then
    y=-huge_dp
    else
    y=base_cont_logpdf(self,x)-log(trunc_norm_cont(self))
    end if
  end function
  real(dp) function cont_dpdf(self,x) result(y)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::x
    real(dp)::h,xl,xr
    if(self%id==D_CUSTOM .and. associated(self%dpdf_cb))then
      y=self%dpdf_cb(x,self%params(1:self%nparams))/trunc_norm_cont(self)
      return
    end if
    h=sqrt(epsilon(1.0_dp))*(1.0_dp+abs(x))
    xl=max(self%lb,x-h)
    xr=min(self%ub,x+h)
    if(xr<=xl)then
      y=0.0_dp
    else
      y=(self%pdf(xr)-self%pdf(xl))/(xr-xl)
    end if
  end function cont_dpdf

  real(dp) function cont_dlogpdf(self,x) result(y)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::x
    real(dp)::f
    f=self%pdf(x)
    if(f<=tiny(1.0_dp))then
      y=0.0_dp
    else
      y=self%dpdf(x)/f
    end if
  end function cont_dlogpdf

  real(dp) function cont_cdf(self,x) result(p)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::x
    real(dp)::a,b
    if(x<=self%lb)then
    p=0.0_dp
    return
    end if
    if(x>=self%ub)then
    p=1.0_dp
    return
    end if
    a=base_cont_cdf(self,self%lb)
    b=base_cont_cdf(self,self%ub)
    p=(base_cont_cdf(self,x)-a)/(b-a)
    p=clamp(p,0.0_dp,1.0_dp)
  end function

  real(dp) function base_quantile(d,p) result(x)
    class(continuous_distribution),intent(in)::d
    real(dp),intent(in)::p
    real(dp)::p1,p2,p3,lo,hi,mid,z
    integer::i
    p1=d%params(1)
    p2=d%params(2)
    p3=d%params(3)
    select case(d%id)
    case(D_NORMAL)
    x=p1+p2*normal_quantile(p)
    case(D_BETA)
    x=beta_quantile(p,p1,p2)
    case(D_CAUCHY)
    x=p1+p2*tan(pi*(p-0.5_dp))
    case(D_CHI)
    x=sqrt(chisq_quantile(p,p1))
    case(D_CHISQ)
    x=chisq_quantile(p,p1)
    case(D_EXP)
    x=-log(1.0_dp-p)/p1
    case(D_F)
    x=f_quantile(p,p1,p2)
    case(D_FRECHET)
    x=p2+p3*(-log(p))**(-1.0_dp/p1)
    case(D_GAMMA)
    x=gamma_quantile(p,p1,p2)
    case(D_GUMBEL)
    x=p1-p2*log(-log(p))
    case(D_LAPLACE)
    if(p<0.5_dp)then
    x=p1+p2*log(2.0_dp*p)
    else
    x=p1-p2*log(2.0_dp*(1.0_dp-p))
    end if
    case(D_LNORM)
    x=exp(p1+p2*normal_quantile(p))
    case(D_LOGIS)
    x=p1+p2*log(p/(1.0_dp-p))
    case(D_LOMAX)
    x=p2*((1.0_dp-p)**(-1.0_dp/p1)-1.0_dp)
    case(D_PARETO)
    x=p1*(1.0_dp-p)**(-1.0_dp/p2)
    case(D_RAYLEIGH)
    x=p1*sqrt(-2.0_dp*log(1.0_dp-p))
    case(D_T)
    x=student_t_quantile(p,p1)
    case(D_TRIANG)
    z=(p2-p1)/(p3-p1)
    if(p<z)then
    x=p1+sqrt(p*(p3-p1)*(p2-p1))
    else
    x=p3-sqrt((1.0_dp-p)*(p3-p1)*(p3-p2))
    end if
    case(D_WEIBULL)
    x=p2*(-log(1.0_dp-p))**(1.0_dp/p1)
    case(D_BURR)
    x=((1.0_dp-p)**(1.0_dp/(1.0_dp-p2))-1.0_dp)**(1.0_dp/p1)
    case default
      lo=d%center-1.0_dp
      hi=d%center+1.0_dp
      if(d%support_lb>-0.5_dp*huge_dp)lo=max(lo,d%support_lb)
      if(d%support_ub<0.5_dp*huge_dp)hi=min(hi,d%support_ub)
      do while(base_cont_cdf(d,lo)>p .and. lo>d%support_lb)
      lo=max(d%support_lb,d%center-2.0_dp*(d%center-lo+1.0_dp))
      end do
      do while(base_cont_cdf(d,hi)<p .and. hi<d%support_ub)
      hi=min(d%support_ub,d%center+2.0_dp*(hi-d%center+1.0_dp))
      end do
      do i=1,120
      mid=0.5_dp*(lo+hi)
      if(base_cont_cdf(d,mid)<p)then
      lo=mid
      else
      hi=mid
      end if
      end do
      x=0.5_dp*(lo+hi)
    end select
  end function base_quantile

  real(dp) function cont_quantile(self,p) result(x)
    class(continuous_distribution),intent(in)::self
    real(dp),intent(in)::p
    real(dp)::a,b,q
    if(p<=0.0_dp)then
    x=self%lb
    return
    end if
    if(p>=1.0_dp)then
    x=self%ub
    return
    end if
    a=base_cont_cdf(self,self%lb)
    b=base_cont_cdf(self,self%ub)
    q=a+p*(b-a)
    x=base_quantile(self,q)
  end function
  real(dp) function sample_gig2(rng,theta,psi,chi) result(x)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::theta,psi,chi
    real(dp)::ymode,l,r,hl,hr,sl,sr,z,la,ra,pick,y,hat,lf,u
    integer::iter
    ! For Y=log(X), log f_Y(y) = theta*y - (chi*exp(-y)+psi*exp(y))/2,
    ! which is strictly concave for psi,chi > 0. Two tangents therefore form
    ! a global integrable rejection envelope.
    if(theta>=0.0_dp)then
      ymode=log((theta+sqrt(theta*theta+psi*chi))/psi)
    else
      ymode=log(chi/(sqrt(theta*theta+psi*chi)-theta))
    end if
    l=ymode-1.0_dp
    r=ymode+1.0_dp
    do while(gig_dlog(l,theta,psi,chi)<=0.0_dp)
      l=l-1.0_dp
    end do
    do while(gig_dlog(r,theta,psi,chi)>=0.0_dp)
      r=r+1.0_dp
    end do
    hl=gig_log(l,theta,psi,chi)
    hr=gig_log(r,theta,psi,chi)
    sl=gig_dlog(l,theta,psi,chi)
    sr=gig_dlog(r,theta,psi,chi)
    z=(hl-hr-sl*l+sr*r)/(sr-sl)
    la=hl+sl*(z-l)-log(sl)
    ra=hr+sr*(z-r)-log(-sr)
    pick=1.0_dp/(1.0_dp+exp(ra-la))
    do iter=1,100000
      if(rng_uniform(rng)<pick)then
        y=z+log(rng_uniform(rng))/sl
        hat=hl+sl*(y-l)
      else
        y=z+log(rng_uniform(rng))/sr
        hat=hr+sr*(y-r)
      end if
      lf=gig_log(y,theta,psi,chi)
      u=log(rng_uniform(rng))
      if(u<=lf-hat)then
        x=exp(y)
        return
      end if
    end do
    error stop 'sample_gig2: rejection sampler failed'
  end function sample_gig2

  pure real(dp) function gig_log(y,theta,psi,chi) result(v)
    real(dp),intent(in)::y,theta,psi,chi
    v=theta*y-0.5_dp*(chi*exp(-y)+psi*exp(y))
  end function gig_log

  pure real(dp) function gig_dlog(y,theta,psi,chi) result(v)
    real(dp),intent(in)::y,theta,psi,chi
    v=theta+0.5_dp*chi*exp(-y)-0.5_dp*psi*exp(y)
  end function gig_dlog

  real(dp) function cont_sample(self,rng) result(x)
    class(continuous_distribution),intent(in)::self
    type(rng_state),intent(inout)::rng
    real(dp)::u,y,w,gam
    integer::k
    logical::full_domain
    full_domain=(self%lb<=self%support_lb .and. self%ub>=self%support_ub)
    if(.not.full_domain)then
      x=self%quantile(rng_uniform(rng))
      return
    end if
    select case(self%id)
    case(D_NORMAL)
    x=self%params(1)+self%params(2)*rng_normal(rng)
    case(D_BETA)
    x=rng_beta(rng,self%params(1),self%params(2))
    case(D_CAUCHY)
    x=rng_cauchy(rng,self%params(1),self%params(2))
    case(D_CHI)
    x=sqrt(rng_chisq(rng,self%params(1)))
    case(D_CHISQ)
    x=rng_chisq(rng,self%params(1))
    case(D_EXP)
    x=rng_exponential(rng,self%params(1))
    case(D_F)
    x=rng_f(rng,self%params(1),self%params(2))
    case(D_FRECHET,D_GUMBEL,D_LOMAX,D_PARETO,D_RAYLEIGH,D_TRIANG,D_WEIBULL,D_BURR)
      x=base_quantile(self,rng_uniform(rng))
    case(D_GAMMA)
    x=rng_gamma(rng,self%params(1),self%params(2))
    case(D_IG)
    x=rng_inverse_gaussian(rng,self%params(1),self%params(2))
    case(D_LAPLACE)
      u=rng_uniform(rng)-0.5_dp
      x=self%params(1)-self%params(2)*sign(1.0_dp,u)*log(1.0_dp-2.0_dp*abs(u))
    case(D_LNORM)
    x=exp(self%params(1)+self%params(2)*rng_normal(rng))
    case(D_LOGIS)
    u=rng_uniform(rng)
    x=self%params(1)+self%params(2)*log(u/(1.0_dp-u))
    case(D_POWEREXP)
      y=rng_gamma(rng,1.0_dp/self%params(1))
      if(rng_uniform(rng)<0.5_dp)then
      x=-y**(1.0_dp/self%params(1))
      else
      x=y**(1.0_dp/self%params(1))
      end if
    case(D_SLASH)
    x=rng_normal(rng)/rng_uniform(rng)
    case(D_T)
    x=rng_student_t(rng,self%params(1))
    case(D_PLANCK)
      ! f(x) propto x^a/(exp(x)-1) is a zeta mixture of Gamma(a+1,rate=k).
      u=rng_uniform(rng)
      y=0.0_dp
      k=0
      do
        k=k+1
        y=y+real(k,dp)**(-(self%params(1)+1.0_dp))/zeta_s(self%params(1)+1.0_dp)
        if(y>=u)exit
      end do
      x=rng_gamma(rng,self%params(1)+1.0_dp,1.0_dp/real(k,dp))
    case(D_GIG)
      x=sample_gig2(rng,self%params(1),self%params(2),self%params(3))
    case(D_GIGA)
      x=sample_gig2(rng,self%params(1),self%params(2)/self%params(3), &
        self%params(2)*self%params(3))
    case(D_GHYP)
      gam=sqrt(self%params(2)**2-self%params(3)**2)
      w=sample_gig2(rng,self%params(1),gam*gam,self%params(4)**2)
      x=self%params(5)+self%params(3)*w+sqrt(w)*rng_normal(rng)
    case(D_HYPERB)
      gam=sqrt(self%params(1)**2-self%params(2)**2)
      w=sample_gig2(rng,1.0_dp,gam*gam,self%params(3)**2)
      x=self%params(4)+self%params(2)*w+sqrt(w)*rng_normal(rng)
    case(D_VG)
      gam=self%params(2)**2-self%params(3)**2
      w=rng_gamma(rng,self%params(1),2.0_dp/gam)
      x=self%params(4)+self%params(3)*w+sqrt(w)*rng_normal(rng)
    case default
      ! Meixner and arbitrary custom laws use numerical inversion.
      x=self%quantile(rng_uniform(rng))
    end select
  end function
  subroutine cont_sample_n(self,rng,x)
    class(continuous_distribution),intent(in)::self
    type(rng_state),intent(inout)::rng
    real(dp),intent(out)::x(:)
    integer::i
    do i=1,size(x)
    x(i)=self%sample(rng)
    end do
  end subroutine

  function make_discr(id,name,p,sl,su,lb,ub) result(d)
    integer,intent(in)::id,sl,su
    character(len=*),intent(in)::name
    real(dp),intent(in)::p(:)
    integer,intent(in),optional::lb,ub
    type(discrete_distribution)::d
    d%id=id
    d%name=name
    d%nparams=min(size(p),8)
    d%params(1:d%nparams)=p(1:d%nparams)
    d%support_lb=sl
    d%support_ub=su
    d%lb=sl
    d%ub=su
    if(present(lb))d%lb=max(sl,lb)
    if(present(ub))d%ub=min(su,ub)
  end function
  function udbinom(size,prob,lb,ub) result(d)
  integer,intent(in)::size
  real(dp),intent(in)::prob
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_BINOM,'binomial',[real(size,dp),prob],0,size,lb,ub)
  end function
  function udgeom(prob,lb,ub) result(d)
  real(dp),intent(in)::prob
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_GEOM,'geometric',[prob],0,huge(1),lb,ub)
  end function
  function udhyper(m,n,k,lb,ub) result(d)
  integer,intent(in)::m,n,k
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_HYPER,'hypergeometric',[real(m,dp),real(n,dp),real(k,dp)],max(0,k-n),min(k,m),lb,ub)
  end function
  function udlogarithmic(shape,lb,ub) result(d)
  real(dp),intent(in)::shape
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_LOGARITHMIC,'logarithmic',[shape],1,huge(1),lb,ub)
  end function
  function udnbinom(size,prob,lb,ub) result(d)
  real(dp),intent(in)::size,prob
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_NBINOM,'negative_binomial',[size,prob],0,huge(1),lb,ub)
  end function
  function udpois(lambda,lb,ub) result(d)
  real(dp),intent(in)::lambda
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_POIS,'poisson',[lambda],0,huge(1),lb,ub)
  end function
  function udzipf(rho,lb,ub) result(d)
  real(dp),intent(in)::rho
  integer,intent(in),optional::lb,ub
  type(discrete_distribution)::d
  d=make_discr(Q_ZIPF,'zipf',[rho],1,huge(1),lb,ub)
  end function

  pure real(dp) function logchoose(n,k) result(v)
  integer,intent(in)::n,k
  if(k<0.or.k>n)then
  v=-huge_dp
  else
  v=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
  end if
  end function
  real(dp) function base_discr_pmf(d,k) result(p)
    class(discrete_distribution),intent(in)::d
    integer,intent(in)::k
    integer::n,m,kk
    real(dp)::a,b
    if(k<d%support_lb.or.k>d%support_ub)then
    p=0.0_dp
    return
    end if
    a=d%params(1)
    b=d%params(2)
    select case(d%id)
    case(Q_CUSTOM)
    if(associated(d%pmf_cb))then
    p=d%pmf_cb(k,d%params(1:d%nparams))
    else
    p=0.0_dp
    end if
    case(Q_BINOM)
    n=nint(a)
    p=exp(logchoose(n,k)+real(k,dp)*log(b)+real(n-k,dp)*log1p_safe(-b))
    case(Q_GEOM)
    p=a*(1.0_dp-a)**k
    case(Q_HYPER)
    m=nint(d%params(1))
    n=nint(d%params(2))
    kk=nint(d%params(3))
    p=exp(logchoose(m,k)+logchoose(n,kk-k)-logchoose(m+n,kk))
    case(Q_LOGARITHMIC)
    p=-a**k/(real(k,dp)*log1p_safe(-a))
    case(Q_NBINOM)
    p=exp(log_gamma(real(k,dp)+a)-log_gamma(a)-log_gamma(real(k+1,dp))+a*log(b)+real(k,dp)*log1p_safe(-b))
    case(Q_POIS)
    if(a<=tiny(1.0_dp).and.k==0)then
    p=1.0_dp
    else
    p=exp(-a+real(k,dp)*log(max(a,tiny(1.0_dp)))-log_gamma(real(k+1,dp)))
    end if
    case(Q_ZIPF)
    p=real(k,dp)**(-a)/zeta_s(a)
    case(Q_TABLE)
    p=d%probvec(k-d%support_lb+1)
    case default
    p=0.0_dp
    end select
  end function
  real(dp) function base_discr_cdf(d,k) result(p)
    class(discrete_distribution),intent(in)::d
    integer,intent(in)::k
    integer::j,hi
    if(k<d%support_lb)then
    p=0.0_dp
    return
    end if
    if(k>=d%support_ub.and.d%support_ub<huge(1))then
    p=1.0_dp
    return
    end if
    if(d%id==Q_CUSTOM.and.associated(d%cdf_cb))then
    p=d%cdf_cb(k,d%params(1:d%nparams))
    return
    end if
    p=0.0_dp
    hi=k
    do j=d%support_lb,hi
    p=p+base_discr_pmf(d,j)
    if(1.0_dp-p<1.0e-15_dp)exit
    end do
    p=min(p,1.0_dp)
  end function
  real(dp) function discr_norm(d) result(z)
  class(discrete_distribution),intent(in)::d
  z=base_discr_cdf(d,d%ub)-base_discr_cdf(d,d%lb-1)
  if(z<=0.0_dp)z=1.0_dp
  end function
  real(dp) function discr_pmf(self,k) result(p)
  class(discrete_distribution),intent(in)::self
  integer,intent(in)::k
  if(k<self%lb.or.k>self%ub)then
  p=0.0_dp
  else
  p=base_discr_pmf(self,k)/discr_norm(self)
  end if
  end function
  real(dp) function discr_cdf(self,k) result(p)
  class(discrete_distribution),intent(in)::self
  integer,intent(in)::k
  real(dp)::a,b
  if(k<self%lb)then
  p=0.0_dp
  else if(k>=self%ub)then
  p=1.0_dp
  else
  a=base_discr_cdf(self,self%lb-1)
  b=base_discr_cdf(self,self%ub)
  p=(base_discr_cdf(self,k)-a)/(b-a)
  end if
  end function
  integer function discr_quantile(self,p) result(k)
  class(discrete_distribution),intent(in)::self
  real(dp),intent(in)::p
  integer::lo,hi,imax2
  imax2=shiftr(huge(hi),1)
  if(p<=0.0_dp)then
  k=self%lb
  return
  end if
  lo=self%lb
  hi=max(lo,1)
  if(self%ub<huge(1))then
  hi=self%ub
  else
  do while(self%cdf(hi)<p .and. hi<imax2)
  hi=min(imax2,max(hi+1,2*hi))
  end do
  end if
  do while(lo<hi)
  k=lo+(hi-lo)/2
  if(self%cdf(k)<p)then
  lo=k+1
  else
  hi=k
  end if
  end do
  k=lo
  end function
  integer function discr_sample(self,rng) result(k)
    class(discrete_distribution),intent(in)::self
    type(rng_state),intent(inout)::rng
    integer::i,m,n,draws,good,bad
    logical::full_domain
    full_domain=(self%lb==self%support_lb .and. self%ub==self%support_ub)
    if(.not.full_domain)then
    k=self%quantile(rng_uniform(rng))
    return
    end if
    select case(self%id)
    case(Q_BINOM)
    k=rng_binomial(rng,nint(self%params(1)),self%params(2))
    case(Q_GEOM)
    k=rng_geometric(rng,self%params(1))
    case(Q_NBINOM)
    k=rng_negative_binomial(rng,self%params(1),self%params(2))
    case(Q_POIS)
    k=rng_poisson(rng,self%params(1))
    case(Q_HYPER)
      m=nint(self%params(1))
      n=nint(self%params(2))
      draws=nint(self%params(3))
      good=m
      bad=n
      k=0
      do i=1,draws
        if(rng_uniform(rng)<real(good,dp)/real(good+bad,dp))then
        k=k+1
        good=good-1
        else
        bad=bad-1
        end if
      end do
    case default
    k=self%quantile(rng_uniform(rng))
    end select
  end function
  subroutine discr_sample_n(self,rng,x)
  class(discrete_distribution),intent(in)::self
  type(rng_state),intent(inout)::rng
  integer,intent(out)::x(:)
  integer::i
  do i=1,size(x)
  x(i)=self%sample(rng)
  end do
  end subroutine
end module runuran_distributions
