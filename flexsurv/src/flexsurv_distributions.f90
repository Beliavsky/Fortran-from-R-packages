! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_distributions
  use flexsurv_kinds, only : dp
  use flexsurv_math, only : fs_pi, normal_pdf, normal_logpdf, normal_cdf, &
    normal_quantile, regularized_gamma_p, regularized_gamma_q, gamma_quantile, &
    regularized_beta, beta_quantile, log_beta, rng_normal, rng_gamma, rng_beta, &
    integrate_gauss_legendre, safe_log, log1p_fs, expm1_fs
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, &
    ieee_quiet_nan, ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dist_exponential = 1
  integer, parameter, public :: dist_weibull = 2
  integer, parameter, public :: dist_weibull_ph = 3
  integer, parameter, public :: dist_gamma = 4
  integer, parameter, public :: dist_lognormal = 5
  integer, parameter, public :: dist_gompertz = 6
  integer, parameter, public :: dist_loglogistic = 7
  integer, parameter, public :: dist_gengamma = 8
  integer, parameter, public :: dist_genf = 9

  public :: dist_name, dist_npar, dist_valid
  public :: dist_logpdf, dist_pdf, dist_cdf, dist_logcdf
  public :: dist_survival, dist_logsurv, dist_hazard, dist_cumhaz
  public :: dist_quantile, dist_random, dist_mean, dist_rmst
  public :: dexp_fs, pexp_fs, qexp_fs, hexp_fs, cumhaz_exp_fs
  public :: dweibull_fs, pweibull_fs, qweibull_fs, hweibull_fs, cumhaz_weibull_fs
  public :: dweibullph_fs, pweibullph_fs, qweibullph_fs
  public :: dgamma_fs, pgamma_fs, qgamma_fs, hgamma_fs, cumhaz_gamma_fs
  public :: dlnorm_fs, plnorm_fs, qlnorm_fs, hlnorm_fs, cumhaz_lnorm_fs
  public :: dgompertz, pgompertz, qgompertz, hgompertz, cumhaz_gompertz, rgompertz
  public :: dllogis, pllogis, qllogis, hllogis, cumhaz_llogis, rllogis
  public :: dgengamma, pgengamma, qgengamma, hgengamma, cumhaz_gengamma, rgengamma
  public :: dgengamma_orig, pgengamma_orig, qgengamma_orig, rgengamma_orig
  public :: dgenf, pgenf, qgenf, hgenf, cumhaz_genf, rgenf
  public :: dgenf_orig, pgenf_orig, qgenf_orig, rgenf_orig

contains

  pure function dist_name(id) result(name)
    integer, intent(in) :: id
    character(len=16) :: name
    select case(id)
    case(dist_exponential); name='exponential'
    case(dist_weibull); name='weibull'
    case(dist_weibull_ph); name='weibullph'
    case(dist_gamma); name='gamma'
    case(dist_lognormal); name='lognormal'
    case(dist_gompertz); name='gompertz'
    case(dist_loglogistic); name='loglogistic'
    case(dist_gengamma); name='gengamma'
    case(dist_genf); name='genf'
    case default; name='unknown'
    end select
  end function dist_name

  pure integer function dist_npar(id) result(n)
    integer,intent(in)::id
    select case(id)
    case(dist_exponential);n=1
    case(dist_weibull,dist_weibull_ph,dist_gamma,dist_lognormal,dist_gompertz,dist_loglogistic);n=2
    case(dist_gengamma);n=3
    case(dist_genf);n=4
    case default;n=0
    end select
  end function dist_npar

  pure logical function dist_valid(id, par) result(ok)
    integer,intent(in)::id
    real(dp),intent(in)::par(:)
    ok=size(par)>=dist_npar(id)
    if(.not.ok)return
    select case(id)
    case(dist_exponential);ok=par(1)>=0.0_dp
    case(dist_weibull);ok=par(1)>0.0_dp.and.par(2)>0.0_dp
    case(dist_weibull_ph);ok=par(1)>0.0_dp.and.par(2)>=0.0_dp
    case(dist_gamma);ok=par(1)>0.0_dp.and.par(2)>0.0_dp
    case(dist_lognormal);ok=par(2)>0.0_dp
    case(dist_gompertz);ok=par(2)>=0.0_dp
    case(dist_loglogistic);ok=par(1)>0.0_dp.and.par(2)>0.0_dp
    case(dist_gengamma);ok=par(2)>0.0_dp
    case(dist_genf);ok=par(2)>0.0_dp.and.par(4)>=0.0_dp
    case default;ok=.false.
    end select
  end function dist_valid

  pure real(dp) function dexp_fs(x,rate,log_pdf) result(y)
    real(dp),intent(in)::x,rate
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(rate<0.0_dp)then;y=nanv();return;end if
    if(x<0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    y=log(rate)-rate*x
    if(.not.lg)y=exp(y)
  end function dexp_fs

  pure real(dp) function pexp_fs(x,rate,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,rate
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(rate<0.0_dp)then;y=nanv();return;end if
    if(x<0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt);else if(lt)then;y=-expm1_fs(-rate*x);else;y=exp(-rate*x);end if
    if(lp)y=log(y)
  end function pexp_fs

  pure real(dp) function qexp_fs(p,rate) result(x)
    real(dp),intent(in)::p,rate
    if(p<=0.0_dp)then;x=0.0_dp;else if(p>=1.0_dp)then;x=posinf();else;x=-log1p_fs(-p)/rate;end if
  end function qexp_fs

  pure real(dp) function hexp_fs(x,rate) result(h)
    real(dp),intent(in)::x,rate;h=merge(rate,0.0_dp,x>=0.0_dp)
  end function hexp_fs
  pure real(dp) function cumhaz_exp_fs(x,rate) result(h)
    real(dp),intent(in)::x,rate;h=merge(rate*x,0.0_dp,x>=0.0_dp)
  end function cumhaz_exp_fs

  pure real(dp) function dweibull_fs(x,shape,scale,log_pdf) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(shape<=0.0_dp.or.scale<=0.0_dp)then;y=nanv();return;end if
    if(x<0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    if(x==0.0_dp)then
      if(shape<1.0_dp)then;y=merge(posinf(),posinf(),lg)
      else if(shape==1.0_dp)then;y=merge(-log(scale),1.0_dp/scale,lg)
      else;y=merge(neginf(),0.0_dp,lg);end if
      return
    end if
    y=log(shape)-log(scale)+(shape-1.0_dp)*log(x/scale)-(x/scale)**shape
    if(.not.lg)y=exp(y)
  end function dweibull_fs

  pure real(dp) function pweibull_fs(x,shape,scale,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    real(dp)::z
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(shape<=0.0_dp.or.scale<=0.0_dp)then;y=nanv();return;end if
    if(x<=0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt)
    else;z=(x/scale)**shape;if(lt)then;y=-expm1_fs(-z);else;y=exp(-z);end if;end if
    if(lp)y=log(y)
  end function pweibull_fs

  pure real(dp) function qweibull_fs(p,shape,scale) result(x)
    real(dp),intent(in)::p,shape,scale
    if(p<=0.0_dp)then;x=0.0_dp;else if(p>=1.0_dp)then;x=posinf();else;x=scale*(-log1p_fs(-p))**(1.0_dp/shape);end if
  end function qweibull_fs
  pure real(dp) function hweibull_fs(x,shape,scale) result(h)
    real(dp),intent(in)::x,shape,scale
    if(x<=0.0_dp)then;h=0.0_dp;else;h=shape/scale*(x/scale)**(shape-1.0_dp);end if
  end function hweibull_fs
  pure real(dp) function cumhaz_weibull_fs(x,shape,scale) result(h)
    real(dp),intent(in)::x,shape,scale
    if(x<=0.0_dp)then;h=0.0_dp;else;h=(x/scale)**shape;end if
  end function cumhaz_weibull_fs

  pure real(dp) function dweibullph_fs(x,shape,scale,log_pdf) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::log_pdf
    y=dweibull_fs(x,shape,scale**(-1.0_dp/shape),log_pdf)
  end function dweibullph_fs
  pure real(dp) function pweibullph_fs(x,shape,scale,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::lower_tail,log_p
    y=pweibull_fs(x,shape,scale**(-1.0_dp/shape),lower_tail,log_p)
  end function pweibullph_fs
  pure real(dp) function qweibullph_fs(p,shape,scale) result(x)
    real(dp),intent(in)::p,shape,scale
    x=qweibull_fs(p,shape,scale**(-1.0_dp/shape))
  end function qweibullph_fs

  pure real(dp) function dgamma_fs(x,shape,rate,log_pdf) result(y)
    real(dp),intent(in)::x,shape,rate
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(shape<=0.0_dp.or.rate<=0.0_dp)then;y=nanv();return;end if
    if(x<0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    if(x==0.0_dp)then
      if(shape<1.0_dp)then;y=posinf()
      else if(shape==1.0_dp)then;y=merge(log(rate),rate,lg)
      else;y=merge(neginf(),0.0_dp,lg);end if
      return
    end if
    y=shape*log(rate)-log_gamma(shape)+(shape-1.0_dp)*log(x)-rate*x
    if(.not.lg)y=exp(y)
  end function dgamma_fs

  pure real(dp) function pgamma_fs(x,shape,rate,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,rate
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt)
    else if(lt)then;y=regularized_gamma_p(shape,rate*x)
    else;y=regularized_gamma_q(shape,rate*x);end if
    if(lp)y=log(y)
  end function pgamma_fs
  real(dp) function qgamma_fs(p,shape,rate) result(x)
    real(dp),intent(in)::p,shape,rate;x=gamma_quantile(p,shape,1.0_dp/rate)
  end function qgamma_fs
  pure real(dp) function hgamma_fs(x,shape,rate) result(h)
    real(dp),intent(in)::x,shape,rate
    h=exp(dgamma_fs(x,shape,rate,.true.)-pgamma_fs(x,shape,rate,.false.,.true.))
  end function hgamma_fs
  pure real(dp) function cumhaz_gamma_fs(x,shape,rate) result(h)
    real(dp),intent(in)::x,shape,rate;h=-pgamma_fs(x,shape,rate,.false.,.true.)
  end function cumhaz_gamma_fs

  pure real(dp) function dlnorm_fs(x,meanlog,sdlog,log_pdf) result(y)
    real(dp),intent(in)::x,meanlog,sdlog
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(sdlog<=0.0_dp)then;y=nanv();return;end if
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    y=normal_logpdf((log(x)-meanlog)/sdlog)-log(sdlog*x)
    if(.not.lg)y=exp(y)
  end function dlnorm_fs
  pure real(dp) function plnorm_fs(x,meanlog,sdlog,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,meanlog,sdlog
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    real(dp)::z
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt)
    else;z=(log(x)-meanlog)/sdlog;y=merge(normal_cdf(z),normal_cdf(-z),lt);end if
    if(lp)y=log(y)
  end function plnorm_fs
  pure real(dp) function qlnorm_fs(p,meanlog,sdlog) result(x)
    real(dp),intent(in)::p,meanlog,sdlog;x=exp(meanlog+sdlog*normal_quantile(p))
  end function qlnorm_fs
  pure real(dp) function hlnorm_fs(x,meanlog,sdlog) result(h)
    real(dp),intent(in)::x,meanlog,sdlog
    h=exp(dlnorm_fs(x,meanlog,sdlog,.true.)-plnorm_fs(x,meanlog,sdlog,.false.,.true.))
  end function hlnorm_fs
  pure real(dp) function cumhaz_lnorm_fs(x,meanlog,sdlog) result(h)
    real(dp),intent(in)::x,meanlog,sdlog;h=-plnorm_fs(x,meanlog,sdlog,.false.,.true.)
  end function cumhaz_lnorm_fs

  pure real(dp) function exprel(z) result(v)
    real(dp),intent(in)::z
    if(abs(z)<1.0e-7_dp)then
      v=1.0_dp+z/2.0_dp+z*z/6.0_dp+z**3/24.0_dp
    else
      v=expm1_fs(z)/z
    end if
  end function exprel

  pure real(dp) function cumhaz_gompertz(x,shape,rate) result(h)
    real(dp),intent(in)::x,shape,rate
    if(x<=0.0_dp)then;h=0.0_dp
    else if(shape==0.0_dp)then;h=rate*x
    else;h=rate*x*exprel(shape*x);end if
  end function cumhaz_gompertz
  pure real(dp) function dgompertz(x,shape,rate,log_pdf) result(y)
    real(dp),intent(in)::x,shape,rate
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(rate<0.0_dp)then;y=nanv();return;end if
    if(x<0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    y=log(rate)+shape*x-cumhaz_gompertz(x,shape,rate)
    if(.not.lg)y=exp(y)
  end function dgompertz
  pure real(dp) function pgompertz(x,shape,rate,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,rate
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    real(dp)::h
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    h=cumhaz_gompertz(x,shape,rate)
    if(lt)then;y=-expm1_fs(-h);else;y=exp(-h);end if
    if(lp)y=log(y)
  end function pgompertz
  pure real(dp) function hgompertz(x,shape,rate) result(h)
    real(dp),intent(in)::x,shape,rate
    if(x<0.0_dp)then;h=0.0_dp;else;h=rate*exp(shape*x);end if
  end function hgompertz
  pure real(dp) function qgompertz(p,shape,rate) result(x)
    real(dp),intent(in)::p,shape,rate
    real(dp)::z,maxp
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(shape<0.0_dp)then
      maxp=1.0_dp-exp(rate/shape)
      if(p>=maxp)then;x=posinf();return;end if
    else if(p>=1.0_dp)then;x=posinf();return;end if
    z=-log1p_fs(-p)
    if(shape==0.0_dp)then;x=z/rate;else;x=log1p_fs(shape*z/rate)/shape;end if
  end function qgompertz
  real(dp) function rgompertz(shape,rate) result(x)
    real(dp),intent(in)::shape,rate;real(dp)::u
    call random_number(u);x=qgompertz(u,shape,rate)
  end function rgompertz

  pure real(dp) function cumhaz_llogis(x,shape,scale) result(h)
    real(dp),intent(in)::x,shape,scale
    if(x<=0.0_dp)then;h=0.0_dp;else;h=log1p_fs((x/scale)**shape);end if
  end function cumhaz_llogis
  pure real(dp) function dllogis(x,shape,scale,log_pdf) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::log_pdf
    logical::lg;real(dp)::z
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(shape<=0.0_dp.or.scale<=0.0_dp)then;y=nanv();return;end if
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    z=shape*(log(x)-log(scale))
    y=log(shape)-log(x)+z-2.0_dp*log1pexp(z)
    if(.not.lg)y=exp(y)
  end function dllogis
  pure real(dp) function pllogis(x,shape,scale,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,scale
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp;real(dp)::z
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt)
    else;z=shape*(log(x)-log(scale));if(lt)then;y=logistic(z);else;y=logistic(-z);end if;end if
    if(lp)y=log(y)
  end function pllogis
  pure real(dp) function qllogis(p,shape,scale) result(x)
    real(dp),intent(in)::p,shape,scale
    if(p<=0.0_dp)then;x=0.0_dp;else if(p>=1.0_dp)then;x=posinf();else;x=scale*(p/(1.0_dp-p))**(1.0_dp/shape);end if
  end function qllogis
  pure real(dp) function hllogis(x,shape,scale) result(h)
    real(dp),intent(in)::x,shape,scale
    if(x<=0.0_dp)then;h=0.0_dp;else;h=shape/x*logistic(shape*(log(x)-log(scale)));end if
  end function hllogis
  real(dp) function rllogis(shape,scale) result(x)
    real(dp),intent(in)::shape,scale;real(dp)::u;call random_number(u);x=qllogis(u,shape,scale)
  end function rllogis

  pure real(dp) function dgengamma(x,mu,sigma,qshape,log_pdf) result(y)
    real(dp),intent(in)::x,mu,sigma,qshape
    logical,intent(in),optional::log_pdf
    logical::lg;real(dp)::w,qi,qw,aq
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(sigma<=0.0_dp)then;y=nanv();return;end if
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    if(qshape==0.0_dp)then;y=dlnorm_fs(x,mu,sigma,lg);return;end if
    w=(log(x)-mu)/sigma;aq=abs(qshape);qi=1.0_dp/(qshape*qshape);qw=qshape*w
    y=-log(sigma*x)+log(aq)*(1.0_dp-2.0_dp*qi)+qi*(qw-exp(qw))-log_gamma(qi)
    if(.not.lg)y=exp(y)
  end function dgengamma

  pure real(dp) function pgengamma(x,mu,sigma,qshape,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,mu,sigma,qshape
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp;real(dp)::w,qq,z,p0
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;y=merge(0.0_dp,1.0_dp,lt);if(lp)y=log(y);return;end if
    if(qshape==0.0_dp)then;y=plnorm_fs(x,mu,sigma,lt,lp);return;end if
    w=(log(x)-mu)/sigma;qq=1.0_dp/(qshape*qshape);z=exp(qshape*w)*qq
    if(qshape>0.0_dp)then;p0=regularized_gamma_p(qq,z);else;p0=regularized_gamma_q(qq,z);end if
    y=merge(p0,1.0_dp-p0,lt);if(lp)y=log(y)
  end function pgengamma
  pure real(dp) function hgengamma(x,mu,sigma,qshape) result(h)
    real(dp),intent(in)::x,mu,sigma,qshape;h=exp(dgengamma(x,mu,sigma,qshape,.true.)-pgengamma(x,mu,sigma,qshape,.false.,.true.))
  end function hgengamma
  pure real(dp) function cumhaz_gengamma(x,mu,sigma,qshape) result(h)
    real(dp),intent(in)::x,mu,sigma,qshape;h=-pgengamma(x,mu,sigma,qshape,.false.,.true.)
  end function cumhaz_gengamma
  real(dp) function qgengamma(p,mu,sigma,qshape) result(x)
    real(dp),intent(in)::p,mu,sigma,qshape
    real(dp)::pp,g
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=posinf();return;end if
    if(qshape==0.0_dp)then;x=qlnorm_fs(p,mu,sigma);return;end if
    pp=merge(p,1.0_dp-p,qshape>0.0_dp)
    g=gamma_quantile(pp,1.0_dp/(qshape*qshape))
    x=exp(mu+sigma*log(qshape*qshape*g)/qshape)
  end function qgengamma
  real(dp) function rgengamma(mu,sigma,qshape) result(x)
    real(dp),intent(in)::mu,sigma,qshape;real(dp)::g
    if(qshape==0.0_dp)then;x=exp(mu+sigma*rng_normal())
    else;g=rng_gamma(1.0_dp/(qshape*qshape));x=exp(mu+sigma*log(qshape*qshape*g)/qshape);end if
  end function rgengamma

  pure real(dp) function dgengamma_orig(x,shape,scale,k,log_pdf) result(y)
    real(dp),intent(in)::x,shape,scale,k
    logical,intent(in),optional::log_pdf
    logical::lg
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    y=log(shape)-log_gamma(k)+(shape*k-1.0_dp)*log(x)-shape*k*log(scale)-(x/scale)**shape
    if(.not.lg)y=exp(y)
  end function dgengamma_orig
  pure real(dp) function pgengamma_orig(x,shape,scale,k,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,shape,scale,k
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp;real(dp)::p0
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;p0=0.0_dp;else;p0=regularized_gamma_p(k,(x/scale)**shape);end if
    y=merge(p0,1.0_dp-p0,lt);if(lp)y=log(y)
  end function pgengamma_orig
  real(dp) function qgengamma_orig(p,shape,scale,k) result(x)
    real(dp),intent(in)::p,shape,scale,k;x=scale*gamma_quantile(p,k)**(1.0_dp/shape)
  end function qgengamma_orig
  real(dp) function rgengamma_orig(shape,scale,k) result(x)
    real(dp),intent(in)::shape,scale,k;x=scale*rng_gamma(k)**(1.0_dp/shape)
  end function rgengamma_orig

  pure subroutine genf_shapes(qshape,p_shape,delta,s1,s2)
    real(dp),intent(in)::qshape,p_shape
    real(dp),intent(out)::delta,s1,s2
    real(dp)::tmp
    tmp=qshape*qshape+2.0_dp*p_shape;delta=sqrt(tmp)
    s1=2.0_dp/(tmp+qshape*delta);s2=2.0_dp/(tmp-qshape*delta)
  end subroutine genf_shapes

  pure real(dp) function dgenf(x,mu,sigma,qshape,p_shape,log_pdf) result(y)
    real(dp),intent(in)::x,mu,sigma,qshape,p_shape
    logical,intent(in),optional::log_pdf
    logical::lg;real(dp)::delta,s1,s2,expw,lx
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(sigma<=0.0_dp.or.p_shape<0.0_dp)then;y=nanv();return;end if
    if(p_shape==0.0_dp)then;y=dgengamma(x,mu,sigma,qshape,lg);return;end if
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    call genf_shapes(qshape,p_shape,delta,s1,s2);lx=log(x)
    expw=exp(delta/sigma*(lx-mu))
    y=log(delta)+(s1*delta/sigma)*(lx-mu)+s1*(log(s1)-log(s2))-log(sigma*x) &
      -(s1+s2)*log1p_fs(s1*expw/s2)-log_beta(s1,s2)
    if(.not.lg)y=exp(y)
  end function dgenf

  pure real(dp) function pgenf(x,mu,sigma,qshape,p_shape,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,mu,sigma,qshape,p_shape
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp;real(dp)::delta,s1,s2,expw,u,p0
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(p_shape==0.0_dp)then;y=pgengamma(x,mu,sigma,qshape,lt,lp);return;end if
    if(x<=0.0_dp)then;p0=0.0_dp
    else
      call genf_shapes(qshape,p_shape,delta,s1,s2)
      expw=exp(delta/sigma*(log(x)-mu))
      u=s1*expw/(s2+s1*expw)
      p0=regularized_beta(u,s1,s2)
    end if
    y=merge(p0,1.0_dp-p0,lt);if(lp)y=log(y)
  end function pgenf
  pure real(dp) function hgenf(x,mu,sigma,qshape,p_shape) result(h)
    real(dp),intent(in)::x,mu,sigma,qshape,p_shape
    h=exp(dgenf(x,mu,sigma,qshape,p_shape,.true.) &
      -pgenf(x,mu,sigma,qshape,p_shape,.false.,.true.))
  end function hgenf
  pure real(dp) function cumhaz_genf(x,mu,sigma,qshape,p_shape) result(h)
    real(dp),intent(in)::x,mu,sigma,qshape,p_shape;h=-pgenf(x,mu,sigma,qshape,p_shape,.false.,.true.)
  end function cumhaz_genf
  real(dp) function qgenf(p,mu,sigma,qshape,p_shape) result(x)
    real(dp),intent(in)::p,mu,sigma,qshape,p_shape
    real(dp)::delta,s1,s2,u,f
    if(p_shape==0.0_dp)then;x=qgengamma(p,mu,sigma,qshape);return;end if
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=posinf();return;end if
    call genf_shapes(qshape,p_shape,delta,s1,s2);u=beta_quantile(p,s1,s2)
    f=(s2/s1)*u/(1.0_dp-u);x=exp(mu+sigma/delta*log(f))
  end function qgenf
  real(dp) function rgenf(mu,sigma,qshape,p_shape) result(x)
    real(dp),intent(in)::mu,sigma,qshape,p_shape
    real(dp)::delta,s1,s2,u,f
    if(p_shape==0.0_dp)then;x=rgengamma(mu,sigma,qshape);return;end if
    call genf_shapes(qshape,p_shape,delta,s1,s2);u=rng_beta(s1,s2);f=(s2/s1)*u/(1.0_dp-u);x=exp(mu+sigma/delta*log(f))
  end function rgenf

  pure real(dp) function dgenf_orig(x,mu,sigma,s1,s2,log_pdf) result(y)
    real(dp),intent(in)::x,mu,sigma,s1,s2
    logical,intent(in),optional::log_pdf
    logical::lg;real(dp)::w,expw
    lg=.false.;if(present(log_pdf))lg=log_pdf
    if(x<=0.0_dp)then;y=merge(neginf(),0.0_dp,lg);return;end if
    w=(log(x)-mu)/sigma;expw=exp(w)
    y=-log(sigma*x)+s1*(log(s1)+w-log(s2))-(s1+s2)*log1p_fs(s1*expw/s2)-log_beta(s1,s2)
    if(.not.lg)y=exp(y)
  end function dgenf_orig
  pure real(dp) function pgenf_orig(x,mu,sigma,s1,s2,lower_tail,log_p) result(y)
    real(dp),intent(in)::x,mu,sigma,s1,s2
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp;real(dp)::u,p0
    lt=.true.;lp=.false.;if(present(lower_tail))lt=lower_tail;if(present(log_p))lp=log_p
    if(x<=0.0_dp)then;p0=0.0_dp;else;u=s1*exp((log(x)-mu)/sigma)/(s2+s1*exp((log(x)-mu)/sigma));p0=regularized_beta(u,s1,s2);end if
    y=merge(p0,1.0_dp-p0,lt);if(lp)y=log(y)
  end function pgenf_orig
  real(dp) function qgenf_orig(p,mu,sigma,s1,s2) result(x)
    real(dp),intent(in)::p,mu,sigma,s1,s2;real(dp)::u,f
    u=beta_quantile(p,s1,s2);f=(s2/s1)*u/(1.0_dp-u);x=exp(mu+sigma*log(f))
  end function qgenf_orig
  real(dp) function rgenf_orig(mu,sigma,s1,s2) result(x)
    real(dp),intent(in)::mu,sigma,s1,s2;real(dp)::u,f
    u=rng_beta(s1,s2);f=(s2/s1)*u/(1.0_dp-u);x=exp(mu+sigma*log(f))
  end function rgenf_orig

  pure real(dp) function dist_logpdf(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:)
    if(.not.dist_valid(id,par))then;v=nanv();return;end if
    select case(id)
    case(dist_exponential);v=dexp_fs(x,par(1),.true.)
    case(dist_weibull);v=dweibull_fs(x,par(1),par(2),.true.)
    case(dist_weibull_ph);v=dweibullph_fs(x,par(1),par(2),.true.)
    case(dist_gamma);v=dgamma_fs(x,par(1),par(2),.true.)
    case(dist_lognormal);v=dlnorm_fs(x,par(1),par(2),.true.)
    case(dist_gompertz);v=dgompertz(x,par(1),par(2),.true.)
    case(dist_loglogistic);v=dllogis(x,par(1),par(2),.true.)
    case(dist_gengamma);v=dgengamma(x,par(1),par(2),par(3),.true.)
    case(dist_genf);v=dgenf(x,par(1),par(2),par(3),par(4),.true.)
    case default;v=nanv()
    end select
  end function dist_logpdf
  pure real(dp) function dist_pdf(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:);v=exp(dist_logpdf(id,x,par))
  end function dist_pdf
  pure real(dp) function dist_cdf(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:)
    select case(id)
    case(dist_exponential);v=pexp_fs(x,par(1))
    case(dist_weibull);v=pweibull_fs(x,par(1),par(2))
    case(dist_weibull_ph);v=pweibullph_fs(x,par(1),par(2))
    case(dist_gamma);v=pgamma_fs(x,par(1),par(2))
    case(dist_lognormal);v=plnorm_fs(x,par(1),par(2))
    case(dist_gompertz);v=pgompertz(x,par(1),par(2))
    case(dist_loglogistic);v=pllogis(x,par(1),par(2))
    case(dist_gengamma);v=pgengamma(x,par(1),par(2),par(3))
    case(dist_genf);v=pgenf(x,par(1),par(2),par(3),par(4))
    case default;v=nanv()
    end select
  end function dist_cdf
  pure real(dp) function dist_logcdf(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:);v=log(dist_cdf(id,x,par))
  end function dist_logcdf
  pure real(dp) function dist_survival(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:)
    select case(id)
    case(dist_exponential);v=pexp_fs(x,par(1),.false.)
    case(dist_weibull);v=pweibull_fs(x,par(1),par(2),.false.)
    case(dist_weibull_ph);v=pweibullph_fs(x,par(1),par(2),.false.)
    case(dist_gamma);v=pgamma_fs(x,par(1),par(2),.false.)
    case(dist_lognormal);v=plnorm_fs(x,par(1),par(2),.false.)
    case(dist_gompertz);v=pgompertz(x,par(1),par(2),.false.)
    case(dist_loglogistic);v=pllogis(x,par(1),par(2),.false.)
    case(dist_gengamma);v=pgengamma(x,par(1),par(2),par(3),.false.)
    case(dist_genf);v=pgenf(x,par(1),par(2),par(3),par(4),.false.)
    case default;v=nanv()
    end select
  end function dist_survival
  pure real(dp) function dist_logsurv(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:);v=log(dist_survival(id,x,par))
  end function dist_logsurv
  pure real(dp) function dist_hazard(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:);v=exp(dist_logpdf(id,x,par)-dist_logsurv(id,x,par))
  end function dist_hazard
  pure real(dp) function dist_cumhaz(id,x,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::x,par(:);v=-dist_logsurv(id,x,par)
  end function dist_cumhaz

  real(dp) function dist_quantile(id,p,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::p,par(:)
    select case(id)
    case(dist_exponential);v=qexp_fs(p,par(1))
    case(dist_weibull);v=qweibull_fs(p,par(1),par(2))
    case(dist_weibull_ph);v=qweibullph_fs(p,par(1),par(2))
    case(dist_gamma);v=qgamma_fs(p,par(1),par(2))
    case(dist_lognormal);v=qlnorm_fs(p,par(1),par(2))
    case(dist_gompertz);v=qgompertz(p,par(1),par(2))
    case(dist_loglogistic);v=qllogis(p,par(1),par(2))
    case(dist_gengamma);v=qgengamma(p,par(1),par(2),par(3))
    case(dist_genf);v=qgenf(p,par(1),par(2),par(3),par(4))
    case default;v=nanv()
    end select
  end function dist_quantile

  real(dp) function dist_random(id,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::par(:);real(dp)::u
    select case(id)
    case(dist_gamma);v=rng_gamma(par(1),1.0_dp/par(2))
    case(dist_lognormal);v=exp(par(1)+par(2)*rng_normal())
    case(dist_gompertz);v=rgompertz(par(1),par(2))
    case(dist_loglogistic);v=rllogis(par(1),par(2))
    case(dist_gengamma);v=rgengamma(par(1),par(2),par(3))
    case(dist_genf);v=rgenf(par(1),par(2),par(3),par(4))
    case default
      call random_number(u);v=dist_quantile(id,u,par)
    end select
  end function dist_random

  real(dp) function dist_mean(id,par) result(v)
    integer,intent(in)::id;real(dp),intent(in)::par(:)
    real(dp)::q,sig,delta,s1,s2,r
    select case(id)
    case(dist_exponential);v=1.0_dp/par(1)
    case(dist_weibull);v=par(2)*gamma(1.0_dp+1.0_dp/par(1))
    case(dist_weibull_ph);v=par(2)**(-1.0_dp/par(1))*gamma(1.0_dp+1.0_dp/par(1))
    case(dist_gamma);v=par(1)/par(2)
    case(dist_lognormal);v=exp(par(1)+0.5_dp*par(2)*par(2))
    case(dist_loglogistic)
      if(par(1)<=1.0_dp)then;v=posinf();else;v=par(2)*fs_pi/par(1)/sin(fs_pi/par(1));end if
    case(dist_gompertz)
      if(par(1)<0.0_dp)then;v=posinf();else;v=dist_rmst(id,posinf(),par);end if
    case(dist_gengamma)
      q=par(3);sig=par(2)
      if(q==0.0_dp)then;v=exp(par(1)+0.5_dp*sig*sig)
      else if(1.0_dp/(q*q)+sig/q<=0.0_dp)then;v=posinf()
      else;v=exp(par(1)+(sig/q)*log(q*q)+log_gamma(1.0_dp/(q*q)+sig/q)-log_gamma(1.0_dp/(q*q)));end if
    case(dist_genf)
      if(par(4)==0.0_dp)then
        q=par(3);sig=par(2)
        if(q==0.0_dp)then;v=exp(par(1)+0.5_dp*sig*sig)
        else if(1.0_dp/(q*q)+sig/q<=0.0_dp)then;v=posinf()
        else;v=exp(par(1)+(sig/q)*log(q*q)+log_gamma(1.0_dp/(q*q)+sig/q)-log_gamma(1.0_dp/(q*q)));end if
      else
        call genf_shapes(par(3),par(4),delta,s1,s2);r=par(2)/delta
        if(s2<=r .or. s1+r<=0.0_dp)then;v=posinf()
        else;v=exp(par(1)+r*log(s2/s1)+log_beta(s1+r,s2-r)-log_beta(s1,s2));end if
      end if
    case default;v=nanv()
    end select
  end function dist_mean

  real(dp) function dist_rmst(id,t,par,start) result(v)
    integer,intent(in)::id
    real(dp),intent(in)::t,par(:)
    real(dp),intent(in),optional::start
    real(dp)::a,b,st,upper,tailprob
    st=0.0_dp;if(present(start))st=max(0.0_dp,start)
    if(t<=st)then;v=0.0_dp;return;end if
    if(.not.ieee_is_finite(t))then
      upper=dist_quantile(id,1.0_dp-1.0e-10_dp,par)
      if(.not.ieee_is_finite(upper))upper=max(100.0_dp,10.0_dp*max(1.0_dp,st))
    else;upper=t;end if
    a=st;b=upper
    v=integrate_gauss_legendre(survfn,a,b,96)
    if(st>0.0_dp)then
      tailprob=dist_survival(id,st,par)
      if(tailprob>0.0_dp)v=v/tailprob
    end if
  contains
    real(dp) function survfn(x) result(s)
      real(dp),intent(in)::x;s=dist_survival(id,x,par)
    end function survfn
  end function dist_rmst

  pure real(dp) function logistic(z) result(y)
    real(dp),intent(in)::z
    if(z>=0.0_dp)then;y=1.0_dp/(1.0_dp+exp(-z));else;y=exp(z)/(1.0_dp+exp(z));end if
  end function logistic
  pure real(dp) function log1pexp(z) result(y)
    real(dp),intent(in)::z
    if(z>35.0_dp)then;y=z;else if(z<-35.0_dp)then;y=exp(z);else;y=log1p_fs(exp(z));end if
  end function log1pexp
  pure real(dp) function posinf() result(x);x=ieee_value(0.0_dp,ieee_positive_inf);end function posinf
  pure real(dp) function neginf() result(x);x=-ieee_value(0.0_dp,ieee_positive_inf);end function neginf
  pure real(dp) function nanv() result(x);x=ieee_value(0.0_dp,ieee_quiet_nan);end function nanv

end module flexsurv_distributions
