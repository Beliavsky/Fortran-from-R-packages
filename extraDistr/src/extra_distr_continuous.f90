! SPDX-License-Identifier: GPL-2.0-only
module extra_distr_continuous
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use extra_distr_kinds, only : dp, pi, sqrt_two_pi
  use extra_distr_math
  use extra_distr_rng
  implicit none
  private

  public :: dbetapr,pbetapr,qbetapr,rbetapr
  public :: dbhatt,pbhatt,rbhatt
  public :: dfatigue,pfatigue,qfatigue,rfatigue
  public :: dfrechet,pfrechet,qfrechet,rfrechet
  public :: dgev,pgev,qgev,rgev,dgompertz,pgompertz,qgompertz,rgompertz
  public :: dgpd,pgpd,qgpd,rgpd,dgumbel,pgumbel,qgumbel,rgumbel
  public :: dhcauchy,phcauchy,qhcauchy,rhcauchy,dhnorm,phnorm,qhnorm,rhnorm
  public :: dht,pht,qht,rht,dhuber,phuber,qhuber,rhuber
  public :: dinvgamma,pinvgamma,qinvgamma,rinvgamma
  public :: dinvchisq,pinvchisq,qinvchisq,rinvchisq
  public :: dkumar,pkumar,qkumar,rkumar
  public :: dlaplace,plaplace,qlaplace,rlaplace
  public :: dlst,plst,qlst,rlst,dlomax,plomax,qlomax,rlomax
  public :: dnsbeta,pnsbeta,qnsbeta,rnsbeta
  public :: dpareto,ppareto,qpareto,rpareto,dpower,ppower,qpower,rpower
  public :: dprop,pprop,qprop,rprop,drayleigh,prayleigh,qrayleigh,rrayleigh
  public :: dsgomp,psgomp,rsgomp,dslash,pslash,rslash
  public :: dtriang,ptriang,qtriang,rtriang
  public :: dtnorm,ptnorm,qtnorm,rtnorm,qtlambda,rtlambda
  public :: dwald,pwald,rwald

contains

  real(dp) function dbetapr(x,shape1,shape2,scale,log_p) result(v)
    real(dp),intent(in)::x,shape1,shape2
    real(dp),intent(in),optional::scale
    logical,intent(in),optional::log_p
    real(dp)::s,z,d
    s=1.0_dp; if(present(scale))s=scale
    if(shape1<=0.0_dp .or. shape2<=0.0_dp .or. s<=0.0_dp) then; d=nan_dp()
    else if(x<=0.0_dp .or. .not.ieee_is_finite(x)) then; d=0.0_dp
    else; z=x/s; d=exp((shape1-1.0_dp)*log(z)-(shape1+shape2)*log(1.0_dp+z)-log_beta(shape1,shape2)-log(s)); end if
    v=apply_density_log(d,log_p)
  end function dbetapr
  real(dp) function pbetapr(q,shape1,shape2,scale,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,shape1,shape2; real(dp),intent(in),optional::scale
    logical,intent(in),optional::lower_tail,log_p; real(dp)::s,p,z
    s=1.0_dp;if(present(scale))s=scale
    if(shape1<=0.0_dp.or.shape2<=0.0_dp.or.s<=0.0_dp)then;p=nan_dp()
    else if(q<=0.0_dp)then;p=0.0_dp
    else if(.not.ieee_is_finite(q))then;p=1.0_dp
    else;z=q/s;p=regularized_beta(z/(1.0_dp+z),shape1,shape2);end if
    v=apply_tail(p,lower_tail,log_p)
  end function pbetapr
  real(dp) function qbetapr(prob,shape1,shape2,scale,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,shape1,shape2;real(dp),intent(in),optional::scale
    logical,intent(in),optional::lower_tail,log_p;real(dp)::p,s,y
    p=decode_probability(prob,lower_tail,log_p);s=1.0_dp;if(present(scale))s=scale
    if(shape1<=0.0_dp.or.shape2<=0.0_dp.or.s<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp()
    else if(p<=0.0_dp)then;x=0.0_dp
    else if(p>=1.0_dp)then;x=pos_inf()
    else;y=beta_quantile(p,shape1,shape2);x=s*y/(1.0_dp-y);end if
  end function qbetapr
  function rbetapr(n,shape1,shape2,scale) result(x)
    integer,intent(in)::n;real(dp),intent(in)::shape1,shape2;real(dp),intent(in),optional::scale
    real(dp),allocatable::x(:);real(dp)::s,y;integer::i
    s=1.0_dp;if(present(scale))s=scale;allocate(x(max(0,n)))
    do i=1,n;y=rbeta_scalar(shape1,shape2);x(i)=s*y/(1.0_dp-y);end do
  end function rbetapr

  pure real(dp) function bhatt_g(z) result(g)
    real(dp),intent(in)::z;g=z*normal_cdf(z)+normal_pdf(z)
  end function bhatt_g
  real(dp) function dbhatt(x,mu,sigma,a,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma,a;logical,intent(in),optional::log_p
    real(dp)::m,s,aa,d,z
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;aa=s;if(present(a))aa=a
    if(s<0.0_dp.or.aa<0.0_dp)then;d=nan_dp()
    else if(s==0.0_dp)then;d=merge(1.0_dp/(2.0_dp*aa),0.0_dp,x>=m-aa.and.x<=m+aa.and.aa>0.0_dp)
    else if(aa==0.0_dp)then;d=normal_pdf((x-m)/s)/s
    else;z=x-m;d=(normal_cdf((z+aa)/s)-normal_cdf((z-aa)/s))/(2.0_dp*aa);end if
    v=apply_density_log(d,log_p)
  end function dbhatt
  real(dp) function pbhatt(q,mu,sigma,a,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma,a;logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,aa,p,z
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;aa=s;if(present(a))aa=a
    if(s<0.0_dp.or.aa<0.0_dp)then;p=nan_dp()
    else if(s==0.0_dp)then;p=clamp_probability((q-(m-aa))/(2.0_dp*aa))
    else if(aa==0.0_dp)then;p=normal_cdf((q-m)/s)
    else;z=q-m;p=s/(2.0_dp*aa)*(bhatt_g((z+aa)/s)-bhatt_g((z-aa)/s));end if
    v=apply_tail(p,lower_tail,log_p)
  end function pbhatt
  function rbhatt(n,mu,sigma,a) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma,a;real(dp),allocatable::x(:)
    real(dp)::m,s,aa;integer::i;m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;aa=s;if(present(a))aa=a
    allocate(x(max(0,n)));do i=1,n;x(i)=m+s*rnorm_std()+aa*(2.0_dp*runif_open()-1.0_dp);end do
  end function rbhatt

  real(dp) function dfatigue(x,alpha,beta,mu,log_p) result(v)
    real(dp),intent(in)::x,alpha;real(dp),intent(in),optional::beta,mu;logical,intent(in),optional::log_p
    real(dp)::b,m,z,zb,bz,d;b=1.0_dp;m=0.0_dp;if(present(beta))b=beta;if(present(mu))m=mu
    if(alpha<=0.0_dp.or.b<=0.0_dp)then;d=nan_dp()
    else if(x<=m.or..not.ieee_is_finite(x))then;d=0.0_dp
    else;z=x-m;zb=sqrt(z/b);bz=sqrt(b/z);d=(zb+bz)/(2.0_dp*alpha*z)*normal_pdf((zb-bz)/alpha);end if
    v=apply_density_log(d,log_p)
  end function dfatigue
  real(dp) function pfatigue(q,alpha,beta,mu,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,alpha;real(dp),intent(in),optional::beta,mu;logical,intent(in),optional::lower_tail,log_p
    real(dp)::b,m,z,p;b=1.0_dp;m=0.0_dp;if(present(beta))b=beta;if(present(mu))m=mu
    if(alpha<=0.0_dp.or.b<=0.0_dp)then;p=nan_dp()
    else if(q<=m)then;p=0.0_dp
    else;z=q-m;p=normal_cdf((sqrt(z/b)-sqrt(b/z))/alpha);end if
    v=apply_tail(p,lower_tail,log_p)
  end function pfatigue
  real(dp) function qfatigue(prob,alpha,beta,mu,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,alpha;real(dp),intent(in),optional::beta,mu;logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,b,m,z,t;p=decode_probability(prob,lower_tail,log_p);b=1.0_dp;m=0.0_dp;if(present(beta))b=beta;if(present(mu))m=mu
    if(alpha<=0.0_dp.or.b<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();return;end if
    z=normal_quantile(p);t=0.5_dp*alpha*z;x=m+b*(t+sqrt(t*t+1.0_dp))**2
  end function qfatigue
  function rfatigue(n,alpha,beta,mu) result(x)
    integer,intent(in)::n;real(dp),intent(in)::alpha;real(dp),intent(in),optional::beta,mu;real(dp),allocatable::x(:)
    real(dp)::b,m,z,t;integer::i;b=1.0_dp;m=0.0_dp;if(present(beta))b=beta;if(present(mu))m=mu;allocate(x(max(0,n)))
    do i=1,n;z=rnorm_std();t=0.5_dp*alpha*z;x(i)=m+b*(t+sqrt(t*t+1.0_dp))**2;end do
  end function rfatigue

  real(dp) function dfrechet(x,lambda,mu,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::lambda,mu,sigma;logical,intent(in),optional::log_p
    real(dp)::la,m,s,z,d;la=1.0_dp;m=0.0_dp;s=1.0_dp;if(present(lambda))la=lambda;if(present(mu))m=mu;if(present(sigma))s=sigma
    if(la<=0.0_dp.or.s<=0.0_dp)then;d=nan_dp();else if(x<=m)then;d=0.0_dp;else;z=(x-m)/s;d=la/s*z**(-1.0_dp-la)*exp(-z**(-la));end if
    v=apply_density_log(d,log_p)
  end function dfrechet
  real(dp) function pfrechet(q,lambda,mu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::lambda,mu,sigma;logical,intent(in),optional::lower_tail,log_p
    real(dp)::la,m,s,p,z;la=1.0_dp;m=0.0_dp;s=1.0_dp;if(present(lambda))la=lambda;if(present(mu))m=mu;if(present(sigma))s=sigma
    if(la<=0.0_dp.or.s<=0.0_dp)then;p=nan_dp();else if(q<=m)then;p=0.0_dp;else;z=(q-m)/s;p=exp(-z**(-la));end if
    v=apply_tail(p,lower_tail,log_p)
  end function pfrechet
  real(dp) function qfrechet(prob,lambda,mu,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::lambda,mu,sigma;logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,la,m,s;p=decode_probability(prob,lower_tail,log_p);la=1.0_dp;m=0.0_dp;s=1.0_dp;if(present(lambda))la=lambda;if(present(mu))m=mu;if(present(sigma))s=sigma
    if(la<=0.0_dp.or.s<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();else if(p==0.0_dp)then;x=m;else if(p==1.0_dp)then;x=pos_inf();else;x=m+s*(-log(p))**(-1.0_dp/la);end if
  end function qfrechet
  function rfrechet(n,lambda,mu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::lambda,mu,sigma;real(dp),allocatable::x(:);real(dp)::la,m,s;integer::i
    la=1.0_dp;m=0.0_dp;s=1.0_dp;if(present(lambda))la=lambda;if(present(mu))m=mu;if(present(sigma))s=sigma;allocate(x(max(0,n)))
    do i=1,n;x(i)=m+s*(-log(runif_open()))**(-1.0_dp/la);end do
  end function rfrechet

  real(dp) function dgev(x,mu,sigma,xi,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::log_p
    real(dp)::m,s,k,z,t,d;m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi
    if(s<=0.0_dp)then;d=nan_dp();else;z=(x-m)/s;t=1.0_dp+k*z;if(k==0.0_dp)then;d=exp(-z-exp(-z))/s;else if(t>0.0_dp)then;d=t**(-1.0_dp-1.0_dp/k)*exp(-t**(-1.0_dp/k))/s;else;d=0.0_dp;end if;end if
    v=apply_density_log(d,log_p)
  end function dgev
  real(dp) function pgev(q,mu,sigma,xi,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,k,z,t,p;m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi
    if(s<=0.0_dp)then;p=nan_dp();else;z=(q-m)/s;t=1.0_dp+k*z;if(k==0.0_dp)then;p=exp(-exp(-z));else if(t>0.0_dp)then;p=exp(-t**(-1.0_dp/k));else;p=merge(1.0_dp,0.0_dp,k<0.0_dp);end if;end if
    v=apply_tail(p,lower_tail,log_p)
  end function pgev
  real(dp) function qgev(prob,mu,sigma,xi,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,m,s,k;p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi
    if(s<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();else if(k==0.0_dp)then;x=m-s*log(-log(p));else;x=m+s*((-log(p))**(-k)-1.0_dp)/k;end if
  end function qgev
  function rgev(n,mu,sigma,xi) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma,xi;real(dp),allocatable::x(:);real(dp)::m,s,k;integer::i
    m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi;allocate(x(max(0,n)))
    do i=1,n;x(i)=qgev(runif_open(),m,s,k);end do
  end function rgev

  real(dp) function dgompertz(x,a,b,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::a,b;logical,intent(in),optional::log_p;real(dp)::aa,bb,d
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b
    if(aa<=0.0_dp.or.bb<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp.or..not.ieee_is_finite(x))then;d=0.0_dp;else;d=aa*exp(bb*x-aa/bb*(exp(bb*x)-1.0_dp));end if
    v=apply_density_log(d,log_p)
  end function dgompertz
  real(dp) function pgompertz(q,a,b,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::aa,bb,p
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b
    if(aa<=0.0_dp.or.bb<=0.0_dp)then;p=nan_dp();else if(q<0.0_dp)then;p=0.0_dp;else;p=1.0_dp-exp(-aa/bb*(exp(bb*q)-1.0_dp));end if;v=apply_tail(p,lower_tail,log_p)
  end function pgompertz
  real(dp) function qgompertz(prob,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,aa,bb
    p=decode_probability(prob,lower_tail,log_p);aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b
    if(aa<=0.0_dp.or.bb<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();else;x=log(1.0_dp-bb/aa*log(1.0_dp-p))/bb;end if
  end function qgompertz
  function rgompertz(n,a,b) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::a,b;real(dp),allocatable::x(:);real(dp)::aa,bb;integer::i
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;allocate(x(max(0,n)));do i=1,n;x(i)=qgompertz(runif_open(),aa,bb);end do
  end function rgompertz

  real(dp) function dgpd(x,mu,sigma,xi,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::log_p;real(dp)::m,s,k,z,t,d
    m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi
    z=(x-m)/s;t=1.0_dp+k*z;if(s<=0.0_dp)then;d=nan_dp();else if(z<=0.0_dp.or.t<=0.0_dp)then;d=0.0_dp;else if(k==0.0_dp)then;d=exp(-z)/s;else;d=t**(-(k+1.0_dp)/k)/s;end if;v=apply_density_log(d,log_p)
  end function dgpd
  real(dp) function pgpd(q,mu,sigma,xi,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,k,z,t,p
    m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi;z=(q-m)/s;t=1.0_dp+k*z
    if(s<=0.0_dp)then;p=nan_dp();else if(z<=0.0_dp)then;p=0.0_dp;else if(t<=0.0_dp)then;p=1.0_dp;else if(k==0.0_dp)then;p=1.0_dp-exp(-z);else;p=1.0_dp-t**(-1.0_dp/k);end if;v=apply_tail(p,lower_tail,log_p)
  end function pgpd
  real(dp) function qgpd(prob,mu,sigma,xi,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mu,sigma,xi;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,m,s,k
    p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi
    if(s<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();else if(k==0.0_dp)then;x=m-s*log(1.0_dp-p);else;x=m+s*((1.0_dp-p)**(-k)-1.0_dp)/k;end if
  end function qgpd
  function rgpd(n,mu,sigma,xi) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma,xi;real(dp),allocatable::x(:);real(dp)::m,s,k;integer::i
    m=0.0_dp;s=1.0_dp;k=0.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(xi))k=xi;allocate(x(max(0,n)));do i=1,n;x(i)=qgpd(runif_open(),m,s,k);end do
  end function rgpd

  real(dp) function dgumbel(x,mu,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::log_p;real(dp)::m,s,z,d
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;z=(x-m)/s
    if(s<=0.0_dp)then;d=nan_dp();else;d=exp(-z-exp(-z))/s;end if;v=apply_density_log(d,log_p)
  end function dgumbel
  real(dp) function pgumbel(q,mu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,p
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(s<=0.0_dp)then;p=nan_dp();else;p=exp(-exp(-(q-m)/s));end if;v=apply_tail(p,lower_tail,log_p)
  end function pgumbel
  real(dp) function qgumbel(prob,mu,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,m,s
    p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(s<=0.0_dp)then;x=nan_dp();else;x=m-s*log(-log(p));end if
  end function qgumbel
  function rgumbel(n,mu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma;real(dp),allocatable::x(:);real(dp)::m,s;integer::i
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=qgumbel(runif_open(),m,s);end do
  end function rgumbel

  real(dp) function dhcauchy(x,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::sigma;logical,intent(in),optional::log_p;real(dp)::s,d
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp)then;d=0.0_dp;else;d=2.0_dp/(pi*s*(1.0_dp+(x/s)**2));end if;v=apply_density_log(d,log_p)
  end function dhcauchy
  real(dp) function phcauchy(q,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::s,p
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;p=nan_dp();else if(q<0.0_dp)then;p=0.0_dp;else;p=2.0_dp/pi*atan(q/s);end if;v=apply_tail(p,lower_tail,log_p)
  end function phcauchy
  real(dp) function qhcauchy(prob,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,s
    p=decode_probability(prob,lower_tail,log_p);s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;x=nan_dp();else;x=s*tan(0.5_dp*pi*p);end if
  end function qhcauchy
  function rhcauchy(n,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::sigma;real(dp),allocatable::x(:);real(dp)::s;integer::i
    s=1.0_dp;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=abs(s*tan(pi*(runif_open()-0.5_dp)));end do
  end function rhcauchy

  real(dp) function dhnorm(x,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::sigma;logical,intent(in),optional::log_p;real(dp)::s,d
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp)then;d=0.0_dp;else;d=2.0_dp*normal_pdf(x/s)/s;end if;v=apply_density_log(d,log_p)
  end function dhnorm
  real(dp) function phnorm(q,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::s,p
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;p=nan_dp();else if(q<0.0_dp)then;p=0.0_dp;else;p=2.0_dp*normal_cdf(q/s)-1.0_dp;end if;v=apply_tail(p,lower_tail,log_p)
  end function phnorm
  real(dp) function qhnorm(prob,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,s
    p=decode_probability(prob,lower_tail,log_p);s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;x=nan_dp();else;x=s*normal_quantile(0.5_dp*(p+1.0_dp));end if
  end function qhnorm
  function rhnorm(n,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::sigma;real(dp),allocatable::x(:);real(dp)::s;integer::i
    s=1.0_dp;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=abs(s*rnorm_std());end do
  end function rhnorm

  real(dp) function dht(x,nu,sigma,log_p) result(v)
    real(dp),intent(in)::x,nu;real(dp),intent(in),optional::sigma;logical,intent(in),optional::log_p;real(dp)::s,d
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp.or.nu<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp)then;d=0.0_dp;else;d=2.0_dp*student_pdf(x/s,nu)/s;end if;v=apply_density_log(d,log_p)
  end function dht
  real(dp) function pht(q,nu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,nu;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::s,p
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp.or.nu<=0.0_dp)then;p=nan_dp();else if(q<0.0_dp)then;p=0.0_dp;else;p=2.0_dp*student_cdf(q/s,nu)-1.0_dp;end if;v=apply_tail(p,lower_tail,log_p)
  end function pht
  real(dp) function qht(prob,nu,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,nu;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,s
    p=decode_probability(prob,lower_tail,log_p);s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp.or.nu<=0.0_dp)then;x=nan_dp();else;x=s*student_quantile(0.5_dp*(p+1.0_dp),nu);end if
  end function qht
  function rht(n,nu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in)::nu;real(dp),intent(in),optional::sigma;real(dp),allocatable::x(:);real(dp)::s;integer::i
    s=1.0_dp;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=abs(s*rt_scalar(nu));end do
  end function rht

  pure real(dp) function huber_norm(c) result(a)
    real(dp),intent(in)::c;a=2.0_dp*(normal_pdf(c)/c-normal_cdf(-c)+0.5_dp)
  end function huber_norm
  real(dp) function dhuber(x,mu,sigma,epsilon,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma,epsilon;logical,intent(in),optional::log_p
    real(dp)::m,s,c,z,rho,d,a;m=0.0_dp;s=1.0_dp;c=1.345_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(epsilon))c=epsilon
    if(s<=0.0_dp.or.c<=0.0_dp)then;d=nan_dp();else;z=abs((x-m)/s);if(z<=c)then;rho=0.5_dp*z*z;else;rho=c*z-0.5_dp*c*c;end if;a=sqrt_two_pi*huber_norm(c);d=exp(-rho)/(a*s);end if;v=apply_density_log(d,log_p)
  end function dhuber
  real(dp) function phuber(q,mu,sigma,epsilon,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma,epsilon;logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,c,z,az,p,a;m=0.0_dp;s=1.0_dp;c=1.345_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(epsilon))c=epsilon
    if(s<=0.0_dp.or.c<=0.0_dp)then;p=nan_dp();else;a=huber_norm(c);z=(q-m)/s;az=-abs(z);if(az<=-c)then;p=exp(0.5_dp*c*c+c*az)/(c*sqrt_two_pi*a);else;p=(normal_pdf(c)/c+normal_cdf(az)-normal_cdf(-c))/a;end if;if(z>0.0_dp)p=1.0_dp-p;end if;v=apply_tail(p,lower_tail,log_p)
  end function phuber
  real(dp) function qhuber(prob,mu,sigma,epsilon,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mu,sigma,epsilon;logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,m,s,c,pm,a,z;p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;c=1.345_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(epsilon))c=epsilon
    if(s<=0.0_dp.or.c<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();return;end if
    a=2.0_dp*sqrt_two_pi*(normal_cdf(c)+normal_pdf(c)/c-0.5_dp);pm=min(p,1.0_dp-p)
    if(pm<=sqrt_two_pi*normal_pdf(c)/(c*a))then;z=log(c*pm*a)/c-0.5_dp*c;else;z=normal_quantile(abs(1.0_dp-normal_cdf(c)+pm*a/sqrt_two_pi-normal_pdf(c)/c));end if
    if(p<0.5_dp)then;x=m+s*z;else;x=m-s*z;end if
  end function qhuber
  function rhuber(n,mu,sigma,epsilon) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma,epsilon;real(dp),allocatable::x(:);real(dp)::m,s,c;integer::i
    m=0.0_dp;s=1.0_dp;c=1.345_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(epsilon))c=epsilon;allocate(x(max(0,n)));do i=1,n;x(i)=qhuber(runif_open(),m,s,c);end do
  end function rhuber

  real(dp) function dinvgamma(x,alpha,beta,log_p) result(v)
    real(dp),intent(in)::x,alpha;real(dp),intent(in),optional::beta;logical,intent(in),optional::log_p;real(dp)::b,d
    b=1.0_dp;if(present(beta))b=beta
    if(alpha<=0.0_dp.or.b<=0.0_dp)then;d=nan_dp();else if(x<=0.0_dp)then;d=0.0_dp;else;d=exp(alpha*log(b)-log_gamma(alpha)-(alpha+1.0_dp)*log(x)-b/x);end if;v=apply_density_log(d,log_p)
  end function dinvgamma
  real(dp) function pinvgamma(q,alpha,beta,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,alpha;real(dp),intent(in),optional::beta;logical,intent(in),optional::lower_tail,log_p;real(dp)::b,p
    b=1.0_dp;if(present(beta))b=beta;if(alpha<=0.0_dp.or.b<=0.0_dp)then;p=nan_dp();else if(q<=0.0_dp)then;p=0.0_dp;else;p=1.0_dp-regularized_gamma_p(alpha,b/q);end if;v=apply_tail(p,lower_tail,log_p)
  end function pinvgamma
  real(dp) function qinvgamma(prob,alpha,beta,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,alpha;real(dp),intent(in),optional::beta;logical,intent(in),optional::lower_tail,log_p;real(dp)::b,p,g
    b=1.0_dp;if(present(beta))b=beta;p=decode_probability(prob,lower_tail,log_p);g=gamma_quantile(1.0_dp-p,alpha,1.0_dp);if(g==0.0_dp)then;x=pos_inf();else;x=b/g;end if
  end function qinvgamma
  function rinvgamma(n,alpha,beta) result(x)
    integer,intent(in)::n;real(dp),intent(in)::alpha;real(dp),intent(in),optional::beta;real(dp),allocatable::x(:);real(dp)::b;integer::i
    b=1.0_dp;if(present(beta))b=beta;allocate(x(max(0,n)));do i=1,n;x(i)=b/rgamma_scalar(alpha,1.0_dp);end do
  end function rinvgamma
  real(dp) function dinvchisq(x,nu,tau,log_p) result(v)
    real(dp),intent(in)::x,nu;real(dp),intent(in),optional::tau;logical,intent(in),optional::log_p;real(dp)::b
    b=0.5_dp;if(present(tau))b=0.5_dp*nu*tau;v=dinvgamma(x,0.5_dp*nu,b,log_p)
  end function dinvchisq
  real(dp) function pinvchisq(q,nu,tau,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,nu;real(dp),intent(in),optional::tau;logical,intent(in),optional::lower_tail,log_p;real(dp)::b
    b=0.5_dp;if(present(tau))b=0.5_dp*nu*tau;v=pinvgamma(q,0.5_dp*nu,b,lower_tail,log_p)
  end function pinvchisq
  real(dp) function qinvchisq(prob,nu,tau,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,nu;real(dp),intent(in),optional::tau;logical,intent(in),optional::lower_tail,log_p;real(dp)::b
    b=0.5_dp;if(present(tau))b=0.5_dp*nu*tau;x=qinvgamma(prob,0.5_dp*nu,b,lower_tail,log_p)
  end function qinvchisq
  function rinvchisq(n,nu,tau) result(x)
    integer,intent(in)::n;real(dp),intent(in)::nu;real(dp),intent(in),optional::tau;real(dp),allocatable::x(:);real(dp)::b;integer::i
    b=0.5_dp;if(present(tau))b=0.5_dp*nu*tau;allocate(x(max(0,n)));do i=1,n;x(i)=b/rgamma_scalar(0.5_dp*nu,1.0_dp);end do
  end function rinvchisq

  real(dp) function dkumar(x,a,b,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::a,b;logical,intent(in),optional::log_p;real(dp)::aa,bb,d
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;if(aa<=0.0_dp.or.bb<=0.0_dp)then;d=nan_dp();else if(x<=0.0_dp.or.x>=1.0_dp)then;d=0.0_dp;else;d=aa*bb*x**(aa-1.0_dp)*(1.0_dp-x**aa)**(bb-1.0_dp);end if;v=apply_density_log(d,log_p)
  end function dkumar
  real(dp) function pkumar(q,a,b,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::aa,bb,p
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;if(aa<=0.0_dp.or.bb<=0.0_dp)then;p=nan_dp();else if(q<=0.0_dp)then;p=0.0_dp;else if(q>=1.0_dp)then;p=1.0_dp;else;p=1.0_dp-(1.0_dp-q**aa)**bb;end if;v=apply_tail(p,lower_tail,log_p)
  end function pkumar
  real(dp) function qkumar(prob,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,aa,bb
    p=decode_probability(prob,lower_tail,log_p);aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;if(aa<=0.0_dp.or.bb<=0.0_dp)then;x=nan_dp();else;x=(1.0_dp-(1.0_dp-p)**(1.0_dp/bb))**(1.0_dp/aa);end if
  end function qkumar
  function rkumar(n,a,b) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::a,b;real(dp),allocatable::x(:);real(dp)::aa,bb;integer::i
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;allocate(x(max(0,n)));do i=1,n;x(i)=qkumar(runif_open(),aa,bb);end do
  end function rkumar

  real(dp) function dlaplace(x,mu,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::log_p;real(dp)::m,s,d
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(s<=0.0_dp)then;d=nan_dp();else;d=0.5_dp/s*exp(-abs(x-m)/s);end if;v=apply_density_log(d,log_p)
  end function dlaplace
  real(dp) function plaplace(q,mu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,z,p
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;z=(q-m)/s;if(s<=0.0_dp)then;p=nan_dp();else if(z<0.0_dp)then;p=0.5_dp*exp(z);else;p=1.0_dp-0.5_dp*exp(-z);end if;v=apply_tail(p,lower_tail,log_p)
  end function plaplace
  real(dp) function qlaplace(prob,mu,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,m,s
    p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(p<0.5_dp)then;x=m+s*log(2.0_dp*p);else;x=m-s*log(2.0_dp*(1.0_dp-p));end if
  end function qlaplace
  function rlaplace(n,mu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma;real(dp),allocatable::x(:);real(dp)::m,s;integer::i
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=qlaplace(runif_open(),m,s);end do
  end function rlaplace

  real(dp) function dlst(x,df,mu,sigma,log_p) result(v)
    real(dp),intent(in)::x,df;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::log_p;real(dp)::m,s,d
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(df<=0.0_dp.or.s<=0.0_dp)then;d=nan_dp();else;d=student_pdf((x-m)/s,df)/s;end if;v=apply_density_log(d,log_p)
  end function dlst
  real(dp) function plst(q,df,mu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,df;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,p
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;if(df<=0.0_dp.or.s<=0.0_dp)then;p=nan_dp();else;p=student_cdf((q-m)/s,df);end if;v=apply_tail(p,lower_tail,log_p)
  end function plst
  real(dp) function qlst(prob,df,mu,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,df;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,m,s
    p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;x=m+s*student_quantile(p,df)
  end function qlst
  function rlst(n,df,mu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in)::df;real(dp),intent(in),optional::mu,sigma;real(dp),allocatable::x(:);real(dp)::m,s;integer::i
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=m+s*rt_scalar(df);end do
  end function rlst

  real(dp) function dlomax(x,lambda,kappa,log_p) result(v)
    real(dp),intent(in)::x,lambda,kappa;logical,intent(in),optional::log_p;real(dp)::d
    if(lambda<=0.0_dp.or.kappa<=0.0_dp)then;d=nan_dp();else if(x<=0.0_dp)then;d=0.0_dp;else;d=lambda*kappa/(1.0_dp+lambda*x)**(kappa+1.0_dp);end if;v=apply_density_log(d,log_p)
  end function dlomax
  real(dp) function plomax(q,lambda,kappa,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,lambda,kappa;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    if(lambda<=0.0_dp.or.kappa<=0.0_dp)then;p=nan_dp();else if(q<=0.0_dp)then;p=0.0_dp;else;p=1.0_dp-(1.0_dp+lambda*q)**(-kappa);end if;v=apply_tail(p,lower_tail,log_p)
  end function plomax
  real(dp) function qlomax(prob,lambda,kappa,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,lambda,kappa;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(prob,lower_tail,log_p);x=((1.0_dp-p)**(-1.0_dp/kappa)-1.0_dp)/lambda
  end function qlomax
  function rlomax(n,lambda,kappa) result(x)
    integer,intent(in)::n;real(dp),intent(in)::lambda,kappa;real(dp),allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=qlomax(runif_open(),lambda,kappa);end do
  end function rlomax

  real(dp) function dnsbeta(x,shape1,shape2,min_value,max_value,log_p) result(v)
    real(dp),intent(in)::x,shape1,shape2;real(dp),intent(in),optional::min_value,max_value;logical,intent(in),optional::log_p;real(dp)::lo,hi,z,d
    lo=0.0_dp;hi=1.0_dp;if(present(min_value))lo=min_value;if(present(max_value))hi=max_value
    if(hi<=lo)then;d=nan_dp();else;z=(x-lo)/(hi-lo);d=beta_pdf(z,shape1,shape2)/(hi-lo);if(x<lo.or.x>hi)d=0.0_dp;end if;v=apply_density_log(d,log_p)
  end function dnsbeta
  real(dp) function pnsbeta(q,shape1,shape2,min_value,max_value,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,shape1,shape2;real(dp),intent(in),optional::min_value,max_value;logical,intent(in),optional::lower_tail,log_p;real(dp)::lo,hi,z,p
    lo=0.0_dp;hi=1.0_dp;if(present(min_value))lo=min_value;if(present(max_value))hi=max_value;z=(q-lo)/(hi-lo);p=regularized_beta(z,shape1,shape2);v=apply_tail(p,lower_tail,log_p)
  end function pnsbeta
  real(dp) function qnsbeta(prob,shape1,shape2,min_value,max_value,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,shape1,shape2;real(dp),intent(in),optional::min_value,max_value;logical,intent(in),optional::lower_tail,log_p;real(dp)::lo,hi,p
    lo=0.0_dp;hi=1.0_dp;if(present(min_value))lo=min_value;if(present(max_value))hi=max_value;p=decode_probability(prob,lower_tail,log_p);x=lo+(hi-lo)*beta_quantile(p,shape1,shape2)
  end function qnsbeta
  function rnsbeta(n,shape1,shape2,min_value,max_value) result(x)
    integer,intent(in)::n;real(dp),intent(in)::shape1,shape2;real(dp),intent(in),optional::min_value,max_value;real(dp),allocatable::x(:);real(dp)::lo,hi;integer::i
    lo=0.0_dp;hi=1.0_dp;if(present(min_value))lo=min_value;if(present(max_value))hi=max_value;allocate(x(max(0,n)));do i=1,n;x(i)=lo+(hi-lo)*rbeta_scalar(shape1,shape2);end do
  end function rnsbeta

  real(dp) function dpareto(x,a,b,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::a,b;logical,intent(in),optional::log_p;real(dp)::aa,bb,d
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;if(aa<=0.0_dp.or.bb<=0.0_dp)then;d=nan_dp();else if(x<bb)then;d=0.0_dp;else;d=aa*bb**aa/x**(aa+1.0_dp);end if;v=apply_density_log(d,log_p)
  end function dpareto
  real(dp) function ppareto(q,a,b,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::aa,bb,p
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;if(q<bb)then;p=0.0_dp;else;p=1.0_dp-(bb/q)**aa;end if;v=apply_tail(p,lower_tail,log_p)
  end function ppareto
  real(dp) function qpareto(prob,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,aa,bb
    p=decode_probability(prob,lower_tail,log_p);aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;x=bb/(1.0_dp-p)**(1.0_dp/aa)
  end function qpareto
  function rpareto(n,a,b) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::a,b;real(dp),allocatable::x(:);real(dp)::aa,bb;integer::i
    aa=1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;allocate(x(max(0,n)));do i=1,n;x(i)=bb/runif_open()**(1.0_dp/aa);end do
  end function rpareto

  real(dp) function dpower(x,alpha,beta,log_p) result(v)
    real(dp),intent(in)::x,alpha,beta;logical,intent(in),optional::log_p;real(dp)::d
    if(alpha<=0.0_dp.or.beta<=0.0_dp)then;d=nan_dp();else if(x<=0.0_dp.or.x>=alpha)then;d=0.0_dp;else;d=beta*x**(beta-1.0_dp)/alpha**beta;end if;v=apply_density_log(d,log_p)
  end function dpower
  real(dp) function ppower(q,alpha,beta,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,alpha,beta;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    if(q<=0.0_dp)then;p=0.0_dp;else if(q>=alpha)then;p=1.0_dp;else;p=(q/alpha)**beta;end if;v=apply_tail(p,lower_tail,log_p)
  end function ppower
  real(dp) function qpower(prob,alpha,beta,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,alpha,beta;logical,intent(in),optional::lower_tail,log_p;x=alpha*decode_probability(prob,lower_tail,log_p)**(1.0_dp/beta)
  end function qpower
  function rpower(n,alpha,beta) result(x)
    integer,intent(in)::n;real(dp),intent(in)::alpha,beta;real(dp),allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=alpha*runif_open()**(1.0_dp/beta);end do
  end function rpower

  real(dp) function dprop(x,size,mean,prior,log_p) result(v)
    real(dp),intent(in)::x,size,mean;real(dp),intent(in),optional::prior;logical,intent(in),optional::log_p;real(dp)::pr,d
    pr=0.0_dp;if(present(prior))pr=prior;d=beta_pdf(x,size*mean+pr,size*(1.0_dp-mean)+pr);v=apply_density_log(d,log_p)
  end function dprop
  real(dp) function pprop(q,size,mean,prior,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,size,mean;real(dp),intent(in),optional::prior;logical,intent(in),optional::lower_tail,log_p;real(dp)::pr,p
    pr=0.0_dp;if(present(prior))pr=prior;p=regularized_beta(q,size*mean+pr,size*(1.0_dp-mean)+pr);v=apply_tail(p,lower_tail,log_p)
  end function pprop
  real(dp) function qprop(prob,size,mean,prior,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,size,mean;real(dp),intent(in),optional::prior;logical,intent(in),optional::lower_tail,log_p;real(dp)::pr,p
    pr=0.0_dp;if(present(prior))pr=prior;p=decode_probability(prob,lower_tail,log_p);x=beta_quantile(p,size*mean+pr,size*(1.0_dp-mean)+pr)
  end function qprop
  function rprop(n,size,mean,prior) result(x)
    integer,intent(in)::n;real(dp),intent(in)::size,mean;real(dp),intent(in),optional::prior;real(dp),allocatable::x(:);real(dp)::pr;integer::i
    pr=0.0_dp;if(present(prior))pr=prior;allocate(x(max(0,n)));do i=1,n;x(i)=rbeta_scalar(size*mean+pr,size*(1.0_dp-mean)+pr);end do
  end function rprop

  real(dp) function drayleigh(x,sigma,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::sigma;logical,intent(in),optional::log_p;real(dp)::s,d
    s=1.0_dp;if(present(sigma))s=sigma;if(s<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp)then;d=0.0_dp;else;d=x/s**2*exp(-x*x/(2.0_dp*s*s));end if;v=apply_density_log(d,log_p)
  end function drayleigh
  real(dp) function prayleigh(q,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::s,p
    s=1.0_dp;if(present(sigma))s=sigma;if(q<0.0_dp)then;p=0.0_dp;else;p=1.0_dp-exp(-q*q/(2.0_dp*s*s));end if;v=apply_tail(p,lower_tail,log_p)
  end function prayleigh
  real(dp) function qrayleigh(prob,sigma,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,s
    p=decode_probability(prob,lower_tail,log_p);s=1.0_dp;if(present(sigma))s=sigma;x=s*sqrt(-2.0_dp*log(1.0_dp-p))
  end function qrayleigh
  function rrayleigh(n,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::sigma;real(dp),allocatable::x(:);real(dp)::s;integer::i
    s=1.0_dp;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=s*sqrt(-2.0_dp*log(runif_open()));end do
  end function rrayleigh

  real(dp) function dsgomp(x,b,eta,log_p) result(v)
    real(dp),intent(in)::x,b,eta;logical,intent(in),optional::log_p;real(dp)::e,d
    if(b<=0.0_dp.or.eta<=0.0_dp)then;d=nan_dp();else if(x<0.0_dp)then;d=0.0_dp;else;e=exp(-b*x);d=b*e*exp(-eta*e)*(1.0_dp+eta*(1.0_dp-e));end if;v=apply_density_log(d,log_p)
  end function dsgomp
  real(dp) function psgomp(q,b,eta,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,b,eta;logical,intent(in),optional::lower_tail,log_p;real(dp)::e,p
    if(q<0.0_dp)then;p=0.0_dp;else;e=exp(-b*q);p=(1.0_dp-e)*exp(-eta*e);end if;v=apply_tail(p,lower_tail,log_p)
  end function psgomp
  function rsgomp(n,b,eta) result(x)
    integer,intent(in)::n;real(dp),intent(in)::b,eta;real(dp),allocatable::x(:);real(dp)::rg,re;integer::i
    allocate(x(max(0,n)));do i=1,n;rg=-log((-log(runif_open()))/eta)/b;re=-log(runif_open())/b;x(i)=max(rg,re);end do
  end function rsgomp

  real(dp) function dslash(x,mu,sigma,log_p,source_compatible) result(v)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma
    logical,intent(in),optional::log_p,source_compatible
    real(dp)::m,s,z,d
    logical::compat
    m=0.0_dp;s=1.0_dp;compat=.true.
    if(present(mu))m=mu
    if(present(sigma))s=sigma
    if(present(source_compatible))compat=source_compatible
    if(s<=0.0_dp)then
      d=nan_dp()
    else
      z=(x-m)/s
      if(z==0.0_dp)then
        d=1.0_dp/(2.0_dp*sqrt_two_pi)
        if(.not.compat)d=d/s
      else
        d=(normal_pdf(0.0_dp)-normal_pdf(z))/(z*z*s)
      end if
    end if
    v=apply_density_log(d,log_p)
  end function dslash
  real(dp) function pslash(q,mu,sigma,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,z,p
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;z=(q-m)/s;if(z==0.0_dp)then;p=0.5_dp;else;p=normal_cdf(z)-(normal_pdf(0.0_dp)-normal_pdf(z))/z;end if;v=apply_tail(p,lower_tail,log_p)
  end function pslash
  function rslash(n,mu,sigma) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma;real(dp),allocatable::x(:);real(dp)::m,s;integer::i
    m=0.0_dp;s=1.0_dp;if(present(mu))m=mu;if(present(sigma))s=sigma;allocate(x(max(0,n)));do i=1,n;x(i)=m+s*rnorm_std()/runif_open();end do
  end function rslash

  real(dp) function dtriang(x,a,b,c,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::a,b,c;logical,intent(in),optional::log_p;real(dp)::aa,bb,cc,d
    aa=-1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;cc=0.5_dp*(aa+bb);if(present(c))cc=c
    if(aa>cc.or.cc>bb.or.aa==bb)then;d=nan_dp();else if(x<aa.or.x>bb)then;d=0.0_dp;else if(x<cc)then;d=2.0_dp*(x-aa)/((bb-aa)*(cc-aa));else if(x>cc)then;d=2.0_dp*(bb-x)/((bb-aa)*(bb-cc));else;d=2.0_dp/(bb-aa);end if;v=apply_density_log(d,log_p)
  end function dtriang
  real(dp) function ptriang(q,a,b,c,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::a,b,c;logical,intent(in),optional::lower_tail,log_p;real(dp)::aa,bb,cc,p
    aa=-1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;cc=0.5_dp*(aa+bb);if(present(c))cc=c
    if(q<aa)then;p=0.0_dp;else if(q>=bb)then;p=1.0_dp;else if(q<=cc)then;p=(q-aa)**2/((bb-aa)*(cc-aa));else;p=1.0_dp-(bb-q)**2/((bb-aa)*(bb-cc));end if;v=apply_tail(p,lower_tail,log_p)
  end function ptriang
  real(dp) function qtriang(prob,a,b,c,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::a,b,c;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,aa,bb,cc,fc
    p=decode_probability(prob,lower_tail,log_p);aa=-1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;cc=0.5_dp*(aa+bb);if(present(c))cc=c;fc=(cc-aa)/(bb-aa)
    if(p<fc)then;x=aa+sqrt(p*(bb-aa)*(cc-aa));else;x=bb-sqrt((1.0_dp-p)*(bb-aa)*(bb-cc));end if
  end function qtriang
  function rtriang(n,a,b,c) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::a,b,c;real(dp),allocatable::x(:);real(dp)::aa,bb,cc;integer::i
    aa=-1.0_dp;bb=1.0_dp;if(present(a))aa=a;if(present(b))bb=b;cc=0.5_dp*(aa+bb);if(present(c))cc=c;allocate(x(max(0,n)));do i=1,n;x(i)=qtriang(runif_open(),aa,bb,cc);end do
  end function rtriang

  real(dp) function dtnorm(x,mean,sd,a,b,log_p) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mean,sd,a,b;logical,intent(in),optional::log_p;real(dp)::m,s,lo,hi,z,d,den
    m=0.0_dp;s=1.0_dp;lo=neg_inf();hi=pos_inf();if(present(mean))m=mean;if(present(sd))s=sd;if(present(a))lo=a;if(present(b))hi=b
    den=normal_cdf((hi-m)/s)-normal_cdf((lo-m)/s);z=(x-m)/s;if(s<=0.0_dp.or.lo>hi)then;d=nan_dp();else if(x<lo.or.x>hi)then;d=0.0_dp;else;d=normal_pdf(z)/(s*den);end if;v=apply_density_log(d,log_p)
  end function dtnorm
  real(dp) function ptnorm(q,mean,sd,a,b,lower_tail,log_p) result(v)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mean,sd,a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::m,s,lo,hi,p,den
    m=0.0_dp;s=1.0_dp;lo=neg_inf();hi=pos_inf();if(present(mean))m=mean;if(present(sd))s=sd;if(present(a))lo=a;if(present(b))hi=b;den=normal_cdf((hi-m)/s)-normal_cdf((lo-m)/s)
    if(q<=lo)then;p=0.0_dp;else if(q>=hi)then;p=1.0_dp;else;p=(normal_cdf((q-m)/s)-normal_cdf((lo-m)/s))/den;end if;v=apply_tail(p,lower_tail,log_p)
  end function ptnorm
  real(dp) function qtnorm(prob,mean,sd,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;real(dp),intent(in),optional::mean,sd,a,b;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,m,s,lo,hi,plo,phi
    p=decode_probability(prob,lower_tail,log_p);m=0.0_dp;s=1.0_dp;lo=neg_inf();hi=pos_inf();if(present(mean))m=mean;if(present(sd))s=sd;if(present(a))lo=a;if(present(b))hi=b;plo=normal_cdf((lo-m)/s);phi=normal_cdf((hi-m)/s);x=m+s*normal_quantile(plo+p*(phi-plo))
  end function qtnorm
  function rtnorm(n,mean,sd,a,b) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mean,sd,a,b;real(dp),allocatable::x(:);real(dp)::m,s,lo,hi;integer::i
    m=0.0_dp;s=1.0_dp;lo=neg_inf();hi=pos_inf();if(present(mean))m=mean;if(present(sd))s=sd;if(present(a))lo=a;if(present(b))hi=b;allocate(x(max(0,n)));do i=1,n;x(i)=qtnorm(runif_open(),m,s,lo,hi);end do
  end function rtnorm

  real(dp) function qtlambda(prob,lambda,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,lambda;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(prob,lower_tail,log_p);if(lambda==0.0_dp)then;x=log(p)-log(1.0_dp-p);else;x=(p**lambda-(1.0_dp-p)**lambda)/lambda;end if
  end function qtlambda
  function rtlambda(n,lambda) result(x)
    integer,intent(in)::n;real(dp),intent(in)::lambda;real(dp),allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=qtlambda(runif_open(),lambda);end do
  end function rtlambda

  real(dp) function dwald(x,mu,lambda,log_p) result(v)
    real(dp),intent(in)::x,mu,lambda;logical,intent(in),optional::log_p;real(dp)::d
    if(mu<=0.0_dp.or.lambda<=0.0_dp)then;d=nan_dp();else if(x<=0.0_dp.or..not.ieee_is_finite(x))then;d=0.0_dp;else;d=sqrt(lambda/(2.0_dp*pi*x**3))*exp(-lambda*(x-mu)**2/(2.0_dp*mu**2*x));end if;v=apply_density_log(d,log_p)
  end function dwald
  real(dp) function pwald(q,mu,lambda,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,mu,lambda;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,t
    if(q<=0.0_dp)then;p=0.0_dp;else;t=sqrt(lambda/q);p=normal_cdf(t*(q/mu-1.0_dp))+exp(2.0_dp*lambda/mu)*normal_cdf(-t*(q/mu+1.0_dp));end if;v=apply_tail(p,lower_tail,log_p)
  end function pwald
  function rwald(n,mu,lambda) result(x)
    integer,intent(in)::n;real(dp),intent(in)::mu,lambda;real(dp),allocatable::x(:);real(dp)::z,y,candidate;integer::i
    allocate(x(max(0,n)));do i=1,n;z=rnorm_std();y=z*z;candidate=mu+mu*mu*y/(2.0_dp*lambda)-mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*y+mu*mu*y*y);if(runif_open()<=mu/(mu+candidate))then;x(i)=candidate;else;x(i)=mu*mu/candidate;end if;end do
  end function rwald

end module extra_distr_continuous
