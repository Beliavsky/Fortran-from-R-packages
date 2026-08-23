
module joker_distributions
  use joker_special
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: dbern,pbern,qbern,rbern, dbeta,pbeta,qbeta,rbeta
  public :: dbinom,pbinom,qbinom,rbinom, dcauchy,pcauchy,qcauchy,rcauchy
  public :: dchisq,pchisq,qchisq,rchisq, dexp_j,pexp_j,qexp_j,rexp_j
  public :: dfisher,pfisher,qfisher,rfisher, dgamma_j,pgamma_j,qgamma_j,rgamma_j
  public :: dgeom,pgeom,qgeom,rgeom, dlaplace,plaplace,qlaplace,rlaplace
  public :: dlnorm,plnorm,qlnorm,rlnorm, dnbinom,pnbinom,qnbinom,rnbinom
  public :: dnorm_j,pnorm_j,qnorm_j,rnorm_j, dpois,ppois,qpois,rpois
  public :: dt_j,pt_j,qt_j,rt_j, dunif_j,punif_j,qunif_j,runif_j
  public :: dweibull,pweibull,qweibull,rweibull
  public :: llbern,llbeta,llbinom,llcauchy,llchisq,llexp,llf,llgamma,llgeom
  public :: lllaplace,lllnorm,llnbinom,llnorm,llpois,llt,llunif,llweibull
contains
  pure real(dp) function dnorm_j(x,mu,sigma,logd) result(v)
    real(dp),intent(in)::x,mu,sigma
    logical,intent(in),optional::logd
    real(dp)::lv
    if(sigma<=0)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    lv=-0.5_dp*log(2*pi)-log(sigma)-0.5_dp*((x-mu)/sigma)**2
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    v=exp(lv)
  end function
  pure real(dp) function pnorm_j(x,mu,sigma) result(v)
    real(dp),intent(in)::x,mu,sigma
    v=normal_cdf((x-mu)/sigma)
  end function
  pure real(dp) function qnorm_j(p,mu,sigma) result(v)
    real(dp),intent(in)::p,mu,sigma
    v=mu+sigma*normal_quantile(p)
  end function
  real(dp) function rnorm_j(mu,sigma) result(v)
    real(dp),intent(in)::mu,sigma
    v=mu+sigma*rng_normal()
  end function

  pure real(dp) function dbern(x,prob,logd) result(v)
    integer,intent(in)::x;real(dp),intent(in)::prob;logical,intent(in),optional::logd
    real(dp)::lv
    if(prob<0.or.prob>1.or.(x/=0.and.x/=1))then
      lv=-huge(1.0_dp)
    else
      if(x==1)then;lv=log(max(prob,tiny(1.0_dp)));else;lv=log(max(1-prob,tiny(1.0_dp)));end if
    end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<=-huge(1.0_dp)/2)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pbern(q,prob) result(v)
    real(dp),intent(in)::q,prob
    if(q<0)then;v=0;else if(q<1)then;v=1-prob;else;v=1;end if
  end function
  pure integer function qbern(p,prob) result(v)
    real(dp),intent(in)::p,prob
    if(p<=1-prob)then;v=0;else;v=1;end if
  end function
  integer function rbern(prob) result(v)
    real(dp),intent(in)::prob;real(dp)::u
    call random_number(u);v=merge(1,0,u<prob)
  end function

  pure real(dp) function dbeta(x,a,b,logd) result(v)
    real(dp),intent(in)::x,a,b;logical,intent(in),optional::logd;real(dp)::lv
    if(x<=0.or.x>=1.or.a<=0.or.b<=0)then;lv=-huge(1.0_dp)
    else;lv=(a-1)*log(x)+(b-1)*log(1-x)-log_beta(a,b);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pbeta(x,a,b) result(v)
    real(dp),intent(in)::x,a,b;v=reg_beta(x,a,b)
  end function
  pure real(dp) function qbeta(p,a,b) result(v)
    real(dp),intent(in)::p,a,b;v=beta_quantile(p,a,b)
  end function
  real(dp) function rbeta(a,b) result(v)
    real(dp),intent(in)::a,b;v=rng_beta(a,b)
  end function

  pure real(dp) function logchoose(n,k) result(v)
    integer,intent(in)::n,k
    if(k<0.or.k>n)then;v=-huge(1.0_dp);else;v=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp));end if
  end function
  pure real(dp) function dbinom(x,n,prob,logd) result(v)
    integer,intent(in)::x,n;real(dp),intent(in)::prob;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.x>n.or.prob<0.or.prob>1)then;lv=-huge(1.0_dp)
    else if(prob<=0)then;lv=merge(0.0_dp,-huge(1.0_dp),x==0)
    else if(prob>=1)then;lv=merge(0.0_dp,-huge(1.0_dp),x==n)
    else;lv=logchoose(n,x)+x*log(prob)+(n-x)*log(1-prob);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pbinom(q,n,prob) result(v)
    integer,intent(in)::q,n;real(dp),intent(in)::prob;integer::k
    if(q<0)then;v=0;return;end if
    if(q>=n)then;v=1;return;end if
    v=0;do k=0,q;v=v+dbinom(k,n,prob);end do;v=min(1.0_dp,v)
  end function
  pure integer function qbinom(p,n,prob) result(v)
    real(dp),intent(in)::p,prob;integer,intent(in)::n;integer::k
    do k=0,n;if(pbinom(k,n,prob)>=p)then;v=k;return;end if;end do;v=n
  end function
  integer function rbinom(n,prob) result(v)
    integer,intent(in)::n;real(dp),intent(in)::prob;integer::i
    v=0;do i=1,n;v=v+rbern(prob);end do
  end function

  pure real(dp) function dcauchy(x,loc,scale,logd) result(v)
    real(dp),intent(in)::x,loc,scale;logical,intent(in),optional::logd;real(dp)::lv
    lv=-log(pi)-log(scale)-log(1+((x-loc)/scale)**2)
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;v=exp(lv)
  end function
  pure real(dp) function pcauchy(x,loc,scale) result(v)
    real(dp),intent(in)::x,loc,scale;v=0.5_dp+atan((x-loc)/scale)/pi
  end function
  pure real(dp) function qcauchy(p,loc,scale) result(v)
    real(dp),intent(in)::p,loc,scale;v=loc+scale*tan(pi*(p-0.5_dp))
  end function
  real(dp) function rcauchy(loc,scale) result(v)
    real(dp),intent(in)::loc,scale;real(dp)::u;call random_number(u);v=qcauchy(u,loc,scale)
  end function

  pure real(dp) function dgamma_j(x,shape,scale,logd) result(v)
    real(dp),intent(in)::x,shape,scale;logical,intent(in),optional::logd;real(dp)::lv
    if(x<=0.or.shape<=0.or.scale<=0)then;lv=-huge(1.0_dp)
    else;lv=(shape-1)*log(x)-x/scale-log_gamma(shape)-shape*log(scale);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pgamma_j(x,shape,scale) result(v)
    real(dp),intent(in)::x,shape,scale
    if(x<=0)then;v=0;else;v=reg_gamma_p(shape,x/scale);end if
  end function
  pure real(dp) function qgamma_j(p,shape,scale) result(v)
    real(dp),intent(in)::p,shape,scale;real(dp)::lo,hi,mid;integer::it
    if(p<=0)then;v=0;return;end if;if(p>=1)then;v=huge(1.0_dp);return;end if
    lo=0;hi=max(scale,shape*scale)
    do while(pgamma_j(hi,shape,scale)<p);hi=2*hi;end do
    do it=1,120;mid=.5_dp*(lo+hi);if(pgamma_j(mid,shape,scale)<p)then;lo=mid;else;hi=mid;end if;end do
    v=.5_dp*(lo+hi)
  end function
  real(dp) function rgamma_j(shape,scale) result(v)
    real(dp),intent(in)::shape,scale;v=rng_gamma(shape,scale)
  end function

  pure real(dp) function dchisq(x,df,logd) result(v)
    real(dp),intent(in)::x,df;logical,intent(in),optional::logd
    if(present(logd))then;v=dgamma_j(x,df/2,2.0_dp,logd);else;v=dgamma_j(x,df/2,2.0_dp);end if
  end function
  pure real(dp) function pchisq(x,df) result(v);real(dp),intent(in)::x,df;v=pgamma_j(x,df/2,2.0_dp);end function
  pure real(dp) function qchisq(p,df) result(v);real(dp),intent(in)::p,df;v=qgamma_j(p,df/2,2.0_dp);end function
  real(dp) function rchisq(df) result(v);real(dp),intent(in)::df;v=rgamma_j(df/2,2.0_dp);end function

  pure real(dp) function dexp_j(x,rate,logd) result(v)
    real(dp),intent(in)::x,rate;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.rate<=0)then;lv=-huge(1.0_dp);else;lv=log(rate)-rate*x;end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pexp_j(x,rate) result(v)
  real(dp),intent(in)::x,rate
  if(x<=0)then
  v=0
  else
  v=1-exp(-rate*x)
  end if
  end function
  pure real(dp) function qexp_j(p,rate) result(v);real(dp),intent(in)::p,rate;v=-log(1-p)/rate;end function
  real(dp) function rexp_j(rate) result(v)
  real(dp),intent(in)::rate
  real(dp)::u
  call random_number(u)
  v=-log(max(u,tiny(1.0_dp)))/rate
  end function

  pure real(dp) function dfisher(x,df1,df2,logd) result(v)
    real(dp),intent(in)::x,df1,df2;logical,intent(in),optional::logd;real(dp)::lv,a,b
    a=df1/2;b=df2/2
    if(x<=0)then;lv=-huge(1.0_dp);else;lv=a*log(df1/df2)+(a-1)*log(x)-(a+b)*log(1+df1*x/df2)-log_beta(a,b);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if
    if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pfisher(x,df1,df2) result(v)
  real(dp),intent(in)::x,df1,df2
  if(x<=0)then
  v=0
  else
  v=reg_beta(df1*x/(df1*x+df2),df1/2,df2/2)
  end if
  end function
  pure real(dp) function qfisher(p,df1,df2) result(v)
  real(dp),intent(in)::p,df1,df2
  real(dp)::z
  z=beta_quantile(p,df1/2,df2/2)
  v=df2*z/(df1*(1-z))
  end function
  real(dp) function rfisher(df1,df2) result(v);real(dp),intent(in)::df1,df2;v=(rchisq(df1)/df1)/(rchisq(df2)/df2);end function

  pure real(dp) function dgeom(x,prob,logd) result(v)
    integer,intent(in)::x;real(dp),intent(in)::prob;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.prob<=0.or.prob>1)then;lv=-huge(1.0_dp);else;lv=log(prob)+x*log(max(1-prob,tiny(1.0_dp)));end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pgeom(q,prob) result(v)
  integer,intent(in)::q
  real(dp),intent(in)::prob
  if(q<0)then
  v=0
  else
  v=1-(1-prob)**(q+1)
  end if
  end function
  pure integer function qgeom(p,prob) result(v)
  real(dp),intent(in)::p,prob
  v=max(0,ceiling(log(1-p)/log(1-prob)-1.0_dp))
  end function
  integer function rgeom(prob) result(v);real(dp),intent(in)::prob;real(dp)::u;call random_number(u);v=qgeom(u,prob);end function

  pure real(dp) function dlaplace(x,mu,sigma,logd) result(v)
    real(dp),intent(in)::x,mu,sigma;logical,intent(in),optional::logd;real(dp)::lv
    lv=-log(2*sigma)-abs(x-mu)/sigma;if(present(logd))then;if(logd)then;v=lv;return;end if;end if;v=exp(lv)
  end function
  pure real(dp) function plaplace(x,mu,sigma) result(v)
    real(dp),intent(in)::x,mu,sigma
    if(x<mu)then;v=.5_dp*exp((x-mu)/sigma);else;v=1-.5_dp*exp(-(x-mu)/sigma);end if
  end function
  pure real(dp) function qlaplace(p,mu,sigma) result(v)
    real(dp),intent(in)::p,mu,sigma
    if(p<.5_dp)then;v=mu+sigma*log(2*p);else;v=mu-sigma*log(2*(1-p));end if
  end function
  real(dp) function rlaplace(mu,sigma) result(v)
  real(dp),intent(in)::mu,sigma
  real(dp)::u
  call random_number(u)
  v=qlaplace(u,mu,sigma)
  end function

  pure real(dp) function dlnorm(x,meanlog,sdlog,logd) result(v)
    real(dp),intent(in)::x,meanlog,sdlog;logical,intent(in),optional::logd;real(dp)::lv
    if(x<=0)then;lv=-huge(1.0_dp);else;lv=dnorm_j(log(x),meanlog,sdlog,.true.)-log(x);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function plnorm(x,meanlog,sdlog) result(v)
  real(dp),intent(in)::x,meanlog,sdlog
  if(x<=0)then
  v=0
  else
  v=pnorm_j(log(x),meanlog,sdlog)
  end if
  end function
  pure real(dp) function qlnorm(p,meanlog,sdlog) result(v)
  real(dp),intent(in)::p,meanlog,sdlog
  v=exp(qnorm_j(p,meanlog,sdlog))
  end function
  real(dp) function rlnorm(meanlog,sdlog) result(v);real(dp),intent(in)::meanlog,sdlog;v=exp(rnorm_j(meanlog,sdlog));end function

  pure real(dp) function dnbinom(x,size,prob,logd) result(v)
    integer,intent(in)::x;real(dp),intent(in)::size,prob;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.size<=0.or.prob<=0.or.prob>1)then;lv=-huge(1.0_dp)
    else;lv=log_gamma(x+size)-log_gamma(size)-log_gamma(real(x+1,dp))+size*log(prob)+x*log(max(1-prob,tiny(1.0_dp)));end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pnbinom(q,size,prob) result(v)
  integer,intent(in)::q
  real(dp),intent(in)::size,prob
  if(q<0)then
  v=0
  else
  v=reg_beta(prob,size,real(q+1,dp))
  end if
  end function
  pure integer function qnbinom(p,size,prob) result(v)
  real(dp),intent(in)::p,size,prob
  integer::k
  k=0
  do while(pnbinom(k,size,prob)<p.and.k<100000)
  k=k+1
  end do
  v=k
  end function
  integer function rnbinom(size,prob) result(v)
  real(dp),intent(in)::size,prob
  real(dp)::lam
  lam=rng_gamma(size,(1-prob)/prob)
  v=rpois(lam)
  end function

  pure real(dp) function dpois(x,lambda,logd) result(v)
    integer,intent(in)::x;real(dp),intent(in)::lambda;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.lambda<0)then;lv=-huge(1.0_dp)
    else if(lambda<=0)then;lv=merge(0.0_dp,-huge(1.0_dp),x==0)
    else;lv=-lambda+x*log(lambda)-log_gamma(real(x+1,dp));end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function ppois(q,lambda) result(v)
  integer,intent(in)::q
  real(dp),intent(in)::lambda
  integer::k
  if(q<0)then
  v=0
  return
  end if
  v=0
  do k=0,q
  v=v+dpois(k,lambda)
  end do
  v=min(1.0_dp,v)
  end function
  pure integer function qpois(p,lambda) result(v)
  real(dp),intent(in)::p,lambda
  integer::k
  k=0
  do while(ppois(k,lambda)<p.and.k<100000)
  k=k+1
  end do
  v=k
  end function
  integer function rpois(lambda) result(v)
    real(dp),intent(in)::lambda;real(dp)::l,p
    if(lambda<30)then;l=exp(-lambda);p=1;v=-1;do;p=p*randu();v=v+1;if(p<=l)exit;end do
    else;v=max(0,nint(lambda+sqrt(lambda)*rng_normal()));end if
  contains
    real(dp) function randu() result(z);call random_number(z);end function
  end function

  pure real(dp) function dt_j(x,df,logd) result(v)
    real(dp),intent(in)::x,df;logical,intent(in),optional::logd;real(dp)::lv
    lv=log_gamma((df+1)/2)-log_gamma(df/2)-.5_dp*log(df*pi)-.5_dp*(df+1)*log(1+x*x/df)
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;v=exp(lv)
  end function
  pure real(dp) function pt_j(x,df) result(v)
    real(dp),intent(in)::x,df;real(dp)::b
    b=reg_beta(df/(df+x*x),df/2,.5_dp)
    if(x>=0)then;v=1-.5_dp*b;else;v=.5_dp*b;end if
  end function
  pure real(dp) function qt_j(p,df) result(v)
    real(dp),intent(in)::p,df;real(dp)::lo,hi,mid;integer::it
    lo=-100;hi=100;do it=1,120;mid=.5_dp*(lo+hi);if(pt_j(mid,df)<p)then;lo=mid;else;hi=mid;end if;end do;v=.5_dp*(lo+hi)
  end function
  real(dp) function rt_j(df) result(v);real(dp),intent(in)::df;v=rng_normal()/sqrt(rchisq(df)/df);end function

  pure real(dp) function dunif_j(x,a,b,logd) result(v)
    real(dp),intent(in)::x,a,b;logical,intent(in),optional::logd;real(dp)::lv
    if(x<a.or.x>b.or.b<=a)then;lv=-huge(1.0_dp);else;lv=-log(b-a);end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function punif_j(x,a,b) result(v);real(dp),intent(in)::x,a,b;v=max(0.0_dp,min(1.0_dp,(x-a)/(b-a)));end function
  pure real(dp) function qunif_j(p,a,b) result(v);real(dp),intent(in)::p,a,b;v=a+p*(b-a);end function
  real(dp) function runif_j(a,b) result(v);real(dp),intent(in)::a,b;real(dp)::u;call random_number(u);v=a+u*(b-a);end function

  pure real(dp) function dweibull(x,shape,scale,logd) result(v)
    real(dp),intent(in)::x,shape,scale;logical,intent(in),optional::logd;real(dp)::lv
    if(x<0.or.shape<=0.or.scale<=0)then
    lv=-huge(1.0_dp)
    else
    lv=log(shape/scale)+(shape-1)*log(max(x/scale,tiny(1.0_dp)))-(x/scale)**shape
    end if
    if(present(logd))then;if(logd)then;v=lv;return;end if;end if;if(lv<-700)then;v=0;else;v=exp(lv);end if
  end function
  pure real(dp) function pweibull(x,shape,scale) result(v)
  real(dp),intent(in)::x,shape,scale
  if(x<=0)then
  v=0
  else
  v=1-exp(-(x/scale)**shape)
  end if
  end function
  pure real(dp) function qweibull(p,shape,scale) result(v)
  real(dp),intent(in)::p,shape,scale
  v=scale*(-log(1-p))**(1/shape)
  end function
  real(dp) function rweibull(shape,scale) result(v)
  real(dp),intent(in)::shape,scale
  real(dp)::u
  call random_number(u)
  v=qweibull(u,shape,scale)
  end function

  pure real(dp) function llbern(x,prob) result(v)
  integer,intent(in)::x(:)
  real(dp),intent(in)::prob
  integer::i
  v=0
  do i=1,size(x)
  v=v+dbern(x(i),prob,.true.)
  end do
  end function
  pure real(dp) function llbeta(x,a,b) result(v)
  real(dp),intent(in)::x(:),a,b
  integer::i
  v=0
  do i=1,size(x)
  v=v+dbeta(x(i),a,b,.true.)
  end do
  end function
  pure real(dp) function llbinom(x,n,prob) result(v)
  integer,intent(in)::x(:),n
  real(dp),intent(in)::prob
  integer::i
  v=0
  do i=1,size(x)
  v=v+dbinom(x(i),n,prob,.true.)
  end do
  end function
  pure real(dp) function llcauchy(x,loc,scale) result(v)
  real(dp),intent(in)::x(:),loc,scale
  integer::i
  v=0
  do i=1,size(x)
  v=v+dcauchy(x(i),loc,scale,.true.)
  end do
  end function
  pure real(dp) function llchisq(x,df) result(v)
  real(dp),intent(in)::x(:),df
  integer::i
  v=0
  do i=1,size(x)
  v=v+dchisq(x(i),df,.true.)
  end do
  end function
  pure real(dp) function llexp(x,rate) result(v)
  real(dp),intent(in)::x(:),rate
  integer::i
  v=0
  do i=1,size(x)
  v=v+dexp_j(x(i),rate,.true.)
  end do
  end function
  pure real(dp) function llf(x,df1,df2) result(v)
  real(dp),intent(in)::x(:),df1,df2
  integer::i
  v=0
  do i=1,size(x)
  v=v+dfisher(x(i),df1,df2,.true.)
  end do
  end function
  pure real(dp) function llgamma(x,shape,scale) result(v)
  real(dp),intent(in)::x(:),shape,scale
  integer::i
  v=0
  do i=1,size(x)
  v=v+dgamma_j(x(i),shape,scale,.true.)
  end do
  end function
  pure real(dp) function llgeom(x,prob) result(v)
  integer,intent(in)::x(:)
  real(dp),intent(in)::prob
  integer::i
  v=0
  do i=1,size(x)
  v=v+dgeom(x(i),prob,.true.)
  end do
  end function
  pure real(dp) function lllaplace(x,mu,sigma) result(v)
  real(dp),intent(in)::x(:),mu,sigma
  integer::i
  v=0
  do i=1,size(x)
  v=v+dlaplace(x(i),mu,sigma,.true.)
  end do
  end function
  pure real(dp) function lllnorm(x,mu,sigma) result(v)
  real(dp),intent(in)::x(:),mu,sigma
  integer::i
  v=0
  do i=1,size(x)
  v=v+dlnorm(x(i),mu,sigma,.true.)
  end do
  end function
  pure real(dp) function llnbinom(x,sizep,prob) result(v)
  integer,intent(in)::x(:)
  real(dp),intent(in)::sizep,prob
  integer::i
  v=0
  do i=1,size(x)
  v=v+dnbinom(x(i),sizep,prob,.true.)
  end do
  end function
  pure real(dp) function llnorm(x,mu,sigma) result(v)
  real(dp),intent(in)::x(:),mu,sigma
  integer::i
  v=0
  do i=1,size(x)
  v=v+dnorm_j(x(i),mu,sigma,.true.)
  end do
  end function
  pure real(dp) function llpois(x,lambda) result(v)
  integer,intent(in)::x(:)
  real(dp),intent(in)::lambda
  integer::i
  v=0
  do i=1,size(x)
  v=v+dpois(x(i),lambda,.true.)
  end do
  end function
  pure real(dp) function llt(x,df) result(v)
  real(dp),intent(in)::x(:),df
  integer::i
  v=0
  do i=1,size(x)
  v=v+dt_j(x(i),df,.true.)
  end do
  end function
  pure real(dp) function llunif(x,a,b) result(v)
  real(dp),intent(in)::x(:),a,b
  integer::i
  v=0
  do i=1,size(x)
  v=v+dunif_j(x(i),a,b,.true.)
  end do
  end function
  pure real(dp) function llweibull(x,shape,scale) result(v)
  real(dp),intent(in)::x(:),shape,scale
  integer::i
  v=0
  do i=1,size(x)
  v=v+dweibull(x(i),shape,scale,.true.)
  end do
  end function
end module joker_distributions
