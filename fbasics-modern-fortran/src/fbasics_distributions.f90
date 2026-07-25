! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_distributions
  use fbasics_kinds, only: dp, pi, clamp
  use fbasics_rng, only: runif_lcg, rnorm_lcg, rt_lcg, random_inverse_gaussian
  use fbasics_special, only: normal_pdf, normal_cdf, normal_quantile, student_pdf, &
    student_cdf, student_quantile, bessel_k_nu, adaptive_simpson
  use fbasics_stats, only: sample_mean, sample_sd, sample_variance, sample_skewness, &
    sample_kurtosis, sample_quantile
  use fbasics_optimize, only: nelder_mead_bounded, numerical_hessian
  use fbasics_linalg, only: matrix_inverse
  implicit none
  private
  type, public :: distribution_fit
    character(len=16) :: family = ''
    real(dp), allocatable :: parameters(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    logical :: converged = .false.
    real(dp), allocatable :: hessian(:,:), covariance(:,:)
  end type
  public :: dnorm_fs, pnorm_fs, qnorm_fs, rnorm_fs
  public :: dt_fs, pt_fs, qt_fs, rt_fs
  public :: dnig, pnig, qnig, rnig, nig_mean, nig_variance, nig_skewness, nig_kurtosis
  public :: dgh, pgh, qgh, rgh, gh_raw_moment, gh_mean, gh_variance, gh_skewness, gh_kurtosis
  public :: standardized_gh_parameters, dsgh, psgh, qsgh, rsgh
  public :: dhyp, phyp, qhyp, rhyp, dght, pght, qght, rght
  public :: qgld_rs, pgld_rs, dgld_rs, rgld_rs, gld_mode
  public :: fit_normal, fit_student, fit_nig, fit_gld_quantiles
  real(dp), allocatable, save :: active_fit_data(:), active_qobs(:), active_probs(:)
contains
  pure elemental real(dp) function dnorm_fs(x,mean,sd,log_density) result(v)
    real(dp),intent(in)::x,mean,sd
    logical,intent(in),optional::log_density
    logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(sd<=0.0_dp)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    v=-log(sd)-0.5_dp*log(2.0_dp*pi)-0.5_dp*((x-mean)/sd)**2
    if(.not.lg)v=exp(v)
  end function
  pure elemental real(dp) function pnorm_fs(x,mean,sd) result(v)
    real(dp),intent(in)::x,mean,sd
    if(sd>0.0_dp)then;v=normal_cdf((x-mean)/sd);else;v=merge(0.0_dp,1.0_dp,x<mean);end if
  end function
  pure elemental real(dp) function qnorm_fs(p,mean,sd) result(v)
    real(dp),intent(in)::p,mean,sd;v=mean+sd*normal_quantile(p)
  end function
  real(dp) function rnorm_fs(mean,sd) result(v)
    real(dp),intent(in)::mean,sd;v=mean+sd*rnorm_lcg()
  end function

  pure elemental real(dp) function dt_fs(x,df,mean,sd,log_density) result(v)
    real(dp),intent(in)::x,df,mean,sd;logical,intent(in),optional::log_density;logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(sd<=0.0_dp.or.df<=0.0_dp)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    v=log(student_pdf((x-mean)/sd,df))-log(sd);if(.not.lg)v=exp(v)
  end function
  real(dp) function pt_fs(x,df,mean,sd) result(v)
    real(dp),intent(in)::x,df,mean,sd;v=student_cdf((x-mean)/sd,df)
  end function
  real(dp) function qt_fs(p,df,mean,sd) result(v)
    real(dp),intent(in)::p,df,mean,sd;v=mean+sd*student_quantile(p,df)
  end function
  real(dp) function rt_fs(df,mean,sd) result(v)
    real(dp),intent(in)::df,mean,sd;v=mean+sd*rt_lcg(df)
  end function

  real(dp) function dnig(x,alpha,beta,delta,mu,log_density) result(v)
    real(dp),intent(in)::x,alpha,beta,delta,mu;logical,intent(in),optional::log_density
    real(dp)::gamma,s,arg,k,lv;logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(alpha<=0.0_dp.or.delta<=0.0_dp.or.abs(beta)>=alpha)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    gamma=sqrt(alpha*alpha-beta*beta);s=sqrt(delta*delta+(x-mu)**2);arg=alpha*s;k=bessel_k_nu(1.0_dp,arg)
    if(k<=0.0_dp)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    lv=delta*gamma+log(delta*alpha/pi)-log(s)+log(k)+beta*(x-mu)
    v=merge(lv,exp(lv),lg)
  end function
  pure real(dp) function nig_mean(alpha,beta,delta,mu) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu;v=mu+delta*beta/sqrt(alpha*alpha-beta*beta)
  end function
  pure real(dp) function nig_variance(alpha,beta,delta) result(v)
    real(dp),intent(in)::alpha,beta,delta;real(dp)::g;g=sqrt(alpha*alpha-beta*beta);v=delta*alpha*alpha/g**3
  end function
  pure real(dp) function nig_skewness(alpha,beta,delta) result(v)
    real(dp),intent(in)::alpha,beta,delta;real(dp)::g;g=sqrt(alpha*alpha-beta*beta);v=3.0_dp*beta/(alpha*sqrt(delta*g))
  end function
  pure real(dp) function nig_kurtosis(alpha,beta,delta) result(v)
    real(dp),intent(in)::alpha,beta,delta;real(dp)::g;g=sqrt(alpha*alpha-beta*beta);v=3.0_dp*(1.0_dp+4.0_dp*beta*beta/(alpha*alpha))/(delta*g)
  end function
  real(dp) function pnig(q,alpha,beta,delta,mu) result(v)
    real(dp),intent(in)::q,alpha,beta,delta,mu
    real(dp)::m,s,lo
    m=nig_mean(alpha,beta,delta,mu)
    s=sqrt(nig_variance(alpha,beta,delta))
    lo=min(q,m-18.0_dp*s-10.0_dp*delta)
    if(q<=lo+epsilon(1.0_dp))then
      v=0.0_dp
    else
      v=integrate_nig_density(lo,q,alpha,beta,delta,mu)
    end if
    v=clamp(v,0.0_dp,1.0_dp)
  end function

  real(dp) function integrate_nig_density(a0,b0,alpha,beta,delta,mu) result(v)
    real(dp),intent(in)::a0,b0,alpha,beta,delta,mu
    integer::n
    real(dp)::oldv,newv
    n=64
    oldv=nig_simpson(a0,b0,alpha,beta,delta,mu,n)
    do while(n<16384)
      n=2*n
      newv=nig_simpson(a0,b0,alpha,beta,delta,mu,n)
      if(abs(newv-oldv)<=1.0e-8_dp*max(1.0_dp,abs(newv)))exit
      oldv=newv
    end do
    v=newv
  end function

  real(dp) function nig_simpson(a0,b0,alpha,beta,delta,mu,n) result(v)
    real(dp),intent(in)::a0,b0,alpha,beta,delta,mu
    integer,intent(in)::n
    real(dp)::h,x0
    integer::i
    h=(b0-a0)/real(n,dp)
    v=dnig(a0,alpha,beta,delta,mu)+dnig(b0,alpha,beta,delta,mu)
    do i=1,n-1
      x0=a0+h*real(i,dp)
      if(mod(i,2)==0)then
        v=v+2.0_dp*dnig(x0,alpha,beta,delta,mu)
      else
        v=v+4.0_dp*dnig(x0,alpha,beta,delta,mu)
      end if
    end do
    v=v*h/3.0_dp
  end function

  real(dp) function qnig(p,alpha,beta,delta,mu) result(v)
    real(dp),intent(in)::p,alpha,beta,delta,mu;real(dp)::lo,hi,mid,s;integer::i
    if(p<=0)then;v=-huge(1.0_dp);return;else if(p>=1)then;v=huge(1.0_dp);return;end if
    s=sqrt(nig_variance(alpha,beta,delta));lo=nig_mean(alpha,beta,delta,mu)-4*s;hi=nig_mean(alpha,beta,delta,mu)+4*s
    do while(pnig(lo,alpha,beta,delta,mu)>p);lo=lo-2*s;end do;do while(pnig(hi,alpha,beta,delta,mu)<p);hi=hi+2*s;end do
    do i=1,80;mid=0.5_dp*(lo+hi);if(pnig(mid,alpha,beta,delta,mu)<p)then;lo=mid;else;hi=mid;end if;end do;v=0.5_dp*(lo+hi)
  end function
  real(dp) function rnig(alpha,beta,delta,mu) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu;real(dp)::g,w
    g=sqrt(alpha*alpha-beta*beta);w=random_inverse_gaussian(delta/g,delta*delta);v=mu+beta*w+sqrt(w)*rnorm_lcg()
  end function

  real(dp) function gig_moment(r,lambda,chi,psi) result(v)
    real(dp),intent(in)::r,lambda,chi,psi;real(dp)::z
    z=sqrt(chi*psi);v=(chi/psi)**(0.5_dp*r)*bessel_k_nu(lambda+r,z)/bessel_k_nu(lambda,z)
  end function
  real(dp) function gh_raw_moment(order,alpha,beta,delta,mu,lambda) result(v)
    integer,intent(in)::order;real(dp),intent(in)::alpha,beta,delta,mu,lambda
    integer::k,j;real(dp)::centered,term,gamma2
    gamma2=alpha*alpha-beta*beta;v=0.0_dp
    do k=0,order
      centered=0.0_dp
      do j=0,k/2
        term=binomial(k,2*j)*beta**(k-2*j)*double_factorial(2*j-1)*gig_moment(real(k-j,dp),lambda,delta*delta,gamma2)
        centered=centered+term
      end do
      v=v+binomial(order,k)*mu**(order-k)*centered
    end do
  contains
    pure real(dp) function binomial(n,k0) result(b);integer,intent(in)::n,k0;integer::i,kk;b=1;kk=min(k0,n-k0);do i=1,kk;b=b*real(n-kk+i,dp)/real(i,dp);end do;end function
    pure real(dp) function double_factorial(n) result(d)
      integer,intent(in)::n
      integer::i
      d=1.0_dp
      if(n>0) then
        do i=n,1,-2
          d=d*real(i,dp)
        end do
      end if
    end function
  end function
  real(dp) function gh_mean(alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;v=gh_raw_moment(1,alpha,beta,delta,mu,lambda)
  end function
  real(dp) function gh_variance(alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;real(dp)::m;m=gh_mean(alpha,beta,delta,mu,lambda);v=gh_raw_moment(2,alpha,beta,delta,mu,lambda)-m*m
  end function
  real(dp) function gh_skewness(alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;real(dp)::m,r2,r3,var
    m=gh_mean(alpha,beta,delta,mu,lambda);r2=gh_raw_moment(2,alpha,beta,delta,mu,lambda);r3=gh_raw_moment(3,alpha,beta,delta,mu,lambda);var=r2-m*m
    v=(r3-3*m*r2+2*m**3)/var**1.5_dp
  end function
  real(dp) function gh_kurtosis(alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;real(dp)::m,r2,r3,r4,var
    m=gh_mean(alpha,beta,delta,mu,lambda);r2=gh_raw_moment(2,alpha,beta,delta,mu,lambda);r3=gh_raw_moment(3,alpha,beta,delta,mu,lambda);r4=gh_raw_moment(4,alpha,beta,delta,mu,lambda);var=r2-m*m
    v=(r4-4*m*r3+6*m*m*r2-3*m**4)/var**2-3.0_dp
  end function
  real(dp) function dgh(x,alpha,beta,delta,mu,lambda,log_density) result(v)
    real(dp),intent(in)::x,alpha,beta,delta,mu,lambda;logical,intent(in),optional::log_density
    real(dp)::g2,z0,s,z1,k0,k1,lv;logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(alpha<=0.or.delta<=0.or.abs(beta)>=alpha)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    g2=alpha*alpha-beta*beta;z0=delta*sqrt(g2);k0=bessel_k_nu(lambda,z0);s=sqrt(delta*delta+(x-mu)**2);z1=alpha*s;k1=bessel_k_nu(lambda-0.5_dp,z1)
    if(k0<=0.or.k1<=0)then;v=merge(-huge(1.0_dp),0.0_dp,lg);return;end if
    lv=0.5_dp*lambda*log(g2)-0.5_dp*log(2.0_dp*pi)-(lambda-0.5_dp)*log(alpha)-lambda*log(delta)-log(k0)+0.5_dp*(lambda-0.5_dp)*log(delta*delta+(x-mu)**2)+log(k1)+beta*(x-mu)
    v=merge(lv,exp(lv),lg)
  end function
  real(dp) function pgh(q,alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::q,alpha,beta,delta,mu,lambda
    real(dp)::m,s,lo
    m=gh_mean(alpha,beta,delta,mu,lambda)
    s=sqrt(max(gh_variance(alpha,beta,delta,mu,lambda),epsilon(1.0_dp)))
    lo=min(q,m-20.0_dp*s-10.0_dp*delta)
    if(q<=lo+epsilon(1.0_dp))then
      v=0.0_dp
    else
      v=integrate_gh_density(lo,q,alpha,beta,delta,mu,lambda)
    end if
    v=clamp(v,0.0_dp,1.0_dp)
  end function

  real(dp) function integrate_gh_density(a0,b0,alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::a0,b0,alpha,beta,delta,mu,lambda
    integer::n
    real(dp)::oldv,newv
    n=64
    oldv=gh_simpson(a0,b0,alpha,beta,delta,mu,lambda,n)
    do while(n<16384)
      n=2*n
      newv=gh_simpson(a0,b0,alpha,beta,delta,mu,lambda,n)
      if(abs(newv-oldv)<=2.0e-8_dp*max(1.0_dp,abs(newv)))exit
      oldv=newv
    end do
    v=newv
  end function

  real(dp) function gh_simpson(a0,b0,alpha,beta,delta,mu,lambda,n) result(v)
    real(dp),intent(in)::a0,b0,alpha,beta,delta,mu,lambda
    integer,intent(in)::n
    real(dp)::h,x0
    integer::i
    h=(b0-a0)/real(n,dp)
    v=dgh(a0,alpha,beta,delta,mu,lambda)+dgh(b0,alpha,beta,delta,mu,lambda)
    do i=1,n-1
      x0=a0+h*real(i,dp)
      if(mod(i,2)==0)then
        v=v+2.0_dp*dgh(x0,alpha,beta,delta,mu,lambda)
      else
        v=v+4.0_dp*dgh(x0,alpha,beta,delta,mu,lambda)
      end if
    end do
    v=v*h/3.0_dp
  end function

  real(dp) function qgh(p,alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::p,alpha,beta,delta,mu,lambda;real(dp)::lo,hi,mid,m,s;integer::i
    if(p<=0)then;v=-huge(1.0_dp);return;else if(p>=1)then;v=huge(1.0_dp);return;end if
    m=gh_mean(alpha,beta,delta,mu,lambda);s=sqrt(max(gh_variance(alpha,beta,delta,mu,lambda),epsilon(1.0_dp)));lo=m-4*s;hi=m+4*s
    do while(pgh(lo,alpha,beta,delta,mu,lambda)>p);lo=lo-2*s;end do;do while(pgh(hi,alpha,beta,delta,mu,lambda)<p);hi=hi+2*s;end do
    do i=1,80;mid=0.5_dp*(lo+hi);if(pgh(mid,alpha,beta,delta,mu,lambda)<p)then;lo=mid;else;hi=mid;end if;end do;v=0.5_dp*(lo+hi)
  end function
  real(dp) function rgig_slice(lambda,chi,psi) result(w)
    real(dp),intent(in)::lambda,chi,psi;real(dp)::y,logy,l,r,prop;integer::it,j
    if(psi<=0.or.chi<=0)then;w=1;return;end if
    w=(lambda+sqrt(lambda*lambda+chi*psi))/psi;if(w<=0)w=sqrt(chi/psi);y=log(w)
    do it=1,25
      logy=log(runif_lcg())+logtarget(y);l=y-1.0_dp*runif_lcg();r=l+1.0_dp
      do j=1,50;if(logtarget(l)<=logy)exit;l=l-1;end do;do j=1,50;if(logtarget(r)<=logy)exit;r=r+1;end do
      do j=1,100;prop=l+(r-l)*runif_lcg();if(logtarget(prop)>=logy)then;y=prop;exit;else if(prop<y)then;l=prop;else;r=prop;end if;end do
    end do;w=exp(y)
  contains
    real(dp) function logtarget(t) result(v);real(dp),intent(in)::t;v=lambda*t-0.5_dp*(chi*exp(-t)+psi*exp(t));end function
  end function
  real(dp) function rgh(alpha,beta,delta,mu,lambda) result(v)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;real(dp)::g2,w
    g2=alpha*alpha-beta*beta
    if(abs(lambda+0.5_dp)<1e-12_dp)then;v=rnig(alpha,beta,delta,mu);return;end if
    w=rgig_slice(lambda,delta*delta,g2);v=mu+beta*w+sqrt(w)*rnorm_lcg()
  end function

  subroutine standardized_gh_parameters(zeta,rho,lambda,alpha,beta,delta,mu)
    real(dp),intent(in)::zeta,rho,lambda;real(dp),intent(out)::alpha,beta,delta,mu
    real(dp)::r2,k0,k1,dk
    r2=1.0_dp-rho*rho;k0=bessel_k_nu(lambda+1.0_dp,zeta)/(bessel_k_nu(lambda,zeta)*zeta);k1=bessel_k_nu(lambda+2.0_dp,zeta)/(bessel_k_nu(lambda+1.0_dp,zeta)*zeta);dk=k1-k0
    alpha=sqrt(zeta*zeta*k0/r2*(1.0_dp+rho*rho*zeta*zeta*dk/r2));beta=alpha*rho;delta=zeta/(alpha*sqrt(r2));mu=-beta*delta*delta*k0
  end subroutine
  real(dp) function dsgh(x,zeta,rho,lambda) result(v);real(dp),intent(in)::x,zeta,rho,lambda;real(dp)::a,b,d,m;call standardized_gh_parameters(zeta,rho,lambda,a,b,d,m);v=dgh(x,a,b,d,m,lambda);end function
  real(dp) function psgh(x,zeta,rho,lambda) result(v);real(dp),intent(in)::x,zeta,rho,lambda;real(dp)::a,b,d,m;call standardized_gh_parameters(zeta,rho,lambda,a,b,d,m);v=pgh(x,a,b,d,m,lambda);end function
  real(dp) function qsgh(p,zeta,rho,lambda) result(v);real(dp),intent(in)::p,zeta,rho,lambda;real(dp)::a,b,d,m;call standardized_gh_parameters(zeta,rho,lambda,a,b,d,m);v=qgh(p,a,b,d,m,lambda);end function
  real(dp) function rsgh(zeta,rho,lambda) result(v);real(dp),intent(in)::zeta,rho,lambda;real(dp)::a,b,d,m;call standardized_gh_parameters(zeta,rho,lambda,a,b,d,m);v=rgh(a,b,d,m,lambda);end function
  real(dp) function dhyp(x,alpha,beta,delta,mu) result(v);real(dp),intent(in)::x,alpha,beta,delta,mu;v=dgh(x,alpha,beta,delta,mu,1.0_dp);end function
  real(dp) function phyp(x,alpha,beta,delta,mu) result(v);real(dp),intent(in)::x,alpha,beta,delta,mu;v=pgh(x,alpha,beta,delta,mu,1.0_dp);end function
  real(dp) function qhyp(p,alpha,beta,delta,mu) result(v);real(dp),intent(in)::p,alpha,beta,delta,mu;v=qgh(p,alpha,beta,delta,mu,1.0_dp);end function
  real(dp) function rhyp(alpha,beta,delta,mu) result(v);real(dp),intent(in)::alpha,beta,delta,mu;v=rgh(alpha,beta,delta,mu,1.0_dp);end function
  real(dp) function dght(x,beta,delta,mu,nu) result(v);real(dp),intent(in)::x,beta,delta,mu,nu;v=dgh(x,abs(beta)+1e-5_dp,beta,delta,mu,-0.5_dp*nu);end function
  real(dp) function pght(x,beta,delta,mu,nu) result(v);real(dp),intent(in)::x,beta,delta,mu,nu;v=pgh(x,abs(beta)+1e-5_dp,beta,delta,mu,-0.5_dp*nu);end function
  real(dp) function qght(p,beta,delta,mu,nu) result(v);real(dp),intent(in)::p,beta,delta,mu,nu;v=qgh(p,abs(beta)+1e-5_dp,beta,delta,mu,-0.5_dp*nu);end function
  real(dp) function rght(beta,delta,mu,nu) result(v);real(dp),intent(in)::beta,delta,mu,nu;v=rgh(abs(beta)+1e-5_dp,beta,delta,mu,-0.5_dp*nu);end function

  pure real(dp) function qgld_rs(p,l1,l2,l3,l4) result(v)
    real(dp),intent(in)::p,l1,l2,l3,l4;real(dp)::u
    u=clamp(p,epsilon(1.0_dp),1.0_dp-epsilon(1.0_dp));v=l1+(u**l3-(1.0_dp-u)**l4)/l2
  end function
  real(dp) function pgld_rs(x,l1,l2,l3,l4) result(v)
    real(dp),intent(in)::x,l1,l2,l3,l4;real(dp)::lo,hi,mid;integer::i
    lo=epsilon(1.0_dp);hi=1.0_dp-epsilon(1.0_dp)
    if(x<=qgld_rs(lo,l1,l2,l3,l4))then;v=0;return;else if(x>=qgld_rs(hi,l1,l2,l3,l4))then;v=1;return;end if
    do i=1,100;mid=0.5_dp*(lo+hi);if(qgld_rs(mid,l1,l2,l3,l4)<x)then;lo=mid;else;hi=mid;end if;end do;v=0.5_dp*(lo+hi)
  end function
  real(dp) function dgld_rs(x,l1,l2,l3,l4) result(v)
    real(dp),intent(in)::x,l1,l2,l3,l4;real(dp)::p,den
    p=pgld_rs(x,l1,l2,l3,l4);if(p<=0.or.p>=1)then;v=0;return;end if;den=l3*p**(l3-1.0_dp)+l4*(1.0_dp-p)**(l4-1.0_dp);v=l2/den
  end function
  real(dp) function rgld_rs(l1,l2,l3,l4) result(v);real(dp),intent(in)::l1,l2,l3,l4;v=qgld_rs(runif_lcg(),l1,l2,l3,l4);end function
  real(dp) function gld_mode(l1,l2,l3,l4) result(v)
    real(dp),intent(in)::l1,l2,l3,l4;real(dp)::p,best,bp,d;integer::i
    best=-1;bp=0.5_dp;do i=1,2000;p=(real(i,dp)-0.5_dp)/2000.0_dp;d=l2/(l3*p**(l3-1.0_dp)+l4*(1-p)**(l4-1.0_dp));if(d>best)then;best=d;bp=p;end if;end do;v=qgld_rs(bp,l1,l2,l3,l4)
  end function

  subroutine fit_normal(x,fit)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;real(dp)::m,s
    m=sample_mean(x);s=sqrt(sum((x-m)**2)/real(size(x),dp));fit%family='normal';allocate(fit%parameters(2));fit%parameters=[m,s];fit%loglik=sum(log(dnorm_fs(x,m,s)));fit%aic=4-2*fit%loglik;fit%bic=2*log(real(size(x),dp))-2*fit%loglik;fit%converged=.true.
  end subroutine
  subroutine fit_student(x,fit)
    real(dp),intent(in)::x(:)
    type(distribution_fit),intent(out)::fit
    real(dp),allocatable::start(:),lower(:),upper(:),best(:)
    real(dp)::fbest
    logical::ok
    integer::n,info
    n=size(x)
    active_fit_data=x
    allocate(start(3),lower(3),upper(3))
    start=[sample_mean(x),max(sample_sd(x),1e-4_dp),8.0_dp]
    lower=[minval(x)-5.0_dp*sample_sd(x),1e-6_dp,2.01_dp]
    upper=[maxval(x)+5.0_dp*sample_sd(x),10.0_dp*sample_sd(x)+1.0_dp,200.0_dp]
    call nelder_mead_bounded(student_fit_objective,start,lower,upper,best,fbest,ok,3000,1e-8_dp)
    fit%family='student';fit%parameters=best;fit%loglik=-fbest
    fit%aic=6.0_dp-2.0_dp*fit%loglik
    fit%bic=3.0_dp*log(real(n,dp))-2.0_dp*fit%loglik
    fit%converged=ok
    call numerical_hessian(student_fit_objective,best,fit%hessian)
    call matrix_inverse(fit%hessian,fit%covariance,info)
    if(allocated(active_fit_data))deallocate(active_fit_data)
  end subroutine

  subroutine fit_nig(x,fit)
    real(dp),intent(in)::x(:)
    type(distribution_fit),intent(out)::fit
    real(dp),allocatable::start(:),lower(:),upper(:),best(:)
    real(dp)::fbest,sdx
    logical::ok
    integer::n,info
    n=size(x);sdx=max(sample_sd(x),1e-3_dp)
    active_fit_data=x
    allocate(start(4),lower(4),upper(4))
    start=[2.0_dp/sdx,0.0_dp,sdx,sample_mean(x)]
    lower=[1e-3_dp,-20.0_dp/sdx,1e-4_dp,minval(x)-5.0_dp*sdx]
    upper=[50.0_dp/sdx,20.0_dp/sdx,20.0_dp*sdx,maxval(x)+5.0_dp*sdx]
    call nelder_mead_bounded(nig_fit_objective,start,lower,upper,best,fbest,ok,1500,1e-7_dp)
    fit%family='nig';fit%parameters=best;fit%loglik=-fbest
    fit%aic=8.0_dp-2.0_dp*fit%loglik
    fit%bic=4.0_dp*log(real(n,dp))-2.0_dp*fit%loglik
    fit%converged=ok
    call numerical_hessian(nig_fit_objective,best,fit%hessian,5e-4_dp)
    call matrix_inverse(fit%hessian,fit%covariance,info)
    if(allocated(active_fit_data))deallocate(active_fit_data)
  end subroutine

  subroutine fit_gld_quantiles(x,fit)
    real(dp),intent(in)::x(:)
    type(distribution_fit),intent(out)::fit
    real(dp),allocatable::start(:),lower(:),upper(:),best(:)
    real(dp)::fbest
    logical::ok
    integer::i,n
    n=size(x)
    active_probs=[.1_dp,.25_dp,.5_dp,.75_dp,.9_dp]
    allocate(active_qobs(5))
    do i=1,5
      active_qobs(i)=sample_quantile(x,active_probs(i))
    end do
    allocate(start(4),lower(4),upper(4))
    start=[active_qobs(3),log(max(sample_sd(x),1e-3_dp)),-2.0_dp,-2.0_dp]
    lower=[minval(x)-5.0_dp*sample_sd(x),-10.0_dp,-8.0_dp,-8.0_dp]
    upper=[maxval(x)+5.0_dp*sample_sd(x),10.0_dp,3.0_dp,3.0_dp]
    call nelder_mead_bounded(gld_fit_objective,start,lower,upper,best,fbest,ok,2000,1e-9_dp)
    fit%family='gld-rs';allocate(fit%parameters(4))
    fit%parameters=[best(1),-exp(best(2)),-exp(best(3)),-exp(best(4))]
    fit%loglik=-0.5_dp*fbest;fit%aic=8.0_dp+fbest
    fit%bic=4.0_dp*log(real(n,dp))+fbest;fit%converged=ok
    if(allocated(active_qobs))deallocate(active_qobs)
    if(allocated(active_probs))deallocate(active_probs)
  end subroutine

  real(dp) function student_fit_objective(p) result(v)
    real(dp),intent(in)::p(:)
    integer::i
    v=0.0_dp
    do i=1,size(active_fit_data)
      v=v-dt_fs(active_fit_data(i),p(3),p(1),p(2),.true.)
    end do
  end function

  real(dp) function nig_fit_objective(p) result(v)
    real(dp),intent(in)::p(:)
    integer::i
    if(abs(p(2))>=p(1))then
      v=huge(1.0_dp)/100.0_dp
      return
    end if
    v=0.0_dp
    do i=1,size(active_fit_data)
      v=v-dnig(active_fit_data(i),p(1),p(2),p(3),p(4),.true.)
    end do
  end function

  real(dp) function gld_fit_objective(p) result(v)
    real(dp),intent(in)::p(:)
    real(dp)::l2,l3,l4
    integer::i
    l2=-exp(p(2));l3=-exp(p(3));l4=-exp(p(4));v=0.0_dp
    do i=1,size(active_probs)
      v=v+(qgld_rs(active_probs(i),p(1),l2,l3,l4)-active_qobs(i))**2
    end do
  end function
end module fbasics_distributions
