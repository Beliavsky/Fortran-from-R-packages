! SPDX-License-Identifier: GPL-3.0-only
module tscopula_timeseries
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tscopula_kinds, only : dp, log_two_pi
  use tscopula_math, only : normal_cdf, normal_quantile, normal_random, &
    optimizer_result, minimize_nelder_mead
  implicit none
  private

  type, public :: arma_copula
    real(dp), allocatable :: ar(:)
    real(dp), allocatable :: ma(:)
  end type arma_copula

  type, public :: sarma_copula
    real(dp), allocatable :: ar(:), ma(:), sar(:), sma(:)
    integer :: period = 4
  end type sarma_copula

  type, public :: arma_filter_result
    real(dp), allocatable :: mean(:), sigma(:), innovations(:), residuals(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
  end type arma_filter_result

  type, public :: arma_fit_result
    type(arma_copula) :: model
    real(dp), allocatable :: parameters(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: convergence = 1
    integer :: iterations = 0
  end type arma_fit_result

  type :: arma_fit_context
    integer :: p = 0, q = 0
    real(dp), allocatable :: u(:)
  end type arma_fit_context

  public :: armacopula, sarmacopula, sarma2arma, expand_ar, expand_ma
  public :: non_stat, non_invert, sigmastarma, arma_autocovariance
  public :: sim_arma_copula, sim_sarma_copula, kfilter, arma_objective
  public :: predict_arma_cdf, predict_arma_quantile, predict_arma_density
  public :: fit_arma_copula, resid_arma_copula, kendall_arma, glag
  public :: acf2pacf, pacf2ar, pacf2acf, kpacf_arma, kpacf_sarma4
  public :: kpacf_sarma12, kpacf_fbn, kpacf_arfima, strank

contains

  function armacopula(ar,ma) result(model)
    real(dp),intent(in),optional :: ar(:),ma(:)
    type(arma_copula) :: model
    if(present(ar))then;allocate(model%ar(size(ar)));model%ar=ar;else;allocate(model%ar(0));end if
    if(present(ma))then;allocate(model%ma(size(ma)));model%ma=ma;else;allocate(model%ma(0));end if
  end function armacopula

  function sarmacopula(ar,ma,sar,sma,period) result(model)
    real(dp),intent(in),optional :: ar(:),ma(:),sar(:),sma(:)
    integer,intent(in),optional :: period
    type(sarma_copula) :: model
    if(present(ar))then;allocate(model%ar(size(ar)));model%ar=ar;else;allocate(model%ar(0));end if
    if(present(ma))then;allocate(model%ma(size(ma)));model%ma=ma;else;allocate(model%ma(0));end if
    if(present(sar))then;allocate(model%sar(size(sar)));model%sar=sar;else;allocate(model%sar(0));end if
    if(present(sma))then;allocate(model%sma(size(sma)));model%sma=sma;else;allocate(model%sma(0));end if
    if(present(period))model%period=period
  end function sarmacopula

  function expand_ar(ar,sar,period) result(phi)
    real(dp),intent(in)::ar(:),sar(:);integer,intent(in)::period
    real(dp),allocatable::a(:),b(:),phi(:),tmp(:)
    integer::i,j
    allocate(a(0:size(ar)));a=0.0_dp;a(0)=1.0_dp
    do i=1,size(ar);a(i)=-ar(i);end do
    allocate(b(0:max(0,period*size(sar))));b=0.0_dp;b(0)=1.0_dp
    do i=1,size(sar);b(i*period)=-sar(i);end do
    allocate(tmp(0:ubound(a,1)+ubound(b,1)));tmp=0.0_dp
    do i=0,ubound(a,1);do j=0,ubound(b,1);tmp(i+j)=tmp(i+j)+a(i)*b(j);end do;end do
    allocate(phi(size(tmp)-1));do i=1,size(phi);phi(i)=-tmp(i);end do
  end function expand_ar

  function expand_ma(ma,sma,period) result(theta)
    real(dp),intent(in)::ma(:),sma(:);integer,intent(in)::period
    real(dp),allocatable::a(:),b(:),theta(:),tmp(:)
    integer::i,j
    allocate(a(0:size(ma)));a=0.0_dp;a(0)=1.0_dp
    do i=1,size(ma);a(i)=ma(i);end do
    allocate(b(0:max(0,period*size(sma))));b=0.0_dp;b(0)=1.0_dp
    do i=1,size(sma);b(i*period)=sma(i);end do
    allocate(tmp(0:ubound(a,1)+ubound(b,1)));tmp=0.0_dp
    do i=0,ubound(a,1);do j=0,ubound(b,1);tmp(i+j)=tmp(i+j)+a(i)*b(j);end do;end do
    allocate(theta(size(tmp)-1));do i=1,size(theta);theta(i)=tmp(i);end do
  end function expand_ma

  function sarma2arma(model) result(out)
    type(sarma_copula),intent(in)::model
    type(arma_copula)::out
    real(dp),allocatable :: a(:), b(:)
    a=expand_ar(model%ar,model%sar,model%period)
    b=expand_ma(model%ma,model%sma,model%period)
    allocate(out%ar(size(a)),out%ma(size(b)))
    out%ar=a;out%ma=b
  end function sarma2arma

  logical function non_stat(ar) result(bad)
    real(dp),intent(in)::ar(:)
    real(dp),allocatable::pacf(:)
    if(size(ar)==0)then;bad=.false.;return;end if
    call ar_to_pacf(ar,pacf,bad)
  end function non_stat

  logical function non_invert(ma) result(bad)
    real(dp),intent(in)::ma(:)
    real(dp),allocatable::ar_like(:)
    if(size(ma)==0)then;bad=.false.;return;end if
    allocate(ar_like(size(ma)));ar_like=-ma
    bad=non_stat(ar_like)
  end function non_invert

  subroutine ar_to_pacf(ar,pacf,bad)
    real(dp),intent(in)::ar(:);real(dp),allocatable,intent(out)::pacf(:);logical,intent(out)::bad
    real(dp),allocatable::work(:),prev(:);real(dp)::den
    integer::m,j
    allocate(work(size(ar)),prev(size(ar)),pacf(size(ar)));work=ar;pacf=0.0_dp;bad=.false.
    do m=size(ar),1,-1
      pacf(m)=work(m)
      if(abs(pacf(m))>=1.0_dp)then;bad=.true.;return;end if
      if(m>1)then
        prev=0.0_dp;den=1.0_dp-pacf(m)*pacf(m)
        do j=1,m-1;prev(j)=(work(j)+pacf(m)*work(m-j))/den;end do
        work(1:m-1)=prev(1:m-1)
      end if
    end do
  end subroutine ar_to_pacf

  function pacf2ar(pacf) result(ar)
    real(dp),intent(in)::pacf(:);real(dp),allocatable::ar(:),old(:)
    integer::m,j
    allocate(ar(size(pacf)),old(size(pacf)));ar=0.0_dp
    do m=1,size(pacf)
      old=ar;ar(m)=pacf(m)
      do j=1,m-1;ar(j)=old(j)-pacf(m)*old(m-j);end do
    end do
  end function pacf2ar

  function arma_psi(model,nmax) result(psi)
    type(arma_copula),intent(in)::model;integer,intent(in)::nmax
    real(dp),allocatable::psi(:);integer::k,j
    allocate(psi(0:nmax));psi=0.0_dp;psi(0)=1.0_dp
    do k=1,nmax
      if(k<=size(model%ma))psi(k)=model%ma(k)
      do j=1,min(k,size(model%ar));psi(k)=psi(k)+model%ar(j)*psi(k-j);end do
    end do
  end function arma_psi

  function arma_autocovariance(model,maxlag,innovation_variance) result(gamma)
    type(arma_copula),intent(in)::model;integer,intent(in)::maxlag
    real(dp),intent(in),optional::innovation_variance
    real(dp),allocatable::gamma(:),psi(:);real(dp)::s2;integer::h,j,m
    s2=1.0_dp;if(present(innovation_variance))s2=innovation_variance
    m=max(2000,maxlag+1000);allocate(psi(0:m+maxlag));psi=arma_psi(model,m+maxlag);allocate(gamma(maxlag+1));gamma=0.0_dp
    do h=0,maxlag;do j=0,m;gamma(h+1)=gamma(h+1)+s2*psi(j)*psi(j+h);end do;end do
  end function arma_autocovariance

  real(dp) function sigmastarma(model) result(sigma)
    class(*),intent(in)::model
    type(arma_copula)::a;real(dp),allocatable::g(:)
    select type(model)
    type is(arma_copula);a=model
    type is(sarma_copula);a=sarma2arma(model)
    class default;sigma=ieee_value(0.0_dp,ieee_quiet_nan);return
    end select
    g=arma_autocovariance(a,0,1.0_dp);sigma=1.0_dp/sqrt(max(g(1),tiny(1.0_dp)))
  end function sigmastarma

  function sim_arma_copula(model,n,burn) result(u)
    type(arma_copula),intent(in)::model;integer,intent(in)::n;integer,intent(in),optional::burn
    real(dp),allocatable::u(:),x(:),eps(:);real(dp)::sigma;integer::b,t,j,nt
    b=max(100,20*(size(model%ar)+size(model%ma)+1));if(present(burn))b=burn;nt=n+b
    allocate(x(nt),eps(nt),u(n));x=0.0_dp;eps=0.0_dp;sigma=sigmastarma(model)
    do t=1,nt
      eps(t)=sigma*normal_random();x(t)=eps(t)
      do j=1,min(size(model%ar),t-1);x(t)=x(t)+model%ar(j)*x(t-j);end do
      do j=1,min(size(model%ma),t-1);x(t)=x(t)+model%ma(j)*eps(t-j);end do
    end do
    do t=1,n;u(t)=normal_cdf(x(b+t));end do
  end function sim_arma_copula

  function sim_sarma_copula(model,n,burn) result(u)
    type(sarma_copula),intent(in)::model;integer,intent(in)::n;integer,intent(in),optional::burn
    real(dp),allocatable::u(:);type(arma_copula)::a
    a=sarma2arma(model);if(present(burn))then;u=sim_arma_copula(a,n,burn);else;u=sim_arma_copula(a,n);end if
  end function sim_sarma_copula

  subroutine state_matrices(model,tmat,qmat,p0)
    type(arma_copula),intent(in)::model
    real(dp),allocatable,intent(out)::tmat(:,:),qmat(:,:),p0(:,:)
    real(dp),allocatable::h(:),pnew(:,:)
    real(dp)::sigma,diff
    integer::m,i,iter
    m=max(1,max(size(model%ar),size(model%ma)+1))
    allocate(tmat(m,m),qmat(m,m),p0(m,m),h(m),pnew(m,m))
    tmat=0.0_dp;h=0.0_dp
    do i=1,size(model%ar);tmat(i,1)=model%ar(i);end do
    do i=2,m;tmat(i-1,i)=1.0_dp;end do
    sigma=sigmastarma(model);h(1)=sigma
    do i=1,min(size(model%ma),m-1);h(i+1)=sigma*model%ma(i);end do
    do i=1,m;qmat(:,i)=h*h(i);end do
    p0=qmat
    do iter=1,10000
      pnew=matmul(tmat,matmul(p0,transpose(tmat)))+qmat
      diff=maxval(abs(pnew-p0));p0=pnew
      if(diff<=1.0e-13_dp*(1.0_dp+maxval(abs(p0))))exit
    end do
  end subroutine state_matrices

  subroutine filter_state(model,u,ans,a_next,p_next)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:)
    type(arma_filter_result),intent(out)::ans
    real(dp),allocatable,intent(out)::a_next(:),p_next(:,:)
    real(dp),allocatable::x(:),tmat(:,:),qmat(:,:),a(:),p(:,:),af(:),pf(:,:),k(:)
    real(dp)::fvar,ll;integer::t,n,m
    n=size(u);call state_matrices(model,tmat,qmat,p);m=size(tmat,1)
    allocate(x(n),ans%mean(n),ans%sigma(n),ans%innovations(n),ans%residuals(n),a(m),af(m),pf(m,m),k(m))
    x=normal_quantile(min(max(u,epsilon(1.0_dp)),1.0_dp-epsilon(1.0_dp)))
    a=0.0_dp;ll=0.0_dp
    do t=1,n
      ans%mean(t)=a(1);fvar=max(p(1,1),tiny(1.0_dp));ans%sigma(t)=sqrt(fvar)
      ans%innovations(t)=x(t)-ans%mean(t);ans%residuals(t)=ans%innovations(t)/ans%sigma(t)
      ll=ll-0.5_dp*(log_two_pi+log(fvar)+ans%residuals(t)**2)
      k=p(:,1)/fvar;af=a+k*ans%innovations(t);pf=p
      pf=pf-spread(k,2,m)*spread(p(1,:),1,m)
      pf=0.5_dp*(pf+transpose(pf));a=matmul(tmat,af);p=matmul(tmat,matmul(pf,transpose(tmat)))+qmat
    end do
    ans%log_likelihood=ll;ans%objective=-ll+sum(-0.5_dp*log_two_pi-0.5_dp*x*x)
    allocate(a_next(m),p_next(m,m));a_next=a;p_next=p
  end subroutine filter_state

  function kfilter(model,u) result(ans)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:)
    type(arma_filter_result)::ans;real(dp),allocatable::a(:),p(:,:)
    call filter_state(model,u,ans,a,p)
  end function kfilter

  real(dp) function arma_objective(model,u) result(value)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:);type(arma_filter_result)::f
    if(non_stat(model%ar).or.non_invert(model%ma))then;value=huge(1.0_dp)/100.0_dp;return;end if
    f=kfilter(model,u);value=f%objective
  end function arma_objective

  function resid_arma_copula(model,u) result(residuals)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:);real(dp),allocatable::residuals(:);type(arma_filter_result)::f
    f=kfilter(model,u);allocate(residuals(size(u)));residuals=normal_cdf(f%residuals)
  end function resid_arma_copula

  subroutine next_state(model,u,mean,sigma)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:);real(dp),intent(out)::mean,sigma
    type(arma_filter_result)::f;real(dp),allocatable::a(:),p(:,:)
    call filter_state(model,u,f,a,p);mean=a(1);sigma=sqrt(max(p(1,1),tiny(1.0_dp)))
  end subroutine next_state

  real(dp) function predict_arma_cdf(model,u,x) result(value)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:),x;real(dp)::m,s
    call next_state(model,u,m,s);value=normal_cdf((normal_quantile(x)-m)/s)
  end function predict_arma_cdf

  real(dp) function predict_arma_quantile(model,u,p) result(value)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:),p;real(dp)::m,s
    call next_state(model,u,m,s);value=normal_cdf(m+s*normal_quantile(p))
  end function predict_arma_quantile

  real(dp) function predict_arma_density(model,u,x) result(value)
    type(arma_copula),intent(in)::model;real(dp),intent(in)::u(:),x;real(dp)::m,s,z,q
    call next_state(model,u,m,s);q=normal_quantile(x);z=(q-m)/s
    value=exp(-0.5_dp*z*z+0.5_dp*q*q)/s
  end function predict_arma_density

  function fit_arma_copula(u,p,q,start,max_iter) result(fit)
    real(dp),intent(in)::u(:);integer,intent(in)::p,q;real(dp),intent(in),optional::start(:);integer,intent(in),optional::max_iter
    type(arma_fit_result)::fit;type(arma_fit_context)::ctx;type(optimizer_result)::opt
    real(dp),allocatable::x0(:),lo(:),hi(:);integer::npar,mit
    npar=p+q;allocate(x0(npar),lo(npar),hi(npar));x0=0.0_dp;if(present(start))x0=start
    lo=-0.98_dp;hi=0.98_dp;ctx%p=p;ctx%q=q;allocate(ctx%u(size(u)));ctx%u=u;mit=1500;if(present(max_iter))mit=max_iter
    call minimize_nelder_mead(arma_fit_objective,x0,lo,hi,ctx,opt,max_iter=mit,tolerance=1.0e-7_dp)
    fit%model=armacopula(opt%par(1:p),opt%par(p+1:npar));allocate(fit%parameters(npar));fit%parameters=opt%par;fit%objective=opt%value
    fit%convergence=opt%convergence;fit%iterations=opt%iterations
  end function fit_arma_copula

  real(dp) function arma_fit_objective(theta,context_any) result(value)
    real(dp),intent(in)::theta(:);class(*),intent(inout)::context_any;type(arma_copula)::model
    select type(context=>context_any)
    type is(arma_fit_context)
      model=armacopula(theta(1:context%p),theta(context%p+1:context%p+context%q));value=arma_objective(model,context%u)
    class default;value=huge(1.0_dp)
    end select
  end function arma_fit_objective

  function acf2pacf(acf) result(pacf)
    real(dp),intent(in)::acf(:);real(dp),allocatable::pacf(:),phi(:,:),v(:);integer::k,j,n
    n=size(acf)-1;allocate(pacf(n),phi(n,n),v(0:n));phi=0.0_dp;v=0.0_dp
    if(n<=0)return;v(0)=acf(1)
    do k=1,n
      if(k==1)then;phi(k,k)=acf(2)/acf(1)
      else
        phi(k,k)=(acf(k+1)-sum(phi(k-1,1:k-1)*acf(k:2:-1)))/v(k-1)
        do j=1,k-1;phi(k,j)=phi(k-1,j)-phi(k,k)*phi(k-1,k-j);end do
      end if
      pacf(k)=phi(k,k);v(k)=v(k-1)*(1.0_dp-phi(k,k)**2)
    end do
  end function acf2pacf

  function pacf2acf(pacf) result(acf)
    real(dp),intent(in)::pacf(:);real(dp),allocatable::acf(:),ar(:);integer::h,j,n
    n=size(pacf);allocate(acf(0:n));acf=0.0_dp;acf(0)=1.0_dp
    do h=1,n
      ar=pacf2ar(pacf(1:h));acf(h)=ar(h)
      do j=1,h-1;acf(h)=acf(h)+ar(j)*acf(h-j);end do
    end do
  end function pacf2acf

  function kpacf_arma(ar,ma,maxlag) result(kpacf)
    real(dp),intent(in)::ar(:),ma(:);integer,intent(in)::maxlag
    type(arma_copula)::m;real(dp),allocatable::g(:),acf(:),pacf(:),kpacf(:)
    m=armacopula(ar,ma);g=arma_autocovariance(m,maxlag);allocate(acf(maxlag+1));acf=g/g(1);pacf=acf2pacf(acf)
    allocate(kpacf(maxlag));kpacf=(2.0_dp/3.1415926535897932384626433832795_dp)*asin(min(max(pacf,-1.0_dp),1.0_dp))
  end function kpacf_arma

  function kpacf_sarma4(ar,ma,sar,sma,maxlag) result(kpacf)
    real(dp),intent(in)::ar(:),ma(:),sar(:),sma(:);integer,intent(in)::maxlag;real(dp),allocatable::kpacf(:);type(sarma_copula)::s;type(arma_copula)::a
    s=sarmacopula(ar,ma,sar,sma,4);a=sarma2arma(s);kpacf=kpacf_arma(a%ar,a%ma,maxlag)
  end function kpacf_sarma4

  function kpacf_sarma12(ar,ma,sar,sma,maxlag) result(kpacf)
    real(dp),intent(in)::ar(:),ma(:),sar(:),sma(:);integer,intent(in)::maxlag;real(dp),allocatable::kpacf(:);type(sarma_copula)::s;type(arma_copula)::a
    s=sarmacopula(ar,ma,sar,sma,12);a=sarma2arma(s);kpacf=kpacf_arma(a%ar,a%ma,maxlag)
  end function kpacf_sarma12

  function kpacf_fbn(hurst,maxlag) result(kpacf)
    real(dp),intent(in)::hurst;integer,intent(in)::maxlag;real(dp),allocatable::acf(:),pacf(:),kpacf(:);integer::h
    allocate(acf(0:maxlag));acf(0)=1.0_dp
    do h=1,maxlag;acf(h)=0.5_dp*((real(h+1,dp))**(2.0_dp*hurst)-2.0_dp*(real(h,dp))**(2.0_dp*hurst)+(real(h-1,dp))**(2.0_dp*hurst));end do
    pacf=acf2pacf(acf);allocate(kpacf(maxlag));kpacf=2.0_dp*asin(min(max(pacf,-1.0_dp),1.0_dp))/acos(-1.0_dp)
  end function kpacf_fbn

  function kpacf_arfima(ar,ma,d,maxlag) result(kpacf)
    real(dp),intent(in)::ar(:),ma(:),d;integer,intent(in)::maxlag
    real(dp),allocatable::frac(:),psi0(:),psi(:),acf(:),pacf(:),kpacf(:);type(arma_copula)::m;integer::j,k,nm
    nm=max(3000,maxlag+1500);allocate(frac(0:nm));frac=0.0_dp;frac(0)=1.0_dp
    do j=1,nm;frac(j)=frac(j-1)*(real(j-1,dp)+d)/real(j,dp);end do
    m=armacopula(ar,ma);allocate(psi0(0:nm));psi0=arma_psi(m,nm);allocate(psi(0:nm));psi=0.0_dp
    do k=0,nm;do j=0,k;psi(k)=psi(k)+psi0(j)*frac(k-j);end do;end do
    allocate(acf(0:maxlag));do k=0,maxlag;acf(k)=sum(psi(0:nm-k)*psi(k:nm));end do;acf=acf/acf(0)
    pacf=acf2pacf(acf);allocate(kpacf(maxlag));kpacf=2.0_dp*asin(min(max(pacf,-1.0_dp),1.0_dp))/acos(-1.0_dp)
  end function kpacf_arfima

  real(dp) function glag(model,lag) result(value)
    class(*),intent(in)::model;integer,intent(in)::lag;real(dp),allocatable::tmp(:),k(:)
    allocate(k(lag));k=0.0_dp
    select type(model)
    type is(arma_copula);tmp=kpacf_arma(model%ar,model%ma,lag);k=tmp
    type is(sarma_copula)
      if(model%period==4)then;tmp=kpacf_sarma4(model%ar,model%ma,model%sar,model%sma,lag)
      else;tmp=kpacf_sarma12(model%ar,model%ma,model%sar,model%sma,lag);end if;k=tmp
    class default;value=0.0_dp;return
    end select
    value=k(lag)
  end function glag

  function kendall_arma(model,maxlag) result(kendall)
    type(arma_copula),intent(in)::model;integer,intent(in)::maxlag;real(dp),allocatable::kendall(:),g(:)
    g=arma_autocovariance(model,maxlag);allocate(kendall(maxlag));kendall=2.0_dp*asin(min(max(g(2:maxlag+1)/g(1),-1.0_dp),1.0_dp))/acos(-1.0_dp)
  end function kendall_arma

  function strank(x) result(u)
    real(dp),intent(in)::x(:);real(dp),allocatable::u(:);integer::i,j,n;real(dp)::less,equal
    n=size(x);allocate(u(n))
    do i=1,n;less=0.0_dp;equal=0.0_dp;do j=1,n;if(x(j)<x(i))less=less+1.0_dp;if(abs(x(j)-x(i))<=8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i))))equal=equal+1.0_dp;end do;u(i)=(less+0.5_dp*equal)/real(n,dp);end do
  end function strank

end module tscopula_timeseries
