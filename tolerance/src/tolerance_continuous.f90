! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of continuous-family tolerance intervals from
! the R package tolerance (D. S. Young and K. Cheng).
module tolerance_continuous
  use tolerance_kinds, only : dp, pi
  use tolerance_types, only : tolerance_interval
  use tolerance_math, only : sample_mean, sample_median, sample_quantile, normal_quantile, &
       chisq_quantile, noncentral_t_quantile, digamma, invert_matrix
  use tolerance_normal, only : k_factor
  use tolerance_optimize, only : nelder_mead, numerical_hessian
  implicit none
  private
  public :: exptol_int, exp2tol_int, uniftol_int, laptol_int, gamtol_int
  public :: logistol_int, cautol_int, paretotol_int, exttol_int

contains

  function exptol_int(x,alpha,p,side,type2) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    logical, intent(in), optional :: type2
    type(tolerance_interval) :: out
    real(dp) :: a,pp,lhat,mx
    integer :: n,r,s,i
    logical :: t2
    a=0.05_dp; if(present(alpha)) a=alpha
    pp=0.99_dp; if(present(p)) pp=p
    s=1; if(present(side)) s=side
    t2=.false.; if(present(type2)) t2=type2
    n=size(x); lhat=sample_mean(x); r=n
    if(t2)then
      mx=maxval(x); r=0
      do i=1,n
        if(x(i)/=mx) r=r+1
      end do
    end if
    out%alpha=a; out%p=pp; out%estimate=lhat
    if(s==2)then
      a=a/2.0_dp
      out%lower=2.0_dp*real(r,dp)*lhat*log(2.0_dp/(1.0_dp+pp))/ &
           chisq_quantile(1.0_dp-a,2.0_dp*real(r,dp))
      out%upper=2.0_dp*real(r,dp)*lhat*log(2.0_dp/(1.0_dp-pp))/ &
           chisq_quantile(a,2.0_dp*real(r,dp))
    else
      out%lower=2.0_dp*real(r,dp)*lhat*log(1.0_dp/pp)/ &
           chisq_quantile(1.0_dp-a,2.0_dp*real(r,dp))
      out%upper=2.0_dp*real(r,dp)*lhat*log(1.0_dp/(1.0_dp-pp))/ &
           chisq_quantile(a,2.0_dp*real(r,dp))
    end if
  end function exptol_int

  function exp2tol_int(x,alpha,p,side,method,type2) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    character(len=*), intent(in), optional :: method
    logical, intent(in), optional :: type2
    type(tolerance_interval) :: out
    real(dp) :: a,pp,t,sx,k1,k2,mm,lam,mx
    integer :: n,r,sd,i
    logical :: t2
    character(len=8) :: meth
    a=0.05_dp; if(present(alpha)) a=alpha
    pp=0.99_dp; if(present(p)) pp=p
    sd=1; if(present(side)) sd=side
    meth='GPU'; if(present(method)) meth=adjustl(method)
    t2=.false.; if(present(type2)) t2=type2
    n=size(x); t=minval(x); sx=sum(x-t); r=n
    if(sd==2)then; a=a/2.0_dp; pp=(pp+1.0_dp)/2.0_dp; end if
    if(t2)then
      mx=maxval(x); r=0
      do i=1,n
        if(x(i)/=mx) r=r+1
      end do
      mm=(1.0_dp+real(n,dp)*log(pp))/(real(r,dp)-2.5_dp)
      k1=(-mm-normal_quantile(1.0_dp-a)*sqrt(mm*mm/real(r,dp)+1.0_dp/real(r*r,dp)))/real(n,dp)
      mm=(1.0_dp+real(n,dp)*log(1.0_dp-pp))/(real(r,dp)-2.5_dp)
      k2=(-mm-normal_quantile(a)*sqrt(mm*mm/real(r,dp)+1.0_dp/real(r*r,dp)))/real(n,dp)
    else
      k1=(1.0_dp-((pp**n)/a)**(1.0_dp/real(n-1,dp)))/real(n,dp)
      if(trim(meth)=='KM')then
        k2=(1.0_dp-(((1.0_dp-pp)**n)/(1.0_dp-a))**(1.0_dp/real(n-1,dp)))/real(n,dp)
      else
        k2=chisq_quantile(pp,2.0_dp)/chisq_quantile(a,2.0_dp*real(n-1,dp))
        if(trim(meth)=='DUN')then
          lam=1.71_dp+1.57_dp*log(log(1.0_dp/a))
          k2=k2-(lam/real(n,dp))**(1.63_dp+0.39_dp*lam)
        end if
      end if
    end if
    out%alpha=merge(2.0_dp*a,a,sd==2)
    out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
    out%estimate=t
    out%lower=t+sx*k1; out%upper=t+sx*k2
  end function exp2tol_int

  function uniftol_int(x,alpha,p,side,known_lower,known_upper) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p,known_lower,known_upper
    integer, intent(in), optional :: side
    type(tolerance_interval) :: out
    real(dp) :: a,pp,x1,xn
    integer :: n,sd
    a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    n=size(x);x1=minval(x);xn=maxval(x)
    if(present(known_lower))x1=known_lower
    if(present(known_upper))xn=known_upper
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
    out%estimate=0.5_dp*(x1+xn)
    out%lower=((1.0_dp-pp)/(1.0_dp-a)**(1.0_dp/real(n,dp)))*(xn-x1)+x1
    out%upper=(pp/a**(1.0_dp/real(n,dp)))*(xn-x1)+x1
  end function uniftol_int

  function laptol_int(x,alpha,p,side) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    type(tolerance_interval) :: out
    real(dp) :: a,pp,kb,mu,beta,k,z
    integer :: n,sd
    a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    n=size(x);mu=sample_median(x);beta=sum(abs(x-mu))/real(n,dp)
    kb=log(2.0_dp*(1.0_dp-pp));z=normal_quantile(1.0_dp-a)
    k=(-real(n,dp)*kb+z*sqrt(real(n,dp)*(1.0_dp+kb*kb)-z*z))/(real(n,dp)-z*z)
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
    out%estimate=mu;out%lower=mu-k*beta;out%upper=mu+k*beta
  end function laptol_int

  function gamtol_int(x,alpha,p,side,method,m,log_gamma_data,shape_hat,scale_hat) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side,m
    character(len=*), intent(in), optional :: method
    logical, intent(in), optional :: log_gamma_data
    real(dp), intent(out), optional :: shape_hat,scale_hat
    type(tolerance_interval) :: out
    real(dp), allocatable :: y(:)
    real(dp) :: a,pp,meanx,varx,shape,scale,mu3,s3,kfac
    integer :: sd,mm,n
    logical :: lg
    character(len=8) :: meth
    a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side
    mm=50;if(present(m))mm=m
    meth='HE';if(present(method))meth=adjustl(method)
    lg=.false.;if(present(log_gamma_data))lg=log_gamma_data
    allocate(y(size(x)));y=x;if(lg)y=log(x)
    n=size(y);meanx=sum(y)/real(n,dp)
    varx=sum((y-meanx)**2)/real(n,dp)
    shape=max(1.0e-8_dp,meanx*meanx/max(varx,1.0e-12_dp))
    call fit_gamma_shape(y,shape)
    scale=meanx/shape
    if(present(shape_hat))shape_hat=shape
    if(present(scale_hat))scale_hat=scale
    mu3=scale**(1.0_dp/3.0_dp)*exp(log_gamma_fn(shape+1.0_dp/3.0_dp)-log_gamma_fn(shape))
    s3=sqrt(max(0.0_dp,scale**(2.0_dp/3.0_dp)* &
         exp(log_gamma_fn(shape+2.0_dp/3.0_dp)-log_gamma_fn(shape))-mu3*mu3))
    kfac=k_factor(n,a,pp,sd,trim(meth),mm)
    out%alpha=a;out%p=pp;out%estimate=meanx
    out%lower=max(0.0_dp,(mu3-s3*kfac)**3)
    out%upper=(mu3+s3*kfac)**3
    if(lg)then;out%lower=exp(out%lower);out%upper=exp(out%upper);out%estimate=exp(mu3);end if
  contains
    subroutine fit_gamma_shape(v,sh)
      real(dp),intent(in)::v(:)
      real(dp),intent(inout)::sh
      real(dp)::target,g,h,newsh
      integer::it
      target=log(sum(v)/real(size(v),dp))-sum(log(v))/real(size(v),dp)
      do it=1,100
        g=log(sh)-digamma(sh)-target
        h=1.0_dp/sh-polygamma_local(1,sh)
        newsh=sh-g/h
        if(newsh<=0.0_dp)newsh=0.5_dp*sh
        if(abs(newsh-sh)<1.0e-10_dp*max(1.0_dp,sh))exit
        sh=newsh
      end do
      sh=max(newsh,1.0e-10_dp)
    end subroutine fit_gamma_shape
    real(dp) function polygamma_local(ord,z) result(v)
      use tolerance_math, only : polygamma
      integer,intent(in)::ord; real(dp),intent(in)::z
      v=polygamma(ord,z)
    end function polygamma_local
    real(dp) function log_gamma_fn(z) result(v)
      real(dp),intent(in)::z
      v=log_gamma(z)
    end function log_gamma_fn
  end function gamtol_int

  function logistol_int(x,alpha,p,side,log_log,location_hat,scale_hat) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    logical, intent(in), optional :: log_log
    real(dp), intent(out), optional :: location_hat,scale_hat
    type(tolerance_interval) :: out
    real(dp),allocatable::y(:)
    real(dp)::a,pp,m0,s0,pars(2),hess(2,2),invh(2,2),kd,z,t1,t2,u,v,disc,kl,ku
    integer::sd,info
    logical::ll
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side;ll=.false.;if(present(log_log))ll=log_log
    allocate(y(size(x)));y=x;if(ll)y=log(x)
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    m0=sum(y)/real(size(y),dp);s0=sqrt(max(1.0e-12_dp,3.0_dp*(sum(y*y)/real(size(y),dp)-m0*m0)))/pi
    pars=[m0,log(max(s0,1.0e-8_dp))]
    call nelder_mead(logistic_nll,pars,step=0.1_dp,tol=1.0e-10_dp,max_iter=2000)
    call numerical_hessian(logistic_nll,pars,hess)
    call invert_matrix(hess,invh,info)
    if(info/=0)then
      out%lower=m0;out%upper=m0;out%estimate=m0
    else
      ! Hessian is in (m, log(s)); convert covariance to (m,s).
      invh(1,2)=invh(1,2)*exp(pars(2));invh(2,1)=invh(1,2)
      invh(2,2)=invh(2,2)*exp(2.0_dp*pars(2))
      kd=sqrt(3.0_dp)/pi*log(pp/(1.0_dp-pp));z=normal_quantile(1.0_dp-a)
      t1=kd-invh(1,2)*z*z;t2=kd+invh(1,2)*z*z
      u=kd*kd-invh(1,1)*z*z;v=1.0_dp-invh(2,2)*z*z
      disc=max(0.0_dp,t1*t1-u*v)
      kl=(t1+sqrt(disc))/v;ku=(t2+sqrt(disc))/v
      out%estimate=pars(1)
      out%lower=pars(1)-kl*exp(pars(2))*pi/sqrt(3.0_dp)
      out%upper=pars(1)+ku*exp(pars(2))*pi/sqrt(3.0_dp)
      if(ll)then;out%lower=exp(out%lower);out%upper=exp(out%upper);end if
    end if
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
    if(present(location_hat))location_hat=pars(1)
    if(present(scale_hat))scale_hat=exp(pars(2))
  contains
    real(dp) function logistic_nll(q) result(val)
      real(dp),intent(in)::q(:)
      real(dp)::sc,zv
      integer::i
      sc=exp(q(2));val=0.0_dp
      do i=1,size(y)
        zv=(y(i)-q(1))/sc
        val=val+q(2)+zv+2.0_dp*log(1.0_dp+exp(-zv))
      end do
    end function logistic_nll
  end function logistol_int

  function cautol_int(x,alpha,p,side,location_hat,scale_hat) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    real(dp), intent(out), optional :: location_hat,scale_hat
    type(tolerance_interval) :: out
    real(dp) :: a,pp,pars(2),sig,qc,cf,k
    integer :: sd,n
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    pars(1)=sample_median(x);pars(2)=log(max(1.0e-8_dp,0.5_dp*(sample_quantile(x,0.75_dp)-sample_quantile(x,0.25_dp))))
    call nelder_mead(cauchy_nll,pars,step=0.1_dp,tol=1.0e-10_dp,max_iter=2000)
    sig=exp(pars(2));n=size(x);qc=tan(pi*((1.0_dp-pp)-0.5_dp))
    cf=2.0_dp+2.0_dp*qc*qc;k=sqrt(cf/real(n,dp))*normal_quantile(1.0_dp-a)-qc
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
    out%estimate=pars(1);out%lower=pars(1)-k*sig;out%upper=pars(1)+k*sig
    if(present(location_hat))location_hat=pars(1);if(present(scale_hat))scale_hat=sig
  contains
    real(dp) function cauchy_nll(q) result(val)
      real(dp),intent(in)::q(:)
      real(dp)::sc
      sc=exp(q(2));val=real(size(x),dp)*(log(pi)+q(2))+sum(log(1.0_dp+((x-q(1))/sc)**2))
    end function cauchy_nll
  end function cautol_int

  function paretotol_int(x,alpha,p,side,method,power_dist) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p
    integer, intent(in), optional :: side
    character(len=*), intent(in), optional :: method
    logical, intent(in), optional :: power_dist
    type(tolerance_interval)::out,tmp
    real(dp),allocatable::y(:)
    real(dp)::a,pp,lo,up
    integer::sd
    logical::pw
    character(len=8)::meth
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side;meth='GPU';if(present(method))meth=adjustl(method)
    pw=.false.;if(present(power_dist))pw=power_dist
    allocate(y(size(x)));if(pw)then;y=log(1.0_dp/x);else;y=log(x);end if
    tmp=exp2tol_int(y,a,pp,sd,meth,.false.);lo=tmp%lower;up=tmp%upper
    out=tmp;out%estimate=0.0_dp
    if(pw)then;out%lower=1.0_dp/exp(up);out%upper=1.0_dp/exp(lo)
    else;out%lower=exp(lo);out%upper=exp(up);end if
  end function paretotol_int

  function exttol_int(x,alpha,p,side,dist,ext,nr_delta,shape1,shape2) result(out)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: alpha,p,nr_delta
    integer, intent(in), optional :: side
    character(len=*), intent(in), optional :: dist,ext
    real(dp), intent(out), optional :: shape1,shape2
    type(tolerance_interval)::out
    real(dp),allocatable::y(:)
    real(dp)::a,pp,tol,mx,delta,xbar,xi,xi0,d0,dnew,f,f1,g,g1,d,d1,lam,lamnew
    real(dp)::lpar,kl,ku,aa,bb,old
    integer::sd,n,it
    logical::scaled,ismin,isweib
    character(len=12)::di,ex
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side;tol=1.0e-8_dp;if(present(nr_delta))tol=nr_delta
    di='Weibull';if(present(dist))di=adjustl(dist);ex='min';if(present(ext))ex=adjustl(ext)
    allocate(y(size(x)));y=x;mx=abs(maxval(y))+1000.0_dp;scaled=any(abs(y)>1000.0_dp)
    if(scaled)y=y/mx
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    isweib=trim(di)=='Weibull';if(isweib)then;y=log(y);ex='min';end if
    ismin=isweib .or. (trim(di)=='Gumbel' .and. trim(ex)=='min')
    n=size(y);xbar=sum(y)/real(n,dp);delta=sqrt(max(1.0e-20_dp,(sum(y*y)/real(n,dp)-xbar*xbar)*6.0_dp/pi**2))
    xi=xbar+digamma(1.0_dp)*merge(-1.0_dp,1.0_dp,ismin)
    xi0=xi;d0=delta
    if(ismin)then
      do it=1,200
        f=sum(y*exp(y/delta));f1=-sum(y*y*exp(y/delta))/(delta*delta)
        g=sum(exp(y/delta));g1=-f/(delta*delta)
        d=delta+xbar-f/g;d1=1.0_dp-(g*f1-f*g1)/(g*g)
        if(abs(d1)<1.0e-14_dp)exit
        dnew=delta-d/d1;if(dnew<=0.0_dp .or. dnew/=dnew)then;xi=xi0;delta=d0;exit;end if
        old=delta;delta=dnew;xi=-delta*log(real(n,dp)/sum(exp(y/delta)))
        if(abs(delta-old)<tol .and. abs(d)<tol)exit
      end do
    else
      lam=1.0_dp/delta
      do it=1,200
        f=sum(y*exp(-lam*y));f1=-sum(y*y*exp(-lam*y));g=sum(exp(-lam*y));g1=-f
        d=1.0_dp/lam-xbar+f/g;d1=f*f/(g*g)+f1/g-1.0_dp/(lam*lam)
        if(abs(d1)<1.0e-14_dp)exit
        lamnew=lam-d/d1;if(lamnew<=0.0_dp .or. lamnew/=lamnew)then;xi=xi0;delta=d0;exit;end if
        old=lam;lam=lamnew;xi=-(1.0_dp/lam)*log(sum(exp(-lam*y))/real(n,dp))
        if(abs(lam-old)<tol .and. abs(d)<tol)exit
      end do
      delta=1.0_dp/lam
    end if
    lpar=log(-log(pp));kl=noncentral_t_quantile(1.0_dp-a,real(n-1,dp),-sqrt(real(n,dp))*lpar)
    lpar=log(-log(1.0_dp-pp));ku=noncentral_t_quantile(a,real(n-1,dp),-sqrt(real(n,dp))*lpar)
    if(ismin)then
      out%lower=xi-delta*kl/sqrt(real(n-1,dp));out%upper=xi-delta*ku/sqrt(real(n-1,dp))
    else
      ! MCMCpack/R source uses reversed tail pairs for maxima.
      lpar=log(-log(1.0_dp-pp));kl=noncentral_t_quantile(a,real(n-1,dp),-sqrt(real(n,dp))*lpar)
      lpar=log(-log(pp));ku=noncentral_t_quantile(1.0_dp-a,real(n-1,dp),-sqrt(real(n,dp))*lpar)
      out%lower=xi+delta*kl/sqrt(real(n-1,dp));out%upper=xi+delta*ku/sqrt(real(n-1,dp))
    end if
    aa=xi;bb=delta
    if(isweib)then
      aa=1.0_dp/delta;bb=exp(xi);out%lower=exp(out%lower);out%upper=exp(out%upper)
    end if
    if(scaled)then;bb=bb*mx;out%lower=out%lower*mx;out%upper=out%upper*mx;end if
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2);out%estimate=aa
    if(present(shape1))shape1=aa;if(present(shape2))shape2=bb
  end function exttol_int

end module tolerance_continuous
