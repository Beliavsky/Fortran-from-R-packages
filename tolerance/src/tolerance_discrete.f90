! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_discrete
  use tolerance_kinds, only : dp
  use tolerance_types, only : discrete_tolerance_interval
  use tolerance_math, only : normal_pdf, normal_cdf, normal_quantile, beta_quantile, chisq_quantile, &
    f_quantile, binom_quantile, poisson_quantile, negbin_quantile, negbin_pmf, hypergeom_cdf, &
    hypergeom_quantile, clamp
  use tolerance_distributions, only : pnhyper, qnhyper
  implicit none
  private
  public :: bintol_int, poistol_int, negbintol_int, hypertol_int, neghypertol_int
  public :: uma_upper, acceptance_sampling

  type, public :: acceptance_plan
    integer :: acceptance_limit=0
    integer :: lot_size=0
    integer :: sample_size=0
    real(dp) :: confidence=0.95_dp
    real(dp) :: p=0.99_dp
    real(dp) :: aql=0.01_dp
    real(dp) :: rql=0.02_dp
    real(dp) :: producer_risk=0.0_dp
    real(dp) :: consumer_risk=0.0_dp
  end type acceptance_plan

contains

  function bintol_int(x,n,m,alpha,p,side,method,a1,a2) result(out)
    integer,intent(in)::x,n
    integer,intent(in),optional::m,side
    real(dp),intent(in),optional::alpha,p,a1,a2
    character(len=*),intent(in),optional::method
    type(discrete_tolerance_interval)::out
    real(dp)::aa,pp,phat,k,xt,nt,pt,lp,up,aa1,aa2,ls,us,vhat,z,mu,gamma
    integer::mm,ss
    character(len=4)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    mm=n;if(present(m))mm=m
    ss=1;if(present(side))ss=side
    meth='LS';if(present(method))meth=adjustl(method)
    aa1=0.5_dp;if(present(a1))aa1=a1
    aa2=0.5_dp;if(present(a2))aa2=a2
    if(ss==2)then;aa=aa/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    phat=real(x,dp)/real(n,dp);k=normal_quantile(1.0_dp-aa)
    xt=real(x,dp)+k*k/2.0_dp;nt=real(n,dp)+k*k;pt=xt/nt
    select case(trim(meth))
    case('LS')
      lp=phat-k*sqrt(max(0.0_dp,phat*(1.0_dp-phat)/real(n,dp)))
      up=phat+k*sqrt(max(0.0_dp,phat*(1.0_dp-phat)/real(n,dp)))
    case('WS')
      lp=pt-(k*sqrt(real(n,dp))/(real(n,dp)+k*k))* &
        sqrt(phat*(1.0_dp-phat)+k*k/(4.0_dp*real(n,dp)))
      up=pt+(k*sqrt(real(n,dp))/(real(n,dp)+k*k))* &
        sqrt(phat*(1.0_dp-phat)+k*k/(4.0_dp*real(n,dp)))
    case('AC')
      lp=pt-k*sqrt(pt*(1.0_dp-pt)/nt);up=pt+k*sqrt(pt*(1.0_dp-pt)/nt)
    case('JF')
      lp=beta_quantile(aa,real(x,dp)+aa1,real(n-x,dp)+aa2)
      up=beta_quantile(1.0_dp-aa,real(x,dp)+aa1,real(n-x,dp)+aa2)
    case('CP')
      if(x==0)then
        lp=0.0_dp
      else
        lp=1.0_dp/(1.0_dp+real(n-x+1,dp)*f_quantile(1.0_dp-aa,2.0_dp*real(n-x+1,dp), &
          2.0_dp*real(x,dp))/real(x,dp))
      end if
      if(x==n)then
        up=1.0_dp
      else
        up=1.0_dp/(1.0_dp+real(n-x,dp)/(real(x+1,dp)*f_quantile(1.0_dp-aa, &
          2.0_dp*real(x+1,dp),2.0_dp*real(n-x,dp))))
      end if
    case('AS')
      pt=(real(x,dp)+3.0_dp/8.0_dp)/(real(n,dp)+3.0_dp/4.0_dp)
      lp=sin(asin(sqrt(pt))-0.5_dp*k/sqrt(real(n,dp)))**2
      up=sin(asin(sqrt(pt))+0.5_dp*k/sqrt(real(n,dp)))**2
    case('LO')
      if(x<=0)then;lp=0.0_dp;up=min(1.0_dp,k/sqrt(real(n,dp)));else if(x>=n)then
        lp=max(0.0_dp,1.0_dp-k/sqrt(real(n,dp)));up=1.0_dp
      else
        z=log(real(x,dp)/real(n-x,dp));vhat=real(n,dp)/real(x*(n-x),dp)
        ls=z-k*sqrt(vhat);us=z+k*sqrt(vhat)
        lp=exp(ls)/(1.0_dp+exp(ls));up=exp(us)/(1.0_dp+exp(us))
      end if
    case('PR')
      if(phat<=0.0_dp)then;lp=0.0_dp;up=0.0_dp
      else if(phat>=1.0_dp)then;lp=1.0_dp;up=1.0_dp
      else
        z=normal_quantile(phat)
        vhat=sqrt(phat*(1.0_dp-phat)/(real(n,dp)*normal_pdf(z)**2))
        lp=normal_cdf(z-k*vhat);up=normal_cdf(z+k*vhat)
      end if
    case('PO')
      if(phat<=0.0_dp)then;lp=0.0_dp;up=0.0_dp
      else
        mu=-log(phat);vhat=k*sqrt((1.0_dp-phat)/(real(n,dp)*phat))
        lp=exp(-(mu+vhat));up=exp(-(mu-vhat))
      end if
    case('CL')
      if(phat<=0.0_dp .or. phat>=1.0_dp)then;lp=phat;up=phat
      else
        mu=-log(phat);gamma=log(mu)
        vhat=k*sqrt((1.0_dp-phat)/(real(n,dp)*phat*mu*mu))
        lp=exp(-exp(gamma+vhat));up=exp(-exp(gamma-vhat))
      end if
    case('CC')
      lp=phat-k*sqrt(phat*(1.0_dp-phat)/real(n,dp))-1.0_dp/(2.0_dp*real(n,dp))
      up=phat+k*sqrt(phat*(1.0_dp-phat)/real(n,dp))+1.0_dp/(2.0_dp*real(n,dp))
    case('CWS')
      lp=(2.0_dp*real(n,dp)*phat+k*k-1.0_dp-k*sqrt(max(0.0_dp,k*k-2.0_dp-1.0_dp/real(n,dp)+ &
        4.0_dp*phat*(real(n,dp)*(1.0_dp-phat)+1.0_dp))))/(2.0_dp*(real(n,dp)+k*k))
      up=(2.0_dp*real(n,dp)*phat+k*k+1.0_dp+k*sqrt(max(0.0_dp,k*k+2.0_dp-1.0_dp/real(n,dp)+ &
        4.0_dp*phat*(real(n,dp)*(1.0_dp-phat)-1.0_dp))))/(2.0_dp*(real(n,dp)+k*k))
    case default
      lp=phat-k*sqrt(phat*(1.0_dp-phat)/real(n,dp));up=phat+k*sqrt(phat*(1.0_dp-phat)/real(n,dp))
    end select
    lp=clamp(lp,0.0_dp,1.0_dp);up=clamp(up,0.0_dp,1.0_dp)
    out%lower=binom_quantile(1.0_dp-pp,mm,lp);out%upper=binom_quantile(pp,mm,up)
    if(ss==2)then;aa=2.0_dp*aa;pp=2.0_dp*pp-1.0_dp;end if
    out%alpha=aa;out%p=pp;out%estimate=phat
  end function bintol_int

  function poistol_int(x,n,m,alpha,p,side,method) result(out)
    integer,intent(in)::x,n
    integer,intent(in),optional::m,side
    real(dp),intent(in),optional::alpha,p
    character(len=*),intent(in),optional::method
    type(discrete_tolerance_interval)::out
    integer::mm,ss
    real(dp)::aa,pp,lam,k,ll,ul,temp,g
    character(len=4)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    mm=n;if(present(m))mm=m;ss=1;if(present(side))ss=side
    meth='TAB';if(present(method))meth=adjustl(method)
    if(ss==2)then;aa=aa/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    lam=real(x,dp)/real(n,dp);k=normal_quantile(1.0_dp-aa)
    select case(trim(meth))
    case('TAB')
      if(x==0)then;ll=0.0_dp;else;ll=0.5_dp*chisq_quantile(aa,2.0_dp*real(x,dp))/real(n,dp);end if
      ul=0.5_dp*chisq_quantile(1.0_dp-aa,2.0_dp*real(x+1,dp))/real(n,dp)
    case('LS')
      ll=lam-k*sqrt(real(x,dp))/real(n,dp);ul=lam+k*sqrt(real(x,dp))/real(n,dp)
    case('SC')
      ll=lam+k*k/(2.0_dp*real(n,dp))-k/sqrt(real(n,dp))*sqrt(lam+k*k/(4.0_dp*real(n,dp)))
      ul=lam+k*k/(2.0_dp*real(n,dp))+k/sqrt(real(n,dp))*sqrt(lam+k*k/(4.0_dp*real(n,dp)))
    case('CC')
      ll=lam-k*sqrt(real(x,dp))/real(n,dp)-0.5_dp/real(n,dp)
      ul=lam+k*sqrt(real(x,dp))/real(n,dp)+0.5_dp/real(n,dp)
    case('VS')
      ll=lam+k*k/(4.0_dp*real(n,dp))-k*sqrt(real(x,dp))/real(n,dp)
      ul=lam+k*k/(4.0_dp*real(n,dp))+k*sqrt(real(x,dp))/real(n,dp)
    case('RVS')
      ll=lam+k*k/(4.0_dp*real(n,dp))-k*sqrt((lam+3.0_dp/8.0_dp)/real(n,dp))
      ul=lam+k*k/(4.0_dp*real(n,dp))+k*sqrt((lam+3.0_dp/8.0_dp)/real(n,dp))
    case('FT')
      temp=sqrt(lam)+sqrt(lam+1.0_dp)-k/sqrt(real(n,dp))
      if(temp>=1.0_dp)then;g=((temp*temp-1.0_dp)/(2.0_dp*temp))**2;ll=g;else;ll=0.0_dp;end if
      temp=sqrt(lam)+sqrt(lam+1.0_dp)+k/sqrt(real(n,dp));ul=((temp*temp-1.0_dp)/(2.0_dp*temp))**2
    case('CSC')
      temp=lam-1.0_dp/(2.0_dp*real(n,dp))+k*k/(2.0_dp*real(n,dp))
      ll=temp-sqrt(max(0.0_dp,temp*temp-lam*lam+lam/real(n,dp)-1.0_dp/(4.0_dp*real(n*n,dp))))
      temp=lam+1.0_dp/(2.0_dp*real(n,dp))+k*k/(2.0_dp*real(n,dp))
      ul=temp+sqrt(max(0.0_dp,temp*temp-lam*lam-lam/real(n,dp)-1.0_dp/(4.0_dp*real(n*n,dp))))
    case default
      ll=lam;ul=lam
    end select
    ll=max(0.0_dp,ll)
    out%lower=poisson_quantile(1.0_dp-pp,real(mm,dp)*ll);out%upper=poisson_quantile(pp,real(mm,dp)*ul)
    if(ss==2)then;aa=2.0_dp*aa;pp=2.0_dp*pp-1.0_dp;end if
    out%alpha=aa;out%p=pp;out%estimate=lam
  end function poistol_int

  function negbintol_int(x,n,m,alpha,p,side,method) result(out)
    integer,intent(in)::x,n
    integer,intent(in),optional::m,side
    real(dp),intent(in),optional::alpha,p
    character(len=*),intent(in),optional::method
    type(discrete_tolerance_interval)::out
    integer::mm,ss,iter
    real(dp)::aa,pp,nu,nut,z,se,lnu,unu,t1,t2,lo,hi,mid,crit
    character(len=4)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    mm=n;if(present(m))mm=m;ss=1;if(present(side))ss=side
    meth='LS';if(present(method))meth=adjustl(method)
    if(ss==2)then;aa=aa/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    nu=real(n,dp)/real(n+x,dp);nut=real(n-1,dp)/real(n+x-1,dp)
    z=normal_quantile(1.0_dp-aa);se=sqrt(nu*nu*(1.0_dp-nu)/real(n,dp))
    select case(trim(meth))
    case('LS')
      lnu=max(1.0e-7_dp,nu-z*se);unu=min(1.0_dp,nu+z*se)
    case('WU')
      se=sqrt(nut*(1.0_dp-nut)/real(n+x-2,dp));lnu=max(1.0e-7_dp,nut-z*se);unu=min(1.0_dp,nut+z*se)
    case('CB')
      lnu=real(n,dp)/(real(n,dp)+real(x+1,dp)*f_quantile(1.0_dp-aa,2.0_dp*real(x+1,dp),2.0_dp*real(n,dp)))
      if(x==0)then;unu=1.0_dp;else
        t1=f_quantile(1.0_dp-aa,2.0_dp*real(n,dp),2.0_dp*real(x,dp))
        unu=real(n,dp)*t1/(real(n,dp)*t1+real(x,dp))
      end if
    case('CS')
      lnu=max(1.0e-7_dp,chisq_quantile(aa,2.0_dp*real(n,dp))/(2.0_dp*real(n+x,dp)))
      unu=min(1.0_dp,chisq_quantile(1.0_dp-aa,2.0_dp*real(n,dp))/(2.0_dp*real(n+x,dp)))
    case('SC')
      t1=2.0_dp*real(n+x,dp)*real(n,dp)-real(n,dp)*z*z
      t2=sqrt(max(0.0_dp,real(n*n,dp)*z**4-4.0_dp*real(n+x,dp)*real(n*n,dp)*z*z+ &
        4.0_dp*real((n+x)*(n+x)*n,dp)*z*z))
      lnu=max(1.0e-7_dp,(t1-t2)/(2.0_dp*real((n+x)*(n+x),dp)))
      unu=min(1.0_dp,(t1+t2)/(2.0_dp*real((n+x)*(n+x),dp)))
    case('LR')
      crit=chisq_quantile(1.0_dp-aa,1.0_dp)
      lo=1.0e-10_dp;hi=nu
      do iter=1,70;mid=0.5_dp*(lo+hi);if(lr_eq(mid)>0.0_dp)then;lo=mid;else;hi=mid;end if;end do
      lnu=0.5_dp*(lo+hi);lo=nu;hi=1.0_dp-1.0e-12_dp
      do iter=1,70;mid=0.5_dp*(lo+hi);if(lr_eq(mid)<0.0_dp)then;lo=mid;else;hi=mid;end if;end do
      unu=0.5_dp*(lo+hi)
    case('SP')
      if(x==0)then;lnu=1.0_dp;unu=1.0_dp
      else
        lo=1.0e-7_dp;hi=nu
        do iter=1,70;mid=0.5_dp*(lo+hi);if(sp_cdf(mid)>aa)then;hi=mid;else;lo=mid;end if;end do
        lnu=0.5_dp*(lo+hi);lo=nu;hi=1.0_dp-1.0e-8_dp
        do iter=1,70;mid=0.5_dp*(lo+hi);if(sp_cdf(mid)>1.0_dp-aa)then;hi=mid;else;lo=mid;end if;end do
        unu=0.5_dp*(lo+hi)
      end if
    case('CC')
      lnu=max(1.0e-7_dp,nu-z*se-0.5_dp/real(n+x,dp));unu=min(1.0_dp,nu+z*se+0.5_dp/real(n+x,dp))
    case default
      lnu=max(1.0e-7_dp,nu-z*se);unu=min(1.0_dp,nu+z*se)
    end select
    out%lower=negbin_quantile(1.0_dp-pp,real(mm,dp),unu)
    if(lnu<=1.0e-7_dp)then;out%upper=huge(out%upper);out%upper_infinite=.true.
    else;out%upper=negbin_quantile(pp,real(mm,dp),lnu);end if
    if(ss==2)then;aa=2.0_dp*aa;pp=2.0_dp*pp-1.0_dp;end if
    out%alpha=aa;out%p=pp;out%estimate=nu
  contains
    real(dp) function lr_eq(prob) result(v)
      real(dp),intent(in)::prob
      v=2.0_dp*(log(max(negbin_pmf(x,real(n,dp),nu),tiny(1.0_dp))) - &
        log(max(negbin_pmf(x,real(n,dp),prob),tiny(1.0_dp))))-crit
    end function lr_eq
    real(dp) function sp_cdf(prob) result(v)
      real(dp),intent(in)::prob
      real(dp)::theta,ky,d1,d2
      theta=log(real(x,dp)/(real(n+x,dp)*(1.0_dp-prob)))
      ky=real(n,dp)*(log(prob)-log(1.0_dp-(1.0_dp-prob)*exp(theta)))
      d1=sign(sqrt(2.0_dp*abs(theta*real(x,dp)-ky)),theta)
      d2=theta*sqrt(real(x,dp)*(1.0_dp+real(x,dp)/real(n,dp)))
      if(abs(d1)<1.0e-10_dp .or. abs(d2)<1.0e-10_dp)then;v=0.5_dp
      else;v=normal_cdf(d1)-normal_pdf(d1)*(1.0_dp/d2-1.0_dp/d1);end if
    end function sp_cdf
  end function negbintol_int

  function hypertol_int(x,n,n_pop,m,alpha,p,side,method) result(out)
    integer,intent(in)::x,n,n_pop
    integer,intent(in),optional::m,side
    real(dp),intent(in),optional::alpha,p
    character(len=*),intent(in),optional::method
    type(discrete_tolerance_interval)::out
    integer::mm,ss,ml,mu,j
    real(dp)::aa,pp,phat,z,fpc,lp,up
    character(len=4)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    mm=n;if(present(m))mm=m;ss=1;if(present(side))ss=side
    meth='EX';if(present(method))meth=adjustl(method)
    if(ss==2)then;aa=aa/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    phat=real(x,dp)/real(n,dp);z=normal_quantile(1.0_dp-aa)
    if(trim(meth)=='EX')then
      ml=x
      do j=x,n_pop
        if(1.0_dp-hypergeom_cdf(x-1,j,n_pop-j,n)>aa)then;ml=j;exit;end if
      end do
      mu=n_pop
      do j=n_pop,x,-1
        if(hypergeom_cdf(x,j,n_pop-j,n)>aa)then;mu=j;exit;end if
      end do
    else
      fpc=sqrt(real(n_pop-n,dp)/real(n_pop-1,dp))
      lp=phat-z*sqrt(phat*(1.0_dp-phat)/real(n,dp))*fpc
      up=phat+z*sqrt(phat*(1.0_dp-phat)/real(n,dp))*fpc
      if(trim(meth)=='CC')then;lp=lp-1.0_dp/(2.0_dp*real(n,dp));up=up+1.0_dp/(2.0_dp*real(n,dp));end if
      ml=max(0,floor(real(n_pop,dp)*clamp(lp,0.0_dp,1.0_dp)))
      mu=min(n_pop,ceiling(real(n_pop,dp)*clamp(up,0.0_dp,1.0_dp)))
    end if
    out%lower=hypergeom_quantile(1.0_dp-pp,ml,n_pop-ml,mm)
    out%upper=hypergeom_quantile(pp,mu,n_pop-mu,mm)
    if(ss==2)then;aa=2.0_dp*aa;pp=2.0_dp*pp-1.0_dp;end if
    out%alpha=aa;out%p=pp;out%estimate=phat
  end function hypertol_int

  function neghypertol_int(x,n,n_pop,m,alpha,p,side,method) result(out)
    integer,intent(in)::x,n,n_pop
    integer,intent(in),optional::m,side
    real(dp),intent(in),optional::alpha,p
    character(len=*),intent(in),optional::method
    type(discrete_tolerance_interval)::out
    integer::mm,ss,ml,mu,j
    real(dp)::aa,pp,nu,z,fpc,se,lp,up
    character(len=4)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    mm=n;if(present(m))mm=m;ss=1;if(present(side))ss=side
    meth='EX';if(present(method))meth=adjustl(method)
    if(ss==2)then;aa=aa/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    nu=real(n,dp)/real(x,dp);z=normal_quantile(1.0_dp-aa)
    if(trim(meth)=='EX')then
      ml=n
      do j=n,n_pop-(x-n)
        if(pnhyper(x,j,n_pop,n)>aa)then;ml=j;exit;end if
      end do
      mu=n_pop-(x-n)
      do j=mu,n,-1
        if(1.0_dp-pnhyper(x-1,j,n_pop,n)>aa)then;mu=j;exit;end if
      end do
    else
      fpc=sqrt(real(n_pop-x,dp)/real(n_pop-1,dp));se=sqrt(nu*nu*(1.0_dp-nu)/real(x,dp))*fpc
      lp=nu-z*se;up=nu+z*se
      if(trim(meth)=='CC')then;lp=lp-1.0_dp/(2.0_dp*real(x,dp));up=up+1.0_dp/(2.0_dp*real(x,dp));end if
      ml=max(mm,floor(real(n_pop,dp)*clamp(lp,1.0e-7_dp,1.0_dp)))
      mu=max(mm,min(n_pop,ceiling(real(n_pop,dp)*clamp(up,1.0e-7_dp,1.0_dp))))
    end if
    out%lower=qnhyper(1.0_dp-pp,mu,n_pop,mm);out%upper=qnhyper(pp,ml,n_pop,mm)
    if(ss==2)then;aa=2.0_dp*aa;pp=2.0_dp*pp-1.0_dp;end if
    out%alpha=aa;out%p=pp;out%estimate=nu
  end function neghypertol_int

  integer function uma_upper(xsum,n,nfuture,dist,alpha,p) result(k)
    integer,intent(in)::xsum,n,nfuture
    character(len=*),intent(in)::dist
    real(dp),intent(in),optional::alpha,p
    real(dp)::aa,pp,r,temp
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    select case(trim(adjustl(dist)))
    case('Bin')
      if(xsum>0)then
        r=max(1.0_dp-beta_quantile(aa,real(nfuture*n-xsum+1,dp),real(xsum,dp)), &
              1.0_dp-beta_quantile(aa,real(nfuture*n-xsum,dp),real(xsum+1,dp)))
      else;r=1.0_dp-aa**(1.0_dp/real(nfuture*n,dp));end if
      k=-1;temp=-1.0_dp
      do while(temp<r);k=k+1;if(k<nfuture)then
        temp=beta_quantile(1.0_dp-pp,real(k+1,dp),real(nfuture-k,dp));else;temp=1.0_dp;end if;end do
    case('NegBin')
      if(xsum>0)then
        r=max(1.0_dp-beta_quantile(aa,real(n*nfuture,dp),real(xsum+1,dp)), &
              1.0_dp-beta_quantile(aa,real(n*nfuture,dp),real(xsum,dp)))
      else;r=1.0_dp-aa**(1.0_dp/real(n*nfuture,dp));end if
      k=-1;temp=1.1_dp
      do while(temp>1.0_dp-r);k=k+1;temp=beta_quantile(pp,real(nfuture,dp),real(k+1,dp));end do
    case('Pois')
      if(xsum>0)then
        r=2.0_dp*max(chisq_quantile(1.0_dp-aa,2.0_dp*real(xsum+1,dp))/(2.0_dp*real(n,dp)), &
          chisq_quantile(1.0_dp-aa,2.0_dp*real(xsum,dp))/(2.0_dp*real(n,dp)))
      else;r=-log(aa)/real(n,dp);end if
      k=-1;temp=-1.0_dp
      do while(temp<r);k=k+1;temp=chisq_quantile(1.0_dp-pp,2.0_dp*real(k+1,dp));end do
    case default
      k=-1
    end select
  end function uma_upper

  function acceptance_sampling(n,n_pop,alpha,p,aql,rql) result(plan)
    integer,intent(in)::n,n_pop
    real(dp),intent(in),optional::alpha,p,aql,rql
    type(acceptance_plan)::plan
    real(dp)::aa,pp,qa,qr,d
    integer::c,mh,nh
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(p))pp=p
    qa=0.01_dp;if(present(aql))qa=aql;qr=0.02_dp;if(present(rql))qr=rql
    d=real(n_pop,dp)*(1.0_dp-pp);mh=nint(d);nh=n_pop-mh;c=0
    do while(c<=ceiling(d))
      if(aa-hypergeom_cdf(c,mh,nh,n)<=0.0_dp)exit
      c=c+1
    end do
    c=max(0,c)
    if(hypergeom_cdf(c,mh,nh,n)>aa)c=max(c-1,0)
    plan%acceptance_limit=c;plan%lot_size=n_pop;plan%sample_size=n
    plan%confidence=1.0_dp-aa;plan%p=pp;plan%aql=qa;plan%rql=qr
    plan%producer_risk=1.0_dp-hypergeom_cdf(c,floor(qa*real(n_pop,dp)),n_pop-floor(qa*real(n_pop,dp)),n)
    plan%consumer_risk=hypergeom_cdf(c,floor(qr*real(n_pop,dp)),n_pop-floor(qr*real(n_pop,dp)),n)
  end function acceptance_sampling

end module tolerance_discrete
