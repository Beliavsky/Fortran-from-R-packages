! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_fiducial
  use tolerance_kinds, only : dp
  use tolerance_types, only : tolerance_interval
  use tolerance_math, only : rng_beta, rng_chisq, rng_binomial, rng_poisson, rng_gamma, &
       sample_quantile, rng_normal, normal_quantile, beta_quantile, polygamma, digamma, &
       noncentral_t_quantile, sample_mean, sample_sd, gamma_quantile
  implicit none
  private
  public :: fidbintol_int, fidpoistol_int, fidnegbintol_int, semiconttol_int

  abstract interface
    real(dp) function two_param_fun(x,y) result(v)
      import dp
      real(dp),intent(in)::x,y
    end function two_param_fun
  end interface
contains

  function fidbintol_int(x1,x2,n1,n2,fun,m1,m2,alpha,p,side,kouter,binner) result(out)
    integer,intent(in)::x1,x2,n1,n2
    procedure(two_param_fun)::fun
    integer,intent(in),optional::m1,m2,kouter,binner
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    type(tolerance_interval)::out
    integer::mm1,mm2,kk,bb,i,j,sd
    real(dp)::a,pp,q1,q2
    real(dp),allocatable::lo(:),up(:),v(:)
    mm1=n1;if(present(m1))mm1=m1;mm2=n2;if(present(m2))mm2=m2
    kk=1000;if(present(kouter))kk=kouter;bb=1000;if(present(binner))bb=binner
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    allocate(lo(kk),up(kk),v(bb));out%estimate=fun(real(x1,dp)/n1,real(x2,dp)/n2)
    do i=1,kk
      q1=rng_beta(real(x1,dp)+0.5_dp,real(n1-x1,dp)+0.5_dp)
      q2=rng_beta(real(x2,dp)+0.5_dp,real(n2-x2,dp)+0.5_dp)
      do j=1,bb;v(j)=fun(real(rng_binomial(mm1,q1),dp)/mm1,real(rng_binomial(mm2,q2),dp)/mm2);end do
      lo(i)=sample_quantile(v,1.0_dp-pp);up(i)=sample_quantile(v,pp)
    end do
    out%lower=sample_quantile(lo,a);out%upper=sample_quantile(up,1.0_dp-a)
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
  end function fidbintol_int

  function fidpoistol_int(x1,x2,n1,n2,fun,m1,m2,alpha,p,side,kouter,binner) result(out)
    integer,intent(in)::x1,x2,n1,n2
    procedure(two_param_fun)::fun
    integer,intent(in),optional::m1,m2,kouter,binner
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    type(tolerance_interval)::out
    integer::mm1,mm2,kk,bb,i,j,sd
    real(dp)::a,pp,q1,q2
    real(dp),allocatable::lo(:),up(:),v(:)
    mm1=n1;if(present(m1))mm1=m1;mm2=n2;if(present(m2))mm2=m2
    kk=1000;if(present(kouter))kk=kouter;bb=1000;if(present(binner))bb=binner
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    allocate(lo(kk),up(kk),v(bb));out%estimate=fun(real(x1,dp)/n1,real(x2,dp)/n2)
    do i=1,kk
      q1=rng_chisq(real(2*x1+1,dp))/(2.0_dp*n1);q2=rng_chisq(real(2*x2+1,dp))/(2.0_dp*n2)
      do j=1,bb;v(j)=fun(real(rng_poisson(mm1*q1),dp)/mm1,real(rng_poisson(mm2*q2),dp)/mm2);end do
      lo(i)=sample_quantile(v,1.0_dp-pp);up(i)=sample_quantile(v,pp)
    end do
    out%lower=sample_quantile(lo,a);out%upper=sample_quantile(up,1.0_dp-a)
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
  end function fidpoistol_int

  function fidnegbintol_int(x1,x2,n1,n2,fun,m1,m2,alpha,p,side,kouter,binner) result(out)
    integer,intent(in)::x1,x2,n1,n2
    procedure(two_param_fun)::fun
    integer,intent(in),optional::m1,m2,kouter,binner
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    type(tolerance_interval)::out
    integer::mm1,mm2,kk,bb,i,j,sd,z1,z2
    real(dp)::a,pp,q1,q2
    real(dp),allocatable::lo(:),up(:),v(:)
    mm1=n1;if(present(m1))mm1=m1;mm2=n2;if(present(m2))mm2=m2
    kk=1000;if(present(kouter))kk=kouter;bb=1000;if(present(binner))bb=binner
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    allocate(lo(kk),up(kk),v(bb));out%estimate=fun(real(n1,dp)/(x1+n1),real(n2,dp)/(x2+n2))
    do i=1,kk
      q1=rng_beta(real(n1,dp),real(x1,dp)+0.5_dp);q2=rng_beta(real(n2,dp),real(x2,dp)+0.5_dp)
      do j=1,bb
        z1=rng_negbin(mm1,q1);z2=rng_negbin(mm2,q2)
        v(j)=fun(real(mm1,dp)/(z1+mm1),real(mm2,dp)/(z2+mm2))
      end do
      lo(i)=sample_quantile(v,1.0_dp-pp);up(i)=sample_quantile(v,pp)
    end do
    out%lower=sample_quantile(lo,a);out%upper=sample_quantile(up,1.0_dp-a)
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2)
  end function fidnegbintol_int

  integer function rng_negbin(size,prob) result(x)
    integer,intent(in)::size
    real(dp),intent(in)::prob
    real(dp)::lambda
    lambda=rng_gamma(real(size,dp),(1.0_dp-prob)/prob);x=rng_poisson(lambda)
  end function rng_negbin

  subroutine semiconttol_int(x,zig_ci,zig_pi,zig_ti,zig_ti_approx,ziln_ci,ziln_pi,ziln_ti,ziln_ti_approx, &
       alpha,content,nsim)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::zig_ci(2),zig_pi,zig_ti,zig_ti_approx,ziln_ci(2),ziln_pi,ziln_ti,ziln_ti_approx
    real(dp),intent(in),optional::alpha,content
    integer,intent(in),optional::nsim
    real(dp),allocatable::xp(:),mstar(:),qp(:),zci(:),zqp(:),zpi1(:),zpi2(:)
    real(dp)::a,pv,xbar,xtilde,sx,sxt,tstat,pistar,u,astar,vstar,bstar,eta,u1,z,etamed,ncp
    integer::n,n0,n1,nn,i,j
    a=0.05_dp;if(present(alpha))a=alpha;pv=0.99_dp;if(present(content))pv=content;nn=1000;if(present(nsim))nn=nsim
    n=size(x);n0=count(x==0.0_dp);n1=n-n0;allocate(xp(n1));j=0
    do i=1,n;if(x(i)>0.0_dp)then;j=j+1;xp(j)=x(i);end if;end do
    xbar=sample_mean(xp);xtilde=sum(log(xp))/real(n1,dp);sx=sample_sd(xp)
    sxt=sqrt(sum((log(xp)-xtilde)**2)/real(max(n1-1,1),dp));tstat=xtilde-log(xbar)
    allocate(mstar(nn),qp(nn),zci(nn),zqp(nn),zpi1(nn),zpi2(nn))
    do i=1,nn
      pistar=gfq_pi(n0,n);call random_number(u);astar=find_alpha(u,n1,tstat)
      vstar=rng_chisq(2.0_dp*real(n1,dp)*astar);bstar=vstar/(2.0_dp*real(n1,dp)*xbar)
      mstar(i)=(1.0_dp-pistar)*astar/bstar;eta=(pv-pistar)/(1.0_dp-pistar)
      if(eta>0.0_dp .and. eta<1.0_dp)then;qp(i)=gamma_quantile(eta,astar,1.0_dp/bstar);else;qp(i)=0.0_dp;end if
      z=rng_normal();u1=sqrt(rng_chisq(real(n1-1,dp))/real(max(n1-1,1),dp))
      zci(i)=exp(log(max(1.0_dp-pistar,tiny(1.0_dp)))+xtilde-(z/u1)*(sxt/sqrt(real(n1,dp)))+ &
           0.5_dp*sxt*sxt/(u1*u1))
      if(eta>0.0_dp .and. eta<1.0_dp)then
        zqp(i)=exp(xtilde+sxt*(z+normal_quantile(eta)*sqrt(real(n1,dp)))/(u1*sqrt(real(n1,dp))))
      else;zqp(i)=0.0_dp;end if
      call random_number(u);if(u<pistar)then;zpi1(i)=0.0_dp;zpi2(i)=0.0_dp
      else;zpi1(i)=rng_gamma(astar,1.0_dp/bstar);zpi2(i)=exp(xtilde+sxt*rng_normal());end if
    end do
    zig_ci=[sample_quantile(mstar,a/2.0_dp),sample_quantile(mstar,1.0_dp-a/2.0_dp)]
    zig_ti=sample_quantile(qp,1.0_dp-a);zig_pi=sample_quantile(zpi1,1.0_dp-a)
    ziln_ci=[sample_quantile(zci,a/2.0_dp),sample_quantile(zci,1.0_dp-a/2.0_dp)]
    ziln_ti=sample_quantile(zqp,1.0_dp-a);ziln_pi=sample_quantile(zpi2,1.0_dp-a)
    pistar=beta_quantile(0.5_dp,real(n0,dp)+0.5_dp,real(n1,dp)+0.5_dp);etamed=(pv-pistar)/(1.0_dp-pistar)
    ncp=normal_quantile(min(max(etamed,1.0e-12_dp),1.0_dp-1.0e-12_dp))*sqrt(real(n1,dp))
    zig_ti_approx=(sample_mean(xp**(1.0_dp/3.0_dp))+noncentral_t_quantile(1.0_dp-a,real(n1-1,dp),ncp)* &
         sample_sd(xp**(1.0_dp/3.0_dp))/sqrt(real(n1,dp)))**3
    ziln_ti_approx=exp(xtilde+noncentral_t_quantile(1.0_dp-a,real(n1-1,dp),ncp)*sxt/sqrt(real(n1,dp)))
    if(sx<0.0_dp)zig_pi=zig_pi
  end subroutine semiconttol_int

  real(dp) function gfq_pi(x,n) result(v)
    integer,intent(in)::x,n
    real(dp)::u
    call random_number(u)
    if(x==0)then
      if(u<0.5_dp)then;v=rng_beta(1.0_dp,real(n-1,dp));else;v=0.0_dp;end if
    else if(x<n)then
      if(u<0.5_dp)then;v=rng_beta(real(x+1,dp),real(n-x,dp));else;v=rng_beta(real(x,dp),real(n-x+1,dp));end if
    else
      if(u<0.5_dp)then;v=1.0_dp;else;v=rng_beta(real(n,dp),1.0_dp);end if
    end if
  end function gfq_pi

  real(dp) function find_alpha(lambda,n,tstat) result(root)
    real(dp),intent(in)::lambda,tstat
    integer,intent(in)::n
    real(dp)::lo,hi,mid,v
    integer::it
    lo=1.0e-6_dp;hi=15.0_dp
    do while(diff_perc(hi,lambda,n,tstat)<0.0_dp .and. hi<1000.0_dp);hi=hi+5.0_dp;end do
    do it=1,100
      mid=0.5_dp*(lo+hi);v=diff_perc(mid,lambda,n,tstat)
      if(v>0.0_dp)then;hi=mid;else;lo=mid;end if
      if(abs(hi-lo)<1.0e-6_dp)exit
    end do
    root=0.5_dp*(lo+hi)
  end function find_alpha

  real(dp) function diff_perc(al,lambda,n,tstat) result(res)
    real(dp),intent(in)::al,lambda,tstat
    integer,intent(in)::n
    real(dp)::k1,k2,k3,k4,k5,d3,d4,d5,z,q
    k1=log(real(n,dp))+digamma(al)-digamma(real(n,dp)*al)
    k2=polygamma(1,al)/real(n,dp)-polygamma(1,real(n,dp)*al)
    k3=polygamma(2,al)/real(n*n,dp)-polygamma(2,real(n,dp)*al)
    k4=polygamma(3,al)/real(n**3,dp)-polygamma(3,real(n,dp)*al)
    k5=polygamma(4,al)/real(n**4,dp)-polygamma(4,real(n,dp)*al)
    d3=k3/k2**1.5_dp;d4=k4/k2**2;d5=k5/k2**2.5_dp;z=normal_quantile(lambda)
    q=z+d3*(z*z-1.0_dp)/6.0_dp+d4*(z**3-3*z)/24.0_dp-d3*d3*(2*z**3-5*z)/36.0_dp+ &
      d5*(z**4-6*z*z+3.0_dp)/120.0_dp-d3*d4*(z**4-5*z*z+2.0_dp)/24.0_dp+ &
      d3**3*(12*z**4-53*z*z+17.0_dp)/324.0_dp
    res=k1+sqrt(max(k2,tiny(1.0_dp)))*q-tstat
  end function diff_perc
end module tolerance_fiducial
