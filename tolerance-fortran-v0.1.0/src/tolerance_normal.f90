! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_normal
  use tolerance_kinds, only : dp, pi
  use tolerance_types, only : tolerance_interval
  use tolerance_math, only : normal_pdf, normal_cdf, normal_quantile, chisq_cdf, chisq_quantile, &
    noncentral_chisq_quantile, noncentral_t_quantile, sample_mean, sample_sd, sample_variance, &
    adaptive_simpson
  implicit none
  private

  public :: k_factor, k_factor_sim
  public :: normtol_int, bayesnormtol_int, simnormtol_int
  public :: diffnormtol_int

contains

  real(dp) function k_factor(n,alpha,p,side,method,f,m) result(k)
    integer,intent(in)::n
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side,f,m
    character(len=*),intent(in),optional::method
    real(dp)::aa,pp,zp,ncp,ta,chi_a,k2,za,dfcut,v,g,r,delta,pnew,diffv
    real(dp)::lo,hi,mid,current_xk,current_delta
    integer::ss,ff,mm,iter
    character(len=8)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    ss=1;if(present(side))ss=side
    ff=n-1;if(present(f))ff=f
    mm=50;if(present(m))mm=m
    meth='HE';if(present(method))meth=adjustl(method)

    if(ss==1)then
      zp=normal_quantile(pp);ncp=sqrt(real(n,dp))*zp
      ta=noncentral_t_quantile(1.0_dp-aa,real(ff,dp),ncp)
      k=ta/sqrt(real(n,dp));return
    end if

    chi_a=chisq_quantile(aa,real(ff,dp))
    k2=sqrt(real(ff,dp)*noncentral_chisq_quantile(pp,1.0_dp,1.0_dp/real(n,dp))/chi_a)
    select case(trim(meth))
    case('HE')
      zp=normal_quantile((1.0_dp+pp)/2.0_dp)
      za=normal_quantile((2.0_dp-aa)/2.0_dp)
      dfcut=real(n*n,dp)*(1.0_dp+1.0_dp/(za*za))
      v=1.0_dp+za*za/real(n,dp)+(3.0_dp-zp*zp)*za**4/(6.0_dp*real(n*n,dp))
      if(real(ff,dp)>dfcut)then
        k=zp*sqrt(v*(1.0_dp+(real(n,dp)*v/(2.0_dp*real(ff,dp)))*(1.0_dp+1.0_dp/(za*za))))
      else
        g=(real(ff-2,dp)-chi_a)/(2.0_dp*real((n+1)*(n+1),dp))
        k=zp*sqrt((real(ff,dp)*(1.0_dp+1.0_dp/real(n,dp))/chi_a)*(1.0_dp+g))
      end if
    case('HE2')
      zp=normal_quantile((1.0_dp+pp)/2.0_dp)
      k=zp*sqrt((1.0_dp+1.0_dp/real(n,dp))*real(ff,dp)/chi_a)
    case('WBE')
      r=0.5_dp
      do iter=1,100
        pnew=normal_cdf(1.0_dp/sqrt(real(n,dp))+r)-normal_cdf(1.0_dp/sqrt(real(n,dp))-r)
        delta=pnew-pp
        diffv=normal_pdf(1.0_dp/sqrt(real(n,dp))+r)+normal_pdf(1.0_dp/sqrt(real(n,dp))-r)
        r=r-delta/diffv
        if(abs(delta)<1.0e-10_dp)exit
      end do
      k=r*sqrt(real(ff,dp)/chi_a)
    case('ELL')
      r=0.5_dp;zp=normal_quantile((1.0_dp+pp)/2.0_dp)
      do iter=1,100
        pnew=normal_cdf(zp/sqrt(real(n,dp))+r)-normal_cdf(zp/sqrt(real(n,dp))-r)
        delta=pnew-pp
        diffv=normal_pdf(zp/sqrt(real(n,dp))+r)+normal_pdf(zp/sqrt(real(n,dp))-r)
        r=r-delta/diffv
        if(abs(delta)<1.0e-10_dp)exit
      end do
      k=r*sqrt(real(ff,dp)/chi_a)
    case('KM')
      k=k2
    case('EXACT')
      lo=1.0e-10_dp;hi=k2+1000.0_dp/real(n,dp)
      do while(exact_eq(hi)<0.0_dp);hi=2.0_dp*hi;end do
      do iter=1,70
        mid=0.5_dp*(lo+hi)
        if(exact_eq(mid)<0.0_dp)then;lo=mid;else;hi=mid;end if
      end do
      k=0.5_dp*(lo+hi)
    case('OCT')
      lo=1.0e-8_dp;hi=max(5.0_dp,k2*2.0_dp)
      do while(oct_eq(hi)<0.0_dp .and. hi<1.0e5_dp);hi=2.0_dp*hi;end do
      do iter=1,70
        mid=0.5_dp*(lo+hi)
        if(oct_eq(mid)<0.0_dp)then;lo=mid;else;hi=mid;end if
      end do
      k=0.5_dp*(lo+hi)
    case default
      k=k2
    end select

  contains
    real(dp) function exact_eq(xk) result(val)
      real(dp),intent(in)::xk
      real(dp)::upper
      current_xk=xk
      upper=12.0_dp/sqrt(real(n,dp))
      val=sqrt(2.0_dp*real(n,dp)/pi)*adaptive_simpson(exact_integrand,0.0_dp,upper,2.0e-8_dp,20) - &
        (1.0_dp-aa)
    end function exact_eq

    function exact_integrand(z) result(y)
      real(dp),intent(in)::z
      real(dp)::y,qnc,arg
      qnc=noncentral_chisq_quantile(pp,1.0_dp,z*z)
      arg=real(ff,dp)*qnc/(current_xk*current_xk)
      y=(1.0_dp-chisq_cdf(arg,real(ff,dp)))*exp(-0.5_dp*real(n,dp)*z*z)
    end function exact_integrand

    real(dp) function oct_eq(xk) result(val)
      real(dp),intent(in)::xk
      real(dp)::lower,upper
      current_xk=xk
      current_delta=sqrt(real(n,dp))*normal_quantile((1.0_dp+pp)/2.0_dp)
      lower=real(ff,dp)*current_delta*current_delta/(xk*xk*real(n,dp))
      upper=chisq_quantile(1.0_dp-1.0e-11_dp,real(ff,dp))
      if(lower>=upper)then
        val=-(1.0_dp-aa);return
      end if
      val=adaptive_simpson(oct_integrand,lower,upper,2.0e-8_dp,20)-(1.0_dp-aa)
    end function oct_eq

    function oct_integrand(z) result(y)
      real(dp),intent(in)::z
      real(dp)::y,cdfv,pdfv
      cdfv=2.0_dp*normal_cdf(-current_delta+current_xk*sqrt(real(n,dp)*z/real(ff,dp)))-1.0_dp
      pdfv=exp((0.5_dp*real(ff,dp)-1.0_dp)*log(z)-0.5_dp*z - &
        0.5_dp*real(ff,dp)*log(2.0_dp)-log_gamma(0.5_dp*real(ff,dp)))
      y=max(0.0_dp,cdfv)*pdfv
    end function oct_integrand
  end function k_factor

  real(dp) function k_factor_sim(n,l,alpha,p,side,method,m) result(k)
    integer,intent(in)::n
    integer,intent(in),optional::l,side,m
    real(dp),intent(in),optional::alpha,p
    character(len=*),intent(in),optional::method
    integer::ll,ss,mm,df,iter
    real(dp)::aa,pp,chi_a,k2,lo,hi,mid,current_xk
    character(len=8)::meth
    ll=1;if(present(l))ll=l
    ss=1;if(present(side))ss=side
    mm=50;if(present(m))mm=m
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    meth='EXACT';if(present(method))meth=adjustl(method)
    df=n*ll-ll
    if(trim(meth)=='BONF')then
      k=k_factor(n,aa/real(ll,dp),pp,ss,'EXACT',df,mm);return
    end if
    chi_a=chisq_quantile(aa,real(df,dp))
    k2=sqrt(real(df,dp)*noncentral_chisq_quantile(pp,1.0_dp,1.0_dp/real(n,dp))/chi_a)
    lo=1.0e-9_dp;hi=k2+100.0_dp
    if(ss==1)then
      do iter=1,70
        mid=0.5_dp*(lo+hi)
        if(one_eq(mid)<0.0_dp)then;lo=mid;else;hi=mid;end if
      end do
    else
      do iter=1,70
        mid=0.5_dp*(lo+hi)
        if(two_eq(mid)<0.0_dp)then;lo=mid;else;hi=mid;end if
      end do
    end if
    k=0.5_dp*(lo+hi)
  contains
    real(dp) function one_eq(xk) result(v)
      real(dp),intent(in)::xk
      real(dp)::up
      current_xk=xk
      up=chisq_quantile(1.0_dp-1.0e-11_dp,real(df,dp))
      v=adaptive_simpson(one_integrand,1.0e-12_dp,up,2.0e-8_dp,20)-(1.0_dp-aa)
    end function one_eq

    function one_integrand(z) result(y)
      real(dp),intent(in)::z
      real(dp)::y,zp,inside,pdf
      zp=normal_quantile(pp)
      inside=sqrt(real(n,dp))*(current_xk*sqrt(z/real(df,dp))-zp)
      pdf=exp((0.5_dp*real(df,dp)-1.0_dp)*log(z)-0.5_dp*z - &
        0.5_dp*real(df,dp)*log(2.0_dp)-log_gamma(0.5_dp*real(df,dp)))
      y=pdf*normal_cdf(inside)**ll
    end function one_integrand

    real(dp) function two_eq(xk) result(v)
      real(dp),intent(in)::xk
      real(dp)::up
      current_xk=xk
      up=8.0_dp
      v=2.0_dp*real(ll,dp)*adaptive_simpson(two_integrand,0.0_dp,up,2.0e-8_dp,20)-(1.0_dp-aa)
    end function two_eq

    function two_integrand(z) result(y)
      real(dp),intent(in)::z
      real(dp)::y,p1,p2,arg,pdf
      arg=real(df,dp)*noncentral_chisq_quantile(pp,1.0_dp,z*z/real(n,dp))/(current_xk*current_xk)
      p1=1.0_dp-chisq_cdf(arg,real(df,dp))
      p2=(2.0_dp*normal_cdf(z)-1.0_dp)**(ll-1)
      pdf=normal_pdf(z)
      y=pdf*p1*p2
    end function two_integrand
  end function k_factor_sim

  function normtol_int(x,alpha,p,side,method,log_norm) result(out)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::log_norm
    type(tolerance_interval)::out
    real(dp),allocatable::y(:)
    real(dp)::aa,pp,xb,s,k
    integer::ss
    logical::ln
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    ss=1;if(present(side))ss=side
    ln=.false.;if(present(log_norm))ln=log_norm
    y=x;if(ln)y=log(y)
    xb=sample_mean(y);s=sample_sd(y)
    if(present(method))then
      k=k_factor(size(y),aa,pp,ss,method)
    else
      k=k_factor(size(y),aa,pp,ss,'HE')
    end if
    out%alpha=aa;out%p=pp;out%estimate=xb;out%lower=xb-s*k;out%upper=xb+s*k
    if(ln)then;out%estimate=exp(out%estimate);out%lower=exp(out%lower);out%upper=exp(out%upper);end if
  end function normtol_int

  function bayesnormtol_int(xbar,s,n,alpha,p,side,method,mu0,sig20,m0,n0) result(out)
    real(dp),intent(in)::xbar,s
    integer,intent(in)::n
    real(dp),intent(in),optional::alpha,p,mu0,sig20,m0,n0
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::method
    type(tolerance_interval)::out
    real(dp)::aa,pp,k,xb,q2
    integer::ss
    logical::prior
    character(len=8)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    ss=1;if(present(side))ss=side
    meth='HE';if(present(method))meth=adjustl(method)
    prior=present(mu0).and.present(sig20).and.present(m0).and.present(n0)
    if(.not.prior)then
      k=k_factor(n,aa,pp,ss,meth);xb=xbar;q2=s*s
    else
      k=k_factor(n+int(n0),aa,pp,ss,meth,int(m0)+n-1)
      xb=(n0*mu0+real(n,dp)*xbar)/(n0+real(n,dp))
      q2=(m0*sig20+real(n-1,dp)*s*s+n0*real(n,dp)*(xbar-mu0)**2/(n0+real(n,dp))) / &
        (m0+real(n-1,dp))
    end if
    out%alpha=aa;out%p=pp;out%estimate=xb;out%lower=xb-sqrt(q2)*k;out%upper=xb+sqrt(q2)*k
  end function bayesnormtol_int

  subroutine simnormtol_int(x,alpha,p,side,method,lower,upper,means,sp)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::method
    real(dp),allocatable,intent(out)::lower(:),upper(:),means(:)
    real(dp),intent(out),optional::sp
    real(dp)::aa,pp,k,s2pool
    integer::ss,n,l,j
    character(len=8)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    ss=1;if(present(side))ss=side
    meth='EXACT';if(present(method))meth=adjustl(method)
    n=size(x,1);l=size(x,2);allocate(lower(l),upper(l),means(l))
    s2pool=0.0_dp
    do j=1,l;means(j)=sample_mean(x(:,j));s2pool=s2pool+real(n-1,dp)*sample_variance(x(:,j));end do
    s2pool=s2pool/real(n*l-l,dp)
    k=k_factor_sim(n,l,aa,pp,ss,meth)
    lower=means-sqrt(s2pool)*k;upper=means+sqrt(s2pool)*k
    if(present(sp))sp=sqrt(s2pool)
  end subroutine simnormtol_int

  function diffnormtol_int(x1,x2,alpha,p,method,var_ratio) result(out)
    real(dp),intent(in)::x1(:),x2(:)
    real(dp),intent(in),optional::alpha,p,var_ratio
    character(len=*),intent(in),optional::method
    type(tolerance_interval)::out
    integer::n1,n2
    real(dp)::aa,pp,m1,m2,s12,s22,zp,q1,q2,f1,f2,nu1,nu2,kl,ku,sdif
    character(len=8)::meth
    aa=0.05_dp;if(present(alpha))aa=alpha
    pp=0.99_dp;if(present(p))pp=p
    meth='HALL';if(present(method))meth=adjustl(method)
    n1=size(x1);n2=size(x2);m1=sample_mean(x1);m2=sample_mean(x2)
    s12=sample_variance(x1);s22=sample_variance(x2);zp=normal_quantile(pp)
    if(present(var_ratio))then
      q1=var_ratio
      nu1=real(n1,dp)*(1.0_dp+q1)/(q1+real(n1,dp)/real(n2,dp))
      sdif=sqrt((1.0_dp+1.0_dp/q1)*(real(n1-1,dp)*s12+real(n2-1,dp)*q1*s22)/real(n1+n2-2,dp))
      kl=noncentral_t_quantile(1.0_dp-aa,real(n1+n2-2,dp),zp*sqrt(nu1))*sdif/sqrt(nu1)
      ku=kl
    else
      if(trim(meth)=='RG')then
        q1=s12/s22
      else
        q1=s12*real(n2-3,dp)/(s22*real(n2-1,dp))
      end if
      f1=real(n1-1,dp)*(q1+1.0_dp)**2/(q1*q1+real(n1-1,dp)/real(n2-1,dp))
      nu1=real(n1,dp)*(1.0_dp+q1)/(q1+real(n1,dp)/real(n2,dp))
      kl=noncentral_t_quantile(1.0_dp-aa,f1,zp*sqrt(nu1))*sqrt((s12+s22)/nu1)
      ku=kl
      if(trim(meth)=='GK')then
        q2=s22*real(n1-3,dp)/(s12*real(n1-1,dp))
        f2=real(n2-1,dp)*(q2+1.0_dp)**2/(q2*q2+real(n2-1,dp)/real(n1-1,dp))
        nu2=real(n2,dp)*(1.0_dp+q2)/(q2+real(n2,dp)/real(n1,dp))
        sdif=noncentral_t_quantile(1.0_dp-aa,f2,zp*sqrt(nu2))*sqrt((s12+s22)/nu2)
        kl=max(kl,sdif);ku=kl
      end if
    end if
    out%alpha=aa;out%p=pp;out%estimate=m1-m2;out%lower=m1-m2-kl;out%upper=m1-m2+ku
  end function diffnormtol_int

end module tolerance_normal
