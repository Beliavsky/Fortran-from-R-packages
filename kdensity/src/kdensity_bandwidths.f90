module kdensity_bandwidths
  use kdensity_kinds, only : dp
  use kdensity_math, only : mean_value, sample_sd, quantile_type7, normal_quantile, &
    adaptive_integral, golden_minimize, kd_pi
  use kdensity_types, only : kd_start, kd_kernel, kd_ok, kd_invalid_input, kd_optimization_failed
  implicit none
  private
  public :: bandwidth_nrd0, bandwidth_nrd, bandwidth_rhe, bandwidth_jh
  public :: bandwidth_hs, bandwidth_ucv, select_bandwidth, standard_bandwidth_name
contains
  function bandwidth_nrd0(x) result(h)
    real(dp),intent(in)::x(:);real(dp)::h,s,iqr
    if(size(x)<2)then;h=0;return;endif
    s=sample_sd(x);iqr=quantile_type7(x,0.75_dp)-quantile_type7(x,0.25_dp)
    if(iqr>0) s=min(s,iqr/1.34_dp)
    if(s<=0)s=abs(x(1));if(s<=0)s=1
    h=0.9_dp*s*real(size(x),dp)**(-0.2_dp)
  end function
  function bandwidth_nrd(x) result(h)
    real(dp),intent(in)::x(:);real(dp)::h
    h=1.06_dp*sample_sd(x)*real(size(x),dp)**(-0.2_dp)
  end function
  function bandwidth_jh(x) result(h)
    real(dp),intent(in)::x(:);real(dp)::h,mu,sigma
    real(dp),allocatable::z(:);integer::i,n
    allocate(z(size(x)));n=0
    do i=1,size(x)
      if(x(i)>0.and.x(i)<1)then;n=n+1;z(n)=normal_quantile(x(i));endif
    enddo
    if(n<2)then;h=0.5_dp;return;endif
    mu=mean_value(z(:n));sigma=sample_sd(z(:n))
    h=min(sigma*(2*mu*mu*sigma*sigma+3*(1-sigma*sigma)**2)**(-0.2_dp)*real(n,dp)**(-0.2_dp),0.5_dp)
  end function
  function bandwidth_rhe(x) result(h)
    real(dp),intent(in)::x(:);real(dp)::h,mu,sigma,z,u,u2,w,d(4);integer::i,n
    n=size(x);mu=mean_value(x);sigma=sample_sd(x);d=0
    if(n<2.or.sigma<=0)then;h=0;return;endif
    do i=1,n
      z=(x(i)-mu)/sigma;u=sqrt(2.0_dp)*z;u2=u*u;w=sqrt(2.0_dp)*exp(-z*z/2)
      d(1)=d(1)+w*(u2-1);d(2)=d(2)+w*u*(u2-3)
      d(3)=d(3)+w*(u2*(u2-6)+3);d(4)=d(4)+w*u*(u2*(u2-10)+15)
    enddo
    d=d/real(n,dp)
    h=(0.25_dp)**0.2_dp*(d(1)**2+d(2)**2+d(3)**2/2+d(4)**2/6)**(-0.2_dp)*sigma*real(n,dp)**(-0.2_dp)
  end function
  function bandwidth_hs(x,status,used_fallback) result(h)
    real(dp),intent(in)::x(:);integer,intent(out),optional::status;logical,intent(out),optional::used_fallback
    real(dp)::h,mu,var,common,a,b,lnum,lden,bvar,bsk,bku,scale,correction
    real(dp),allocatable::z(:);integer::n
    if(present(status))status=kd_ok;if(present(used_fallback))used_fallback=.false.
    z=pack(x,x>0.and.x<1);n=size(z)
    if(n==0)then;h=0;if(present(status))status=kd_invalid_input;return;endif
    mu=mean_value(z);var=sum((z-mu)**2)/real(n,dp)
    if(var<=0)then;h=0;if(present(status))status=kd_invalid_input;return;endif
    common=mu*(1-mu)/var-1;a=mu*common;b=(1-mu)*common
    if(a>1.5_dp.and.b>1.5_dp.and.a+b>3.and.(a-1)*(b-1)>0.and.6-4*b+a*(3*b-4)>0)then
      lnum=log(2*a+2*b-5)+log(2*a+2*b-3)+log_gamma(2*a+2*b-6)+log_gamma(a)+log_gamma(b)+log_gamma(a-0.5_dp)+log_gamma(b-0.5_dp)
      lden=log((a-1)*(b-1))+log(6-4*b+a*(3*b-4))+log_gamma(2*a-3)+log_gamma(2*b-3)+log_gamma(a+b)+log_gamma(a+b-1)
      h=exp((2.0_dp/5)*(lnum-lden-log(2.0_dp)-log(real(n,dp))-0.5_dp*log(kd_pi)))
    else
      if(present(used_fallback))used_fallback=.true.
      if(a<=0.or.b<=0)then;h=bandwidth_nrd0(z);return;endif
      bvar=a*b/((a+b)**2*(a+b+1));bsk=2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b))
      bku=6*((a-b)**2*(a+b+1)-a*b*(a+b+2))/(a*b*(a+b+2)*(a+b+3))
      scale=sqrt(bvar);correction=1+abs(bsk)+abs(bku);h=scale/correction*real(n,dp)**(-0.4_dp)
    endif
  end function

  function standard_bandwidth_name(kernel,start) result(name)
    type(kd_kernel),intent(in)::kernel;type(kd_start),intent(in)::start;character(len=32)::name
    logical::uniform_start
    uniform_start=trim(start%name)=='uniform'.or.trim(start%name)=='constant'.or.trim(start%name)=='unif'
    if(trim(kernel%name)=='gcopula'.and.uniform_start)then;name='JH'
    else if(trim(kernel%name)=='beta'.and.uniform_start)then;name='HS'
    else if(.not.uniform_start)then;if(kernel%has_sd)then;name='RHE';else;name='ucv';endif
    else;if(kernel%has_sd)then;name='nrd0';else;name='ucv';endif
    endif
  end function

  function select_bandwidth(x,kernel,start,support,name,status) result(h)
    real(dp),intent(in)::x(:),support(2);type(kd_kernel),intent(in)::kernel;type(kd_start),intent(in)::start
    character(len=*),intent(in)::name;integer,intent(out),optional::status;real(dp)::h
    character(len=:),allocatable::key;integer::istat
    key=lowercase(trim(name));istat=kd_ok
    select case(key)
    case('nrd0');h=bandwidth_nrd0(x)
    case('nrd');h=bandwidth_nrd(x)
    case('rhe');h=bandwidth_rhe(x)
    case('jh');h=bandwidth_jh(x)
    case('hs');h=bandwidth_hs(x,istat)
    case('bcv');h=1.02_dp*bandwidth_nrd(x)
    case('sj');h=bandwidth_nrd(x)
    case('ucv');h=bandwidth_ucv(x,kernel,start,support,istat)
    case default;h=0;istat=kd_invalid_input
    end select
    if(present(status))status=istat
  end function

  function bandwidth_ucv(x,kernel,start,support,status) result(hbest)
    real(dp),intent(in)::x(:),support(2)
    type(kd_kernel),intent(in)::kernel
    type(kd_start),intent(in)::start
    integer,intent(out),optional::status
    real(dp)::hbest,using,lower,upper
    integer::istat
    if(trim(start%name)=='uniform'.and.kernel%has_sd)then
      using=bandwidth_nrd0(x);lower=0.2_dp*using;upper=5*using
    else if(trim(kernel%name)=='gcopula'.or.index(trim(kernel%name),'beta')==1)then
      using=bandwidth_jh(x);lower=0.25_dp*using;upper=0.25_dp-1e-10_dp
    else if(kernel%has_sd)then
      using=bandwidth_rhe(x);lower=0.2_dp*using;upper=5*using
    else
      using=bandwidth_nrd0(x);lower=0.2_dp*using;upper=5*using
    endif
    if(lower<=0.or.upper<=lower)then
      hbest=using
      if(present(status))status=kd_optimization_failed
      return
    endif
    hbest=golden_minimize(objective,lower,upper,1e-4_dp,istat)
    if(present(status))status=merge(kd_ok,kd_optimization_failed,istat==0)
  contains
    function objective(h) result(obj)
      real(dp),intent(in)::h
      real(dp)::obj,term1,term2,fi,di,lo,hi,step,y,raw,den,norm
      real(dp),allocatable::pfull(:),ploo(:),xd(:),raw_grid(:)
      integer::i,j,k,s,n_grid
      call start%estimator(x,pfull,s)
      if(s/=0)then;obj=huge(1.0_dp);return;endif
      lo=max(support(1),minval(x)-6.0_dp*h)
      hi=min(support(2),maxval(x)+6.0_dp*h)
      if(abs(lo)>0.1_dp*huge(1.0_dp))lo=minval(x)-6.0_dp*h
      if(abs(hi)>0.1_dp*huge(1.0_dp))hi=maxval(x)+6.0_dp*h
      n_grid=301;step=(hi-lo)/real(n_grid-1,dp);allocate(raw_grid(n_grid))
      do k=1,n_grid
        y=lo+real(k-1,dp)*step;raw=0.0_dp
        do j=1,size(x)
          den=start%density(x(j),pfull)
          if(den>tiny(1.0_dp))raw=raw+kernel%evaluate(y,x(j),h)/den
        enddo
        raw_grid(k)=raw/real(size(x),dp)/h*start%density(y,pfull)
      enddo
      norm=step*(0.5_dp*raw_grid(1)+sum(raw_grid(2:n_grid-1))+0.5_dp*raw_grid(n_grid))
      if(norm<=tiny(1.0_dp))then;obj=huge(1.0_dp);return;endif
      term1=step*(0.5_dp*raw_grid(1)**2+sum(raw_grid(2:n_grid-1)**2)+0.5_dp*raw_grid(n_grid)**2)/(norm*norm)
      term2=0.0_dp
      do i=1,size(x)
        allocate(xd(size(x)-1));xd=[x(:i-1),x(i+1:)]
        call start%estimator(xd,ploo,s)
        if(s/=0)then;obj=huge(1.0_dp);return;endif
        fi=0.0_dp;di=start%density(x(i),ploo)
        do j=1,size(xd)
          den=start%density(xd(j),ploo)
          if(den>tiny(1.0_dp))fi=fi+kernel%evaluate(x(i),xd(j),h)/den
        enddo
        term2=term2+fi/real(size(xd),dp)/h*di
        deallocate(xd,ploo)
      enddo
      obj=term1-2.0_dp*term2/real(size(x),dp)
    end function objective
  end function bandwidth_ucv

  pure function lowercase(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,c
    do i=1,len(s);c=iachar(s(i:i));out(i:i)=s(i:i);if(c>=65.and.c<=90)out(i:i)=achar(c+32);enddo
  end function
end module kdensity_bandwidths
