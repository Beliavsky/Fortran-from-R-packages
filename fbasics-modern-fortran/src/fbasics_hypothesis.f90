! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_hypothesis
  use fbasics_kinds, only: dp, pi, clamp
  use fbasics_stats, only: sample_mean, sample_variance, sample_sd, sample_skewness, sample_kurtosis
  use fbasics_special, only: normal_cdf, student_cdf, chi_square_cdf, f_cdf
  implicit none
  private
  type, public :: test_result
    character(len=64) :: method=''
    real(dp) :: statistic=0.0_dp
    real(dp) :: p_value=1.0_dp
    real(dp) :: estimate=0.0_dp
    integer :: df=0
  end type
  public :: pearson_test, spearman_test, kendall_test
  public :: location_t_test, variance_f_test, jarque_bera_test
  public :: dagostino_skewness_test, dagostino_kurtosis_test, dagostino_omnibus_test
  public :: ks_two_sample_test, cvm_normal_test, anderson_darling_normal_test
  public :: lilliefors_normal_test, shapiro_francia_test
contains
  function pearson_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r;real(dp)::mx,my,sx,sy,corr,t;integer::n
    n=min(size(x),size(y));mx=sample_mean(x(1:n));my=sample_mean(y(1:n));sx=sample_sd(x(1:n));sy=sample_sd(y(1:n))
    corr=sum((x(1:n)-mx)*(y(1:n)-my))/real(n-1,dp)/(sx*sy);corr=clamp(corr,-1.0_dp,1.0_dp)
    if(n>2.and.abs(corr)<1)then;t=corr*sqrt(real(n-2,dp)/(1-corr*corr));r%p_value=2*(1-student_cdf(abs(t),real(n-2,dp)));else;r%p_value=0;end if
    r%method='Pearson correlation test';r%statistic=corr;r%estimate=corr;r%df=n-2
  end function
  function spearman_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r;real(dp),allocatable::rx(:),ry(:)
    rx=ranks(x);ry=ranks(y);r=pearson_test(rx,ry);r%method='Spearman rank correlation test'
  end function
  function kendall_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r;integer::i,j,n,c,d,ties;real(dp)::tau,z,var0
    n=min(size(x),size(y));c=0;d=0;ties=0
    do i=1,n-1;do j=i+1,n
      if((x(i)-x(j))*(y(i)-y(j))>0)then;c=c+1;else if((x(i)-x(j))*(y(i)-y(j))<0)then;d=d+1;else;ties=ties+1;end if
    end do;end do
    tau=real(c-d,dp)/sqrt(real((c+d+ties)*(c+d+ties),dp));var0=2.0_dp*(2*n+5)/real(9*n*(n-1),dp);z=tau/sqrt(var0)
    r%method='Kendall rank correlation test';r%statistic=tau;r%estimate=tau;r%p_value=2*(1-normal_cdf(abs(z)))
  end function
  function location_t_test(x,y,equal_variance) result(r)
    real(dp),intent(in)::x(:),y(:);logical,intent(in),optional::equal_variance;type(test_result)::r
    real(dp)::mx,my,vx,vy,se,t,dfv,sp;integer::nx,ny;logical::eq
    nx=size(x);ny=size(y);mx=sample_mean(x);my=sample_mean(y);vx=sample_variance(x);vy=sample_variance(y);eq=.false.;if(present(equal_variance))eq=equal_variance
    if(eq)then;sp=((nx-1)*vx+(ny-1)*vy)/real(nx+ny-2,dp);se=sqrt(sp*(1.0_dp/nx+1.0_dp/ny));dfv=nx+ny-2
    else;se=sqrt(vx/nx+vy/ny);dfv=(vx/nx+vy/ny)**2/((vx/nx)**2/(nx-1)+(vy/ny)**2/(ny-1));end if
    t=(mx-my)/se;r%method='Two-sample t test';r%statistic=t;r%estimate=mx-my;r%df=nint(dfv);r%p_value=2*(1-student_cdf(abs(t),dfv))
  end function
  function variance_f_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r;real(dp)::f,p;integer::d1,d2
    f=sample_variance(x)/sample_variance(y);d1=size(x)-1;d2=size(y)-1;p=f_cdf(f,real(d1,dp),real(d2,dp));r%method='Variance F test';r%statistic=f;r%estimate=f;r%df=d1;r%p_value=2*min(p,1-p)
  end function
  function jarque_bera_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp)::m,m2,m3,m4,b1,b2;integer::n
    n=size(x);m=sample_mean(x);m2=sum((x-m)**2)/n;m3=sum((x-m)**3)/n;m4=sum((x-m)**4)/n;b1=(m3/m2**1.5_dp)**2;b2=m4/m2**2
    r%method='Jarque-Bera normality test';r%statistic=n*b1/6.0_dp+n*(b2-3)**2/24.0_dp;r%p_value=1-chi_square_cdf(r%statistic,2.0_dp);r%df=2
  end function
  function dagostino_skewness_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;integer::n;real(dp)::m,s,a3,sd3,u3,b,w2,delta,a,z
    n=size(x);m=sample_mean(x);s=sqrt(sum((x-m)**2)/n);a3=sum((x-m)**3)/n/s**3;sd3=sqrt(6.0_dp*(n-2)/real((n+1)*(n+3),dp));u3=a3/sd3
    b=3.0_dp*(n*n+27*n-70)*(n+1)*(n+3)/real((n-2)*(n+5)*(n+7)*(n+9),dp);w2=sqrt(2.0_dp*(b-1))-1;delta=1/sqrt(log(sqrt(w2)));a=sqrt(2/(w2-1));z=delta*asinh(u3/a)
    r%method="D'Agostino skewness test";r%statistic=z;r%p_value=2*(1-normal_cdf(abs(z)))
  end function
  function dagostino_kurtosis_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;integer::n;real(dp)::m,s,a4,sd4,u4,b,a,jm,pos0,pos,z
    n=size(x);m=sample_mean(x);s=sqrt(sum((x-m)**2)/n);a4=sum((x-m)**4)/n/s**4;sd4=sqrt(24.0_dp*(n-2)*(n-3)*n/real((n+1)**2*(n+3)*(n+5),dp));u4=(a4-3+6.0_dp/(n+1))/sd4
    b=6.0_dp*(n*n-5*n+2)/real((n+7)*(n+9),dp)*sqrt(6.0_dp*(n+3)*(n+5)/real(n*(n-2)*(n-3),dp));a=6+8/b*(2/b+sqrt(1+4/b**2));jm=sqrt(2/(9*a));pos0=(1-2/a)/(1+u4*sqrt(2/(a-4)));pos=sign(abs(pos0)**(1.0_dp/3),pos0);z=(1-2/(9*a)-pos)/jm
    r%method="D'Agostino kurtosis test";r%statistic=z;r%p_value=2*(1-normal_cdf(abs(z)))
  end function
  function dagostino_omnibus_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;type(test_result)::s,k
    s=dagostino_skewness_test(x);k=dagostino_kurtosis_test(x);r%method="D'Agostino omnibus test";r%statistic=s%statistic**2+k%statistic**2;r%p_value=1-chi_square_cdf(r%statistic,2.0_dp);r%df=2
  end function
  function ks_two_sample_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:)
    type(test_result)::r
    real(dp),allocatable::sx(:),sy(:)
    integer::i,j,nx,ny
    real(dp)::d,fx,fy,en,p,nextv
    sx=x;sy=y;call sort_values(sx);call sort_values(sy)
    nx=size(x);ny=size(y);i=1;j=1;d=0.0_dp
    do while(i<=nx .or. j<=ny)
      if(i>nx) then
        nextv=sy(j)
      else if(j>ny) then
        nextv=sx(i)
      else
        nextv=min(sx(i),sy(j))
      end if
      do while(i<=nx)
        if(sx(i)>nextv) exit
        i=i+1
      end do
      do while(j<=ny)
        if(sy(j)>nextv) exit
        j=j+1
      end do
      fx=real(i-1,dp)/real(nx,dp)
      fy=real(j-1,dp)/real(ny,dp)
      d=max(d,abs(fx-fy))
    end do
    en=sqrt(real(nx*ny,dp)/real(nx+ny,dp))
    p=kolmogorov_survival((en+0.12_dp+0.11_dp/en)*d)
    r%method='Two-sample Kolmogorov-Smirnov test';r%statistic=d;r%p_value=p
  end function
  function cvm_normal_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp),allocatable::z(:);real(dp)::m,s,w2;integer::i,n
    n=size(x);m=sample_mean(x);s=sample_sd(x);allocate(z(n));z=(x-m)/s;call sort_values(z);w2=1.0_dp/(12*n)
    do i=1,n;w2=w2+(normal_cdf(z(i))-(2.0_dp*i-1)/(2.0_dp*n))**2;end do
    r%method='Cramer-von Mises normality test';r%statistic=w2;r%p_value=exp(-max(0.0_dp,1.2337141_dp/w2)*(1.0_dp+0.026_dp/w2))
    r%p_value=clamp(r%p_value,0.0_dp,1.0_dp)
  end function
  function anderson_darling_normal_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp),allocatable::z(:);real(dp)::m,s,a2,p,fi,fj;integer::i,n
    n=size(x);m=sample_mean(x);s=sample_sd(x);allocate(z(n));z=(x-m)/s;call sort_values(z);a2=-n
    do i=1,n;fi=clamp(normal_cdf(z(i)),1e-15_dp,1-1e-15_dp);fj=clamp(normal_cdf(z(n+1-i)),1e-15_dp,1-1e-15_dp);a2=a2-(2*i-1.0_dp)/n*(log(fi)+log(1-fj));end do
    a2=a2*(1+0.75_dp/n+2.25_dp/n**2)
    if(a2<0.2)then;p=1-exp(-13.436+101.14*a2-223.73*a2*a2);else if(a2<0.34)then;p=1-exp(-8.318+42.796*a2-59.938*a2*a2);else if(a2<0.6)then;p=exp(0.9177-4.279*a2-1.38*a2*a2);else;p=exp(1.2937-5.709*a2+0.0186*a2*a2);end if
    r%method='Anderson-Darling normality test';r%statistic=a2;r%p_value=clamp(p,0.0_dp,1.0_dp)
  end function
  function lilliefors_normal_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp),allocatable::z(:);real(dp)::m,s,d,fn,p;integer::i,n
    n=size(x);m=sample_mean(x);s=sample_sd(x);allocate(z(n));z=(x-m)/s;call sort_values(z);d=0
    do i=1,n;fn=normal_cdf(z(i));d=max(d,max(abs(fn-real(i,dp)/n),abs(fn-real(i-1,dp)/n)));end do
    p=exp(-7.01256_dp*d*d*(n+2.78019_dp)+2.99587_dp*d*sqrt(n+2.78019_dp)-0.122119_dp+0.974598_dp/sqrt(real(n,dp))+1.67997_dp/real(n,dp))
    r%method='Lilliefors normality test';r%statistic=d;r%p_value=clamp(p,0.0_dp,1.0_dp)
  end function
  function shapiro_francia_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp),allocatable::sx(:),m(:);real(dp)::xb,s2,w,y,mu,sigma,z;integer::i,n
    n=size(x);sx=x;call sort_values(sx);allocate(m(n));do i=1,n;m(i)=inverse_normal((i-0.375_dp)/(n+0.25_dp));end do;m=m/sqrt(sum(m*m));xb=sample_mean(x);s2=sum((x-xb)**2);w=(dot_product(m,sx))**2/s2;y=log(1-w);mu=-1.2725_dp+1.0521_dp*log(log(real(n,dp)));sigma=1.0308_dp-0.26758_dp*log(log(real(n,dp)));z=(y-mu)/sigma
    r%method='Shapiro-Francia normality test';r%statistic=w;r%p_value=1-normal_cdf(z)
  contains
    real(dp) function inverse_normal(p) result(v);use fbasics_special,only:normal_quantile;real(dp),intent(in)::p;v=normal_quantile(p);end function
  end function
  function ranks(x) result(r)
    real(dp),intent(in)::x(:);real(dp),allocatable::r(:);integer,allocatable::idx(:);integer::i,j,k,n;real(dp)::avg,key;integer::it
    n=size(x);allocate(r(n),idx(n));idx=[(i,i=1,n)]
    do i=2,n;it=idx(i);key=x(it);j=i-1;do while(j>=1);if(x(idx(j))<=key)exit;idx(j+1)=idx(j);j=j-1;end do;idx(j+1)=it;end do
    i=1
    do while(i<=n)
      j=i
      do while(j<n)
        if(x(idx(j+1))/=x(idx(i))) exit
        j=j+1
      end do
      avg=0.5_dp*real(i+j,dp)
      do k=i,j
        r(idx(k))=avg
      end do
      i=j+1
    end do
  end function
  subroutine sort_values(a)
    real(dp),intent(inout)::a(:);integer::i,j;real(dp)::key
    do i=2,size(a);key=a(i);j=i-1;do while(j>=1);if(a(j)<=key)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=key;end do
  end subroutine
  real(dp) function kolmogorov_survival(x) result(p)
    real(dp),intent(in)::x;integer::k;real(dp)::term
    p=0;do k=1,100;term=2*(-1.0_dp)**(k-1)*exp(-2.0_dp*k*k*x*x);p=p+term;if(abs(term)<1e-14_dp)exit;end do;p=clamp(p,0.0_dp,1.0_dp)
  end function
end module fbasics_hypothesis
