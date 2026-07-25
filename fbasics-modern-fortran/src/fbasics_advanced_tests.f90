! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_advanced_tests
  use fbasics_kinds, only: dp, clamp
  use fbasics_stats, only: sample_mean, sample_variance, sample_median
  use fbasics_special, only: normal_cdf, normal_quantile, chi_square_cdf
  use fbasics_hypothesis, only: test_result
  implicit none
  private
  public :: shapiro_wilk_test, ansari_bradley_test, mood_scale_test
  public :: bartlett_two_sample_test, fligner_killeen_test
  public :: wilcoxon_rank_sum_test, kruskal_wallis_two_sample_test
  public :: ks_normal_test, pearson_chi_square_normal_test
  public :: adjusted_jarque_bera_test
contains
  function shapiro_wilk_test(x) result(r)
    real(dp),intent(in)::x(:)
    type(test_result)::r
    real(dp),allocatable::sx(:),m(:),a(:)
    real(dp)::ssm,ssx,u,a1,a2,phi,w,y,mu,sigma,z,lnn,gamma
    integer::n,i
    n=size(x);r%method='Shapiro-Wilk normality test'
    if(n<3)then;r%statistic=0.0_dp;r%p_value=0.0_dp;return;end if
    sx=x;call sort_values(sx);allocate(m(n),a(n));
    do i=1,n;m(i)=normal_quantile((real(i,dp)-0.375_dp)/(real(n,dp)+0.25_dp));end do
    ssm=sum(m*m);u=1.0_dp/sqrt(real(n,dp))
    a1=-2.706056_dp*u**5+4.434685_dp*u**4-2.071190_dp*u**3-0.147981_dp*u**2+0.221157_dp*u+m(n)/sqrt(ssm)
    if(n>5)then
      a2=-3.582633_dp*u**5+5.682633_dp*u**4-1.752461_dp*u**3-0.293762_dp*u**2+0.042981_dp*u+m(n-1)/sqrt(ssm)
      phi=(ssm-2.0_dp*m(n)**2-2.0_dp*m(n-1)**2)/(1.0_dp-2.0_dp*a1*a1-2.0_dp*a2*a2)
      a=m/sqrt(phi);a(1)=-a1;a(n)=a1;a(2)=-a2;a(n-1)=a2
    else
      phi=(ssm-2.0_dp*m(n)**2)/(1.0_dp-2.0_dp*a1*a1)
      a=m/sqrt(phi);a(1)=-a1;a(n)=a1
    end if
    ssx=sum((x-sample_mean(x))**2)
    w=min(1.0_dp,max(0.0_dp,dot_product(a,sx)**2/max(ssx,tiny(1.0_dp))))
    if(n<=11)then
      gamma=-2.273_dp+0.459_dp*real(n,dp)
      if(gamma<=-log(max(1.0e-16_dp,1.0_dp-w)))then;r%p_value=1.0e-16_dp;r%statistic=w;return;end if
      y=-log(gamma-log(max(1.0e-16_dp,1.0_dp-w)))
      mu=0.5440_dp-0.39978_dp*n+0.025054_dp*n*n-0.0006714_dp*n*n*n
      sigma=exp(1.3822_dp-0.77857_dp*n+0.062767_dp*n*n-0.0020322_dp*n*n*n)
    else
      lnn=log(real(n,dp));y=log(max(1.0e-16_dp,1.0_dp-w))
      mu=-1.5861_dp-0.31082_dp*lnn-0.083751_dp*lnn**2+0.0038915_dp*lnn**3
      sigma=exp(-0.4803_dp-0.082676_dp*lnn+0.0030302_dp*lnn**2)
    end if
    z=(y-mu)/sigma;r%statistic=w;r%p_value=clamp(1.0_dp-normal_cdf(z),0.0_dp,1.0_dp)
  end function shapiro_wilk_test

  function ansari_bradley_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r
    real(dp),allocatable::all(:),score(:);integer,allocatable::group(:)
    real(dp)::stat,mean0,var0,z,ms,vs;integer::nx,ny,n,i
    nx=size(x);ny=size(y);n=nx+ny;allocate(all(n),group(n),score(n));all(1:nx)=x;all(nx+1:n)=y;group(1:nx)=1;group(nx+1:n)=2
    call sort_with_group(all,group);do i=1,n;score(i)=real(min(i,n+1-i),dp);end do;call average_tied_scores(all,score)
    stat=0.0_dp;do i=1,n;if(group(i)==1)stat=stat+score(i);end do
    ms=sum(score)/n;vs=sum((score-ms)**2)/n;mean0=nx*ms;var0=real(nx*ny,dp)*vs/real(max(1,n-1),dp);z=(stat-mean0)/sqrt(max(var0,tiny(1.0_dp)))
    r%method='Ansari-Bradley scale test';r%statistic=stat;r%estimate=z;r%p_value=2.0_dp*(1.0_dp-normal_cdf(abs(z)))
  end function ansari_bradley_test

  function mood_scale_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r
    real(dp),allocatable::all(:),rk(:);integer,allocatable::group(:)
    real(dp)::stat,mean0,var0,z;integer::nx,ny,n,i
    nx=size(x);ny=size(y);n=nx+ny;allocate(all(n),group(n),rk(n));all(1:nx)=x;all(nx+1:n)=y;group(1:nx)=1;group(nx+1:n)=2
    call sort_with_group(all,group);do i=1,n;rk(i)=real(i,dp);end do;call average_tied_scores(all,rk)
    stat=0.0_dp;do i=1,n;if(group(i)==1)stat=stat+(rk(i)-0.5_dp*(n+1))**2;end do
    mean0=real(nx,dp)*(real(n*n,dp)-1.0_dp)/12.0_dp
    var0=real(nx*ny,dp)*(real(n+1,dp))*(real(n*n,dp)-4.0_dp)/180.0_dp
    z=(stat-mean0)/sqrt(max(var0,tiny(1.0_dp)));r%method='Mood two-sample scale test';r%statistic=stat;r%estimate=z;r%p_value=2.0_dp*(1.0_dp-normal_cdf(abs(z)))
  end function mood_scale_test

  function bartlett_two_sample_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r
    real(dp)::v1,v2,sp,num,c,stat;integer::n1,n2,n
    n1=size(x);n2=size(y);n=n1+n2;v1=sample_variance(x);v2=sample_variance(y);sp=((n1-1)*v1+(n2-1)*v2)/real(n-2,dp)
    num=(n-2)*log(sp)-(n1-1)*log(v1)-(n2-1)*log(v2);c=1.0_dp+(1.0_dp/(n1-1)+1.0_dp/(n2-1)-1.0_dp/(n-2))/3.0_dp;stat=num/c
    r%method='Bartlett two-sample variance test';r%statistic=stat;r%df=1;r%p_value=1.0_dp-chi_square_cdf(stat,1.0_dp)
  end function bartlett_two_sample_test

  function fligner_killeen_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r
    real(dp),allocatable::dev(:),rk(:),score(:);integer,allocatable::group(:)
    real(dp)::mx,my,overall,between,within,stat;integer::nx,ny,n,i
    nx=size(x);ny=size(y);n=nx+ny;allocate(dev(n),rk(n),score(n),group(n));dev(1:nx)=abs(x-sample_median(x));dev(nx+1:n)=abs(y-sample_median(y));group(1:nx)=1;group(nx+1:n)=2
    call sort_with_group(dev,group);do i=1,n;rk(i)=i;end do;call average_tied_scores(dev,rk);do i=1,n;score(i)=normal_quantile(0.5_dp*(1.0_dp+(rk(i)-0.5_dp)/real(n,dp)));end do
    overall=sum(score)/n;mx=0.0_dp;my=0.0_dp;do i=1,n;if(group(i)==1)then;mx=mx+score(i);else;my=my+score(i);end if;end do;mx=mx/nx;my=my/ny
    between=nx*(mx-overall)**2+ny*(my-overall)**2;within=sum((score-overall)**2)/real(n-1,dp);stat=between/max(within,tiny(1.0_dp))
    r%method='Fligner-Killeen scale test';r%statistic=stat;r%df=1;r%p_value=1.0_dp-chi_square_cdf(stat,1.0_dp)
  end function fligner_killeen_test

  function wilcoxon_rank_sum_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r
    real(dp),allocatable::all(:),rk(:);integer,allocatable::group(:)
    real(dp)::w,mean0,var0,z,tiecorr;integer::nx,ny,n,i,j
    nx=size(x);ny=size(y);n=nx+ny;allocate(all(n),rk(n),group(n));all(1:nx)=x;all(nx+1:n)=y;group(1:nx)=1;group(nx+1:n)=2
    call sort_with_group(all,group);do i=1,n;rk(i)=i;end do;call average_tied_scores(all,rk);w=0.0_dp;do i=1,n;if(group(i)==1)w=w+rk(i);end do
    tiecorr=0.0_dp;i=1;do while(i<=n);j=i;do while(j<n);if(all(j+1)/=all(i))exit;j=j+1;end do;tiecorr=tiecorr+real((j-i+1)**3-(j-i+1),dp);i=j+1;end do
    mean0=0.5_dp*nx*(n+1);var0=real(nx*ny,dp)/12.0_dp*((n+1.0_dp)-tiecorr/real(n*(n-1),dp));z=(w-mean0)/sqrt(max(var0,tiny(1.0_dp)))
    r%method='Wilcoxon rank-sum test';r%statistic=w;r%estimate=z;r%p_value=2.0_dp*(1.0_dp-normal_cdf(abs(z)))
  end function wilcoxon_rank_sum_test

  function kruskal_wallis_two_sample_test(x,y) result(r)
    real(dp),intent(in)::x(:),y(:);type(test_result)::r;type(test_result)::w
    real(dp)::h;integer::nx,ny,n
    w=wilcoxon_rank_sum_test(x,y);nx=size(x);ny=size(y);n=nx+ny;h=12.0_dp/real(n*(n+1),dp)*(w%statistic**2/nx+(real(n*(n+1)/2,dp)-w%statistic)**2/ny)-3.0_dp*(n+1)
    r%method='Kruskal-Wallis two-sample test';r%statistic=h;r%df=1;r%p_value=1.0_dp-chi_square_cdf(h,1.0_dp)
  end function kruskal_wallis_two_sample_test

  function adjusted_jarque_bera_test(x) result(r)
    real(dp),intent(in)::x(:)
    type(test_result)::r
    integer::n
    real(dp)::m,m2,m3,m4,b1,b2,eb2,varb1,varb2,stat
    n=size(x);r%method='Urzua adjusted Jarque-Bera normality test';r%df=2
    if(n<8)then;r%statistic=0.0_dp;r%p_value=0.0_dp;return;end if
    m=sample_mean(x);m2=sum((x-m)**2)/real(n,dp)
    if(m2<=tiny(1.0_dp))then;r%statistic=0.0_dp;r%p_value=0.0_dp;return;end if
    m3=sum((x-m)**3)/real(n,dp);m4=sum((x-m)**4)/real(n,dp)
    b1=(m3/m2**1.5_dp)**2;b2=m4/(m2*m2)
    varb1=6.0_dp*real(n-2,dp)/real((n+1)*(n+3),dp)
    eb2=3.0_dp*real(n-1,dp)/real(n+1,dp)
    varb2=24.0_dp*real(n*(n-2)*(n-3),dp)/real((n+1)*(n+1)*(n+3)*(n+5),dp)
    stat=b1/max(varb1,tiny(1.0_dp))+(b2-eb2)**2/max(varb2,tiny(1.0_dp))
    r%statistic=stat;r%p_value=1.0_dp-chi_square_cdf(stat,2.0_dp)
  end function adjusted_jarque_bera_test

  function ks_normal_test(x) result(r)
    real(dp),intent(in)::x(:);type(test_result)::r;real(dp),allocatable::sx(:);real(dp)::d,p,f;integer::i,n
    n=size(x);sx=x;call sort_values(sx);d=0.0_dp;do i=1,n;f=normal_cdf(sx(i));d=max(d,max(abs(f-real(i,dp)/n),abs(f-real(i-1,dp)/n)));end do;p=kolmogorov_survival((sqrt(real(n,dp))+0.12_dp+0.11_dp/sqrt(real(n,dp)))*d)
    r%method='One-sample normal Kolmogorov-Smirnov test';r%statistic=d;r%p_value=p
  end function ks_normal_test

  function pearson_chi_square_normal_test(x,bins) result(r)
    real(dp),intent(in)::x(:);integer,intent(in),optional::bins;type(test_result)::r
    integer::k,i,j,n,count;real(dp)::m,s,lo,hi,e,stat
    n=size(x);k=max(4,nint(2.0_dp*real(n,dp)**0.4_dp));if(present(bins))k=max(3,bins);m=sample_mean(x);s=sqrt(sample_variance(x));stat=0.0_dp
    do j=1,k
      if(j==1)then;lo=-huge(1.0_dp);else;lo=m+s*normal_quantile(real(j-1,dp)/k);end if
      if(j==k)then;hi=huge(1.0_dp);else;hi=m+s*normal_quantile(real(j,dp)/k);end if
      count=0;do i=1,n;if(x(i)>lo.and.x(i)<=hi)count=count+1;end do;e=real(n,dp)/k;stat=stat+(count-e)**2/e
    end do
    r%method='Pearson chi-square normality test';r%statistic=stat;r%df=max(1,k-3);r%p_value=1.0_dp-chi_square_cdf(stat,real(r%df,dp))
  end function pearson_chi_square_normal_test

  subroutine sort_with_group(a,g)
    real(dp),intent(inout)::a(:);integer,intent(inout)::g(:);integer::i,j,kg;real(dp)::key
    do i=2,size(a);key=a(i);kg=g(i);j=i-1;do while(j>=1);if(a(j)<=key)exit;a(j+1)=a(j);g(j+1)=g(j);j=j-1;end do;a(j+1)=key;g(j+1)=kg;end do
  end subroutine
  subroutine sort_values(a)
    real(dp),intent(inout)::a(:);integer::i,j;real(dp)::key
    do i=2,size(a);key=a(i);j=i-1;do while(j>=1);if(a(j)<=key)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=key;end do
  end subroutine
  subroutine average_tied_scores(a,s)
    real(dp),intent(in)::a(:);real(dp),intent(inout)::s(:);integer::i,j;i=1;do while(i<=size(a));j=i;do while(j<size(a));if(a(j+1)/=a(i))exit;j=j+1;end do;s(i:j)=sum(s(i:j))/real(j-i+1,dp);i=j+1;end do
  end subroutine
  real(dp) function kolmogorov_survival(x) result(p)
    real(dp),intent(in)::x;integer::k;real(dp)::term;p=0.0_dp;do k=1,100;term=2.0_dp*(-1.0_dp)**(k-1)*exp(-2.0_dp*k*k*x*x);p=p+term;if(abs(term)<1.0e-14_dp)exit;end do;p=clamp(p,0.0_dp,1.0_dp)
  end function
end module fbasics_advanced_tests
