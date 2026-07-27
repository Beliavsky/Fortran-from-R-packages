! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use fbasics_kinds, only: dp
  use fbasics_special, only: student_quantile
  implicit none
  private
  type, public :: basic_stats_result
    integer :: nobs = 0
    integer :: missing = 0
    real(dp) :: minimum=0.0_dp, maximum=0.0_dp
    real(dp) :: q1=0.0_dp, q3=0.0_dp, mean=0.0_dp, median=0.0_dp
    real(dp) :: sum=0.0_dp, se_mean=0.0_dp, lcl_mean=0.0_dp, ucl_mean=0.0_dp
    real(dp) :: variance=0.0_dp, stdev=0.0_dp, skewness=0.0_dp, kurtosis=0.0_dp
  end type
  public :: sample_mean, sample_variance, sample_sd, sample_quantile
  public :: sample_skewness, sample_kurtosis, basic_stats
  public :: row_means, row_sds, row_vars, row_skewness, row_kurtosis
  public :: row_maxs, row_mins, row_prods, row_quantiles
  public :: sample_median, sample_iqr, sample_robust_skewness, sample_robust_kurtosis
  public :: sample_lmoments, covariance_matrix, correlation_matrix
contains
  pure real(dp) function sample_mean(x) result(v)
    real(dp),intent(in)::x(:); integer::i,n
    v=0.0_dp;n=0
    do i=1,size(x); if(.not.ieee_is_nan(x(i)))then;v=v+x(i);n=n+1;end if;end do
    if(n>0)v=v/real(n,dp)
  end function
  pure real(dp) function sample_variance(x) result(v)
    real(dp),intent(in)::x(:); real(dp)::m,d; integer::i,n
    m=sample_mean(x); v=0.0_dp;n=0
    do i=1,size(x);if(.not.ieee_is_nan(x(i)))then;d=x(i)-m;v=v+d*d;n=n+1;end if;end do
    if(n>1)then;v=v/real(n-1,dp);else;v=0.0_dp;end if
  end function
  pure real(dp) function sample_sd(x) result(v)
    real(dp),intent(in)::x(:);v=sqrt(max(sample_variance(x),0.0_dp))
  end function
  real(dp) function sample_quantile(x,p) result(q)
    real(dp),intent(in)::x(:),p; real(dp),allocatable::y(:); real(dp)::h
    integer::i,n,j
    n=count(.not.ieee_is_nan(x)); if(n==0)then;q=0.0_dp;return;end if
    allocate(y(n));j=0;do i=1,size(x);if(.not.ieee_is_nan(x(i)))then;j=j+1;y(j)=x(i);end if;end do
    call sort_real(y)
    if(n==1)then;q=y(1);return;end if
    h=1.0_dp+(real(n-1,dp))*min(max(p,0.0_dp),1.0_dp);i=int(floor(h));j=min(i+1,n)
    q=y(i)+(h-real(i,dp))*(y(j)-y(i))
  end function
  real(dp) function sample_median(x) result(v)
    real(dp),intent(in)::x(:);v=sample_quantile(x,0.5_dp)
  end function
  real(dp) function sample_iqr(x) result(v)
    real(dp),intent(in)::x(:);v=sample_quantile(x,0.75_dp)-sample_quantile(x,0.25_dp)
  end function
  pure real(dp) function sample_skewness(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::m,s,d,s3;integer::i,n
    m=sample_mean(x);s=sample_sd(x);s3=0.0_dp;n=0
    if(s<=0.0_dp)then;v=0.0_dp;return;end if
    do i=1,size(x);if(.not.ieee_is_nan(x(i)))then;d=(x(i)-m)/s;s3=s3+d**3;n=n+1;end if;end do
    if(n>2)then;v=real(n,dp)*s3/real((n-1)*(n-2),dp);else;v=0.0_dp;end if
  end function
  pure real(dp) function sample_kurtosis(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::m,s,d,s4;integer::i,n
    m=sample_mean(x);s=sample_sd(x);s4=0.0_dp;n=0
    if(s<=0.0_dp)then;v=0.0_dp;return;end if
    do i=1,size(x);if(.not.ieee_is_nan(x(i)))then;d=(x(i)-m)/s;s4=s4+d**4;n=n+1;end if;end do
    if(n>3)then
      v=(real(n*(n+1),dp)*s4/real(n-1,dp)-3.0_dp*real((n-1)*(n-1),dp))/real((n-2)*(n-3),dp)
    else;v=0.0_dp;end if
  end function
  function basic_stats(x,ci) result(r)
    real(dp),intent(in)::x(:);real(dp),intent(in),optional::ci;type(basic_stats_result)::r
    real(dp)::c,tcrit;integer::n
    c=0.95_dp;if(present(ci))c=ci;n=count(.not.ieee_is_nan(x))
    r%nobs=size(x);r%missing=size(x)-n
    if(n==0)return
    r%minimum=minval(x,mask=.not.ieee_is_nan(x));r%maximum=maxval(x,mask=.not.ieee_is_nan(x))
    r%q1=sample_quantile(x,0.25_dp);r%q3=sample_quantile(x,0.75_dp);r%mean=sample_mean(x)
    r%median=sample_median(x);r%sum=r%mean*real(n,dp);r%variance=sample_variance(x);r%stdev=sqrt(r%variance)
    r%se_mean=r%stdev/sqrt(real(n,dp));r%skewness=sample_skewness(x);r%kurtosis=sample_kurtosis(x)
    if(n>1)then;tcrit=student_quantile(0.5_dp+0.5_dp*c,real(n-1,dp));r%lcl_mean=r%mean-tcrit*r%se_mean;r%ucl_mean=r%mean+tcrit*r%se_mean;end if
  end function
  function row_means(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_mean(x(i,:));end do
  end function
  function row_sds(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_sd(x(i,:));end do
  end function
  function row_vars(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_variance(x(i,:));end do
  end function
  function row_skewness(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_skewness(x(i,:));end do
  end function
  function row_kurtosis(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_kurtosis(x(i,:));end do
  end function
  function row_maxs(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=maxval(x(i,:),mask=.not.ieee_is_nan(x(i,:)));end do
  end function
  function row_mins(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=minval(x(i,:),mask=.not.ieee_is_nan(x(i,:)));end do
  end function
  function row_prods(x) result(v)
    real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i,j
    v=1.0_dp;do i=1,size(x,1);do j=1,size(x,2);if(.not.ieee_is_nan(x(i,j)))v(i)=v(i)*x(i,j);end do;end do
  end function
  function row_quantiles(x,p) result(v)
    real(dp),intent(in)::x(:,:),p;real(dp)::v(size(x,1));integer::i
    do i=1,size(x,1);v(i)=sample_quantile(x(i,:),p);end do
  end function
  real(dp) function sample_robust_skewness(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::q1,q2,q3
    q1=sample_quantile(x,0.25_dp);q2=sample_quantile(x,0.5_dp);q3=sample_quantile(x,0.75_dp)
    if(q3>q1)then;v=(q3+q1-2.0_dp*q2)/(q3-q1);else;v=0.0_dp;end if
  end function
  real(dp) function sample_robust_kurtosis(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::q025,q25,q75,q975
    q025=sample_quantile(x,0.025_dp);q25=sample_quantile(x,0.25_dp);q75=sample_quantile(x,0.75_dp);q975=sample_quantile(x,0.975_dp)
    if(q75>q25)then;v=(q975-q025)/(q75-q25);else;v=0.0_dp;end if
  end function
  subroutine sample_lmoments(x,l1,l2,t3,t4)
    real(dp),intent(in)::x(:);real(dp),intent(out)::l1,l2,t3,t4
    real(dp),allocatable::y(:);real(dp)::b0,b1,b2,b3;integer::n,i
    n=size(x);allocate(y(n));y=x;call sort_real(y)
    b0=sum(y)/real(n,dp);b1=0.0_dp;b2=0.0_dp;b3=0.0_dp
    if(n>1) then
      do i=1,n
        b1=b1+real(i-1,dp)/real(n-1,dp)*y(i)
      end do
    end if
    if(n>2) then
      do i=1,n
        b2=b2+real((i-1)*(i-2),dp)/real((n-1)*(n-2),dp)*y(i)
      end do
    end if
    if(n>3) then
      do i=1,n
        b3=b3+real((i-1)*(i-2)*(i-3),dp)/real((n-1)*(n-2)*(n-3),dp)*y(i)
      end do
    end if
    b1=b1/real(n,dp);b2=b2/real(n,dp);b3=b3/real(n,dp)
    l1=b0;l2=2.0_dp*b1-b0
    if(abs(l2)>epsilon(1.0_dp))then;t3=(6.0_dp*b2-6.0_dp*b1+b0)/l2;t4=(20.0_dp*b3-30.0_dp*b2+12.0_dp*b1-b0)/l2;else;t3=0.0_dp;t4=0.0_dp;end if
  end subroutine
  function covariance_matrix(x) result(c)
    real(dp),intent(in)::x(:,:);real(dp)::c(size(x,2),size(x,2));real(dp)::mi,mj
    integer::i,j,k,n
    c=0.0_dp
    do j=1,size(x,2);do i=1,j;mi=sample_mean(x(:,i));mj=sample_mean(x(:,j));n=0
      do k=1,size(x,1);if(.not.ieee_is_nan(x(k,i)).and..not.ieee_is_nan(x(k,j)))then;c(i,j)=c(i,j)+(x(k,i)-mi)*(x(k,j)-mj);n=n+1;end if;end do
      if(n>1)c(i,j)=c(i,j)/real(n-1,dp);c(j,i)=c(i,j)
    end do;end do
  end function
  function correlation_matrix(x) result(r)
    real(dp),intent(in)::x(:,:);real(dp)::r(size(x,2),size(x,2)),c(size(x,2),size(x,2));integer::i,j
    c=covariance_matrix(x);do j=1,size(x,2);do i=1,size(x,2);if(c(i,i)>0.and.c(j,j)>0)then;r(i,j)=c(i,j)/sqrt(c(i,i)*c(j,j));else;r(i,j)=0;end if;end do;end do
  end function
  subroutine sort_real(a)
    real(dp),intent(inout)::a(:);integer::i,j;real(dp)::key
    do i=2,size(a);key=a(i);j=i-1;do while(j>=1);if(a(j)<=key)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=key;end do
  end subroutine
end module fbasics_stats
