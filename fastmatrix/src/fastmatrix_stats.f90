module fastmatrix_stats
  use fastmatrix_base, only: dp, pi, gamma_p, normal_cdf, normal_rand, chol_lower, inverse_matrix
  implicit none
  private
  type, public :: moments_result
    real(dp) :: mean=0.0_dp, second=0.0_dp, third=0.0_dp, fourth=0.0_dp, skewness=0.0_dp, kurtosis=0.0_dp
  end type
  public :: moments, geomean, skewness, kurtosis, mahalanobis, cov_weighted, cov_mssd, mediancenter, ccc
  public :: dchi, pchi, qchi, rchi, jarque_bera, wilson_hilferty
contains
  pure function moments(x) result(r)
    real(dp),intent(in)::x(:)
    type(moments_result)::r
    real(dp)::z(size(x))
    integer::n
    n=size(x)
    r%mean=sum(x)/real(n,dp)
    z=x-r%mean
    r%second=sum(z**2)/n
    r%third=sum(z**3)/n
    r%fourth=sum(z**4)/n
    if(r%second>0)then
    r%skewness=r%third/r%second**1.5_dp
    r%kurtosis=r%fourth/r%second**2
    end if
  end function
  pure function skewness(x) result(v)
  real(dp),intent(in)::x(:)
  real(dp)::v
  type(moments_result)::r
  r=moments(x)
  v=r%skewness
  end function
  pure function kurtosis(x) result(v)
  real(dp),intent(in)::x(:)
  real(dp)::v
  type(moments_result)::r
  r=moments(x)
  v=r%kurtosis
  end function
  pure function geomean(x) result(v)
  real(dp),intent(in)::x(:)
  real(dp)::v
  if(any(x<=0))then
  v=0
  else
  v=exp(sum(log(x))/size(x))
  end if
  end function

  function mahalanobis(x,center,cov) result(d2)
    real(dp),intent(in)::x(:,:),center(:),cov(:,:)
    real(dp)::d2(size(x,1)),ci(size(cov,1),size(cov,2)),z(size(center))
    integer::i,info
    call inverse_matrix(cov,ci,info)
    do i=1,size(x,1)
    z=x(i,:)-center
    d2(i)=dot_product(z,matmul(ci,z))
    end do
  end function

  subroutine cov_weighted(x,w,center,cov)
    real(dp),intent(in)::x(:,:),w(:)
    real(dp),intent(out)::center(:),cov(:,:)
    real(dp)::sw,z(size(x,2))
    integer::i,j
    sw=sum(w)
    center=matmul(transpose(x),w)/sw
    cov=0
    do i=1,size(x,1)
    z=x(i,:)-center
    do j=1,size(z)
    cov(j,:)=cov(j,:)+w(i)*z(j)*z
    end do
    end do
    cov=cov/sw
  end subroutine

  subroutine cov_mssd(x,cov)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::cov(:,:)
    real(dp)::z(size(x,2))
    integer::i,j,n
    n=size(x,1)
    cov=0
    do i=1,n-1
    z=x(i+1,:)-x(i,:)
    do j=1,size(z)
    cov(j,:)=cov(j,:)+z(j)*z
    end do
    end do
    cov=cov/(2.0_dp*real(n-1,dp))
  end subroutine

  pure function median1(x) result(m)
    real(dp),intent(in)::x(:)
    real(dp)::m,y(size(x)),t
    integer::i,j,n
    y=x
    n=size(y)
    do i=2,n
    t=y(i)
    j=i-1
    do while(j>=1.and.y(j)>t)
    y(j+1)=y(j)
    j=j-1
    end do
    y(j+1)=t
    end do
    if(mod(n,2)==1)then
    m=y((n+1)/2)
    else
    m=0.5_dp*(y(n/2)+y(n/2+1))
    end if
  end function
  subroutine mediancenter(x,centered,center)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::centered(:,:),center(:)
    integer::j
    do j=1,size(x,2)
    center(j)=median1(x(:,j))
    centered(:,j)=x(:,j)-center(j)
    end do
  end subroutine

  pure function ccc(x,y) result(v)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::v,mx,my,vx,vy,cxy
    mx=sum(x)/size(x)
    my=sum(y)/size(y)
    vx=sum((x-mx)**2)/size(x)
    vy=sum((y-my)**2)/size(y)
    cxy=sum((x-mx)*(y-my))/size(x)
    v=2*cxy/(vx+vy+(mx-my)**2)
  end function

  pure function dchi(x,df,log_density) result(v)
    real(dp),intent(in)::x,df
    logical,intent(in),optional::log_density
    real(dp)::v,lv
    logical::lg
    lg=.false.
    if(present(log_density))lg=log_density
    if(x<0.or.df<=0)then
    lv=-huge(1.0_dp)
    else
    lv=(1-df/2)*log(2.0_dp)-log_gamma(df/2)+(df-1)*log(max(x,tiny(1.0_dp)))-x*x/2
    end if
    if(lg)then
    v=lv
    else
    v=exp(lv)
    end if
  end function
  pure function pchi(x,df) result(v)
  real(dp),intent(in)::x,df
  real(dp)::v
  if(x<=0)then
  v=0
  else
  v=gamma_p(df/2,x*x/2)
  end if
  end function
  function qchi(p,df) result(x)
    real(dp),intent(in)::p,df
    real(dp)::x,lo,hi,mid
    integer::i
    if(p<=0)then
    x=0
    return
    else if(p>=1)then
    x=huge(1.0_dp)
    return
    end if
    lo=0
    hi=max(1.0_dp,sqrt(df))
    do while(pchi(hi,df)<p)
    hi=2*hi
    end do
    do i=1,100
    mid=.5_dp*(lo+hi)
    if(pchi(mid,df)<p)then
    lo=mid
    else
    hi=mid
    end if
    end do
    x=.5_dp*(lo+hi)
  end function
  function rchi(df) result(x)
  real(dp),intent(in)::df
  real(dp)::x
  integer::k
  real(dp)::s,z
    ! exact for integer df, Wilson-Hilferty fallback otherwise
    if(abs(df-nint(df))<1e-12_dp.and.df<=1000)then
    s=0
    do k=1,nint(df)
    z=normal_rand()
    s=s+z*z
    end do
    x=sqrt(s)
    else
    call random_number(s)
    x=qchi(s,df)
    end if
  end function

  subroutine jarque_bera(x,stat,pvalue)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::stat,pvalue
    type(moments_result)::r
    real(dp)::n
    r=moments(x)
    n=real(size(x),dp)
    stat=n/6.0_dp*(r%skewness**2+(r%kurtosis-3.0_dp)**2/4.0_dp)
    pvalue=exp(-stat/2.0_dp)
  end subroutine
  pure function wilson_hilferty(x,df) result(z)
  real(dp),intent(in)::x,df
  real(dp)::z
  z=((x/df)**(1.0_dp/3.0_dp)-(1-2/(9*df)))/sqrt(2/(9*df))
  end function
end module
