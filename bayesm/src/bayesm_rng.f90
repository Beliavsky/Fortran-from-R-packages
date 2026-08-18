module bayesm_rng
  use bayesm_kinds, only: dp, pi
  use bayesm_math, only: normal_cdf, normal_quantile
  use bayesm_linalg, only: chol_upper, inverse_upper
  implicit none
  private
  public :: rng_seed, randn, rand_gamma, rand_chisq, rand_beta, rand_dirichlet
  public :: rand_categorical, rand_poisson, rand_binomial, rand_mvn, rand_mvst
  public :: rand_truncnorm, rwishart_draw, rand_inverse_gamma, rand_uniform
contains
  subroutine rng_seed(seed)
    integer, intent(in) :: seed
    integer :: n,i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=mod(abs(seed)+104729*i+7919*i*i,huge(1)-1)+1
    end do
    call random_seed(put=put)
  end subroutine rng_seed

  real(dp) function rand_uniform() result(u)
    call random_number(u)
    if (u<=0.0_dp) u=tiny(1.0_dp)
    if (u>=1.0_dp) u=1.0_dp-epsilon(1.0_dp)
  end function rand_uniform

  real(dp) function randn() result(z)
    real(dp) :: u1,u2
    u1=rand_uniform(); u2=rand_uniform()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function randn

  recursive real(dp) function rand_gamma(shape,scale) result(x)
    real(dp), intent(in) :: shape,scale
    real(dp) :: d,c,z,u
    if (shape<=0.0_dp .or. scale<=0.0_dp) then
      x=0.0_dp
      return
    end if
    if (shape<1.0_dp) then
      u=rand_uniform()
      x=rand_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp
    c=1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z=randn()
        if (1.0_dp+c*z>0.0_dp) exit
      end do
      x=(1.0_dp+c*z)**3
      u=rand_uniform()
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-x+log(x))) exit
    end do
    x=d*x*scale
  end function rand_gamma

  real(dp) function rand_chisq(nu) result(x)
    real(dp), intent(in) :: nu
    x=rand_gamma(0.5_dp*nu,2.0_dp)
  end function rand_chisq

  real(dp) function rand_inverse_gamma(shape,scale) result(x)
    real(dp), intent(in) :: shape,scale
    x=1.0_dp/rand_gamma(shape,1.0_dp/scale)
  end function rand_inverse_gamma

  real(dp) function rand_beta(a,b) result(x)
    real(dp), intent(in) :: a,b
    real(dp) :: g1,g2
    g1=rand_gamma(a,1.0_dp); g2=rand_gamma(b,1.0_dp)
    x=g1/(g1+g2)
  end function rand_beta

  subroutine rand_dirichlet(alpha,x)
    real(dp), intent(in) :: alpha(:)
    real(dp), intent(out) :: x(size(alpha))
    integer :: i
    real(dp) :: s
    do i=1,size(alpha)
      x(i)=rand_gamma(alpha(i),1.0_dp)
    end do
    s=sum(x)
    if (s>0.0_dp) then
      x=x/s
    else
      x=1.0_dp/real(size(alpha),dp)
    end if
  end subroutine rand_dirichlet

  integer function rand_categorical(prob) result(k)
    real(dp), intent(in) :: prob(:)
    real(dp) :: u,s,total
    integer :: i
    total=sum(max(prob,0.0_dp))
    if (total<=0.0_dp) then
      k=1
      return
    end if
    u=rand_uniform()*total
    s=0.0_dp
    do i=1,size(prob)
      s=s+max(0.0_dp,prob(i))
      if (u<=s) then
        k=i
        return
      end if
    end do
    k=size(prob)
  end function rand_categorical

  integer function rand_binomial(n,p) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: i
    real(dp) :: pp
    pp=min(1.0_dp,max(0.0_dp,p))
    x=0
    if (n<=0) return
    if (n<100 .or. n*min(pp,1.0_dp-pp)<20.0_dp) then
      do i=1,n
        if (rand_uniform()<pp) x=x+1
      end do
    else
      x=nint(real(n,dp)*pp+sqrt(real(n,dp)*pp*(1.0_dp-pp))*randn())
      x=min(n,max(0,x))
    end if
  end function rand_binomial

  integer function rand_poisson(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: l,p,u,z
    integer :: i
    if (lambda<=0.0_dp) then
      k=0
    else if (lambda<30.0_dp) then
      l=exp(-lambda); p=1.0_dp; i=0
      do
        i=i+1; u=rand_uniform(); p=p*u
        if (p<=l) exit
      end do
      k=i-1
    else
      do
        z=lambda+sqrt(lambda)*randn()
        if (z>=-0.5_dp) exit
      end do
      k=max(0,nint(z))
    end if
  end function rand_poisson

  subroutine rand_mvn(mu,cov,x,info)
    real(dp), intent(in) :: mu(:),cov(:,:)
    real(dp), intent(out) :: x(size(mu))
    integer, intent(out) :: info
    real(dp) :: r(size(mu),size(mu)),z(size(mu))
    integer :: i
    call chol_upper(cov,r,info)
    if (info/=0) then
      x=mu
      return
    end if
    do i=1,size(mu); z(i)=randn(); end do
    x=mu+matmul(transpose(r),z)
  end subroutine rand_mvn

  subroutine rand_mvst(nu,mu,root,x)
    real(dp), intent(in) :: nu,mu(:),root(:,:)
    real(dp), intent(out) :: x(size(mu))
    real(dp) :: z(size(mu)),s
    integer :: i
    do i=1,size(mu); z(i)=randn(); end do
    s=sqrt(rand_chisq(nu)/nu)
    x=mu+matmul(transpose(root),z)/s
  end subroutine rand_mvst

  real(dp) function rand_truncnorm(mu,sigma,a,b) result(x)
    real(dp), intent(in) :: mu,sigma,a,b
    real(dp) :: lo,hi,u
    if (sigma<=0.0_dp) then
      x=mu
      return
    end if
    if (a <= -0.5_dp*huge(1.0_dp)) then
      lo=0.0_dp
    else
      lo=normal_cdf((a-mu)/sigma)
    end if
    if (b >= 0.5_dp*huge(1.0_dp)) then
      hi=1.0_dp
    else
      hi=normal_cdf((b-mu)/sigma)
    end if
    if (hi<=lo+10.0_dp*epsilon(1.0_dp)) then
      ! Rejection fallback is more reliable in extremely narrow tails.
      do
        x=mu+sigma*randn()
        if (x>=a .and. x<=b) return
      end do
    end if
    u=lo+(hi-lo)*rand_uniform()
    x=mu+sigma*normal_quantile(u)
    x=min(b,max(a,x))
  end function rand_truncnorm

  subroutine rwishart_draw(nu,v,w,iw,c,ci,info)
    real(dp), intent(in) :: nu,v(:,:)
    real(dp), intent(out) :: w(size(v,1),size(v,2)),iw(size(v,1),size(v,2))
    real(dp), intent(out) :: c(size(v,1),size(v,2)),ci(size(v,1),size(v,2))
    integer, intent(out) :: info
    integer :: m,i,j
    real(dp) :: t(size(v,1),size(v,2)),rv(size(v,1),size(v,2))
    m=size(v,1)
    t=0.0_dp
    call chol_upper(v,rv,info)
    if (info/=0) return
    do i=1,m
      if (nu-real(i-1,dp)<=0.0_dp) then
        info=-2
        return
      end if
      t(i,i)=sqrt(rand_chisq(nu-real(i-1,dp)))
      do j=1,i-1
        t(i,j)=randn()
      end do
    end do
    c=matmul(transpose(t),rv)
    call inverse_upper(c,ci,info)
    if (info/=0) return
    w=matmul(transpose(c),c)
    iw=matmul(ci,transpose(ci))
  end subroutine rwishart_draw
end module bayesm_rng
