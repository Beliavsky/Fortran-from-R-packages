! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_sun
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use iso_fortran_env, only : int64
  use sn_kinds, only : dp, tiny_dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_dimension_mismatch, &
                        sn_not_positive_definite
  use sn_math, only : clamp_probability
  use sn_rng, only : sn_rng_state
  use sn_linalg, only : covariance_to_correlation, inverse_spd, mvn_logpdf, &
                        rmvn, outer_product, identity_matrix
  use sn_mvn, only : mvn_cdf, rtmvn_lower
  use sn_multivariate, only : sn_mv_params, delta_etc_mv
  implicit none
  private

  type, public :: sun_params
    real(dp), allocatable :: xi(:)
    real(dp), allocatable :: omega(:,:)
    real(dp), allocatable :: delta(:,:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: gamma(:,:)
  contains
    procedure :: dimension => sun_dimension
    procedure :: hidden_dimension => sun_hidden_dimension
    procedure :: validate => sun_validate
  end type sun_params

  type, public :: sun_moments_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    integer :: status = sn_ok
    integer :: accepted_samples = 0
  end type sun_moments_result

  public :: dsun, psun, rsun, sun_moments
  public :: marginal_sun, affine_transform_sun, convolution_sun, join_sun
  public :: conditional_sun_equal, conditional_sun_greater
  public :: sn_to_sun, csn_to_sun

contains

  integer function sun_dimension(self) result(d)
    class(sun_params), intent(in) :: self
    if (allocated(self%xi)) then
      d = size(self%xi)
    else
      d = 0
    end if
  end function sun_dimension

  integer function sun_hidden_dimension(self) result(m)
    class(sun_params), intent(in) :: self
    if (allocated(self%tau)) then
      m = size(self%tau)
    else
      m = 0
    end if
  end function sun_hidden_dimension

  integer function sun_validate(self) result(status)
    class(sun_params), intent(in) :: self
    real(dp), allocatable :: inv(:,:), sd(:), cor(:,:), big(:,:)
    integer :: d,m,info
    d=self%dimension()
    m=self%hidden_dimension()
    if (d<1 .or. m<1 .or. .not. allocated(self%omega) .or. &
        .not. allocated(self%delta) .or. .not. allocated(self%gamma)) then
      status=sn_invalid_argument
      return
    end if
    if (size(self%omega,1)/=d .or. size(self%omega,2)/=d .or. &
        size(self%delta,1)/=d .or. size(self%delta,2)/=m .or. &
        size(self%gamma,1)/=m .or. size(self%gamma,2)/=m) then
      status=sn_dimension_mismatch
      return
    end if
    call inverse_spd(self%omega,inv,info)
    if (info/=sn_ok) then
      status=info
      return
    end if
    call inverse_spd(self%gamma,inv,info)
    if (info/=sn_ok) then
      status=info
      return
    end if
    call covariance_to_correlation(self%omega,cor,sd,info)
    allocate(big(d+m,d+m))
    big(1:d,1:d)=cor
    big(1:d,d+1:d+m)=self%delta
    big(d+1:d+m,1:d)=transpose(self%delta)
    big(d+1:d+m,d+1:d+m)=self%gamma
    call inverse_spd(big,inv,info)
    if (info/=sn_ok) then
      status=sn_not_positive_definite
    else
      status=sn_ok
    end if
  end function sun_validate

  real(dp) function dsun(x,params,log_pdf,info,samples) result(value)
    real(dp), intent(in) :: x(:)
    type(sun_params), intent(in) :: params
    logical, intent(in), optional :: log_pdf
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: cor(:,:),sd(:),invcor(:,:),z(:),mean0(:),cond_cov(:,:),upper(:)
    real(dp) :: p1,p2,lp
    integer :: d,m,ierr,ns
    logical :: give_log

    give_log=.false.
    if (present(log_pdf)) give_log=log_pdf
    ierr=params%validate()
    d=params%dimension()
    m=params%hidden_dimension()
    if (ierr/=sn_ok .or. size(x)/=d) then
      value=ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info=merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    call inverse_spd(cor,invcor,ierr)
    allocate(z(d),mean0(m),cond_cov(m,m),upper(m))
    z=(x-params%xi)/sd
    mean0=matmul(transpose(params%delta),matmul(invcor,z))
    cond_cov=params%gamma-matmul(transpose(params%delta),matmul(invcor,params%delta))
    upper=params%tau+mean0
    ns=32768
    if (present(samples)) ns=samples
    p1=mvn_cdf(upper,[ (0.0_dp,ierr=1,m) ],cond_cov,ierr,ns)
    p2=mvn_cdf(params%tau,[ (0.0_dp,ierr=1,m) ],params%gamma,ierr,ns)
    lp=mvn_logpdf(x,params%xi,params%omega,ierr)+log(max(p1,tiny_dp))-log(max(p2,tiny_dp))
    if (give_log) then
      value=lp
    else
      value=exp(lp)
    end if
    if (present(info)) info=ierr
  end function dsun

  real(dp) function psun(x,params,info,samples) result(value)
    real(dp), intent(in) :: x(:)
    type(sun_params), intent(in) :: params
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: cor(:,:),sd(:),big(:,:),upper(:),zero(:)
    real(dp) :: p1,p2
    integer :: d,m,ierr,ns

    ierr=params%validate()
    d=params%dimension()
    m=params%hidden_dimension()
    if (ierr/=sn_ok .or. size(x)/=d) then
      value=ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info=merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    allocate(big(d+m,d+m),upper(d+m),zero(d+m))
    big(1:d,1:d)=cor
    big(1:d,d+1:d+m)=-params%delta
    big(d+1:d+m,1:d)=-transpose(params%delta)
    big(d+1:d+m,d+1:d+m)=params%gamma
    upper(1:d)=(x-params%xi)/sd
    upper(d+1:d+m)=params%tau
    zero=0.0_dp
    ns=65536
    if (present(samples)) ns=samples
    p1=mvn_cdf(upper,zero,big,ierr,ns)
    p2=mvn_cdf(params%tau,zero(1:m),params%gamma,ierr,ns)
    value=clamp_probability(p1/max(p2,tiny_dp))
    if (present(info)) info=ierr
  end function psun

  subroutine rsun(rng,n,params,x,info)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    type(sun_params), intent(in) :: params
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: cor(:,:),sd(:),invgamma(:,:),a(:,:),psi(:,:),u0(:,:),u1(:,:),zero_d(:),lower(:)
    integer :: d,m,i

    info=params%validate()
    d=params%dimension()
    m=params%hidden_dimension()
    if (info/=sn_ok .or. n<0) then
      allocate(x(0,0))
      return
    end if
    call covariance_to_correlation(params%omega,cor,sd,info)
    call inverse_spd(params%gamma,invgamma,info)
    a=matmul(params%delta,invgamma)
    psi=cor-matmul(a,transpose(params%delta))
    allocate(zero_d(d),lower(m))
    zero_d=0.0_dp
    lower=-params%tau
    call rmvn(rng,n,zero_d,psi,u0,info)
    if (info/=sn_ok) then
      allocate(x(0,0))
      return
    end if
    call rtmvn_lower(rng,n,[ (0.0_dp,i=1,m) ],params%gamma,lower,u1,info,max_attempts=max(100000,10000*n))
    if (info/=sn_ok) then
      allocate(x(0,0))
      return
    end if
    allocate(x(n,d))
    do i=1,n
      x(i,:)=params%xi+sd*(u0(i,:)+matmul(a,u1(i,:)))
    end do
    info=sn_ok
  end subroutine rsun

  function sun_moments(params,samples) result(out)
    type(sun_params), intent(in) :: params
    integer, intent(in), optional :: samples
    type(sun_moments_result) :: out
    type(sn_rng_state) :: rng
    real(dp), allocatable :: cor(:,:),sd(:),invgamma(:,:),a(:,:),psi(:,:),u(:,:),lower(:), &
                            eu(:),vu(:,:),z(:)
    integer :: d,m,n,info,i

    out%status=params%validate()
    if (out%status/=sn_ok) return
    d=params%dimension()
    m=params%hidden_dimension()
    n=100000
    if (present(samples)) n=max(1000,samples)
    call covariance_to_correlation(params%omega,cor,sd,info)
    call inverse_spd(params%gamma,invgamma,info)
    a=matmul(params%delta,invgamma)
    psi=cor-matmul(a,transpose(params%delta))
    allocate(lower(m))
    lower=-params%tau
    call rng%seed(78193085935_int64)
    call rtmvn_lower(rng,n,[ (0.0_dp,i=1,m) ],params%gamma,lower,u,info,max_attempts=max(1000000,100*n))
    if (info/=sn_ok) then
      out%status=info
      return
    end if
    allocate(eu(m),vu(m,m),z(m),out%mean(d),out%covariance(d,d))
    eu=sum(u,dim=1)/real(n,dp)
    vu=0.0_dp
    do i=1,n
      z=u(i,:)-eu
      vu=vu+outer_product(z,z)
    end do
    vu=vu/real(max(1,n-1),dp)
    out%mean=params%xi+sd*matmul(a,eu)
    out%covariance=diag_left_right(sd,psi+matmul(a,matmul(vu,transpose(a))))
    out%accepted_samples=n
    out%status=sn_ok
  end function sun_moments

  subroutine marginal_sun(params,components,out,info)
    type(sun_params), intent(in) :: params
    integer, intent(in) :: components(:)
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    integer :: d,m,k,i,j
    info=params%validate()
    d=params%dimension()
    m=params%hidden_dimension()
    k=size(components)
    if (info/=sn_ok .or. k<1 .or. any(components<1) .or. any(components>d)) then
      info=sn_invalid_argument
      return
    end if
    allocate(out%xi(k),out%omega(k,k),out%delta(k,m),out%tau(m),out%gamma(m,m))
    do i=1,k
      out%xi(i)=params%xi(components(i))
      out%delta(i,:)=params%delta(components(i),:)
      do j=1,k
        out%omega(i,j)=params%omega(components(i),components(j))
      end do
    end do
    out%tau=params%tau
    out%gamma=params%gamma
    info=out%validate()
  end subroutine marginal_sun

  subroutine affine_transform_sun(params,a,transform,out,info)
    type(sun_params), intent(in) :: params
    real(dp), intent(in) :: a(:),transform(:,:)
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: sd(:),cor(:,:),newsd(:),newcor(:,:),cross(:,:)
    integer :: d,m,k,i
    info=params%validate()
    d=params%dimension()
    m=params%hidden_dimension()
    k=size(transform,1)
    if (info/=sn_ok .or. size(transform,2)/=d .or. size(a)/=k) then
      info=sn_dimension_mismatch
      return
    end if
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(out%xi(k),out%omega(k,k),out%delta(k,m),out%tau(m),out%gamma(m,m),cross(d,m))
    out%xi=a+matmul(transform,params%xi)
    out%omega=matmul(transform,matmul(params%omega,transpose(transform)))
    call covariance_to_correlation(out%omega,newcor,newsd,info)
    do i=1,d
      cross(i,:)=sd(i)*params%delta(i,:)
    end do
    out%delta=matmul(transform,cross)
    do i=1,k
      out%delta(i,:)=out%delta(i,:)/newsd(i)
    end do
    out%tau=params%tau
    out%gamma=params%gamma
    info=out%validate()
  end subroutine affine_transform_sun

  subroutine convolution_sun(a,b,out,info)
    type(sun_params), intent(in) :: a,b
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: sda(:),sdb(:),sd(:),cora(:,:),corb(:,:),cor(:,:)
    integer :: d,m1,m2,i
    info=a%validate()
    if (info/=sn_ok) return
    info=b%validate()
    if (info/=sn_ok) return
    d=a%dimension()
    if (b%dimension()/=d) then
      info=sn_dimension_mismatch
      return
    end if
    m1=a%hidden_dimension()
    m2=b%hidden_dimension()
    call covariance_to_correlation(a%omega,cora,sda,info)
    call covariance_to_correlation(b%omega,corb,sdb,info)
    allocate(out%xi(d),out%omega(d,d),out%delta(d,m1+m2),out%tau(m1+m2),out%gamma(m1+m2,m1+m2))
    out%xi=a%xi+b%xi
    out%omega=a%omega+b%omega
    call covariance_to_correlation(out%omega,cor,sd,info)
    do i=1,d
      out%delta(i,1:m1)=sda(i)/sd(i)*a%delta(i,:)
      out%delta(i,m1+1:m1+m2)=sdb(i)/sd(i)*b%delta(i,:)
    end do
    out%tau=[a%tau,b%tau]
    out%gamma=0.0_dp
    out%gamma(1:m1,1:m1)=a%gamma
    out%gamma(m1+1:m1+m2,m1+1:m1+m2)=b%gamma
    info=out%validate()
  end subroutine convolution_sun

  subroutine join_sun(a,b,out,info)
    type(sun_params), intent(in) :: a,b
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    integer :: d1,d2,m1,m2
    info=a%validate()
    if (info/=sn_ok) return
    info=b%validate()
    if (info/=sn_ok) return
    d1=a%dimension(); d2=b%dimension(); m1=a%hidden_dimension(); m2=b%hidden_dimension()
    allocate(out%xi(d1+d2),out%omega(d1+d2,d1+d2),out%delta(d1+d2,m1+m2), &
             out%tau(m1+m2),out%gamma(m1+m2,m1+m2))
    out%xi=[a%xi,b%xi]
    out%omega=0.0_dp
    out%omega(1:d1,1:d1)=a%omega
    out%omega(d1+1:d1+d2,d1+1:d1+d2)=b%omega
    out%delta=0.0_dp
    out%delta(1:d1,1:m1)=a%delta
    out%delta(d1+1:d1+d2,m1+1:m1+m2)=b%delta
    out%tau=[a%tau,b%tau]
    out%gamma=0.0_dp
    out%gamma(1:m1,1:m1)=a%gamma
    out%gamma(m1+1:m1+m2,m1+1:m1+m2)=b%gamma
    info=out%validate()
  end subroutine join_sun

  subroutine conditional_sun_equal(params,components,values,out,info)
    type(sun_params), intent(in) :: params
    integer, intent(in) :: components(:)
    real(dp), intent(in) :: values(:)
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    integer, allocatable :: freec(:)
    real(dp), allocatable :: sd(:),cor(:,:),o11(:,:),o12(:,:),o21(:,:),o22(:,:),io11(:,:), &
                            r11(:,:),ir11(:,:),delta1(:,:),delta2(:,:),tmp(:,:),gamma_c(:,:), &
                            s(:),reg(:,:),v0(:)
    integer :: d,m,h,nf,i,j,k

    info=params%validate()
    d=params%dimension(); m=params%hidden_dimension(); h=size(components)
    if (info/=sn_ok .or. size(values)/=h .or. h<1 .or. h>=d .or. &
        any(components<1) .or. any(components>d)) then
      info=sn_invalid_argument
      return
    end if
    nf=d-h
    allocate(freec(nf)); k=0
    do i=1,d
      if (.not. any(components==i)) then
        k=k+1; freec(k)=i
      end if
    end do
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(o11(h,h),o12(h,nf),o21(nf,h),o22(nf,nf),r11(h,h),delta1(h,m),delta2(nf,m),v0(h))
    do i=1,h
      v0(i)=values(i)-params%xi(components(i))
      delta1(i,:)=params%delta(components(i),:)
      do j=1,h
        o11(i,j)=params%omega(components(i),components(j))
        r11(i,j)=cor(components(i),components(j))
      end do
      do j=1,nf
        o12(i,j)=params%omega(components(i),freec(j))
      end do
    end do
    do i=1,nf
      delta2(i,:)=params%delta(freec(i),:)
      do j=1,h
        o21(i,j)=params%omega(freec(i),components(j))
      end do
      do j=1,nf
        o22(i,j)=params%omega(freec(i),freec(j))
      end do
    end do
    call inverse_spd(o11,io11,info)
    call inverse_spd(r11,ir11,info)
    reg=matmul(o21,io11)
    allocate(out%xi(nf),out%omega(nf,nf))
    out%xi=params%xi(freec)+matmul(reg,v0)
    out%omega=o22-matmul(reg,o12)
    tmp=delta2-matmul(cor(freec,components),matmul(ir11,delta1))
    gamma_c=params%gamma-matmul(transpose(delta1),matmul(ir11,delta1))
    allocate(s(m),out%delta(nf,m),out%tau(m),out%gamma(m,m))
    do i=1,m
      s(i)=sqrt(gamma_c(i,i))
    end do
    out%delta=tmp
    do j=1,m
      out%delta(:,j)=out%delta(:,j)/s(j)
    end do
    out%tau=params%tau+matmul(transpose(delta1),matmul(ir11,v0/sd(components)))
    out%tau=out%tau/s
    out%gamma=gamma_c
    do j=1,m
      do i=1,m
        out%gamma(i,j)=out%gamma(i,j)/(s(i)*s(j))
      end do
    end do
    info=out%validate()
  end subroutine conditional_sun_equal

  subroutine conditional_sun_greater(params,components,values,out,info)
    type(sun_params), intent(in) :: params
    integer, intent(in) :: components(:)
    real(dp), intent(in) :: values(:)
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    integer, allocatable :: freec(:)
    real(dp), allocatable :: sd(:),cor(:,:)
    integer :: d,m,h,nf,i,j,k

    info=params%validate()
    d=params%dimension(); m=params%hidden_dimension(); h=size(components)
    if (info/=sn_ok .or. size(values)/=h .or. h<1 .or. h>=d .or. &
        any(components<1) .or. any(components>d)) then
      info=sn_invalid_argument
      return
    end if
    nf=d-h
    allocate(freec(nf)); k=0
    do i=1,d
      if (.not. any(components==i)) then
        k=k+1; freec(k)=i
      end if
    end do
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(out%xi(nf),out%omega(nf,nf),out%delta(nf,m+h),out%tau(m+h),out%gamma(m+h,m+h))
    out%xi=params%xi(freec)
    do i=1,nf
      do j=1,nf
        out%omega(i,j)=params%omega(freec(i),freec(j))
      end do
      out%delta(i,1:h)=cor(freec(i),components)
      out%delta(i,h+1:h+m)=params%delta(freec(i),:)
    end do
    out%tau(1:h)=(params%xi(components)-values)/sd(components)
    out%tau(h+1:h+m)=params%tau
    out%gamma(1:h,1:h)=cor(components,components)
    out%gamma(1:h,h+1:h+m)=params%delta(components,:)
    out%gamma(h+1:h+m,1:h)=transpose(params%delta(components,:))
    out%gamma(h+1:h+m,h+1:h+m)=params%gamma
    info=out%validate()
  end subroutine conditional_sun_greater

  subroutine sn_to_sun(params,out,info)
    type(sn_mv_params), intent(in) :: params
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: delta(:),cor(:,:)
    real(dp) :: ds,as
    integer :: d
    info=params%validate()
    if (info/=sn_ok) return
    d=params%dimension()
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    allocate(out%xi(d),out%omega(d,d),out%delta(d,1),out%tau(1),out%gamma(1,1))
    out%xi=params%xi
    out%omega=params%omega
    out%delta(:,1)=delta
    out%tau(1)=params%tau
    out%gamma(1,1)=1.0_dp
    info=out%validate()
  end subroutine sn_to_sun

  subroutine csn_to_sun(mu,sigma,dmat,nu,delta_csn,out,info)
    real(dp), intent(in) :: mu(:),sigma(:,:),dmat(:,:),nu(:),delta_csn(:,:)
    type(sun_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: ds(:,:),mstar(:,:),cor(:,:),sd(:)
    integer :: p,q
    p=size(mu); q=size(nu)
    if (size(sigma,1)/=p .or. size(sigma,2)/=p .or. size(dmat,1)/=q .or. &
        size(dmat,2)/=p .or. size(delta_csn,1)/=q .or. size(delta_csn,2)/=q) then
      info=sn_dimension_mismatch
      return
    end if
    ds=matmul(dmat,sigma)
    allocate(mstar(p+q,p+q))
    mstar(1:p,1:p)=sigma
    mstar(1:p,p+1:p+q)=transpose(ds)
    mstar(p+1:p+q,1:p)=ds
    mstar(p+1:p+q,p+1:p+q)=delta_csn+matmul(ds,transpose(dmat))
    call covariance_to_correlation(mstar,cor,sd,info)
    if (info/=sn_ok) return
    allocate(out%xi(p),out%omega(p,p),out%delta(p,q),out%tau(q),out%gamma(q,q))
    out%xi=mu
    out%omega=sigma
    out%delta=cor(1:p,p+1:p+q)
    out%tau=-nu
    out%gamma=cor(p+1:p+q,p+1:p+q)
    info=out%validate()
  end subroutine csn_to_sun

  pure function diag_left_right(sd,a) result(out)
    real(dp), intent(in) :: sd(:),a(:,:)
    real(dp) :: out(size(sd),size(sd))
    integer :: i,j
    do j=1,size(sd)
      do i=1,size(sd)
        out(i,j)=sd(i)*a(i,j)*sd(j)
      end do
    end do
  end function diag_left_right

end module sn_sun
