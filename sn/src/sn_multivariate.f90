! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_multivariate
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use sn_kinds, only : dp, pi, tiny_dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_dimension_mismatch, &
                        sn_not_positive_definite
  use sn_math, only : normal_cdf, normal_logcdf, student_t_cdf, clamp_probability
  use sn_rng, only : sn_rng_state
  use sn_linalg, only : covariance_to_correlation, inverse_spd, logdet_spd, &
                        mvn_logpdf, rmvn, outer_product, quadratic_form
  use sn_mvn, only : mvn_cdf, mvt_cdf
  use sn_univariate, only : dsn, psn, rsn, dst, pst, rst, dsc, psc, rsc, &
                            delta_from_alpha, alpha_from_delta, zeta, b_nu
  implicit none
  private

  type, public :: sn_mv_params
    real(dp), allocatable :: xi(:)
    real(dp), allocatable :: omega(:,:)
    real(dp), allocatable :: alpha(:)
    real(dp) :: tau = 0.0_dp
  contains
    procedure :: dimension => sn_mv_dimension
    procedure :: validate => sn_mv_validate
  end type sn_mv_params

  type, public :: st_mv_params
    real(dp), allocatable :: xi(:)
    real(dp), allocatable :: omega(:,:)
    real(dp), allocatable :: alpha(:)
    real(dp) :: nu = 10.0_dp
  contains
    procedure :: dimension => st_mv_dimension
    procedure :: validate => st_mv_validate
  end type st_mv_params

  type, public :: mv_moments
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: marginal_skewness(:)
    real(dp) :: mardia_skewness = 0.0_dp
    real(dp) :: mardia_kurtosis = 0.0_dp
    integer :: status = sn_ok
  end type mv_moments

  public :: dmsn, pmsn, rmsn
  public :: dmst, pmst, rmst
  public :: dmsc, pmsc, rmsc
  public :: msn_moments, mst_moments
  public :: marginal_sn, affine_transform_sn, conditional_sn
  public :: delta_etc_mv, alpha_from_delta_mv

contains

  integer function sn_mv_dimension(self) result(d)
    class(sn_mv_params), intent(in) :: self
    if (allocated(self%alpha)) then
      d = size(self%alpha)
    else
      d = 0
    end if
  end function sn_mv_dimension

  integer function sn_mv_validate(self) result(status)
    class(sn_mv_params), intent(in) :: self
    real(dp), allocatable :: inv(:,:)
    integer :: d, ierr
    d = self%dimension()
    if (d < 1 .or. .not. allocated(self%xi) .or. .not. allocated(self%omega)) then
      status = sn_invalid_argument
      return
    end if
    if (size(self%xi) /= d .or. size(self%omega,1) /= d .or. size(self%omega,2) /= d) then
      status = sn_dimension_mismatch
      return
    end if
    call inverse_spd(self%omega,inv,ierr)
    status = ierr
  end function sn_mv_validate

  integer function st_mv_dimension(self) result(d)
    class(st_mv_params), intent(in) :: self
    if (allocated(self%alpha)) then
      d = size(self%alpha)
    else
      d = 0
    end if
  end function st_mv_dimension

  integer function st_mv_validate(self) result(status)
    class(st_mv_params), intent(in) :: self
    real(dp), allocatable :: inv(:,:)
    integer :: d, ierr
    d = self%dimension()
    if (d < 1 .or. self%nu <= 0.0_dp .or. .not. allocated(self%xi) .or. &
        .not. allocated(self%omega)) then
      status = sn_invalid_argument
      return
    end if
    if (size(self%xi) /= d .or. size(self%omega,1) /= d .or. size(self%omega,2) /= d) then
      status = sn_dimension_mismatch
      return
    end if
    call inverse_spd(self%omega,inv,ierr)
    status = ierr
  end function st_mv_validate

  subroutine delta_etc_mv(alpha, omega, delta, delta_star, alpha_star, cor, info)
    real(dp), intent(in) :: alpha(:), omega(:,:)
    real(dp), allocatable, intent(out) :: delta(:), cor(:,:)
    real(dp), intent(out) :: delta_star, alpha_star
    integer, intent(out) :: info
    real(dp), allocatable :: sd(:), ca(:)
    real(dp) :: asq

    if (size(omega,1) /= size(alpha) .or. size(omega,2) /= size(alpha)) then
      allocate(delta(0),cor(0,0))
      delta_star = 0.0_dp
      alpha_star = 0.0_dp
      info = sn_dimension_mismatch
      return
    end if
    call covariance_to_correlation(omega,cor,sd,info)
    if (info /= sn_ok) then
      allocate(delta(0))
      delta_star = 0.0_dp
      alpha_star = 0.0_dp
      return
    end if
    allocate(ca(size(alpha)),delta(size(alpha)))
    ca = matmul(cor,alpha)
    asq = dot_product(alpha,ca)
    if (asq < 0.0_dp) then
      info = sn_not_positive_definite
      delta = 0.0_dp
      delta_star = 0.0_dp
      alpha_star = 0.0_dp
      return
    end if
    delta = ca/sqrt(1.0_dp+asq)
    alpha_star = sqrt(asq)
    delta_star = sqrt(asq/(1.0_dp+asq))
    info = sn_ok
  end subroutine delta_etc_mv

  subroutine alpha_from_delta_mv(delta, cor, alpha, info)
    real(dp), intent(in) :: delta(:), cor(:,:)
    real(dp), allocatable, intent(out) :: alpha(:)
    integer, intent(out) :: info
    real(dp), allocatable :: inv(:,:), tmp(:)
    real(dp) :: q
    if (size(cor,1) /= size(delta) .or. size(cor,2) /= size(delta)) then
      allocate(alpha(0))
      info = sn_dimension_mismatch
      return
    end if
    call inverse_spd(cor,inv,info)
    if (info /= sn_ok) then
      allocate(alpha(0))
      return
    end if
    allocate(tmp(size(delta)),alpha(size(delta)))
    tmp = matmul(inv,delta)
    q = dot_product(delta,tmp)
    if (q >= 1.0_dp) then
      alpha = sign(huge(1.0_dp),tmp)
    else
      alpha = tmp/sqrt(max(tiny_dp,1.0_dp-q))
    end if
    info = sn_ok
  end subroutine alpha_from_delta_mv

  real(dp) function dmsn(x, params, log_pdf, info) result(value)
    real(dp), intent(in) :: x(:)
    type(sn_mv_params), intent(in) :: params
    logical, intent(in), optional :: log_pdf
    integer, intent(out), optional :: info
    real(dp), allocatable :: cor(:,:), sd(:), z(:), ca(:)
    real(dp) :: lp, alpha0, linear
    integer :: d, ierr
    logical :: give_log

    give_log = .false.
    if (present(log_pdf)) give_log = log_pdf
    ierr = params%validate()
    d = params%dimension()
    if (ierr /= sn_ok .or. size(x) /= d) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    allocate(z(d),ca(d))
    z = (x-params%xi)/sd
    ca = matmul(cor,params%alpha)
    alpha0 = params%tau*sqrt(1.0_dp+dot_product(params%alpha,ca))
    linear = alpha0+dot_product(params%alpha,z)
    lp = mvn_logpdf(x,params%xi,params%omega,ierr)+normal_logcdf(linear)-normal_logcdf(params%tau)
    if (give_log) then
      value = lp
    else
      value = exp(lp)
    end if
    if (present(info)) info = ierr
  end function dmsn

  real(dp) function pmsn(x, params, info, samples) result(value)
    real(dp), intent(in) :: x(:)
    type(sn_mv_params), intent(in) :: params
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: delta(:), cor(:,:), sd(:), big(:,:), upper(:), zero(:)
    real(dp) :: ds, as, den
    integer :: d, ierr, ns

    ierr = params%validate()
    d = params%dimension()
    if (ierr /= sn_ok .or. size(x) /= d) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,ierr)
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    allocate(big(d+1,d+1),upper(d+1),zero(d+1))
    big = 0.0_dp
    big(1,1) = 1.0_dp
    big(1,2:d+1) = -delta
    big(2:d+1,1) = -delta
    big(2:d+1,2:d+1) = cor
    upper(1) = params%tau
    upper(2:d+1) = (x-params%xi)/sd
    zero = 0.0_dp
    ns = 32768
    if (present(samples)) ns = samples
    den = normal_cdf(params%tau)
    value = clamp_probability(mvn_cdf(upper,zero,big,ierr,ns)/max(den,tiny_dp))
    if (present(info)) info = ierr
  end function pmsn

  subroutine rmsn(rng, n, params, x, info)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    type(sn_mv_params), intent(in) :: params
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: delta(:), cor(:,:), sd(:), residual(:,:), u(:,:), z(:)
    real(dp) :: ds, as, u0
    integer :: d, ierr, i

    info = params%validate()
    d = params%dimension()
    if (info /= sn_ok .or. n < 0) then
      allocate(x(0,0))
      return
    end if
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,ierr)
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    residual = cor-outer_product(delta,delta)
    call rmvn(rng,n,[ (0.0_dp,i=1,d) ],residual,u,ierr)
    if (ierr /= sn_ok) then
      allocate(x(0,0))
      info = ierr
      return
    end if
    allocate(x(n,d),z(d))
    do i=1,n
      if (abs(params%tau) <= tiny_dp) then
        u0 = abs(rng%normal())
      else
        u0 = rng%truncated_normal_lower(-params%tau)
      end if
      z = delta*u0+u(i,:)
      x(i,:) = params%xi+sd*z
    end do
    info = sn_ok
  end subroutine rmsn

  real(dp) function dmst(x, params, log_pdf, info) result(value)
    real(dp), intent(in) :: x(:)
    type(st_mv_params), intent(in) :: params
    logical, intent(in), optional :: log_pdf
    integer, intent(out), optional :: info
    real(dp), allocatable :: inv(:,:), sd(:), cor(:,:), z(:), diff(:)
    real(dp) :: q, linear, logdet, lp
    integer :: d, ierr
    logical :: give_log

    give_log = .false.
    if (present(log_pdf)) give_log = log_pdf
    ierr = params%validate()
    d = params%dimension()
    if (ierr /= sn_ok .or. size(x) /= d) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    if (params%nu > 1.0e8_dp) then
      sn_block: block
        type(sn_mv_params) :: p
        allocate(p%xi(d),p%omega(d,d),p%alpha(d))
        p%xi=params%xi; p%omega=params%omega; p%alpha=params%alpha; p%tau=0.0_dp
        value=dmsn(x,p,give_log,ierr)
      end block sn_block
      if (present(info)) info=ierr
      return
    end if
    call inverse_spd(params%omega,inv,ierr)
    call logdet_spd(params%omega,logdet,ierr)
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    allocate(diff(d),z(d))
    diff=x-params%xi
    z=diff/sd
    q=quadratic_form(diff,inv)
    linear=dot_product(params%alpha,z)*sqrt((params%nu+real(d,dp))/(params%nu+q))
    lp=log(2.0_dp)+log_gamma(0.5_dp*(params%nu+real(d,dp)))-log_gamma(0.5_dp*params%nu) &
       -0.5_dp*(real(d,dp)*log(params%nu*pi)+logdet) &
       -0.5_dp*(params%nu+real(d,dp))*log(1.0_dp+q/params%nu) &
       +log(max(student_t_cdf(linear,params%nu+real(d,dp)),tiny_dp))
    if (give_log) then
      value=lp
    else
      value=exp(lp)
    end if
    if (present(info)) info=ierr
  end function dmst

  real(dp) function pmst(x, params, info, samples) result(value)
    real(dp), intent(in) :: x(:)
    type(st_mv_params), intent(in) :: params
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: delta(:), cor(:,:), sd(:), big(:,:), upper(:), zero(:)
    real(dp) :: ds, as
    integer :: d, ierr, ns

    ierr=params%validate()
    d=params%dimension()
    if (ierr/=sn_ok .or. size(x)/=d) then
      value=ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info=merge(ierr,sn_dimension_mismatch,ierr/=sn_ok)
      return
    end if
    if (d==1) then
      value=pst(x(1),params%xi(1),sqrt(params%omega(1,1)),params%alpha(1),params%nu)
      if (present(info)) info=sn_ok
      return
    end if
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,ierr)
    call covariance_to_correlation(params%omega,cor,sd,ierr)
    allocate(big(d+1,d+1),upper(d+1),zero(d+1))
    big=0.0_dp
    big(1,1)=1.0_dp
    big(1,2:d+1)=-delta
    big(2:d+1,1)=-delta
    big(2:d+1,2:d+1)=cor
    upper(1)=0.0_dp
    upper(2:d+1)=(x-params%xi)/sd
    zero=0.0_dp
    ns=65536
    if (present(samples)) ns=samples
    value=clamp_probability(2.0_dp*mvt_cdf(upper,zero,big,params%nu,ierr,ns))
    if (present(info)) info=ierr
  end function pmst

  subroutine rmst(rng,n,params,x,info)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    type(st_mv_params), intent(in) :: params
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    type(sn_mv_params) :: p
    real(dp), allocatable :: z(:,:)
    real(dp) :: v
    integer :: d,i

    info=params%validate()
    d=params%dimension()
    if (info/=sn_ok .or. n<0) then
      allocate(x(0,0))
      return
    end if
    allocate(p%xi(d),p%omega(d,d),p%alpha(d))
    p%xi=0.0_dp
    p%omega=params%omega
    p%alpha=params%alpha
    p%tau=0.0_dp
    call rmsn(rng,n,p,z,info)
    if (info/=sn_ok) then
      allocate(x(0,0))
      return
    end if
    allocate(x(n,d))
    do i=1,n
      if (params%nu<1.0e8_dp) then
        v=rng%chi_square(params%nu)/params%nu
      else
        v=1.0_dp
      end if
      x(i,:)=params%xi+z(i,:)/sqrt(v)
    end do
    info=sn_ok
  end subroutine rmst

  real(dp) function dmsc(x,params,log_pdf,info) result(value)
    real(dp), intent(in) :: x(:)
    type(st_mv_params), intent(in) :: params
    logical, intent(in), optional :: log_pdf
    integer, intent(out), optional :: info
    type(st_mv_params) :: p
    p=params
    p%nu=1.0_dp
    value=dmst(x,p,log_pdf,info)
  end function dmsc

  real(dp) function pmsc(x,params,info,samples) result(value)
    real(dp), intent(in) :: x(:)
    type(st_mv_params), intent(in) :: params
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    type(st_mv_params) :: p
    p=params
    p%nu=1.0_dp
    value=pmst(x,p,info,samples)
  end function pmsc

  subroutine rmsc(rng,n,params,x,info)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    type(st_mv_params), intent(in) :: params
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    type(st_mv_params) :: p
    p=params
    p%nu=1.0_dp
    call rmst(rng,n,p,x,info)
  end subroutine rmsc

  function msn_moments(params) result(out)
    type(sn_mv_params), intent(in) :: params
    type(mv_moments) :: out
    real(dp), allocatable :: delta(:),cor(:,:),sd(:),shift(:),sdz(:)
    real(dp) :: ds,as,ratio2
    integer :: d,i,info

    out%status=params%validate()
    if (out%status/=sn_ok) return
    d=params%dimension()
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(out%mean(d),out%covariance(d,d),out%marginal_skewness(d),shift(d),sdz(d))
    shift=sd*delta
    out%mean=params%xi+zeta(1,params%tau)*shift
    out%covariance=params%omega+zeta(2,params%tau)*outer_product(shift,shift)
    do i=1,d
      sdz(i)=sqrt(1.0_dp+zeta(2,params%tau)*delta(i)**2)
      out%marginal_skewness(i)=zeta(3,params%tau)*(delta(i)/sdz(i))**3
    end do
    ratio2=ds*ds/(1.0_dp+zeta(2,params%tau)*ds*ds)
    out%mardia_skewness=zeta(3,params%tau)**2*ratio2**3
    out%mardia_kurtosis=zeta(4,params%tau)*ratio2**2
    out%status=sn_ok
  end function msn_moments

  function mst_moments(params) result(out)
    type(st_mv_params), intent(in) :: params
    type(mv_moments) :: out
    real(dp), allocatable :: delta(:),cor(:,:),sd(:),shift(:)
    real(dp) :: ds,as,bv
    integer :: d,info

    out%status=params%validate()
    if (out%status/=sn_ok) return
    d=params%dimension()
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(out%mean(d),out%covariance(d,d),out%marginal_skewness(d),shift(d))
    shift=sd*delta
    if (params%nu>1.0_dp) then
      bv=b_nu(params%nu)
      out%mean=params%xi+bv*shift
    else
      out%mean=ieee_value(0.0_dp,ieee_quiet_nan)
      bv=ieee_value(0.0_dp,ieee_quiet_nan)
    end if
    if (params%nu>2.0_dp) then
      out%covariance=params%omega*params%nu/(params%nu-2.0_dp)-outer_product(bv*shift,bv*shift)
    else
      out%covariance=ieee_value(0.0_dp,ieee_quiet_nan)
    end if
    out%marginal_skewness=ieee_value(0.0_dp,ieee_quiet_nan)
    out%status=sn_ok
  end function mst_moments

  subroutine marginal_sn(params,components,out,info)
    type(sn_mv_params), intent(in) :: params
    integer, intent(in) :: components(:)
    type(sn_mv_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: delta(:),cor(:,:),delta1(:),cor1(:,:),alpha1(:)
    real(dp) :: ds,as
    integer :: d,k,i,j

    info=params%validate()
    d=params%dimension()
    k=size(components)
    if (info/=sn_ok .or. k<1 .or. any(components<1) .or. any(components>d)) then
      info=sn_invalid_argument
      return
    end if
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    allocate(out%xi(k),out%omega(k,k),delta1(k),cor1(k,k))
    do i=1,k
      out%xi(i)=params%xi(components(i))
      delta1(i)=delta(components(i))
      do j=1,k
        out%omega(i,j)=params%omega(components(i),components(j))
        cor1(i,j)=cor(components(i),components(j))
      end do
    end do
    call alpha_from_delta_mv(delta1,cor1,alpha1,info)
    allocate(out%alpha(k))
    out%alpha=alpha1
    out%tau=params%tau
  end subroutine marginal_sn

  subroutine affine_transform_sn(params,a,transform,out,info)
    type(sn_mv_params), intent(in) :: params
    real(dp), intent(in) :: a(:),transform(:,:)
    type(sn_mv_params), intent(out) :: out
    integer, intent(out) :: info
    real(dp), allocatable :: delta(:),cor(:,:),sd(:),cross(:),newcor(:,:),newsd(:),newdelta(:),newalpha(:)
    real(dp) :: ds,as
    integer :: d,k

    info=params%validate()
    d=params%dimension()
    k=size(transform,1)
    if (info/=sn_ok .or. size(transform,2)/=d .or. size(a)/=k) then
      info=sn_dimension_mismatch
      return
    end if
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(cross(d))
    cross=sd*delta
    allocate(out%xi(k),out%omega(k,k))
    out%xi=a+matmul(transform,params%xi)
    out%omega=matmul(transform,matmul(params%omega,transpose(transform)))
    call covariance_to_correlation(out%omega,newcor,newsd,info)
    allocate(newdelta(k))
    newdelta=matmul(transform,cross)/newsd
    call alpha_from_delta_mv(newdelta,newcor,newalpha,info)
    allocate(out%alpha(k))
    out%alpha=newalpha
    out%tau=params%tau
  end subroutine affine_transform_sn

  subroutine conditional_sn(params,fixed_components,fixed_values,out,info)
    type(sn_mv_params), intent(in) :: params
    integer, intent(in) :: fixed_components(:)
    real(dp), intent(in) :: fixed_values(:)
    type(sn_mv_params), intent(out) :: out
    integer, intent(out) :: info
    integer, allocatable :: free_components(:)
    real(dp), allocatable :: sd(:),cor(:,:),r11(:,:),r12(:,:),r21(:,:),r22(:,:), &
                            ir11(:,:),alpha1(:),alpha2(:),tmp(:),rcond(:,:), &
                            o11(:,:),o12(:,:),o21(:,:),o22(:,:),io11(:,:),reg(:,:), &
                            new_sd(:),alpha_cond(:)
    real(dp) :: asum,tau_cond
    integer :: d,h,k,i,j,nfree

    info=params%validate()
    d=params%dimension()
    h=size(fixed_components)
    if (info/=sn_ok .or. size(fixed_values)/=h .or. h<1 .or. h>=d .or. &
        any(fixed_components<1) .or. any(fixed_components>d)) then
      info=sn_invalid_argument
      return
    end if
    allocate(free_components(d-h))
    k=0
    do i=1,d
      if (.not. any(fixed_components==i)) then
        k=k+1
        free_components(k)=i
      end if
    end do
    nfree=d-h
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(r11(h,h),r12(h,nfree),r21(nfree,h),r22(nfree,nfree),alpha1(h),alpha2(nfree))
    do i=1,h
      alpha1(i)=params%alpha(fixed_components(i))
      do j=1,h
        r11(i,j)=cor(fixed_components(i),fixed_components(j))
      end do
      do j=1,nfree
        r12(i,j)=cor(fixed_components(i),free_components(j))
      end do
    end do
    do i=1,nfree
      alpha2(i)=params%alpha(free_components(i))
      do j=1,h
        r21(i,j)=cor(free_components(i),fixed_components(j))
      end do
      do j=1,nfree
        r22(i,j)=cor(free_components(i),free_components(j))
      end do
    end do
    call inverse_spd(r11,ir11,info)
    if (info/=sn_ok) return
    rcond=r22-matmul(r21,matmul(ir11,r12))
    asum=dot_product(alpha2,matmul(rcond,alpha2))
    allocate(tmp(h))
    tmp=alpha1+matmul(ir11,matmul(r12,alpha2))
    tmp=tmp/sqrt(1.0_dp+asum)
    tau_cond=params%tau*sqrt(1.0_dp+dot_product(tmp,matmul(ir11,tmp)))+ &
             dot_product(tmp,(fixed_values-params%xi(fixed_components))/sd(fixed_components))

    allocate(o11(h,h),o12(h,nfree),o21(nfree,h),o22(nfree,nfree))
    do i=1,h
      do j=1,h
        o11(i,j)=params%omega(fixed_components(i),fixed_components(j))
      end do
      do j=1,nfree
        o12(i,j)=params%omega(fixed_components(i),free_components(j))
      end do
    end do
    do i=1,nfree
      do j=1,h
        o21(i,j)=params%omega(free_components(i),fixed_components(j))
      end do
      do j=1,nfree
        o22(i,j)=params%omega(free_components(i),free_components(j))
      end do
    end do
    call inverse_spd(o11,io11,info)
    if (info/=sn_ok) return
    reg=matmul(o21,io11)
    allocate(out%xi(nfree),out%omega(nfree,nfree),new_sd(nfree),alpha_cond(nfree),out%alpha(nfree))
    out%xi=params%xi(free_components)+matmul(reg,fixed_values-params%xi(fixed_components))
    out%omega=o22-matmul(reg,o12)
    do i=1,nfree
      new_sd(i)=sqrt(out%omega(i,i))
      alpha_cond(i)=new_sd(i)/sd(free_components(i))*alpha2(i)
    end do
    out%alpha=alpha_cond
    out%tau=tau_cond
    info=sn_ok
  end subroutine conditional_sn

end module sn_multivariate
