! SPDX-License-Identifier: GPL-3.0-only
module matrix_normal_distribution
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
  use mvtnorm_kinds, only : dp, log_two_pi
  use mvtnorm_types, only : probability_control, probability_result
  use mvtnorm_linalg, only : cholesky_lower, logdet_cholesky, solve_spd
  use mvtnorm_distributions, only : rmvnorm
  use mvtnorm_probabilities, only : pmvnorm
  use matrix_normal_utils, only : vec, kronecker_product, is_positive_definite
  implicit none
  private
  public :: dmatnorm, pmatnorm, rmatnorm, check_matnorm

  interface pmatnorm
    module procedure pmatnorm_bounds
    module procedure pmatnorm_upper
    module procedure pmatnorm_full
  end interface pmatnorm

  interface rmatnorm
    module procedure rmatnorm_one
    module procedure rmatnorm_many
  end interface rmatnorm

contains

  subroutine check_matnorm(m,u,v,ok,message,tol)
    real(dp), intent(in) :: m(:,:),u(:,:),v(:,:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: tol
    real(dp) :: eps
    eps=sqrt(epsilon(1.0_dp))
    if(present(tol)) eps=tol
    ok=.false.
    message=''
    if(size(u,1)/=size(u,2)) then
    message='U must be square'
    return
    end if
    if(size(v,1)/=size(v,2)) then
    message='V must be square'
    return
    end if
    if(size(m,1)/=size(u,1)) then
    message='rows of M must match order of U'
    return
    end if
    if(size(m,2)/=size(v,1)) then
    message='columns of M must match order of V'
    return
    end if
    if(.not.is_positive_definite(u,eps)) then
    message='U is not positive definite'
    return
    end if
    if(.not.is_positive_definite(v,eps)) then
    message='V is not positive definite'
    return
    end if
    ok=.true.
  end subroutine check_matnorm

  real(dp) function dmatnorm(a,m,u,v,log_density,ok,message,tol) result(value)
    real(dp), intent(in) :: a(:,:),m(:,:),u(:,:),v(:,:)
    logical, intent(in), optional :: log_density
    logical, intent(out), optional :: ok
    character(len=*), intent(out), optional :: message
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: e(:,:),uinve(:,:),vinvet(:,:),lu(:,:),lv(:,:)
    real(dp) :: q,ldu,ldv,logd
    logical :: valid,solve_ok,llog
    character(len=256) :: msg
    integer :: n,p
    llog=.true.
    if(present(log_density)) llog=log_density
    valid=.false.
    msg=''
    if(any(shape(a)/=shape(m))) then
      msg='A and M must have identical dimensions'
      goto 900
    end if
    call check_matnorm(m,u,v,valid,msg,tol)
    if(.not.valid) goto 900
    n=size(a,1)
    p=size(a,2)
    e=a-m
    call solve_spd(u,e,uinve,solve_ok,msg)
    if(.not.solve_ok) then
    valid=.false.
    goto 900
    end if
    call solve_spd(v,transpose(e),vinvet,solve_ok,msg)
    if(.not.solve_ok) then
    valid=.false.
    goto 900
    end if
    call cholesky_lower(u,lu,solve_ok,msg)
    if(.not.solve_ok) then
    valid=.false.
    goto 900
    end if
    call cholesky_lower(v,lv,solve_ok,msg)
    if(.not.solve_ok) then
    valid=.false.
    goto 900
    end if
    q=sum(uinve*transpose(vinvet))
    ldu=logdet_cholesky(lu)
    ldv=logdet_cholesky(lv)
    logd=-0.5_dp*real(n*p,dp)*log_two_pi-0.5_dp*real(p,dp)*ldu &
      -0.5_dp*real(n,dp)*ldv-0.5_dp*q
    if(llog) then
    value=logd
    else
    value=exp(logd)
    end if
    valid=.true.
900 continue
    if(.not.valid) value=merge(-huge(1.0_dp),0.0_dp,llog)
    if(present(ok)) ok=valid
    if(present(message)) message=trim(msg)
  end function dmatnorm

  function pmatnorm_bounds(lower,upper,m,u,v,control,legacy_covariance_order) result(res)
    real(dp), intent(in) :: lower(:,:),upper(:,:),m(:,:),u(:,:),v(:,:)
    type(probability_control), intent(in), optional :: control
    logical, intent(in), optional :: legacy_covariance_order
    type(probability_result) :: res
    real(dp), allocatable :: sigma(:,:)
    logical :: ok,legacy
    character(len=256) :: msg
    if(any(shape(lower)/=shape(m)) .or. any(shape(upper)/=shape(m))) then
      res%inform=2
      res%message='bounds and M must have identical dimensions'
      return
    end if
    call check_matnorm(m,u,v,ok,msg)
    if(.not.ok) then
    res%inform=3
    res%message=msg
    return
    end if
    legacy=.false.
    if(present(legacy_covariance_order)) legacy=legacy_covariance_order
    if(legacy) then
      sigma=kronecker_product(u,v)
    else
      sigma=kronecker_product(v,u)
    end if
    res=pmvnorm(vec(lower),vec(upper),vec(m),sigma,control)
  end function pmatnorm_bounds

  function pmatnorm_upper(upper,m,u,v,control,legacy_covariance_order) result(res)
    real(dp), intent(in) :: upper(:,:),m(:,:),u(:,:),v(:,:)
    type(probability_control), intent(in), optional :: control
    logical, intent(in), optional :: legacy_covariance_order
    type(probability_result) :: res
    real(dp), allocatable :: lower(:,:)
    allocate(lower(size(m,1),size(m,2)))
    lower=ieee_value(0.0_dp,ieee_negative_inf)
    res=pmatnorm_bounds(lower,upper,m,u,v,control,legacy_covariance_order)
  end function pmatnorm_upper

  function pmatnorm_full(m,u,v,control,legacy_covariance_order) result(res)
    real(dp), intent(in) :: m(:,:),u(:,:),v(:,:)
    type(probability_control), intent(in), optional :: control
    logical, intent(in), optional :: legacy_covariance_order
    type(probability_result) :: res
    real(dp), allocatable :: lower(:,:),upper(:,:)
    allocate(lower(size(m,1),size(m,2)),upper(size(m,1),size(m,2)))
    lower=ieee_value(0.0_dp,ieee_negative_inf)
    upper=ieee_value(0.0_dp,ieee_positive_inf)
    res=pmatnorm_bounds(lower,upper,m,u,v,control,legacy_covariance_order)
  end function pmatnorm_full

  function rmatnorm_one(m,u,v,seed) result(x)
    real(dp), intent(in) :: m(:,:),u(:,:),v(:,:)
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: sigma(:,:),draws(:,:)
    logical :: ok
    character(len=256) :: msg
    call check_matnorm(m,u,v,ok,msg)
    allocate(x(size(m,1),size(m,2)))
    x=0.0_dp
    if(.not.ok) return
    sigma=kronecker_product(v,u)
    allocate(draws(1,size(m)))
    draws=0.0_dp
    if(present(seed)) then
      draws=rmvnorm(1,vec(m),sigma,seed)
    else
      draws=rmvnorm(1,vec(m),sigma)
    end if
    x=reshape(draws(1,:),shape(m))
  end function rmatnorm_one

  function rmatnorm_many(s,m,u,v,seed) result(x)
    integer, intent(in) :: s
    real(dp), intent(in) :: m(:,:),u(:,:),v(:,:)
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:,:,:)
    real(dp), allocatable :: sigma(:,:),draws(:,:)
    logical :: ok
    character(len=256) :: msg
    integer :: k
    allocate(x(size(m,1),size(m,2),max(0,s)))
    x=0.0_dp
    if(s<=0) return
    call check_matnorm(m,u,v,ok,msg)
    if(.not.ok) return
    sigma=kronecker_product(v,u)
    allocate(draws(s,size(m)))
    draws=0.0_dp
    if(present(seed)) then
      draws=rmvnorm(s,vec(m),sigma,seed)
    else
      draws=rmvnorm(s,vec(m),sigma)
    end if
    do k=1,s
      x(:,:,k)=reshape(draws(k,:),shape(m))
    end do
  end function rmatnorm_many

end module matrix_normal_distribution
