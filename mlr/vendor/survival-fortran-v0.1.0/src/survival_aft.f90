! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_aft
  use survival_kinds, only : dp, pi
  use survival_types, only : aft_result
  use survival_linalg, only : solve_sym, invert_matrix
  implicit none
  private
  public :: survreg_fit, survreg_loglik
contains

  subroutine survreg_fit(time,status,x,dist,result,weights,maxiter,eps)
    real(dp), intent(in) :: time(:), x(:,:)
    integer, intent(in) :: status(:)
    character(len=*), intent(in) :: dist
    type(aft_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: eps

    integer :: n,p,q,iter,max_it
    real(dp) :: tolerance,ll,newll,alpha
    real(dp), allocatable :: theta(:),grad(:),hess(:,:),step(:),w(:),vinv(:,:)
    real(dp), allocatable :: ty(:)
    logical :: ok,fixed_scale

    n=size(time)
    p=size(x,2)
    fixed_scale=(trim(dist)=='exponential' .or. trim(dist)=='rayleigh')
    if(fixed_scale) then
      q=p
    else
      q=p+1
    end if
    max_it=30
    if(present(maxiter)) max_it=maxiter
    tolerance=1.0e-8_dp
    if(present(eps)) tolerance=eps
    allocate(theta(q),grad(q),hess(q,q),step(q),w(n),vinv(q,q),ty(n))
    if(present(weights)) then
      w=weights
    else
      w=1.0_dp
    end if
    ty=transform_time(time,dist)
    theta=0.0_dp
    if(p>=1) theta(1)=sum(w*ty)/sum(w)
    if(.not.fixed_scale) then
      theta(q)=log(max(0.1_dp,sqrt(sum(w*(ty-sum(w*ty)/sum(w))**2)/sum(w))))
    end if

    ll=survreg_loglik(time,status,x,theta,dist,w)
    result%converged=.false.
    iter=0
    do iter=1,max_it
      call numeric_grad_hess(time,status,x,theta,dist,w,grad,hess)
      call solve_sym(-hess,grad,step,ok)
      if(.not.ok) exit
      alpha=1.0_dp
      do
        newll=survreg_loglik(time,status,x,theta+alpha*step,dist,w)
        if(newll>=ll .or. alpha<1.0e-7_dp) exit
        alpha=alpha/2.0_dp
      end do
      theta=theta+alpha*step
      if(abs(newll-ll)<=tolerance*(1.0_dp+abs(newll))) then
        result%converged=.true.
        ll=newll
        exit
      end if
      ll=newll
    end do

    call numeric_grad_hess(time,status,x,theta,dist,w,grad,hess)
    call invert_matrix(-hess,vinv,ok)
    if(.not.ok) vinv=0.0_dp
    allocate(result%coef(p),result%var(q,q))
    result%coef=theta(1:p)
    result%var=vinv
    result%loglik=ll
    result%iterations=min(iter,max_it)
    if(fixed_scale) then
      if(trim(dist)=='rayleigh') then
        result%scale=0.5_dp
      else
        result%scale=1.0_dp
      end if
    else
      result%scale=exp(theta(q))
    end if
  end subroutine survreg_fit

  real(dp) function survreg_loglik(time,status,x,theta,dist,weights) result(ll)
    real(dp), intent(in) :: time(:),x(:,:),theta(:),weights(:)
    integer, intent(in) :: status(:)
    character(len=*), intent(in) :: dist
    integer :: i,p
    real(dp) :: scale,z,lp,pdf,surv,jac
    logical :: fixed_scale

    p=size(x,2)
    fixed_scale=(trim(dist)=='exponential' .or. trim(dist)=='rayleigh')
    if(fixed_scale) then
      if(trim(dist)=='rayleigh') then
        scale=0.5_dp
      else
        scale=1.0_dp
      end if
    else
      scale=exp(theta(p+1))
    end if
    ll=0.0_dp
    do i=1,size(time)
      lp=dot_product(x(i,:),theta(1:p))
      z=(transform_scalar(time(i),dist)-lp)/scale
      call base_density(z,dist,pdf,surv)
      jac=transform_jacobian(time(i),dist)
      if(status(i)/=0) then
        ll=ll+weights(i)*log(max(tiny(1.0_dp),pdf*jac/scale))
      else
        ll=ll+weights(i)*log(max(tiny(1.0_dp),surv))
      end if
    end do
  end function survreg_loglik

  subroutine base_density(z,dist,pdf,surv)
    real(dp), intent(in) :: z
    character(len=*), intent(in) :: dist
    real(dp), intent(out) :: pdf,surv
    real(dp) :: ez
    select case(trim(base_dist(dist)))
    case('extreme')
      ez=exp(min(z,700.0_dp))
      surv=exp(-ez)
      pdf=ez*surv
    case('logistic')
      if(z>=0.0_dp) then
        ez=exp(-z)
        surv=ez/(1.0_dp+ez)
        pdf=ez/(1.0_dp+ez)**2
      else
        ez=exp(z)
        surv=1.0_dp/(1.0_dp+ez)
        pdf=ez/(1.0_dp+ez)**2
      end if
    case default
      pdf=exp(-0.5_dp*z*z)/sqrt(2.0_dp*pi)
      surv=0.5_dp*erfc(z/sqrt(2.0_dp))
    end select
  end subroutine base_density

  pure function base_dist(dist) result(base)
    character(len=*), intent(in) :: dist
    character(len=12) :: base
    select case(trim(dist))
    case('weibull','exponential','rayleigh','extreme')
      base='extreme'
    case('loglogistic','logistic')
      base='logistic'
    case default
      base='gaussian'
    end select
  end function base_dist

  pure real(dp) function transform_scalar(y,dist) result(value)
    real(dp), intent(in) :: y
    character(len=*), intent(in) :: dist
    if(is_log_time_distribution(dist)) then
      value=log(y)
    else
      value=y
    end if
  end function transform_scalar

  pure real(dp) function transform_jacobian(y,dist) result(value)
    real(dp), intent(in) :: y
    character(len=*), intent(in) :: dist
    if(is_log_time_distribution(dist)) then
      value=1.0_dp/y
    else
      value=1.0_dp
    end if
  end function transform_jacobian

  pure function transform_time(y,dist) result(value)
    real(dp), intent(in) :: y(:)
    character(len=*), intent(in) :: dist
    real(dp) :: value(size(y))
    integer :: i
    do i=1,size(y)
      value(i)=transform_scalar(y(i),dist)
    end do
  end function transform_time

  pure logical function is_log_time_distribution(dist) result(answer)
    character(len=*), intent(in) :: dist
    answer = trim(dist)=='weibull' .or. trim(dist)=='exponential' .or. &
             trim(dist)=='rayleigh' .or. trim(dist)=='lognormal' .or. &
             trim(dist)=='loggaussian' .or. trim(dist)=='loglogistic'
  end function is_log_time_distribution

  subroutine numeric_grad_hess(time,status,x,theta,dist,w,grad,hess)
    real(dp), intent(in) :: time(:),x(:,:),theta(:),w(:)
    integer, intent(in) :: status(:)
    character(len=*), intent(in) :: dist
    real(dp), intent(out) :: grad(:),hess(:,:)
    real(dp) :: hi,hj,f0,fpp,fpm,fmp,fmm
    real(dp), allocatable :: trial(:)
    integer :: i,j,q

    q=size(theta)
    allocate(trial(q))
    f0=survreg_loglik(time,status,x,theta,dist,w)
    do i=1,q
      hi=1.0e-5_dp*max(1.0_dp,abs(theta(i)))
      trial=theta
      trial(i)=trial(i)+hi
      fpp=survreg_loglik(time,status,x,trial,dist,w)
      trial(i)=theta(i)-hi
      fpm=survreg_loglik(time,status,x,trial,dist,w)
      grad(i)=(fpp-fpm)/(2.0_dp*hi)
      hess(i,i)=(fpp-2.0_dp*f0+fpm)/(hi*hi)
    end do
    do i=1,q
      hi=1.0e-5_dp*max(1.0_dp,abs(theta(i)))
      do j=i+1,q
        hj=1.0e-5_dp*max(1.0_dp,abs(theta(j)))
        trial=theta
        trial(i)=trial(i)+hi; trial(j)=trial(j)+hj
        fpp=survreg_loglik(time,status,x,trial,dist,w)
        trial=theta
        trial(i)=trial(i)+hi; trial(j)=trial(j)-hj
        fpm=survreg_loglik(time,status,x,trial,dist,w)
        trial=theta
        trial(i)=trial(i)-hi; trial(j)=trial(j)+hj
        fmp=survreg_loglik(time,status,x,trial,dist,w)
        trial=theta
        trial(i)=trial(i)-hi; trial(j)=trial(j)-hj
        fmm=survreg_loglik(time,status,x,trial,dist,w)
        hess(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numeric_grad_hess
end module survival_aft
