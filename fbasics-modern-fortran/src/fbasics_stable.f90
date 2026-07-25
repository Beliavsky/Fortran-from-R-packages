! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_stable
  use fbasics_kinds, only: dp, pi, clamp
  use fbasics_rng, only: runif_lcg, rnorm_lcg
  use fbasics_special, only: normal_pdf, normal_cdf, normal_quantile
  use fbasics_stats, only: sample_quantile, sample_sd
  use fbasics_optimize, only: nelder_mead_bounded, numerical_hessian
  use fbasics_linalg, only: matrix_inverse
  implicit none
  private
  type, public :: stable_fit_result
    real(dp) :: alpha = 2.0_dp
    real(dp) :: beta = 0.0_dp
    real(dp) :: gamma = 1.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    character(len=12) :: method = 'ecf'
    logical :: converged = .false.
    real(dp), allocatable :: hessian(:,:), covariance(:,:)
  end type stable_fit_result
  public :: dstable_s1, pstable_s1, qstable_s1, rstable_s1
  public :: stable_characteristic, fit_stable, fit_stable_ecf, fit_stable_mle
  real(dp), allocatable, save :: fit_data(:), fit_freq(:)
contains
  pure complex(dp) function stable_characteristic(t,alpha,beta,gamma,delta) result(phi)
    real(dp), intent(in) :: t,alpha,beta,gamma,delta
    real(dp) :: at, real_part, imag_part, tangent
    if (gamma<=0.0_dp .or. alpha<=0.0_dp .or. alpha>2.0_dp .or. abs(beta)>1.0_dp) then
      phi=(0.0_dp,0.0_dp)
      return
    end if
    if (t==0.0_dp) then
      phi=(1.0_dp,0.0_dp)
      return
    end if
    at=abs(gamma*t)
    if (abs(alpha-1.0_dp)>1.0e-10_dp) then
      tangent=tan(0.5_dp*pi*alpha)
      real_part=-at**alpha
      imag_part=delta*t+beta*sign(1.0_dp,t)*at**alpha*tangent
    else
      real_part=-at
      imag_part=delta*t-beta*(2.0_dp/pi)*gamma*t*log(abs(t))
    end if
    phi=exp(cmplx(real_part,imag_part,kind=dp))
  end function stable_characteristic

  real(dp) function dstable_s1(x,alpha,beta,gamma,delta,log_density) result(v)
    real(dp), intent(in) :: x,alpha,beta,gamma,delta
    logical, intent(in), optional :: log_density
    logical :: lg
    integer :: i,n
    real(dp) :: h,t,tmax,sumv,term
    complex(dp) :: phi
    lg=.false.; if (present(log_density)) lg=log_density
    if (gamma<=0.0_dp .or. alpha<=0.0_dp .or. alpha>2.0_dp .or. abs(beta)>1.0_dp) then
      v=merge(-huge(1.0_dp),0.0_dp,lg)
      return
    end if
    if (abs(alpha-2.0_dp)<1.0e-12_dp) then
      v=normal_pdf((x-delta)/(sqrt(2.0_dp)*gamma))/(sqrt(2.0_dp)*gamma)
    else if (abs(alpha-1.0_dp)<1.0e-12_dp .and. abs(beta)<1.0e-12_dp) then
      v=gamma/(pi*((x-delta)**2+gamma*gamma))
    else
      tmax=(-log(1.0e-13_dp))**(1.0_dp/alpha)/gamma
      n=max(800,2*int(90.0_dp*tmax/pi)+2)
      if (mod(n,2)/=0) n=n+1
      n=min(n,12000)
      h=tmax/real(n,dp)
      sumv=0.0_dp
      do i=0,n
        t=h*real(i,dp)
        phi=stable_characteristic(t,alpha,beta,gamma,delta)
        term=real(phi*exp(cmplx(0.0_dp,-t*x,kind=dp)),dp)
        if (i==0 .or. i==n) then
          sumv=sumv+term
        else if (mod(i,2)==0) then
          sumv=sumv+2.0_dp*term
        else
          sumv=sumv+4.0_dp*term
        end if
      end do
      v=max(0.0_dp,h*sumv/(3.0_dp*pi))
    end if
    if (lg) then
      if (v>0.0_dp) then
        v=log(v)
      else
        v=-huge(1.0_dp)
      end if
    end if
  end function dstable_s1

  real(dp) function pstable_s1(x,alpha,beta,gamma,delta) result(v)
    real(dp), intent(in) :: x,alpha,beta,gamma,delta
    integer :: i,n
    real(dp) :: h,t,tmax,sumv,integrand
    complex(dp) :: phi,z
    if (gamma<=0.0_dp .or. alpha<=0.0_dp .or. alpha>2.0_dp .or. abs(beta)>1.0_dp) then
      v=0.0_dp
      return
    end if
    if (abs(alpha-2.0_dp)<1.0e-12_dp) then
      v=normal_cdf((x-delta)/(sqrt(2.0_dp)*gamma))
      return
    else if (abs(alpha-1.0_dp)<1.0e-12_dp .and. abs(beta)<1.0e-12_dp) then
      v=0.5_dp+atan((x-delta)/gamma)/pi
      return
    end if
    tmax=(-log(1.0e-13_dp))**(1.0_dp/alpha)/gamma
    n=max(1200,int(120.0_dp*tmax/pi)+1)
    n=min(n,20000)
    h=tmax/real(n,dp)
    sumv=0.0_dp
    do i=1,n
      t=(real(i,dp)-0.5_dp)*h
      phi=stable_characteristic(t,alpha,beta,gamma,delta)
      z=phi*exp(cmplx(0.0_dp,-t*x,kind=dp))
      integrand=aimag(z)/t
      sumv=sumv+integrand
    end do
    v=clamp(0.5_dp-h*sumv/pi,0.0_dp,1.0_dp)
  end function pstable_s1

  real(dp) function qstable_s1(p,alpha,beta,gamma,delta) result(v)
    real(dp), intent(in) :: p,alpha,beta,gamma,delta
    real(dp) :: lo,hi,mid,width
    integer :: i
    if (p<=0.0_dp) then
      v=-huge(1.0_dp); return
    else if (p>=1.0_dp) then
      v=huge(1.0_dp); return
    end if
    if (abs(alpha-2.0_dp)<1.0e-12_dp) then
      v=delta+sqrt(2.0_dp)*gamma*normal_quantile(p)
      return
    else if (abs(alpha-1.0_dp)<1.0e-12_dp .and. abs(beta)<1.0e-12_dp) then
      v=delta+gamma*tan(pi*(p-0.5_dp))
      return
    end if
    width=4.0_dp*gamma
    lo=delta-width; hi=delta+width
    do while (pstable_s1(lo,alpha,beta,gamma,delta)>p)
      width=2.0_dp*width; lo=delta-width
      if (width>1.0e8_dp*gamma) exit
    end do
    do while (pstable_s1(hi,alpha,beta,gamma,delta)<p)
      width=2.0_dp*width; hi=delta+width
      if (width>1.0e8_dp*gamma) exit
    end do
    do i=1,70
      mid=0.5_dp*(lo+hi)
      if (pstable_s1(mid,alpha,beta,gamma,delta)<p) then
        lo=mid
      else
        hi=mid
      end if
    end do
    v=0.5_dp*(lo+hi)
  end function qstable_s1

  real(dp) function rstable_s1(alpha,beta,gamma,delta) result(v)
    real(dp), intent(in) :: alpha,beta,gamma,delta
    real(dp) :: u,w,bangle,scale,x,tangent,den
    u=pi*(runif_lcg()-0.5_dp)
    w=-log(max(runif_lcg(),tiny(1.0_dp)))
    if (abs(alpha-2.0_dp)<1.0e-12_dp) then
      v=delta+sqrt(2.0_dp)*gamma*rnorm_lcg()
      return
    end if
    if (abs(alpha-1.0_dp)>1.0e-10_dp) then
      tangent=tan(0.5_dp*pi*alpha)
      bangle=atan(beta*tangent)/alpha
      scale=(1.0_dp+(beta*tangent)**2)**(0.5_dp/alpha)
      den=max(cos(u),tiny(1.0_dp))
      x=scale*sin(alpha*(u+bangle))/den**(1.0_dp/alpha) * &
        (cos(u-alpha*(u+bangle))/w)**((1.0_dp-alpha)/alpha)
      v=delta+gamma*x
    else
      den=0.5_dp*pi+beta*u
      x=(2.0_dp/pi)*(den*tan(u)-beta*log((0.5_dp*pi*w*cos(u))/den))
      v=delta+gamma*x+(2.0_dp/pi)*beta*gamma*log(gamma)
    end if
  end function rstable_s1

  subroutine fit_stable_ecf(x,fit,max_iter)
    real(dp), intent(in) :: x(:)
    type(stable_fit_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: start(:),lower(:),upper(:),best(:),h(:,:),cov(:,:)
    real(dp) :: fbest,sdx,q25,q75,scale0
    logical :: conv
    integer :: info
    integer :: imax
    if (size(x)<8) then
      fit%converged=.false.; return
    end if
    if (allocated(fit_data)) deallocate(fit_data)
    if (allocated(fit_freq)) deallocate(fit_freq)
    allocate(fit_data(size(x)),fit_freq(6))
    fit_data=x
    sdx=max(sample_sd(x),1.0e-6_dp)
    q25=sample_quantile(x,0.25_dp); q75=sample_quantile(x,0.75_dp)
    scale0=max((q75-q25)/1.35_dp,0.25_dp*sdx)
    fit_freq=[0.20_dp,0.35_dp,0.55_dp,0.80_dp,1.10_dp,1.45_dp]/scale0
    allocate(start(4),lower(4),upper(4))
    start=[1.7_dp,0.0_dp,scale0,sample_quantile(x,0.5_dp)]
    lower=[0.35_dp,-0.99_dp,1.0e-5_dp,minval(x)-5.0_dp*sdx]
    upper=[2.0_dp,0.99_dp,20.0_dp*sdx,maxval(x)+5.0_dp*sdx]
    imax=1200; if (present(max_iter)) imax=max_iter
    call nelder_mead_bounded(stable_ecf_objective,start,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-9_dp)
    fit%alpha=best(1); fit%beta=best(2); fit%gamma=best(3); fit%delta=best(4)
    fit%objective=fbest; fit%converged=conv; fit%method='ecf'
    call numerical_hessian(stable_ecf_objective,best,h)
    call matrix_inverse(h,cov,info)
    fit%hessian=h
    if (info==0) fit%covariance=cov
  end subroutine fit_stable_ecf


  subroutine fit_stable(x,method,fit,max_iter)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: method
    type(stable_fit_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    select case(trim(adjustl(method)))
    case('mle','MLE')
      call fit_stable_mle(x,fit,max_iter)
    case default
      call fit_stable_ecf(x,fit,max_iter)
    end select
  end subroutine fit_stable

  subroutine fit_stable_mle(x,fit,max_iter)
    real(dp), intent(in) :: x(:)
    type(stable_fit_result), intent(out) :: fit
    integer, intent(in), optional :: max_iter
    type(stable_fit_result) :: initial
    real(dp), allocatable :: start(:),lower(:),upper(:),best(:),h(:,:),cov(:,:)
    real(dp) :: fbest,sdx
    logical :: conv
    integer :: info,imax
    if(size(x)<8)then
      fit%converged=.false.
      fit%method='mle'
      return
    end if
    imax=600; if(present(max_iter))imax=max_iter
    call fit_stable_ecf(x,initial,max_iter=min(240,max(40,imax/2)))
    if(allocated(fit_data))deallocate(fit_data)
    allocate(fit_data(size(x)));fit_data=x
    sdx=max(sample_sd(x),1.0e-6_dp)
    allocate(start(4),lower(4),upper(4))
    start=[initial%alpha,initial%beta,max(initial%gamma,1.0e-5_dp),initial%delta]
    lower=[0.35_dp,-0.99_dp,1.0e-5_dp,minval(x)-5.0_dp*sdx]
    upper=[2.0_dp,0.99_dp,20.0_dp*sdx,maxval(x)+5.0_dp*sdx]
    call nelder_mead_bounded(stable_mle_objective,start,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-8_dp)
    fit%alpha=best(1);fit%beta=best(2);fit%gamma=best(3);fit%delta=best(4)
    fit%objective=fbest;fit%loglik=-fbest;fit%aic=2.0_dp*4.0_dp-2.0_dp*fit%loglik
    fit%bic=log(real(size(x),dp))*4.0_dp-2.0_dp*fit%loglik
    fit%converged=conv;fit%method='mle'
    call numerical_hessian(stable_mle_objective,best,h)
    call matrix_inverse(h,cov,info)
    fit%hessian=h
    if(info==0)fit%covariance=cov
  end subroutine fit_stable_mle

  real(dp) function stable_mle_objective(p) result(v)
    real(dp), intent(in) :: p(:)
    real(dp) :: ld
    integer :: i
    if(size(p)/=4 .or. p(1)<=0.0_dp .or. p(1)>2.0_dp .or. abs(p(2))>=1.0_dp .or. p(3)<=0.0_dp)then
      v=huge(1.0_dp);return
    end if
    v=0.0_dp
    do i=1,size(fit_data)
      ld=dstable_s1(fit_data(i),p(1),p(2),p(3),p(4),.true.)
      if(.not.(ld>-huge(1.0_dp)/2.0_dp) .or. ld/=ld)then
        v=huge(1.0_dp);return
      end if
      v=v-ld
    end do
  end function stable_mle_objective

  real(dp) function stable_ecf_objective(p) result(v)
    real(dp), intent(in) :: p(:)
    integer :: i,j,n
    complex(dp) :: empirical,model,z
    real(dp) :: t
    if (size(p)/=4 .or. p(1)<=0.0_dp .or. p(1)>2.0_dp .or. abs(p(2))>=1.0_dp .or. p(3)<=0.0_dp) then
      v=huge(1.0_dp); return
    end if
    n=size(fit_data); v=0.0_dp
    do j=1,size(fit_freq)
      t=fit_freq(j); empirical=(0.0_dp,0.0_dp)
      do i=1,n
        z=exp(cmplx(0.0_dp,t*fit_data(i),kind=dp))
        empirical=empirical+z
      end do
      empirical=empirical/real(n,dp)
      model=stable_characteristic(t,p(1),p(2),p(3),p(4))
      v=v+abs(empirical-model)**2/(1.0_dp+t*t)
    end do
  end function stable_ecf_objective
end module fbasics_stable
