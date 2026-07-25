! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_fit
  use fextremes_kinds, only: dp, euler_gamma, huge_penalty, pi
  use fextremes_stats, only: mean_value, variance_value, sort_real, nan_value
  use fextremes_distributions, only: gev_logpdf, gpd_logpdf
  use fextremes_optimize, only: nelder_mead_bounded, numerical_hessian, invert_matrix
  implicit none
  private
  public :: gev_fit_result, gpd_fit_result, fit_gev, fit_gumbel, fit_gpd
  public :: gev_nll, gpd_nll, gev_pwm, gumbel_pwm, gpd_pwm

  type :: gev_fit_result
    real(dp) :: xi=0.0_dp, mu=0.0_dp, beta=1.0_dp
    real(dp) :: se(3)=0.0_dp, covariance(3,3)=0.0_dp
    real(dp) :: nll=0.0_dp
    logical :: converged=.false., gumbel=.false.
    integer :: iterations=0, evaluations=0
    character(len=8) :: method='mle'
    real(dp), allocatable :: residuals(:)
  end type gev_fit_result

  type :: gpd_fit_result
    real(dp) :: xi=0.0_dp, beta=1.0_dp, threshold=0.0_dp, exceedance_probability=0.0_dp
    real(dp) :: se(2)=0.0_dp, covariance(2,2)=0.0_dp
    real(dp) :: nll=0.0_dp
    logical :: converged=.false.
    integer :: iterations=0, evaluations=0, n_exceed=0
    character(len=8) :: method='mle'
    real(dp), allocatable :: residuals(:), exceedances(:)
  end type gpd_fit_result

  real(dp), allocatable, save :: fit_data(:), fit_excess(:)
contains
  real(dp) function gev_nll(x, xi, mu, beta) result(v)
    real(dp), intent(in) :: x(:),xi,mu,beta
    integer :: i
    v=0.0_dp
    if(beta<=0.0_dp) then; v=huge_penalty; return; end if
    do i=1,size(x)
      if(gev_logpdf(x(i),xi,mu,beta)<=-0.5_dp*huge(1.0_dp)) then; v=huge_penalty; return; end if
      v=v-gev_logpdf(x(i),xi,mu,beta)
    end do
  end function gev_nll

  real(dp) function gpd_nll(excess, xi, beta) result(v)
    real(dp), intent(in) :: excess(:),xi,beta
    integer :: i
    v=0.0_dp
    if(beta<=0.0_dp) then; v=huge_penalty; return; end if
    do i=1,size(excess)
      if(gpd_logpdf(excess(i),xi,0.0_dp,beta)<=-0.5_dp*huge(1.0_dp)) then; v=huge_penalty; return; end if
      v=v-gpd_logpdf(excess(i),xi,0.0_dp,beta)
    end do
  end function gpd_nll

  subroutine gumbel_pwm(x,mu,beta,ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: mu,beta
    logical, intent(out) :: ok
    real(dp), allocatable :: s(:)
    real(dp) :: b0,b1
    integer :: i,n
    n=size(x); ok=.false.; mu=nan_value(); beta=nan_value()
    if(n<3) return
    allocate(s(n)); s=x; call sort_real(s)
    b0=mean_value(s); b1=0.0_dp
    do i=1,n
      b1=b1+real(i-1,dp)/real(n-1,dp)*s(i)/real(n,dp)
    end do
    beta=(2.0_dp*b1-b0)/log(2.0_dp)
    mu=b0-euler_gamma*beta
    ok=beta>0.0_dp
  end subroutine gumbel_pwm

  subroutine gev_pwm(x,xi,mu,beta,ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: xi,mu,beta
    logical, intent(out) :: ok
    real(dp), allocatable :: s(:)
    real(dp) :: b0,b1,b2,target,lo,hi,mid,flo,fhi,fmid,den
    integer :: i,n,it
    n=size(x); ok=.false.; xi=0.0_dp; mu=0.0_dp; beta=1.0_dp
    if(n<5) return
    allocate(s(n)); s=x; call sort_real(s)
    b0=mean_value(s); b1=0.0_dp; b2=0.0_dp
    do i=1,n
      b1=b1+real(i-1,dp)/real(n-1,dp)*s(i)/real(n,dp)
      b2=b2+real((i-1)*(i-2),dp)/real((n-1)*(n-2),dp)*s(i)/real(n,dp)
    end do
    den=2.0_dp*b1-b0
    if(abs(den)<1.0e-12_dp) return
    target=(3.0_dp*b2-b0)/den
    lo=-0.95_dp; hi=0.95_dp; flo=pwm_eq(lo,target); fhi=pwm_eq(hi,target)
    if(flo*fhi>0.0_dp) then
      xi=0.0_dp; beta=den/log(2.0_dp); mu=b0-euler_gamma*beta; ok=beta>0.0_dp; return
    end if
    do it=1,100
      mid=0.5_dp*(lo+hi); fmid=pwm_eq(mid,target)
      if(abs(fmid)<1.0e-10_dp) exit
      if(flo*fmid<=0.0_dp) then; hi=mid; fhi=fmid; else; lo=mid; flo=fmid; end if
    end do
    xi=mid
    if(abs(xi)<1.0e-7_dp) then
      beta=den/log(2.0_dp); mu=b0-euler_gamma*beta
    else
      beta=den*xi/(gamma(1.0_dp-xi)*(2.0_dp**xi-1.0_dp))
      mu=b0+beta*(1.0_dp-gamma(1.0_dp-xi))/xi
    end if
    ok=beta>0.0_dp
  contains
    pure real(dp) function pwm_eq(a,t) result(f)
      real(dp),intent(in)::a,t
      if(abs(a)<1.0e-8_dp) then
        f=log(3.0_dp)/log(2.0_dp)-t
      else
        f=(3.0_dp**a-1.0_dp)/(2.0_dp**a-1.0_dp)-t
      end if
    end function pwm_eq
  end subroutine gev_pwm

  subroutine gpd_pwm(x,threshold,xi,beta,varcov,ok)
    real(dp),intent(in)::x(:),threshold
    real(dp),intent(out)::xi,beta,varcov(2,2)
    logical,intent(out)::ok
    real(dp),allocatable::e(:)
    real(dp)::a0,a1,denom,one,two,cov
    integer::i,j,n
    n=count(x>threshold); ok=.false.; xi=nan_value(); beta=nan_value(); varcov=nan_value()
    if(n<3) return
    allocate(e(n)); j=0
    do i=1,size(x); if(x(i)>threshold) then; j=j+1; e(j)=x(i)-threshold; end if; end do
    call sort_real(e)
    a0=mean_value(e); a1=0.0_dp
    do i=1,n; a1=a1+e(i)*(1.0_dp-(real(i,dp)-0.35_dp)/real(n,dp))/real(n,dp); end do
    if(abs(a0-2.0_dp*a1)<1.0e-12_dp) return
    xi=2.0_dp-a0/(a0-2.0_dp*a1); beta=2.0_dp*a0*a1/(a0-2.0_dp*a1)
    if(beta<=0.0_dp) return
    varcov=nan_value()
    if(xi<0.5_dp) then
      denom=real(n,dp)*(1.0_dp-2.0_dp*xi)*(3.0_dp-2.0_dp*xi)
      one=(1.0_dp-xi)*(1.0_dp-xi+2.0_dp*xi*xi)*(2.0_dp-xi)**2
      two=(7.0_dp-18.0_dp*xi+11.0_dp*xi*xi-2.0_dp*xi**3)*beta*beta
      cov=beta*(2.0_dp-xi)*(2.0_dp-6.0_dp*xi+7.0_dp*xi*xi-2.0_dp*xi**3)
      varcov=reshape([one,cov,cov,two],[2,2])/denom
    end if
    ok=.true.
  end subroutine gpd_pwm

  subroutine fit_gev(x,result,method)
    real(dp),intent(in)::x(:)
    type(gev_fit_result),intent(out)::result
    character(len=*),intent(in),optional::method
    character(len=8)::meth
    real(dp)::xi0,mu0,beta0,start(3),lo(3),hi(3),best(3),fbest,h(3,3),covt(3,3),jac(3,3)
    logical::ok,invok
    integer::i
    meth='mle'; if(present(method)) meth=adjustl(method)
    call gev_pwm(x,xi0,mu0,beta0,ok)
    if (.not. ok) then
      beta0 = sqrt(max(variance_value(x),1.0e-8_dp))*sqrt(6.0_dp)/pi
      mu0 = mean_value(x)-euler_gamma*beta0
      xi0 = 0.1_dp
    end if
    result%method=meth; result%gumbel=.false.; result%xi=xi0; result%mu=mu0; result%beta=beta0
    if(trim(meth)=='pwm') then
      result%converged=ok; result%nll=gev_nll(x,xi0,mu0,beta0); result%se=nan_value(); result%covariance=nan_value()
    else
      fit_data=x; start=[xi0,mu0,log(max(beta0,1.0e-6_dp))]
      lo=[-2.0_dp,minval(x)-5.0_dp*max(1.0_dp,sqrt(max(variance_value(x),0.0_dp))),log(1.0e-8_dp)]
      hi = [2.0_dp, &
        maxval(x)+5.0_dp*max(1.0_dp,sqrt(max(variance_value(x),0.0_dp))), &
        log(max(1.0_dp,100.0_dp*(maxval(x)-minval(x)+1.0_dp)))]
      call nelder_mead_bounded(gev_obj,start,lo,hi,best,fbest,result%converged,result%iterations,result%evaluations)
      result%xi=best(1); result%mu=best(2); result%beta=exp(best(3)); result%nll=fbest
      call numerical_hessian(gev_obj,best,h); call invert_matrix(h,covt,invok)
      if(invok) then
        jac=0.0_dp; jac(1,1)=1.0_dp; jac(2,2)=1.0_dp; jac(3,3)=result%beta
        result%covariance=matmul(jac,matmul(covt,transpose(jac)))
        do i=1,3; result%se(i)=sqrt(max(0.0_dp,result%covariance(i,i))); end do
      else; result%covariance=nan_value(); result%se=nan_value(); end if
      if(allocated(fit_data)) deallocate(fit_data)
    end if
    allocate(result%residuals(size(x)))
    if(abs(result%xi)<1.0e-10_dp) then
      result%residuals=exp(-(x-result%mu)/result%beta)
    else
      result%residuals=(1.0_dp+result%xi*(x-result%mu)/result%beta)**(-1.0_dp/result%xi)
    end if
  end subroutine fit_gev

  subroutine fit_gumbel(x,result,method)
    real(dp),intent(in)::x(:)
    type(gev_fit_result),intent(out)::result
    character(len=*),intent(in),optional::method
    character(len=8)::meth
    real(dp)::mu0,beta0,start(2),lo(2),hi(2),best(2),fbest,h(2,2),covt(2,2),jac2(2,2)
    logical::ok,invok
    meth='mle'; if(present(method)) meth=adjustl(method)
    call gumbel_pwm(x,mu0,beta0,ok)
    if (.not. ok) then
      beta0 = sqrt(max(variance_value(x),1.0e-8_dp))*sqrt(6.0_dp)/pi
      mu0 = mean_value(x)-euler_gamma*beta0
    end if
    result%gumbel=.true.; result%method=meth; result%xi=0.0_dp; result%mu=mu0; result%beta=beta0
    if(trim(meth)=='pwm') then
      result%converged=ok; result%nll=gev_nll(x,0.0_dp,mu0,beta0); result%covariance=nan_value(); result%se=nan_value()
    else
      fit_data=x; start=[mu0,log(max(beta0,1.0e-6_dp))]
      lo=[minval(x)-5.0_dp*max(1.0_dp,sqrt(max(variance_value(x),0.0_dp))),log(1.0e-8_dp)]
      hi = [maxval(x)+5.0_dp*max(1.0_dp,sqrt(max(variance_value(x),0.0_dp))), &
        log(max(1.0_dp,100.0_dp*(maxval(x)-minval(x)+1.0_dp)))]
      call nelder_mead_bounded(gum_obj,start,lo,hi,best,fbest,result%converged,result%iterations,result%evaluations)
      result%mu=best(1); result%beta=exp(best(2)); result%nll=fbest
      call numerical_hessian(gum_obj,best,h); call invert_matrix(h,covt,invok)
      result%covariance=0.0_dp; result%se=0.0_dp
      if(invok) then
        jac2=0.0_dp; jac2(1,1)=1.0_dp; jac2(2,2)=result%beta
        covt=matmul(jac2,matmul(covt,transpose(jac2)))
        result%covariance(2:3,2:3)=covt
        result%se(2)=sqrt(max(0.0_dp,covt(1,1))); result%se(3)=sqrt(max(0.0_dp,covt(2,2)))
      else; result%covariance=nan_value(); result%se=nan_value(); end if
      if(allocated(fit_data)) deallocate(fit_data)
    end if
    allocate(result%residuals(size(x))); result%residuals=exp(-(x-result%mu)/result%beta)
  end subroutine fit_gumbel

  subroutine fit_gpd(x,threshold,result,method,information)
    real(dp),intent(in)::x(:),threshold
    type(gpd_fit_result),intent(out)::result
    character(len=*),intent(in),optional::method,information
    character(len=8)::meth,info
    real(dp)::xi0,beta0,vc0(2,2),start(2),lo(2),hi(2),best(2),fbest,h(2,2),covt(2,2),jac(2,2)
    logical::ok,invok
    integer::i,j,n
    meth='mle'; info='observed'
    if (present(method)) meth=adjustl(method)
    if (present(information)) info=adjustl(information)
    n=count(x>threshold)
    result%n_exceed=n
    result%threshold=threshold
    result%exceedance_probability=real(n,dp)/real(size(x),dp)
    result%method=meth
    allocate(result%exceedances(n),result%residuals(n)); j=0
    do i=1,size(x); if(x(i)>threshold) then; j=j+1; result%exceedances(j)=x(i); end if; end do
    call gpd_pwm(x,threshold,xi0,beta0,vc0,ok)
    if(.not.ok) then; xi0=0.1_dp; beta0=max(mean_value(result%exceedances-threshold),1.0e-6_dp); end if
    result%xi=xi0; result%beta=beta0
    if(trim(meth)=='pwm') then
      result%converged=ok; result%covariance=vc0; result%nll=gpd_nll(result%exceedances-threshold,xi0,beta0)
    else
      fit_excess=result%exceedances-threshold; start=[xi0,log(max(beta0,1.0e-8_dp))]
      lo=[-2.0_dp,log(1.0e-10_dp)]; hi=[5.0_dp,log(max(1.0_dp,100.0_dp*(maxval(fit_excess)+1.0_dp)))]
      call nelder_mead_bounded(gpd_obj,start,lo,hi,best,fbest,result%converged,result%iterations,result%evaluations)
      result%xi=best(1); result%beta=exp(best(2)); result%nll=fbest
      if(trim(info)=='expected' .and. result%xi>-1.0_dp) then
        result%covariance(1,1)=(1.0_dp+result%xi)**2/real(n,dp)
        result%covariance(2,2)=2.0_dp*(1.0_dp+result%xi)*result%beta**2/real(n,dp)
        result%covariance(1,2)=-(1.0_dp+result%xi)*result%beta/real(n,dp); result%covariance(2,1)=result%covariance(1,2)
      else
        call numerical_hessian(gpd_obj,best,h); call invert_matrix(h,covt,invok)
        if(invok) then
          jac=0.0_dp; jac(1,1)=1.0_dp; jac(2,2)=result%beta
          result%covariance=matmul(jac,matmul(covt,transpose(jac)))
        else; result%covariance=nan_value(); end if
      end if
      if(allocated(fit_excess)) deallocate(fit_excess)
    end if
    do i=1,2; result%se(i)=sqrt(max(0.0_dp,result%covariance(i,i))); end do
    if(abs(result%xi)<1.0e-10_dp) then; result%residuals=(result%exceedances-threshold)/result%beta
    else; result%residuals=log(1.0_dp+result%xi*(result%exceedances-threshold)/result%beta)/result%xi; end if
  end subroutine fit_gpd

  real(dp) function gev_obj(p) result(v)
    real(dp),intent(in)::p(:)
    v=gev_nll(fit_data,p(1),p(2),exp(p(3)))
  end function gev_obj
  real(dp) function gum_obj(p) result(v)
    real(dp),intent(in)::p(:)
    v=gev_nll(fit_data,0.0_dp,p(1),exp(p(2)))
  end function gum_obj
  real(dp) function gpd_obj(p) result(v)
    real(dp),intent(in)::p(:)
    v=gpd_nll(fit_excess,p(1),exp(p(2)))
  end function gpd_obj
end module fextremes_fit
