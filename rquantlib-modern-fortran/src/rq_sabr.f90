! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_sabr
  use rq_kinds, only: dp
  implicit none
  private
  public :: sabr_result, sabr_lognormal_volatility, calibrate_sabr

  type :: sabr_result
    real(dp) :: alpha=0.2_dp
    real(dp) :: beta=0.5_dp
    real(dp) :: rho=0.0_dp
    real(dp) :: nu=0.5_dp
    real(dp) :: rmse=0.0_dp
    integer :: iterations=0
    integer :: status=1
    real(dp),allocatable :: fitted_vols(:)
  end type sabr_result
contains
  pure elemental real(dp) function sabr_lognormal_volatility(forward,strike,t, &
      alpha,beta,rho,nu) result(vol)
    real(dp),intent(in)::forward,strike,t,alpha,beta,rho,nu
    real(dp)::one_minus_beta,fk,logfk,z,xz,a,b,c,denom
    if(forward<=0.0_dp.or.strike<=0.0_dp.or.alpha<=0.0_dp) then
      vol=0.0_dp
      return
    end if
    one_minus_beta=1.0_dp-beta
    logfk=log(forward/strike)
    fk=(forward*strike)**(0.5_dp*one_minus_beta)
    a=alpha/(fk*(1.0_dp+one_minus_beta**2*logfk**2/24.0_dp+ &
      one_minus_beta**4*logfk**4/1920.0_dp))
    b=(one_minus_beta**2*alpha**2/(24.0_dp*fk**2)+ &
      rho*beta*nu*alpha/(4.0_dp*fk)+ &
      (2.0_dp-3.0_dp*rho**2)*nu**2/24.0_dp)*t
    if(abs(logfk)<1.0e-10_dp) then
      vol=alpha/forward**one_minus_beta*(1.0_dp+b)
      return
    end if
    z=nu*fk*logfk/alpha
    denom=sqrt(max(1.0_dp-2.0_dp*rho*z+z*z,0.0_dp))+z-rho
    c=(1.0_dp-rho)
    if(abs(denom-c)<1.0e-14_dp) then
      xz=1.0_dp
    else
      xz=log(denom/c)
    end if
    if(abs(xz)<1.0e-12_dp) then
      vol=a*(1.0_dp+b)
    else
      vol=a*z/xz*(1.0_dp+b)
    end if
  end function sabr_lognormal_volatility

  subroutine calibrate_sabr(forward,strikes,t,market_vols,result,beta_fixed, &
                            max_iter,tolerance)
    real(dp),intent(in)::forward,strikes(:),t,market_vols(:)
    type(sabr_result),intent(out)::result
    real(dp),intent(in),optional::beta_fixed,tolerance
    integer,intent(in),optional::max_iter
    real(dp)::p(3),trial(3),steps(3),best,score,tol,beta
    integer::iter,j,nmax
    beta=0.5_dp
    if(present(beta_fixed)) beta=min(max(beta_fixed,0.0_dp),1.0_dp)
    nmax=500
    if(present(max_iter)) nmax=max_iter
    tol=1.0e-8_dp
    if(present(tolerance)) tol=tolerance
    p=[max(0.01_dp,market_vols(minloc(abs(strikes-forward),dim=1))* &
       forward**(1.0_dp-beta)),0.0_dp,0.5_dp]
    steps=[0.25_dp*max(p(1),0.05_dp),0.25_dp,0.25_dp]
    best=sabr_sse(p)
    do iter=1,nmax
      do j=1,3
        trial=p
        trial(j)=bounded(j,p(j)+steps(j))
        score=sabr_sse(trial)
        if(score<best) then
          p=trial
          best=score
          cycle
        end if
        trial=p
        trial(j)=bounded(j,p(j)-steps(j))
        score=sabr_sse(trial)
        if(score<best) then
          p=trial
          best=score
        else
          steps(j)=0.7_dp*steps(j)
        end if
      end do
      if(maxval(steps)<tol) exit
    end do
    result%alpha=p(1)
    result%beta=beta
    result%rho=p(2)
    result%nu=p(3)
    result%iterations=iter
    allocate(result%fitted_vols(size(strikes)))
    result%fitted_vols=sabr_lognormal_volatility(forward,strikes,t, &
      result%alpha,beta,result%rho,result%nu)
    result%rmse=sqrt(sum((result%fitted_vols-market_vols)**2)/ &
      real(size(strikes),dp))
    result%status=0
  contains
    real(dp) function sabr_sse(x) result(s)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::v(:)
      allocate(v(size(strikes)))
      v=sabr_lognormal_volatility(forward,strikes,t,x(1),beta,x(2),x(3))
      s=sum((v-market_vols)**2)
    end function sabr_sse
    pure real(dp) function bounded(index,value) result(v)
      integer,intent(in)::index
      real(dp),intent(in)::value
      select case(index)
      case(1)
        v=min(max(value,1.0e-5_dp),5.0_dp)
      case(2)
        v=min(max(value,-0.999_dp),0.999_dp)
      case default
        v=min(max(value,1.0e-5_dp),5.0_dp)
      end select
    end function bounded
  end subroutine calibrate_sabr
end module rq_sabr
