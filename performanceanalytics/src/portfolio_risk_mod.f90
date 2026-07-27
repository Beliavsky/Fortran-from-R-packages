! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module portfolio_risk_mod
  use kinds_mod, only: dp
  use probability_mod, only: normal_pdf, normal_quantile
  use comoments_mod, only: portfolio_mean, portfolio_variance, portfolio_skewness, portfolio_kurtosis
  implicit none
  private
  public :: portfolio_gaussian_var, portfolio_modified_var
  public :: portfolio_gaussian_es, portfolio_modified_es
  public :: portfolio_risk_contributions, herfindahl_index
contains
  pure real(dp) function portfolio_gaussian_var(w,mu,cov,p) result(v)
    real(dp),intent(in)::w(:),mu(:),cov(:,:),p
    v=-portfolio_mean(w,mu)-normal_quantile(1.0_dp-p)*sqrt(max(0.0_dp,portfolio_variance(w,cov)))
  end function portfolio_gaussian_var

  pure real(dp) function portfolio_gaussian_es(w,mu,cov,p) result(v)
    real(dp),intent(in)::w(:),mu(:),cov(:,:),p
    real(dp)::z
    z=normal_quantile(1.0_dp-p)
    v=-portfolio_mean(w,mu)+normal_pdf(z)*sqrt(max(0.0_dp,portfolio_variance(w,cov)))/(1.0_dp-p)
  end function portfolio_gaussian_es

  pure real(dp) function portfolio_modified_var(w,mu,cov,m3,m4,p) result(v)
    real(dp),intent(in)::w(:),mu(:),cov(:,:),m3(:,:),m4(:,:),p
    real(dp)::z,s,k,zcf,sigma
    sigma=sqrt(max(0.0_dp,portfolio_variance(w,cov)))
    s=portfolio_skewness(w,cov,m3);k=portfolio_kurtosis(w,cov,m4)-3.0_dp
    z=normal_quantile(1.0_dp-p)
    zcf=z+(z*z-1.0_dp)*s/6.0_dp+(z**3-3.0_dp*z)*k/24.0_dp- &
      (2.0_dp*z**3-5.0_dp*z)*s*s/36.0_dp
    v=-portfolio_mean(w,mu)-zcf*sigma
  end function portfolio_modified_var

  pure real(dp) function portfolio_modified_es(w,mu,cov,m3,m4,p) result(v)
    real(dp),intent(in)::w(:),mu(:),cov(:,:),m3(:,:),m4(:,:),p
    real(dp)::z,h,s,k,e,sigma
    sigma=sqrt(max(0.0_dp,portfolio_variance(w,cov)))
    s=portfolio_skewness(w,cov,m3);k=portfolio_kurtosis(w,cov,m4)-3.0_dp
    z=normal_quantile(1.0_dp-p)
    h=z+(z*z-1.0_dp)*s/6.0_dp+(z**3-3.0_dp*z)*k/24.0_dp- &
      (2.0_dp*z**3-5.0_dp*z)*s*s/36.0_dp
    e=normal_pdf(h)*(1.0_dp+h**3*s/6.0_dp+ &
      (h**6-9.0_dp*h**4+9.0_dp*h*h+3.0_dp)*s*s/72.0_dp+ &
      (h**4-2.0_dp*h*h-1.0_dp)*k/24.0_dp)/(1.0_dp-p)
    v=-portfolio_mean(w,mu)-sigma*min(-e,h)
  end function portfolio_modified_es

  subroutine portfolio_risk_contributions(w,mu,cov,m3,m4,p,method,components,total)
    real(dp),intent(in)::w(:),mu(:),cov(:,:),m3(:,:),m4(:,:),p
    character(len=*),intent(in)::method
    real(dp),intent(out)::components(:),total
    real(dp),allocatable::wp(:),wm(:)
    real(dp)::h,vp,vm
    integer::i,n
    n=size(w);allocate(wp(n),wm(n));h=1.0e-5_dp
    select case(trim(adjustl(method)))
    case('gaussian_var');total=portfolio_gaussian_var(w,mu,cov,p)
    case('gaussian_es');total=portfolio_gaussian_es(w,mu,cov,p)
    case('modified_es');total=portfolio_modified_es(w,mu,cov,m3,m4,p)
    case default;total=portfolio_modified_var(w,mu,cov,m3,m4,p)
    end select
    do i=1,min(n,size(components))
      wp=w;wm=w;wp(i)=wp(i)+h;wm(i)=wm(i)-h
      select case(trim(adjustl(method)))
      case('gaussian_var');vp=portfolio_gaussian_var(wp,mu,cov,p);vm=portfolio_gaussian_var(wm,mu,cov,p)
      case('gaussian_es');vp=portfolio_gaussian_es(wp,mu,cov,p);vm=portfolio_gaussian_es(wm,mu,cov,p)
      case('modified_es');vp=portfolio_modified_es(wp,mu,cov,m3,m4,p);vm=portfolio_modified_es(wm,mu,cov,m3,m4,p)
      case default;vp=portfolio_modified_var(wp,mu,cov,m3,m4,p);vm=portfolio_modified_var(wm,mu,cov,m3,m4,p)
      end select
      components(i)=w(i)*(vp-vm)/(2.0_dp*h)
    end do
  end subroutine portfolio_risk_contributions

  pure real(dp) function herfindahl_index(weights,normalize) result(v)
    real(dp),intent(in)::weights(:)
    logical,intent(in),optional::normalize
    logical::norm
    integer::n
    real(dp)::h
    n=size(weights);h=sum(weights*weights);norm=.false.;if(present(normalize))norm=normalize
    if(norm .and. n>1)then;v=(h-1.0_dp/real(n,dp))/(1.0_dp-1.0_dp/real(n,dp));else;v=h;end if
  end function herfindahl_index
end module portfolio_risk_mod
