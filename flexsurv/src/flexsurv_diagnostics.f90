! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_diagnostics
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, flexsurv_spec, flexsurv_result, &
    predict_cumhaz, predict_hazard
  use flexsurv_math, only : near_positive_definite, rng_normal
  implicit none
  private
  public :: coxsnell_residuals, hazard_ratio, normboot_parameters
contains

  subroutine coxsnell_residuals(data,spec,residual_theta,resid)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::residual_theta(:)
    real(dp),intent(out)::resid(size(data%lower))
    integer::i
    do i=1,size(resid)
      resid(i)=predict_cumhaz(spec,residual_theta,i,data%lower(i))
      if(data%bhazard(i)>0.0_dp) &
        resid(i)=resid(i)-log(max(data%bcondsurv(i),tiny(1.0_dp)))
    end do
  end subroutine coxsnell_residuals

  real(dp) function hazard_ratio(spec,theta,row1,row0,t) result(hr)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    integer,intent(in)::row1,row0
    real(dp)::h0
    h0=predict_hazard(spec,theta,row0,t)
    if(h0<=0.0_dp)then;hr=huge(1.0_dp);else;hr=predict_hazard(spec,theta,row1,t)/h0;end if
  end function hazard_ratio

  subroutine normboot_parameters(res,nsim,draws,seed)
    type(flexsurv_result),intent(in)::res
    integer,intent(in)::nsim
    real(dp),allocatable,intent(out)::draws(:,:)
    integer,intent(in),optional::seed
    real(dp),allocatable::cov(:,:),covpd(:,:),l(:,:),z(:)
    integer::p,i,k,st
    if(present(seed))call set_seed(seed)
    p=size(res%theta);allocate(draws(p,nsim),cov(p,p),covpd(p,p),l(p,p),z(p))
    cov=res%covariance;call near_positive_definite(cov,covpd,1.0e-10_dp)
    call chol_lower(covpd,l,st)
    if(st/=0)then
      draws=spread(res%theta,2,nsim);return
    end if
    do k=1,nsim
      do i=1,p;z(i)=rng_normal();end do
      draws(:,k)=res%theta+matmul(l,z)
    end do
  end subroutine normboot_parameters

  subroutine chol_lower(a,l,status)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(size(a,1),size(a,2))
    integer,intent(out)::status
    real(dp)::s
    integer::i,j,k,n
    n=size(a,1);l=0.0_dp;status=0
    do i=1,n
      do j=1,i
        s=a(i,j);do k=1,j-1;s=s-l(i,k)*l(j,k);end do
        if(i==j)then
          if(s<=0.0_dp)then;status=1;return;end if
          l(i,j)=sqrt(s)
        else;l(i,j)=s/l(j,j);end if
      end do
    end do
  end subroutine chol_lower

  subroutine set_seed(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n));do i=1,n;put(i)=mod(abs(seed)+104729*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed
end module flexsurv_diagnostics
