! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_parameters
  use sn_kinds, only : dp, tiny_dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_dimension_mismatch
  use sn_linalg, only : covariance_to_correlation, inverse_spd
  use sn_univariate, only : sn_uv_params, st_uv_params, delta_from_alpha, alpha_from_delta, &
                            st_cumulants
  use sn_multivariate, only : sn_mv_params, st_mv_params, delta_etc_mv
  implicit none
  private

  type, public :: uv_operational_params
    real(dp) :: xi = 0.0_dp
    real(dp) :: psi = 1.0_dp
    real(dp) :: lambda = 0.0_dp
    real(dp) :: extra = 0.0_dp
  end type uv_operational_params

  type, public :: mv_operational_params
    real(dp), allocatable :: xi(:)
    real(dp), allocatable :: psi(:,:)
    real(dp), allocatable :: lambda(:)
    real(dp) :: extra = 0.0_dp
  end type mv_operational_params

  public :: dp_to_op_uv, op_to_dp_uv, dp_to_op_mv, op_to_dp_mv
  public :: dp_to_cp_st

contains

  pure function dp_to_op_uv(params) result(op)
    type(sn_uv_params), intent(in) :: params
    type(uv_operational_params) :: op
    real(dp) :: delta
    delta = delta_from_alpha(params%alpha)
    op%xi = params%xi
    op%psi = params%omega*sqrt(max(tiny_dp,1.0_dp-delta*delta))
    op%lambda = params%alpha
    op%extra = params%tau
  end function dp_to_op_uv

  pure function op_to_dp_uv(op) result(params)
    type(uv_operational_params), intent(in) :: op
    type(sn_uv_params) :: params
    real(dp) :: delta
    delta = delta_from_alpha(op%lambda)
    params%xi = op%xi
    params%omega = op%psi/sqrt(max(tiny_dp,1.0_dp-delta*delta))
    params%alpha = op%lambda
    params%tau = op%extra
  end function op_to_dp_uv

  subroutine dp_to_op_mv(params,op,info)
    type(sn_mv_params), intent(in) :: params
    type(mv_operational_params), intent(out) :: op
    integer, intent(out) :: info
    real(dp), allocatable :: delta(:),cor(:,:),sd(:)
    real(dp) :: ds,as
    integer :: d,i,j
    info=params%validate()
    if(info/=sn_ok) return
    d=params%dimension()
    call delta_etc_mv(params%alpha,params%omega,delta,ds,as,cor,info)
    call covariance_to_correlation(params%omega,cor,sd,info)
    allocate(op%xi(d),op%psi(d,d),op%lambda(d))
    op%xi=params%xi
    do j=1,d
      do i=1,d
        op%psi(i,j)=params%omega(i,j)-sd(i)*delta(i)*sd(j)*delta(j)
      end do
    end do
    op%lambda=delta/sqrt(max(tiny_dp,1.0_dp-delta*delta))
    op%extra=params%tau
  end subroutine dp_to_op_mv

  subroutine op_to_dp_mv(op,params,info)
    type(mv_operational_params), intent(in) :: op
    type(sn_mv_params), intent(out) :: params
    integer, intent(out) :: info
    real(dp), allocatable :: psibar(:,:),psi_sd(:),inv(:,:),tmp(:),delta(:)
    real(dp) :: denom
    integer :: d,i,j
    if(.not.allocated(op%xi).or..not.allocated(op%psi).or..not.allocated(op%lambda)) then
      info=sn_invalid_argument; return
    end if
    d=size(op%xi)
    if(size(op%psi,1)/=d.or.size(op%psi,2)/=d.or.size(op%lambda)/=d) then
      info=sn_dimension_mismatch; return
    end if
    call covariance_to_correlation(op%psi,psibar,psi_sd,info)
    if(info/=sn_ok) return
    call inverse_spd(psibar,inv,info)
    if(info/=sn_ok) return
    allocate(tmp(d),delta(d),params%xi(d),params%omega(d,d),params%alpha(d))
    delta=op%lambda/sqrt(1.0_dp+op%lambda*op%lambda)
    do j=1,d
      do i=1,d
        params%omega(i,j)=op%psi(i,j)+psi_sd(i)*op%lambda(i)*psi_sd(j)*op%lambda(j)
      end do
    end do
    tmp=matmul(inv,op%lambda)
    denom=sqrt(1.0_dp+dot_product(op%lambda,tmp))
    params%alpha=(tmp/sqrt(max(tiny_dp,1.0_dp-delta*delta)))/denom
    params%xi=op%xi
    params%tau=op%extra
    info=params%validate()
  end subroutine op_to_dp_mv

  subroutine dp_to_cp_st(params,mean,sd,skewness,kurtosis,info)
    type(st_uv_params), intent(in) :: params
    real(dp), intent(out) :: mean,sd,skewness,kurtosis
    integer, intent(out) :: info
    real(dp),allocatable :: k(:)
    call st_cumulants(params,4,k,info)
    if(info/=sn_ok.or.size(k)<4.or.params%nu<=4.0_dp) then
      mean=0.0_dp; sd=0.0_dp; skewness=0.0_dp; kurtosis=0.0_dp
      info=sn_invalid_argument
      return
    end if
    mean=k(1); sd=sqrt(k(2)); skewness=k(3)/sd**3; kurtosis=k(4)/sd**4
  end subroutine dp_to_cp_st

end module sn_parameters
