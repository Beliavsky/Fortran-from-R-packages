! SPDX-License-Identifier: MIT
module bekks_forecast
  use bekks_kinds, only: dp
  use bekks_types
  use bekks_model
  use bekks_linalg, only: symmetric_sqrt, outer_product
  use bekks_matrix, only: vech_lower
  use bekks_math, only: normal_quantile
  implicit none
  private
  public :: forecast_bekk, virf_bekk

contains

  subroutine forecast_path(par,last_h,last_return,n_ahead,signs,expected_indicator,h,status)
    type(bekk_parameters), intent(in) :: par
    real(dp), intent(in) :: last_h(:,:),last_return(:),signs(:),expected_indicator
    integer, intent(in) :: n_ahead
    real(dp), allocatable, intent(out) :: h(:,:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: previous(:,:)
    integer :: k,n,ind
    n=size(last_h,1);allocate(h(n,n,n_ahead),previous(n,n));previous=last_h
    do k=1,n_ahead
      if(k==1)then
        ind=merge(indicator_function(last_return,signs),0,par%asymmetric)
        call covariance_step_public(par,previous,last_return,ind,h(:,:,k))
      else
        call expected_covariance_step(par,previous,expected_indicator,h(:,:,k))
      end if
      previous=h(:,:,k)
    end do
    status=bekk_ok
  contains
    subroutine covariance_step_public(p,hprev,rprev,indicator,hnew)
      type(bekk_parameters), intent(in) :: p
      real(dp), intent(in) :: hprev(:,:),rprev(:)
      integer, intent(in) :: indicator
      real(dp), intent(out) :: hnew(:,:)
      real(dp) :: coeff
      hnew=matmul(p%c,transpose(p%c))
      select case(p%model_type)
      case(bekk_scalar)
        coeff=p%a_scalar
        if(p%asymmetric)coeff=coeff+real(indicator,dp)*p%b_scalar
        hnew=hnew+coeff*outer_product(rprev)+p%g_scalar*hprev
      case default
        hnew=hnew+matmul(transpose(p%a),matmul(outer_product(rprev),p%a)) &
          +matmul(transpose(p%g),matmul(hprev,p%g))
        if(p%asymmetric .and. indicator==1) &
          hnew=hnew+matmul(transpose(p%b),matmul(outer_product(rprev),p%b))
      end select
      hnew=0.5_dp*(hnew+transpose(hnew))
    end subroutine covariance_step_public
  end subroutine forecast_path

  subroutine expected_covariance_step(par,hprev,e,hnew)
    type(bekk_parameters), intent(in) :: par
    real(dp), intent(in) :: hprev(:,:),e
    real(dp), intent(out) :: hnew(:,:)
    real(dp) :: coeff
    hnew=matmul(par%c,transpose(par%c))
    select case(par%model_type)
    case(bekk_scalar)
      coeff=par%a_scalar+par%g_scalar
      if(par%asymmetric)coeff=coeff+e*par%b_scalar
      hnew=hnew+coeff*hprev
    case default
      hnew=hnew+matmul(transpose(par%a),matmul(hprev,par%a)) &
        +matmul(transpose(par%g),matmul(hprev,par%g))
      if(par%asymmetric)hnew=hnew+e*matmul(transpose(par%b),matmul(hprev,par%b))
    end select
    hnew=0.5_dp*(hnew+transpose(hnew))
  end subroutine expected_covariance_step

  subroutine forecast_bekk(fit,n_ahead,result,confidence_level)
    type(bekk_fit_result), intent(in) :: fit
    integer, intent(in) :: n_ahead
    type(bekk_forecast_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence_level
    type(bekk_parameters) :: pl,pu
    real(dp), allocatable :: lower_theta(:),upper_theta(:),hl(:,:,:),hu(:,:,:),el(:,:),eu(:,:)
    real(dp) :: ci,z
    integer :: n,t,st
    if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)then;result%status=fit%status;return;end if
    if(n_ahead<1)then;result%status=bekk_invalid_input;return;end if
    n=size(fit%data,2);t=size(fit%data,1)
    call forecast_path(fit%parameters,fit%h(:,:,t),fit%data(t,:),n_ahead,fit%signs,fit%expected_indicator,result%h,st)
    if(st/=bekk_ok)then;result%status=st;return;end if
    call covariance_to_volatility(result%h,result%standard_deviation,result%correlation)
    ci=0.95_dp;if(present(confidence_level))ci=confidence_level
    if(allocated(fit%standard_error))then
      z=normal_quantile(0.5_dp+0.5_dp*ci)
      lower_theta=fit%theta-z*fit%standard_error;upper_theta=fit%theta+z*fit%standard_error
      call unpack_parameters(lower_theta,n,fit%spec%model_type,fit%spec%asymmetric,pl,st)
      if(st==bekk_ok)then
        if(valid_parameters(pl,fit%expected_indicator))then
          call filter_bekk(lower_theta,fit%data,fit%spec%model_type, &
            fit%spec%asymmetric,fit%signs,hl,el,st)
          if(st==bekk_ok)then
            call forecast_path(pl,hl(:,:,t),fit%data(t,:),n_ahead,fit%signs, &
              fit%expected_indicator,result%covariance_lower,st)
          end if
        end if
      end if
      call unpack_parameters(upper_theta,n,fit%spec%model_type,fit%spec%asymmetric,pu,st)
      if(st==bekk_ok)then
        if(valid_parameters(pu,fit%expected_indicator))then
          call filter_bekk(upper_theta,fit%data,fit%spec%model_type, &
            fit%spec%asymmetric,fit%signs,hu,eu,st)
          if(st==bekk_ok)then
            call forecast_path(pu,hu(:,:,t),fit%data(t,:),n_ahead,fit%signs, &
              fit%expected_indicator,result%covariance_upper,st)
          end if
        end if
      end if
    end if
    if(.not.allocated(result%covariance_lower))then
      allocate(result%covariance_lower(n,n,n_ahead));result%covariance_lower=result%h
    end if
    if(.not.allocated(result%covariance_upper))then
      allocate(result%covariance_upper(n,n,n_ahead));result%covariance_upper=result%h
    end if
    result%status=bekk_ok
  end subroutine forecast_bekk

  subroutine virf_bekk(fit,h0,shock,periods,result,confidence_level)
    type(bekk_fit_result), intent(in) :: fit
    real(dp), intent(in) :: h0(:,:),shock(:)
    integer, intent(in) :: periods
    type(bekk_virf_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence_level
    real(dp), allocatable :: tp(:),tm(:),rp(:,:),rm(:,:),jacobian(:,:),parameter_covariance(:,:)
    real(dp), allocatable :: variance(:),se(:)
    real(dp) :: ci,z,step
    integer :: n,m,p,j,k,i,sp,sm,flat

    n=size(h0,1);m=n*(n+1)/2
    if(size(shock)/=n .or. periods<1)then
      result%status=bekk_invalid_input
      return
    end if
    call response_for_theta(fit%theta,result%response,result%status)
    if(result%status/=bekk_ok)return

    ci=0.90_dp
    if(present(confidence_level))ci=confidence_level
    z=normal_quantile(0.5_dp+0.5_dp*ci)
    allocate(result%lower(periods,m),result%upper(periods,m))
    if(.not.allocated(fit%covariance) .and. .not.allocated(fit%robust_covariance))then
      result%lower=result%response
      result%upper=result%response
      result%status=bekk_ok
      return
    end if

    if(allocated(fit%robust_covariance))then
      parameter_covariance=fit%robust_covariance
    else
      parameter_covariance=fit%covariance
    end if
    p=size(fit%theta)
    allocate(tp(p),tm(p),jacobian(periods*m,p),variance(periods*m),se(periods*m))
    jacobian=0.0_dp
    do j=1,p
      step=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(fit%theta(j)))
      tp=fit%theta;tm=fit%theta
      tp(j)=tp(j)+step;tm(j)=tm(j)-step
      call response_for_theta(tp,rp,sp)
      call response_for_theta(tm,rm,sm)
      if(sp==bekk_ok .and. sm==bekk_ok)then
        do k=1,periods
          do i=1,m
            flat=(k-1)*m+i
            jacobian(flat,j)=(rp(k,i)-rm(k,i))/(2.0_dp*step)
          end do
        end do
      else if(sp==bekk_ok)then
        do k=1,periods
          do i=1,m
            flat=(k-1)*m+i
            jacobian(flat,j)=(rp(k,i)-result%response(k,i))/step
          end do
        end do
      else if(sm==bekk_ok)then
        do k=1,periods
          do i=1,m
            flat=(k-1)*m+i
            jacobian(flat,j)=(result%response(k,i)-rm(k,i))/step
          end do
        end do
      end if
    end do
    variance=0.0_dp
    do flat=1,periods*m
      variance(flat)=dot_product(jacobian(flat,:),matmul(parameter_covariance,jacobian(flat,:)))
    end do
    se=sqrt(max(variance,0.0_dp))
    do k=1,periods
      do i=1,m
        flat=(k-1)*m+i
        result%lower(k,i)=result%response(k,i)-z*se(flat)
        result%upper(k,i)=result%response(k,i)+z*se(flat)
      end do
    end do
    result%status=bekk_ok

  contains

    subroutine response_for_theta(theta,response,status)
      real(dp), intent(in) :: theta(:)
      real(dp), allocatable, intent(out) :: response(:,:)
      integer, intent(out) :: status
      type(bekk_parameters) :: parameters
      real(dp), allocatable :: root(:,:),innovation_delta(:,:),delta(:,:),next_delta(:,:)
      integer :: info,indicator,period

      call unpack_parameters(theta,n,fit%spec%model_type,fit%spec%asymmetric,parameters,status)
      if(status/=bekk_ok)return
      allocate(root(n,n),innovation_delta(n,n),delta(n,n),next_delta(n,n),response(periods,m))
      call symmetric_sqrt(h0,root,info)
      if(info/=0)then
        status=bekk_linalg_failure
        return
      end if
      innovation_delta=matmul(root,matmul(outer_product(shock)-identity(n),root))
      select case(parameters%model_type)
      case(bekk_scalar)
        delta=parameters%a_scalar*innovation_delta
        if(parameters%asymmetric)then
          indicator=indicator_function(matmul(root,shock),fit%signs)
          delta=delta+real(indicator,dp)*parameters%b_scalar*innovation_delta
        end if
      case default
        delta=matmul(transpose(parameters%a),matmul(innovation_delta,parameters%a))
        if(parameters%asymmetric)then
          indicator=indicator_function(matmul(root,shock),fit%signs)
          if(indicator==1)delta=delta+matmul(transpose(parameters%b), &
            matmul(innovation_delta,parameters%b))
        end if
      end select
      do period=1,periods
        response(period,:)=vech_lower(delta)
        select case(parameters%model_type)
        case(bekk_scalar)
          next_delta=(parameters%a_scalar+parameters%g_scalar+ &
            merge(fit%expected_indicator*parameters%b_scalar,0.0_dp,parameters%asymmetric))*delta
        case default
          next_delta=matmul(transpose(parameters%a),matmul(delta,parameters%a))+ &
            matmul(transpose(parameters%g),matmul(delta,parameters%g))
          if(parameters%asymmetric)next_delta=next_delta+fit%expected_indicator* &
            matmul(transpose(parameters%b),matmul(delta,parameters%b))
        end select
        delta=0.5_dp*(next_delta+transpose(next_delta))
      end do
      status=bekk_ok
    end subroutine response_for_theta

    function identity(nn) result(a)
      integer, intent(in) :: nn
      real(dp) :: a(nn,nn)
      integer :: ii
      a=0.0_dp
      do ii=1,nn
        a(ii,ii)=1.0_dp
      end do
    end function identity

  end subroutine virf_bekk

end module bekks_forecast
