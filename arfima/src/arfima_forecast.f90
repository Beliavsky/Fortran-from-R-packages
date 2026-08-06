module arfima_forecast_mod
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_singular
  use arfima_types, only : arfima_spec, arfima_parameters, arfima_forecast_result, arfima_error, set_error
  use arfima_autocov, only : tacvf_arfima
  use arfima_polynomial, only : difference_series, differencing_polynomial
  use arfima_linalg, only : solve_linear
  implicit none
  private
  public :: arfima_forecast, stationary_forecast
contains

  subroutine stationary_forecast(y,r,mean_value,h,forecast,error)
    real(dp),intent(in)::y(:),r(:),mean_value
    integer,intent(in)::h
    type(arfima_forecast_result),intent(out)::forecast
    type(arfima_error),intent(out)::error
    real(dp),allocatable::tmat(:,:),kmat(:,:),centered(:),weights(:),cond(:,:)
    integer::n,i,j,info

    call set_error(error,arfima_ok,'')
    n=size(y)
    if(n<1 .or. h<1 .or. size(r)<n+h) then
      call set_error(error,arfima_invalid_input,'r must contain n+h autocovariances'); forecast%error=error; return
    end if
    allocate(tmat(n,n),kmat(h,n),centered(n),forecast%mean(h),forecast%covariance(h,h),forecast%standard_error(h))
    do i=1,n
      do j=1,n
        tmat(i,j)=r(abs(i-j)+1)
      end do
    end do
    do i=1,h
      do j=1,n
        kmat(i,j)=r(n+i-j+1)
      end do
    end do
    centered=y-mean_value
    allocate(cond(h,h)); cond=0.0_dp
    do i=1,h
      call solve_linear(tmat,kmat(i,:),weights,info)
      if(info/=0) then
        call set_error(error,arfima_singular,'past covariance matrix is singular'); forecast%error=error; return
      end if
      forecast%mean(i)=mean_value+dot_product(weights,centered)
      do j=1,h
        cond(i,j)=dot_product(weights,kmat(j,:))
      end do
    end do
    do i=1,h
      do j=1,h
        forecast%covariance(i,j)=r(abs(i-j)+1)-cond(i,j)
      end do
      forecast%standard_error(i)=sqrt(max(0.0_dp,forecast%covariance(i,i)))
    end do
    forecast%error=error
  end subroutine stationary_forecast

  subroutine arfima_forecast(spec,params,z,h,forecast,error,xreg_future)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    real(dp),intent(in)::z(:)
    integer,intent(in)::h
    type(arfima_forecast_result),intent(out)::forecast
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::xreg_future(:,:)
    real(dp),allocatable::base(:),y(:),r(:),a(:),zf(:),w(:,:),covz(:,:)
    type(arfima_forecast_result)::sf
    type(arfima_error)::err
    integer::n,lag,i,j,u,idx

    call set_error(error,arfima_ok,'')
    if(size(z)<1 .or. h<1) then
      call set_error(error,arfima_invalid_input,'z and horizon must be nonempty'); forecast%error=error; return
    end if
    allocate(base(size(z))); base=z
    if(spec%use_transfer) then
      call set_error(error,arfima_invalid_input,'forecasting dynamic transfer functions is not implemented'); forecast%error=error; return
    end if
    if(spec%use_regression) then
      if(.not.allocated(spec%xreg) .or. .not.allocated(params%beta) .or. size(spec%xreg,1)/=size(z) .or. &
         size(spec%xreg,2)/=size(params%beta)) then
        call set_error(error,arfima_invalid_input,'historical xreg is incompatible'); forecast%error=error; return
      end if
      if(.not.present(xreg_future)) then
        call set_error(error,arfima_invalid_input,'xreg_future is required'); forecast%error=error; return
      end if
      if(size(xreg_future,1)/=h .or. size(xreg_future,2)/=size(params%beta)) then
        call set_error(error,arfima_invalid_input,'xreg_future has incompatible dimensions'); forecast%error=error; return
      end if
      base=base-matmul(spec%xreg,params%beta)
    end if
    lag=spec%dint+spec%dseas*spec%period
    if(lag>0) then
      call difference_series(base,spec%dint,spec%dseas,spec%period,y,err)
      if(err%code/=arfima_ok) then; error=err; forecast%error=error; return; end if
    else
      allocate(y(size(base))); y=base
    end if
    n=size(y)
    call tacvf_arfima(spec,params,n+h-1,1.0_dp,r,err)
    if(err%code/=arfima_ok) then; error=err; forecast%error=error; return; end if
    call stationary_forecast(y,r,params%mean,h,sf,err)
    if(err%code/=arfima_ok) then; error=err; forecast%error=error; return; end if
    if(lag==0) then
      forecast=sf
      if(spec%use_regression) forecast%mean=forecast%mean+matmul(xreg_future,params%beta)
      error=forecast%error
      return
    end if
    a=differencing_polynomial(spec%dint,spec%dseas,spec%period)
    allocate(zf(h),w(h,h)); zf=0.0_dp; w=0.0_dp
    do i=1,h
      zf(i)=sf%mean(i)
      w(i,i)=1.0_dp
      do j=1,lag
        idx=size(base)+i-j
        if(idx>size(base)) then
          zf(i)=zf(i)-a(j+1)*zf(idx-size(base))
        else
          zf(i)=zf(i)-a(j+1)*base(idx)
        end if
        if(i-j>=1) then
          do u=1,h
            w(i,u)=w(i,u)-a(j+1)*w(i-j,u)
          end do
        end if
      end do
    end do
    covz=matmul(w,matmul(sf%covariance,transpose(w)))
    allocate(forecast%mean(h),forecast%covariance(h,h),forecast%standard_error(h))
    forecast%mean=zf
    if(spec%use_regression) forecast%mean=forecast%mean+matmul(xreg_future,params%beta)
    forecast%covariance=covz
    do i=1,h; forecast%standard_error(i)=sqrt(max(0.0_dp,covz(i,i))); end do
    forecast%error=error
  end subroutine arfima_forecast
end module arfima_forecast_mod
