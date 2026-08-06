module arfima_durbin
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_not_positive_definite
  use arfima_types, only : dl_result, arfima_error, set_error
  implicit none
  private
  public :: durbin_levinson, dl_loglikelihood, dl_residuals, dl_simulate

contains

  subroutine durbin_levinson(r, y, result)
    real(dp), intent(in) :: r(:), y(:)
    type(dl_result), intent(out) :: result
    integer :: n, k, j
    real(dp) :: reflection, vprev, vnew, s, prediction
    real(dp), allocatable :: phi(:), old(:)

    call set_error(result%error,arfima_ok,'')
    n=size(y)
    if(n<1 .or. size(r)<n) then
      allocate(result%residuals(0),result%innovation_variance(0))
      call set_error(result%error,arfima_invalid_input,'r must contain at least n autocovariances')
      return
    end if
    allocate(result%residuals(n),result%innovation_variance(n))
    if(r(1)<=epsilon(1.0_dp)) then
      call set_error(result%error,arfima_not_positive_definite,'lag-zero autocovariance is not positive')
      result%loglik=-huge(1.0_dp)
      return
    end if
    result%residuals(1)=y(1)
    result%innovation_variance(1)=r(1)
    vprev=r(1)
    allocate(phi(0))
    do k=1,n-1
      s=0.0_dp
      do j=1,k-1
        s=s+phi(j)*r(k-j+1)
      end do
      reflection=(r(k+1)-s)/vprev
      if(abs(reflection)>=1.0_dp-100.0_dp*epsilon(1.0_dp)) then
        call set_error(result%error,arfima_not_positive_definite,'autocovariance sequence is not positive definite')
        result%loglik=-huge(1.0_dp)
        return
      end if
      allocate(old(k-1))
      if(k>1) old=phi
      deallocate(phi)
      allocate(phi(k))
      do j=1,k-1
        phi(j)=old(j)-reflection*old(k-j)
      end do
      phi(k)=reflection
      if(allocated(old)) deallocate(old)
      prediction=0.0_dp
      do j=1,k
        prediction=prediction+phi(j)*y(k-j+1)
      end do
      result%residuals(k+1)=y(k+1)-prediction
      vnew=vprev*(1.0_dp-reflection*reflection)
      if(vnew<=epsilon(1.0_dp)*r(1)) then
        call set_error(result%error,arfima_not_positive_definite,'innovation variance became nonpositive')
        result%loglik=-huge(1.0_dp)
        return
      end if
      result%innovation_variance(k+1)=vnew
      vprev=vnew
    end do
    s=sum(result%residuals**2/result%innovation_variance)
    result%sigma2_mle=s/real(n,dp)
    if(result%sigma2_mle<=0.0_dp) then
      result%loglik=huge(1.0_dp)
    else
      result%loglik=-0.5_dp*real(n,dp)*log(result%sigma2_mle)-0.5_dp*sum(log(result%innovation_variance))
    end if
  end subroutine durbin_levinson

  real(dp) function dl_loglikelihood(r,y,error,sigma2_mle) result(loglik)
    real(dp),intent(in)::r(:),y(:)
    type(arfima_error),intent(out),optional::error
    real(dp),intent(out),optional::sigma2_mle
    type(dl_result)::res
    call durbin_levinson(r,y,res)
    loglik=res%loglik
    if(present(error)) error=res%error
    if(present(sigma2_mle)) sigma2_mle=res%sigma2_mle
  end function dl_loglikelihood

  subroutine dl_residuals(r,y,residuals,innovation_variance,error)
    real(dp),intent(in)::r(:),y(:)
    real(dp),allocatable,intent(out)::residuals(:)
    real(dp),allocatable,intent(out),optional::innovation_variance(:)
    type(arfima_error),intent(out)::error
    type(dl_result)::res
    call durbin_levinson(r,y,res)
    error=res%error
    residuals=res%residuals
    if(present(innovation_variance)) innovation_variance=res%innovation_variance
  end subroutine dl_residuals

  subroutine dl_simulate(r,innov,z,error)
    real(dp),intent(in)::r(:),innov(:)
    real(dp),allocatable,intent(out)::z(:)
    type(arfima_error),intent(out)::error
    integer::n,k,j
    real(dp)::reflection,vprev,vnew,s,prediction
    real(dp),allocatable::phi(:),old(:)
    call set_error(error,arfima_ok,'')
    n=size(innov)
    if(n<1 .or. size(r)<n .or. r(1)<=epsilon(1.0_dp)) then
      allocate(z(0)); call set_error(error,arfima_invalid_input,'invalid simulation inputs'); return
    end if
    allocate(z(n)); z(1)=sqrt(r(1))*innov(1)
    vprev=r(1); allocate(phi(0))
    do k=1,n-1
      s=0.0_dp
      do j=1,k-1; s=s+phi(j)*r(k-j+1); end do
      reflection=(r(k+1)-s)/vprev
      if(abs(reflection)>=1.0_dp) then
        call set_error(error,arfima_not_positive_definite,'autocovariance sequence is not positive definite'); return
      end if
      allocate(old(k-1)); if(k>1) old=phi
      deallocate(phi); allocate(phi(k))
      do j=1,k-1; phi(j)=old(j)-reflection*old(k-j); end do
      phi(k)=reflection
      if(allocated(old)) deallocate(old)
      vnew=vprev*(1.0_dp-reflection*reflection)
      if(vnew<=epsilon(1.0_dp)*r(1)) then
        call set_error(error,arfima_not_positive_definite,'innovation variance became nonpositive'); return
      end if
      prediction=0.0_dp
      do j=1,k; prediction=prediction+phi(j)*z(k-j+1); end do
      z(k+1)=prediction+sqrt(vnew)*innov(k+1)
      vprev=vnew
    end do
  end subroutine dl_simulate

end module arfima_durbin
