module nfcp_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use nfcp_types, only : dp, nfcp_model_t, nfcp_ok, nfcp_invalid_input
  implicit none
  private
  public :: tsfit_volatility

contains

  subroutine tsfit_volatility(model, futures, futures_ttm, dt, theoretical, empirical, status)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: futures(:,:), futures_ttm(:), dt
    real(dp), allocatable, intent(out) :: theoretical(:), empirical(:)
    integer, intent(out), optional :: status
    integer :: j, i, k, nret
    real(dp) :: variance, corr, mean_return
    real(dp), allocatable :: returns(:)
    real(dp) :: nanv

    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    if(present(status)) status=nfcp_invalid_input
    if(size(futures,2)/=size(futures_ttm) .or. dt<=0.0_dp .or. any(futures_ttm<0.0_dp)) then
      allocate(theoretical(0),empirical(0)); return
    end if
    allocate(theoretical(size(futures_ttm)),empirical(size(futures_ttm)))
    theoretical=nanv; empirical=nanv
    do j=1,size(futures_ttm)
      variance=0.0_dp
      do i=1,model%n_factors
        do k=1,model%n_factors
          corr=model%rho(i,k)
          variance=variance+model%sigma(i)*model%sigma(k)*corr* &
            exp(-(model%kappa(i)+model%kappa(k))*futures_ttm(j))
        end do
      end do
      theoretical(j)=sqrt(max(0.0_dp,variance))
      allocate(returns(max(0,size(futures,1)-1)))
      nret=0
      do i=2,size(futures,1)
        if(ieee_is_finite(futures(i,j)) .and. ieee_is_finite(futures(i-1,j)) .and. &
           futures(i,j)>0.0_dp .and. futures(i-1,j)>0.0_dp) then
          nret=nret+1
          returns(nret)=log(futures(i,j)/futures(i-1,j))
        end if
      end do
      if(nret>=2) then
        mean_return=sum(returns(:nret))/real(nret,dp)
        empirical(j)=sqrt(sum((returns(:nret)-mean_return)**2)/(real(nret+1,dp)*dt))
      end if
      deallocate(returns)
    end do
    if(present(status)) status=nfcp_ok
  end subroutine tsfit_volatility

end module nfcp_analysis
