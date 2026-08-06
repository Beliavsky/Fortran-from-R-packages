! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_forecast
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input
  use tvmvp_linalg, only : sample_mean
  implicit none
  private
  public :: arma_candidate_fit, comp_expected_returns

  type, public :: arma_candidate_fit
    integer :: p = 0
    integer :: q = 0
    real(dp) :: mean = 0.0_dp
    real(dp) :: ar = 0.0_dp
    real(dp) :: ma = 0.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: aic = huge(1.0_dp)
    real(dp), allocatable :: residuals(:)
  end type arma_candidate_fit
contains
  real(dp) function css_objective(y,p,q,par,residuals)
    real(dp), intent(in) :: y(:),par(:)
    integer, intent(in) :: p,q
    real(dp), intent(out), optional :: residuals(:)
    real(dp) :: mu,phi,theta,e_prev,pred,e
    integer :: t
    mu=par(1); phi=0.0_dp; theta=0.0_dp
    if (p==1) phi=par(2)
    if (q==1) theta=par(2+p)
    e_prev=y(1)-mu
    css_objective=e_prev*e_prev
    if (present(residuals)) residuals(1)=e_prev
    do t=2,size(y)
      pred=mu
      if (p==1) pred=pred+phi*(y(t-1)-mu)
      if (q==1) pred=pred+theta*e_prev
      e=y(t)-pred
      if (present(residuals)) residuals(t)=e
      css_objective=css_objective+e*e
      e_prev=e
    end do
  end function css_objective

  subroutine fit_candidate(y,p,q,fit)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: p,q
    type(arma_candidate_fit), intent(out) :: fit
    real(dp), allocatable :: par(:),trial(:),steps(:),res(:)
    real(dp) :: best,val
    integer :: k,j,iter,npar,ncoef
    npar=1+p+q; allocate(par(npar),trial(npar),steps(npar),res(size(y)))
    par=0.0_dp; par(1)=sample_mean(y)
    steps(1)=max(0.1_dp*sqrt(sum((y-par(1))**2)/real(max(1,size(y)-1),dp)),1.0e-5_dp)
    if (p==1) steps(2)=0.2_dp
    if (q==1) steps(2+p)=0.2_dp
    best=css_objective(y,p,q,par)
    do iter=1,120
      do k=1,npar
        do j=-1,1,2
          trial=par; trial(k)=par(k)+real(j,dp)*steps(k)
          if (k>1) trial(k)=max(-0.98_dp,min(0.98_dp,trial(k)))
          val=css_objective(y,p,q,trial)
          if (val<best) then
            best=val; par=trial
          end if
        end do
      end do
      steps=steps*0.85_dp
      if (maxval(steps)<1.0e-7_dp) exit
    end do
    val=css_objective(y,p,q,par,res)
    fit%p=p; fit%q=q; fit%mean=par(1)
    if (p==1) fit%ar=par(2)
    if (q==1) fit%ma=par(2+p)
    fit%variance=max(val/real(size(y),dp),tiny(1.0_dp))
    ncoef=1+p+q+1
    fit%aic=real(size(y),dp)*log(fit%variance)+2.0_dp*real(ncoef,dp)
    allocate(fit%residuals(size(y))); fit%residuals=res
  end subroutine fit_candidate

  real(dp) function mean_forecast(y,fit,horizon)
    real(dp), intent(in) :: y(:)
    type(arma_candidate_fit), intent(in) :: fit
    integer, intent(in) :: horizon
    real(dp) :: previous,forecast,last_e
    integer :: h
    previous=y(size(y)); last_e=fit%residuals(size(y)); mean_forecast=0.0_dp
    do h=1,horizon
      forecast=fit%mean
      if (fit%p==1) forecast=forecast+fit%ar*(previous-fit%mean)
      if (fit%q==1 .and. h==1) forecast=forecast+fit%ma*last_e
      mean_forecast=mean_forecast+forecast
      previous=forecast
    end do
    mean_forecast=mean_forecast/real(horizon,dp)
  end function mean_forecast

  subroutine comp_expected_returns(returns,horizon,expected,err,selected_orders)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: horizon
    real(dp), allocatable, intent(out) :: expected(:)
    type(tvmvp_error), intent(out) :: err
    integer, allocatable, intent(out), optional :: selected_orders(:,:)
    type(arma_candidate_fit) :: fits(4)
    integer, parameter :: ps(4)=[0,1,0,1], qs(4)=[0,0,1,1]
    integer :: j,k,best
    call clear_error(err)
    if (size(returns,1)<3 .or. size(returns,2)<1 .or. horizon<1) then
      allocate(expected(0)); call set_error(err,tvmvp_invalid_input,'invalid expected-return forecast input'); return
    end if
    allocate(expected(size(returns,2)))
    if (present(selected_orders)) allocate(selected_orders(size(returns,2),2))
    do j=1,size(returns,2)
      do k=1,4
        call fit_candidate(returns(:,j),ps(k),qs(k),fits(k))
      end do
      best=1
      do k=2,4
        if (fits(k)%aic<fits(best)%aic) best=k
      end do
      expected(j)=mean_forecast(returns(:,j),fits(best),horizon)
      if (present(selected_orders)) selected_orders(j,:)=[fits(best)%p,fits(best)%q]
    end do
  end subroutine comp_expected_returns
end module tvmvp_forecast
