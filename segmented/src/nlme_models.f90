! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_models
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_DIMENSION_ERROR
  implicit none
  private
  public :: exponential_decay_model, logistic3_model, michaelis_menten_model
contains
  subroutine exponential_decay_model(theta,x,yhat,status)
    real(dp), intent(in) :: theta(:),x(:,:)
    real(dp), intent(out) :: yhat(:)
    integer, intent(out) :: status
    if(size(theta)/=3 .or. size(x,2)<1 .or. size(yhat)/=size(x,1))then
      status=NLME_DIMENSION_ERROR;return
    end if
    yhat=theta(1)*exp(-theta(2)*x(:,1))+theta(3)
    status=NLME_SUCCESS
  end subroutine exponential_decay_model

  subroutine logistic3_model(theta,x,yhat,status)
    real(dp), intent(in) :: theta(:),x(:,:)
    real(dp), intent(out) :: yhat(:)
    integer, intent(out) :: status
    if(size(theta)/=3 .or. size(x,2)<1 .or. size(yhat)/=size(x,1) .or. abs(theta(3))<=tiny(1.0_dp))then
      status=NLME_DIMENSION_ERROR;return
    end if
    yhat=theta(1)/(1.0_dp+exp((theta(2)-x(:,1))/theta(3)))
    status=NLME_SUCCESS
  end subroutine logistic3_model

  subroutine michaelis_menten_model(theta,x,yhat,status)
    real(dp), intent(in) :: theta(:),x(:,:)
    real(dp), intent(out) :: yhat(:)
    integer, intent(out) :: status
    if(size(theta)/=2 .or. size(x,2)<1 .or. size(yhat)/=size(x,1) .or. &
       any(abs(theta(2)+x(:,1))<=tiny(1.0_dp)))then
      status=NLME_DIMENSION_ERROR;return
    end if
    yhat=theta(1)*x(:,1)/(theta(2)+x(:,1))
    status=NLME_SUCCESS
  end subroutine michaelis_menten_model
end module nlme_models
