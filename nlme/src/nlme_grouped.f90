! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_grouped
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_DIMENSION_ERROR
  use nlme_types, only : nonlinear_model, nonlinear_result, nlme_control
  use nlme_linalg, only : unique_integers, find_group_indices, solve_least_squares
  use nlme_nonlinear, only : fit_gnls
  implicit none
  private
  public :: is_balanced, group_summary, fit_lm_list, fit_nls_list
contains
  logical function is_balanced(group,covariate) result(balanced)
    integer, intent(in) :: group(:)
    real(dp), intent(in), optional :: covariate(:)
    integer, allocatable :: levels(:),idx(:),first_idx(:)
    integer :: status,g,n0
    balanced=.false.
    if(size(group)==0)return
    if(present(covariate))then
      if(size(covariate)/=size(group))return
    end if
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)return
    call find_group_indices(group,levels(1),first_idx);n0=size(first_idx)
    do g=2,size(levels)
      call find_group_indices(group,levels(g),idx)
      if(size(idx)/=n0)return
      if(present(covariate))then
        if(maxval(abs(covariate(idx)-covariate(first_idx)))>100.0_dp*epsilon(1.0_dp))return
      end if
    end do
    balanced=.true.
  end function is_balanced

  subroutine group_summary(values,group,levels,means,sds,counts,status)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: group(:)
    integer, allocatable, intent(out) :: levels(:),counts(:)
    real(dp), allocatable, intent(out) :: means(:),sds(:)
    integer, intent(out) :: status
    integer, allocatable :: idx(:)
    integer :: g
    if(size(values)/=size(group))then
      allocate(levels(0),counts(0),means(0),sds(0));status=NLME_DIMENSION_ERROR;return
    end if
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)return
    allocate(counts(size(levels)),means(size(levels)),sds(size(levels)))
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx);counts(g)=size(idx)
      means(g)=sum(values(idx))/real(counts(g),dp)
      if(counts(g)>1)then
        sds(g)=sqrt(sum((values(idx)-means(g))**2)/real(counts(g)-1,dp))
      else
        sds(g)=0.0_dp
      end if
    end do
    status=NLME_SUCCESS
  end subroutine group_summary

  subroutine fit_lm_list(y,x,group,levels,coefficients,sigma,status)
    real(dp), intent(in) :: y(:),x(:,:)
    integer, intent(in) :: group(:)
    integer, allocatable, intent(out) :: levels(:)
    real(dp), allocatable, intent(out) :: coefficients(:,:),sigma(:)
    integer, intent(out) :: status
    integer, allocatable :: idx(:)
    real(dp), allocatable :: beta(:),cov(:,:)
    real(dp) :: rss
    integer :: g,p
    if(size(x,1)/=size(y) .or. size(group)/=size(y))then
      allocate(levels(0),coefficients(0,0),sigma(0));status=NLME_DIMENSION_ERROR;return
    end if
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)return
    p=size(x,2);allocate(coefficients(size(levels),p),sigma(size(levels)))
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx)
      call solve_least_squares(x(idx,:),y(idx),beta,cov,rss,status)
      if(status/=NLME_SUCCESS)return
      coefficients(g,:)=beta;sigma(g)=sqrt(rss/real(max(1,size(idx)-p),dp))
    end do
  end subroutine fit_lm_list

  subroutine fit_nls_list(model,y,xdata,group,theta0,levels,parameters,sigma,status,control)
    procedure(nonlinear_model) :: model
    real(dp), intent(in) :: y(:),xdata(:,:),theta0(:)
    integer, intent(in) :: group(:)
    integer, allocatable, intent(out) :: levels(:)
    real(dp), allocatable, intent(out) :: parameters(:,:),sigma(:)
    integer, intent(out) :: status
    type(nlme_control), intent(in), optional :: control
    integer, allocatable :: idx(:)
    type(nonlinear_result) :: fit
    integer :: g
    if(size(xdata,1)/=size(y) .or. size(group)/=size(y))then
      allocate(levels(0),parameters(0,0),sigma(0));status=NLME_DIMENSION_ERROR;return
    end if
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)return
    allocate(parameters(size(levels),size(theta0)),sigma(size(levels)))
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx)
      call fit_gnls(model,y(idx),xdata(idx,:),theta0,fit,control=control)
      if(fit%status/=NLME_SUCCESS)then;status=fit%status;return;end if
      parameters(g,:)=fit%parameters;sigma(g)=fit%sigma
    end do
    status=NLME_SUCCESS
  end subroutine fit_nls_list
end module nlme_grouped
