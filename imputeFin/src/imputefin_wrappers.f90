! SPDX-License-Identifier: GPL-3.0-only
module imputefin_wrappers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  use imputefin_kinds, only : dp
  use imputefin_types, only : ar1_options, impute_ok, impute_invalid_input
  use imputefin_ar1_gaussian, only : impute_rolling_ar1_gaussian
  implicit none
  private
  public :: impute_ohlc, impute_vol
contains
  elemental function safe_log_value(x) result(v)
    real(dp),intent(in)::x
    real(dp)::v
    if(ieee_is_nan(x))then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=log(x)
    end if
  end function safe_log_value

  elemental function safe_exp_value(x) result(v)
    real(dp),intent(in)::x
    real(dp)::v
    if(ieee_is_nan(x))then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=exp(x)
    end if
  end function safe_exp_value

  elemental function safe_difference(x,y) result(v)
    real(dp),intent(in)::x,y
    real(dp)::v
    if(ieee_is_nan(x).or.ieee_is_nan(y))then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=x-y
    end if
  end function safe_difference

  elemental function safe_sum(x,y) result(v)
    real(dp),intent(in)::x,y
    real(dp)::v
    if(ieee_is_nan(x).or.ieee_is_nan(y))then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=x+y
    end if
  end function safe_sum

  logical function has_nonpositive_observation(x)
    real(dp),intent(in)::x(:,:)
    integer::i,j
    has_nonpositive_observation=.false.
    do j=1,size(x,2)
      do i=1,size(x,1)
        if(.not.ieee_is_nan(x(i,j)))then
          if(x(i,j)<=0.0_dp)then
            has_nonpositive_observation=.true.
            return
          end if
        end if
      end do
    end do
  end function has_nonpositive_observation

  logical function vector_has_nonpositive_observation(x)
    real(dp),intent(in)::x(:)
    integer::i
    vector_has_nonpositive_observation=.false.
    do i=1,size(x)
      if(.not.ieee_is_nan(x(i)))then
        if(x(i)<=0.0_dp)then
          vector_has_nonpositive_observation=.true.
          return
        end if
      end if
    end do
  end function vector_has_nonpositive_observation

  subroutine impute_ohlc(y_ohlc,y_imputed,rolling_window,remove_outliers,outlier_prob_th,tol,maxiter,seed,status)
    real(dp),intent(in)::y_ohlc(:,:)
    real(dp),allocatable,intent(out)::y_imputed(:,:)
    integer,intent(in),optional::rolling_window,maxiter
    logical,intent(in),optional::remove_outliers
    real(dp),intent(in),optional::outlier_prob_th,tol
    integer(kind=8),intent(in),optional::seed
    integer,intent(out),optional::status
    real(dp),allocatable::ly(:,:),close_imp(:),spreadv(:),work(:)
    type(ar1_options)::opt
    integer::rw,st,i
    if(size(y_ohlc,2)/=4)then
      allocate(y_imputed(0,0));if(present(status))status=impute_invalid_input;return
    end if
    if(has_nonpositive_observation(y_ohlc))then
      allocate(y_imputed(0,0));if(present(status))status=impute_invalid_input;return
    end if
    rw=252;if(present(rolling_window))rw=rolling_window
    opt=ar1_options();if(present(remove_outliers))opt%remove_outliers=remove_outliers
    if(present(outlier_prob_th))opt%outlier_prob_th=outlier_prob_th
    if(present(tol))opt%tol=tol;if(present(maxiter))opt%maxiter=maxiter
    ly=safe_log_value(y_ohlc);allocate(y_imputed(size(y_ohlc,1),4),work(size(y_ohlc,1)))
    opt%random_walk=.true.;opt%zero_mean=.false.
    call impute_rolling_ar1_gaussian(ly(:,4),close_imp,opt,rw,seed,st)
    ly(:,4)=close_imp
    opt%random_walk=.false.;opt%zero_mean=.true.
    work=safe_difference(ly(:,1),ly(:,4))
    call impute_rolling_ar1_gaussian(work,spreadv,opt,rw,seed,st)
    ly(:,1)=safe_sum(ly(:,4),spreadv)
    opt%zero_mean=.false.
    work=safe_difference(ly(:,2),ly(:,4))
    call impute_rolling_ar1_gaussian(work,spreadv,opt,rw,seed,st)
    do i=1,size(spreadv)
      if(.not.ieee_is_nan(spreadv(i)))spreadv(i)=max(0.0_dp,spreadv(i))
    end do
    ly(:,2)=safe_sum(ly(:,4),spreadv)
    work=safe_difference(ly(:,3),ly(:,4))
    call impute_rolling_ar1_gaussian(work,spreadv,opt,rw,seed,st)
    do i=1,size(spreadv)
      if(.not.ieee_is_nan(spreadv(i)))spreadv(i)=min(0.0_dp,spreadv(i))
    end do
    ly(:,3)=safe_sum(ly(:,4),spreadv)
    y_imputed=safe_exp_value(ly)
    if(present(status))status=impute_ok
  end subroutine impute_ohlc

  subroutine impute_vol(y_vol,y_imputed,rolling_window,remove_outliers,outlier_prob_th,tol,maxiter,seed,status)
    real(dp),intent(in)::y_vol(:)
    real(dp),allocatable,intent(out)::y_imputed(:)
    integer,intent(in),optional::rolling_window,maxiter
    logical,intent(in),optional::remove_outliers
    real(dp),intent(in),optional::outlier_prob_th,tol
    integer(kind=8),intent(in),optional::seed
    integer,intent(out),optional::status
    real(dp),allocatable::ly(:),ly_imp(:)
    type(ar1_options)::opt
    integer::rw,st
    if(vector_has_nonpositive_observation(y_vol))then
      allocate(y_imputed(0));if(present(status))status=impute_invalid_input;return
    end if
    rw=252;if(present(rolling_window))rw=rolling_window
    opt=ar1_options();if(present(remove_outliers))opt%remove_outliers=remove_outliers
    if(present(outlier_prob_th))opt%outlier_prob_th=outlier_prob_th
    if(present(tol))opt%tol=tol;if(present(maxiter))opt%maxiter=maxiter
    ly=safe_log_value(y_vol)
    call impute_rolling_ar1_gaussian(ly,ly_imp,opt,rw,seed,st)
    y_imputed=safe_exp_value(ly_imp)
    if(present(status))status=st
  end subroutine impute_vol
end module imputefin_wrappers
