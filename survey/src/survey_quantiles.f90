! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_quantiles
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use r_quantiles, only : r_weighted_quantile_survey
  use r_distributions, only : r_qnorm, r_qt
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t, quantile_result_t
  use survey_types, only : QRULE_MATH, QRULE_HF9
  use survey_estimators, only : svy_mean
  use survey_replicates, only : svr_var
  use survey_design, only : design_degf
  implicit none
  private
  public :: weighted_quantile, svy_quantile, rep_quantile, normal_quantile, student_t_quantile
contains

  real(dp) function weighted_quantile(x,w,p,rule) result(q)
    real(dp), intent(in) :: x(:),w(:),p
    integer, intent(in), optional :: rule
    integer :: r
    if(size(x)/=size(w)) error stop 'weighted_quantile: size mismatch'
    if(p<0 .or. p>1) error stop 'weighted_quantile: p outside [0,1]'
    r=QRULE_MATH; if(present(rule)) r=rule
    if(r<QRULE_MATH .or. r>QRULE_HF9) error stop 'weighted_quantile: unknown rule'
    if(any(w<0.0_dp)) error stop 'weighted_quantile: negative weight'
    q=r_weighted_quantile_survey(x,w,p,r)
    if(ieee_is_nan(q)) error stop 'weighted_quantile: invalid data or weights'
  end function weighted_quantile

  function svy_quantile(x,probs,design,rule,alpha) result(ans)
    real(dp), intent(in) :: x(:),probs(:)
    type(survey_design_t), intent(in) :: design
    integer, intent(in), optional :: rule
    real(dp), intent(in), optional :: alpha
    type(quantile_result_t) :: ans
    real(dp) :: a,qhat,phat,sep,crit,pl,pu
    integer :: j,r,df
    ! local workaround avoided: directly use svy_mean result below
    if(size(x)/=design%n) error stop 'svy_quantile: size mismatch'
    r=QRULE_MATH; if(present(rule)) r=rule; a=0.05_dp; if(present(alpha)) a=alpha
    allocate(ans%quantile(size(probs)),ans%se(size(probs)),ans%lower(size(probs)),ans%upper(size(probs)))
    df=design_degf(design); crit=student_t_quantile(1-a/2,max(df,1))
    do j=1,size(probs)
      qhat=weighted_quantile(x,design%weight,probs(j),r); ans%quantile(j)=qhat
      call cdf_point_uncertainty(x,qhat,design,phat,sep)
      pl=max(0.0_dp,phat-crit*sep); pu=min(1.0_dp,phat+crit*sep)
      ans%lower(j)=weighted_quantile(x,design%weight,pl,r); ans%upper(j)=weighted_quantile(x,design%weight,pu,r)
      ans%se(j)=(ans%upper(j)-ans%lower(j))/(2*crit)
    end do
  end function svy_quantile

  function rep_quantile(x,probs,design,rule,alpha) result(ans)
    real(dp), intent(in) :: x(:),probs(:)
    type(rep_design_t), intent(in) :: design
    integer, intent(in), optional :: rule
    real(dp), intent(in), optional :: alpha
    type(quantile_result_t) :: ans
    real(dp), allocatable :: theta(:,:), coef(:), vv(:,:)
    real(dp) :: a,crit
    integer :: j,r,qr
    if(size(x)/=design%n) error stop 'rep_quantile: size mismatch'
    qr=QRULE_MATH; if(present(rule)) qr=rule; a=0.05_dp; if(present(alpha)) a=alpha; crit=normal_quantile(1-a/2)
    allocate(ans%quantile(size(probs)),ans%se(size(probs)),ans%lower(size(probs)),ans%upper(size(probs)), &
             theta(design%r,size(probs)),coef(size(probs)),vv(size(probs),size(probs)))
    do j=1,size(probs); coef(j)=weighted_quantile(x,design%weight,probs(j),qr); end do
    do r=1,design%r; do j=1,size(probs); theta(r,j)=weighted_quantile(x,design%repweights(:,r),probs(j),qr); end do; end do
    vv=svr_var(theta,design%scale,design%rscales,design%mse,coef)
    ans%quantile=coef
    do j=1,size(probs); ans%se(j)=sqrt(max(0.0_dp,vv(j,j))); ans%lower(j)=coef(j)-crit*ans%se(j); ans%upper(j)=coef(j)+crit*ans%se(j); end do
  end function rep_quantile

  subroutine cdf_point_uncertainty(x,q,design,p,se)
    use survey_types, only : svystat_t
    real(dp),intent(in)::x(:),q; type(survey_design_t),intent(in)::design; real(dp),intent(out)::p,se
    real(dp),allocatable::z(:,:); type(svystat_t)::m
    allocate(z(size(x),1)); z(:,1)=merge(1.0_dp,0.0_dp,x<=q); m=svy_mean(z,design); p=m%estimate(1); se=sqrt(max(0.0_dp,m%variance(1,1)))
  end subroutine cdf_point_uncertainty

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    if(p<=0) then; x=-huge(1.0_dp); return; else if(p>=1) then; x=huge(1.0_dp); return; end if
    x=r_qnorm(p)
  end function normal_quantile

  pure real(dp) function student_t_quantile(p,df) result(t)
    real(dp),intent(in)::p;integer,intent(in)::df
    if(df<=0) then; t=normal_quantile(p); return; end if
    if(p<=0.0_dp .or. p>=1.0_dp) then; t=normal_quantile(p); return; end if
    t=r_qt(p,real(df,dp))
  end function student_t_quantile

end module survey_quantiles
