! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_multivariate
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, pca_result_t
  use survey_estimators, only : svy_covariance
  use survey_linalg, only : symmetric_eigen
  implicit none
  private
  public :: svy_cralpha, svy_pca, weighted_correlation
contains

  real(dp) function svy_cralpha(items,design) result(alpha)
    real(dp),intent(in)::items(:,:);type(survey_design_t),intent(in)::design
    real(dp),allocatable::covitems(:,:),tot(:,:),covtot(:,:);integer::k,i
    if(size(items,1)/=design%n)error stop 'svy_cralpha: row mismatch';k=size(items,2);if(k<2)error stop 'svy_cralpha: need at least two items'
    allocate(covitems(k,k),tot(design%n,1),covtot(1,1));call svy_covariance(items,design,covitems);tot(:,1)=sum(items,dim=2);call svy_covariance(tot,design,covtot)
    if(covtot(1,1)<=0)then;alpha=0;return;end if
    alpha=real(k,dp)/real(k-1,dp)*(1-sum([(covitems(i,i),i=1,k)])/covtot(1,1))
  end function svy_cralpha

  subroutine weighted_correlation(x,design,corr)
    real(dp),intent(in)::x(:,:);type(survey_design_t),intent(in)::design;real(dp),intent(out)::corr(:,:)
    real(dp),allocatable::cov(:,:),sd(:);integer::p,i,j
    p=size(x,2);if(any(shape(corr)/=[p,p]))error stop 'weighted_correlation: shape mismatch';allocate(cov(p,p),sd(p));call svy_covariance(x,design,cov)
    do i=1,p;sd(i)=sqrt(max(0.0_dp,cov(i,i)));end do
    do i=1,p;do j=1,p;if(sd(i)>0.and.sd(j)>0)then;corr(i,j)=cov(i,j)/(sd(i)*sd(j));else;corr(i,j)=0;end if;end do;end do
  end subroutine weighted_correlation

  subroutine svy_pca(x,design,result,center,scale)
    real(dp),intent(in)::x(:,:);type(survey_design_t),intent(in)::design;type(pca_result_t),intent(out)::result
    logical,intent(in),optional::center,scale
    real(dp),allocatable::z(:,:),cov(:,:),eval(:),evec(:,:);real(dp)::sw;integer::n,p,i,info;logical::doc,dos
    n=size(x,1);p=size(x,2);if(n/=design%n)error stop 'svy_pca: row mismatch';doc=.true.;if(present(center))doc=center;dos=.false.;if(present(scale))dos=scale
    allocate(z(n,p),cov(p,p),eval(p),evec(p,p),result%center(p),result%scale(p),result%rotation(p,p),result%sdev(p));z=x;sw=sum(design%weight)
    if(doc)then;do i=1,p;result%center(i)=dot_product(design%weight,x(:,i))/sw;z(:,i)=z(:,i)-result%center(i);end do;else;result%center=0;end if
    call svy_covariance(z,design,cov)
    result%scale=1
    if(dos)then;do i=1,p;result%scale(i)=sqrt(max(cov(i,i),tiny(1.0_dp)));z(:,i)=z(:,i)/result%scale(i);end do;call svy_covariance(z,design,cov);end if
    call symmetric_eigen(cov,eval,evec,info);result%rotation=evec;result%sdev=sqrt(max(eval,0.0_dp))
  end subroutine svy_pca
end module survey_multivariate
