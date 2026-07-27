! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_pca_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use yc_kinds, only : dp
  use yc_types, only : pca_result_t
  use yc_linalg, only : symmetric_eigen, sort_eigen_descending
  implicit none
  private
  public :: yc_pca

contains

  function yc_pca(curves_matrix,n_components,scale_data) result(out)
    real(dp),intent(in)::curves_matrix(:,:)
    integer,intent(in),optional::n_components
    logical,intent(in),optional::scale_data
    type(pca_result_t)::out
    real(dp),allocatable::x(:,:),cov(:,:),eig(:),vec(:,:)
    real(dp)::total
    integer::nobs,nvar,nc,i,j
    logical::do_scale,ok
    nobs=size(curves_matrix,1);nvar=size(curves_matrix,2)
    if(nobs<3.or.nvar<2)then;out%ok=.false.;out%message='PCA needs at least 3 rows and 2 columns.';return;end if
    do i=1,nobs;do j=1,nvar
      if(.not.ieee_is_finite(curves_matrix(i,j)))then;out%ok=.false.;out%message='PCA data must be finite.';return;end if
    end do;end do
    nc=min(3,min(nvar,nobs-1));if(present(n_components))nc=min(max(1,n_components),min(nvar,nobs-1))
    do_scale=.false.;if(present(scale_data))do_scale=scale_data
    allocate(x(nobs,nvar),out%center(nvar),out%scale(nvar));x=curves_matrix
    do j=1,nvar
      out%center(j)=sum(x(:,j))/real(nobs,dp)
      x(:,j)=x(:,j)-out%center(j)
      out%scale(j)=sqrt(sum(x(:,j)**2)/real(nobs-1,dp))
      if(do_scale)then
        if(out%scale(j)<=sqrt(epsilon(1.0_dp)))then;out%ok=.false.;out%message='Cannot scale a constant tenor.';return;end if
        x(:,j)=x(:,j)/out%scale(j)
      else
        out%scale(j)=1.0_dp
      end if
    end do
    allocate(cov(nvar,nvar));cov=matmul(transpose(x),x)/real(nobs-1,dp)
    call symmetric_eigen(cov,eig,vec,ok)
    if(.not.ok)then;out%ok=.false.;out%message='PCA eigendecomposition failed.';return;end if
    call sort_eigen_descending(eig,vec)
    eig=max(eig,0.0_dp);total=sum(eig)
    if(total<=0.0_dp)then;out%ok=.false.;out%message='PCA total variance is zero.';return;end if
    do j=1,nvar
      i=maxloc(abs(vec(:,j)),dim=1)
      if(vec(i,j)<0.0_dp)vec(:,j)=-vec(:,j)
    end do
    out%n_components=nc
    allocate(out%loadings(nvar,nc),out%scores(nobs,nc),out%variance_explained(nc), &
      out%cumulative_variance(nc),out%sdev(nc))
    out%loadings=vec(:,1:nc);out%scores=matmul(x,out%loadings)
    out%sdev=sqrt(eig(1:nc));out%variance_explained=eig(1:nc)/total
    out%cumulative_variance(1)=out%variance_explained(1)
    do i=2,nc;out%cumulative_variance(i)=out%cumulative_variance(i-1)+out%variance_explained(i);end do
  end function yc_pca

end module yc_pca_mod
