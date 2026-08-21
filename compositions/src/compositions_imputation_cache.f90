! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_imputation_cache
  !! Source-oriented translation of the experimental gsiCImpAcomp cache and
  !! gsiCImpAcompClrExpectation machinery.  The upstream C routine returns an
  !! ALR-coordinate covariance while centering only the mean; both that
  !! source-compatible covariance and a properly CLR-projected covariance are
  !! returned here.
  use compositions_kinds, only: dp
  use compositions_imputation, only: mt_observed, mt_bdl, missing_pattern_indices, conditional_alr_moments
  implicit none
  private

  type, public :: imputation_pattern_cache
    integer :: nmissing = 0
    integer, allocatable :: order(:)
    real(dp), allocatable :: lambda(:,:)
    real(dp), allocatable :: residual_covariance(:,:)
    real(dp), allocatable :: residual_cholesky(:,:)
    logical :: cholesky_ok = .false.
  end type imputation_pattern_cache

  type, public :: clr_expectation_result
    real(dp), allocatable :: mean_clr(:)
    real(dp), allocatable :: covariance_source(:,:)
    real(dp), allocatable :: covariance_clr(:,:)
    real(dp), allocatable :: xlr_imputed(:,:)
    integer, allocatable :: used(:)
    integer, allocatable :: pattern_id(:)
    type(imputation_pattern_cache), allocatable :: cache(:)
    logical :: ok = .false.
  end type clr_expectation_result

  public :: build_imputation_cache, acomp_clr_expectation

contains

  subroutine build_imputation_cache(mt, clr_cov, cache, pattern_id)
    integer, intent(in) :: mt(:,:)
    real(dp), intent(in) :: clr_cov(:,:)
    type(imputation_pattern_cache), allocatable, intent(out) :: cache(:)
    integer, allocatable, intent(out) :: pattern_id(:)
    integer, allocatable :: representative(:), order(:,:), nmissing(:)
    integer :: k, np, info

    if (size(clr_cov,1)/=size(mt,2) .or. size(clr_cov,2)/=size(mt,2)) &
      error stop 'build_imputation_cache: covariance mismatch'
    call missing_pattern_indices(mt,pattern_id,representative,order,nmissing)
    np=size(representative)
    allocate(cache(np))
    do k=1,np
      cache(k)%nmissing=nmissing(k)
      cache(k)%order=order(k,:)
      call conditional_alr_moments(clr_cov,cache(k)%order,cache(k)%nmissing, &
        cache(k)%lambda,cache(k)%residual_covariance)
      if(cache(k)%nmissing>0) then
        call safe_cholesky(cache(k)%residual_covariance,cache(k)%residual_cholesky,info)
        cache(k)%cholesky_ok=(info==0)
      else
        allocate(cache(k)%residual_cholesky(0,0))
        cache(k)%cholesky_ok=.true.
      end if
    end do
  end subroutine build_imputation_cache

  function acomp_clr_expectation(comp,pred_clr,mt,clr_cov,norm_residuals,dl,source_dl_last) result(res)
    real(dp), intent(in) :: comp(:,:),pred_clr(:,:),clr_cov(:,:),norm_residuals(:,:)
    integer, intent(in) :: mt(:,:)
    real(dp), intent(in), optional :: dl(:,:)
    logical, intent(in), optional :: source_dl_last
    type(clr_expectation_result) :: res
    integer :: n,d,i,j,k,t,c,nm,ng,ref,nacc,truej
    real(dp), allocatable :: xlr(:),obs(:),pg(:),pm(:),mis(:),draw(:),sumv(:),sum2(:,:),row2(:,:)
    real(dp), allocatable :: z(:),h(:,:),tmp(:,:),limit(:)
    logical :: has_bdl,accept,use_last
    real(dp) :: lref,pref,dlref

    n=size(comp,1); d=size(comp,2)
    if(any(shape(pred_clr)/=[n,d]) .or. any(shape(mt)/=[n,d])) &
      error stop 'acomp_clr_expectation: shape mismatch'
    if(any(shape(clr_cov)/=[d,d])) error stop 'acomp_clr_expectation: covariance mismatch'
    if(present(dl)) then
      if(any(shape(dl)/=[n,d])) error stop 'acomp_clr_expectation: detection-limit mismatch'
    end if
    use_last=.false.; if(present(source_dl_last)) use_last=source_dl_last

    call build_imputation_cache(mt,clr_cov,res%cache,res%pattern_id)
    allocate(res%xlr_imputed(n,d),res%used(n),sumv(d),sum2(d,d),xlr(d),row2(d,d))
    res%xlr_imputed=0.0_dp; res%used=-1; sumv=0.0_dp; sum2=0.0_dp

    do i=1,n
      c=res%pattern_id(i); nm=res%cache(c)%nmissing; ng=d-nm; ref=res%cache(c)%order(d)
      xlr=0.0_dp; row2=0.0_dp

      if(nm==0) then
        lref=log(comp(i,ref))
        do j=1,d
          xlr(j)=log(comp(i,j))-lref
        end do
        res%xlr_imputed(i,:)=xlr
        sumv=sumv+xlr; sum2=sum2+outer_product(xlr,xlr)
        cycle
      end if

      if(ng==0) then
        pref=pred_clr(i,ref)
        do j=1,d; xlr(j)=pred_clr(i,j)-pref; end do
        res%xlr_imputed(i,:)=xlr
        row2=outer_product(xlr,xlr)
        do j=1,nm
          do k=1,nm
            row2(res%cache(c)%order(j),res%cache(c)%order(k)) = &
              row2(res%cache(c)%order(j),res%cache(c)%order(k)) + &
              res%cache(c)%residual_covariance(j,k)
          end do
        end do
        sumv=sumv+xlr; sum2=sum2+row2
        cycle
      end if

      lref=log(comp(i,ref)); pref=pred_clr(i,ref)
      allocate(obs(ng),pg(ng),pm(nm),mis(nm),draw(nm),limit(nm))
      do j=1,ng
        truej=res%cache(c)%order(nm+j)
        obs(j)=log(comp(i,truej))-lref
        pg(j)=pred_clr(i,truej)-pref
        xlr(truej)=obs(j)
      end do
      do j=1,nm
        truej=res%cache(c)%order(j)
        pm(j)=pred_clr(i,truej)-pref
      end do
      mis=pm
      if(size(res%cache(c)%lambda,2)>0) mis=pm+matmul(res%cache(c)%lambda,obs-pg)
      do j=1,nm; xlr(res%cache(c)%order(j))=mis(j); end do

      has_bdl=.false.; limit=huge(1.0_dp)
      if(present(dl)) then
        if(use_last .and. comp(i,d)>0.0_dp) then
          dlref=log(comp(i,d))
        else
          dlref=lref
        end if
        do j=1,nm
          truej=res%cache(c)%order(j)
          if(mt(i,truej)==mt_bdl .and. dl(i,truej)>0.0_dp) then
            has_bdl=.true.; limit(j)=log(dl(i,truej))-dlref
          end if
        end do
      end if

      if(.not.has_bdl) then
        res%xlr_imputed(i,:)=xlr
        row2=outer_product(xlr,xlr)
        do j=1,nm
          do k=1,nm
            row2(res%cache(c)%order(j),res%cache(c)%order(k)) = &
              row2(res%cache(c)%order(j),res%cache(c)%order(k)) + &
              res%cache(c)%residual_covariance(j,k)
          end do
        end do
        sumv=sumv+xlr; sum2=sum2+row2
      else
        if(size(norm_residuals,2)<nm) error stop 'acomp_clr_expectation: Norm has too few columns'
        nacc=0; draw=0.0_dp; row2=0.0_dp
        if(res%cache(c)%cholesky_ok) then
          allocate(z(nm))
          do t=1,size(norm_residuals,1)
            z=mis+matmul(res%cache(c)%residual_cholesky,norm_residuals(t,1:nm))
            accept=.true.
            do j=1,nm
              truej=res%cache(c)%order(j)
              if(mt(i,truej)==mt_bdl .and. z(j)>=limit(j)) then
                accept=.false.; exit
              end if
            end do
            if(accept) then
              nacc=nacc+1; draw=draw+z
              tmp=xlr_as_matrix(xlr,res%cache(c)%order,z,nm)
              row2=row2+outer_product(tmp(:,1),tmp(:,1))
            end if
          end do
          deallocate(z)
        end if
        res%used(i)=nacc
        if(nacc>0) then
          mis=draw/real(nacc,dp)
          do j=1,nm; xlr(res%cache(c)%order(j))=mis(j); end do
          res%xlr_imputed(i,:)=xlr
          sumv=sumv+xlr; sum2=sum2+row2/real(nacc,dp)
        else
          ! The upstream fallback has an out-of-bounds second-moment expression.
          ! Use its intended bounded mean rule and a defined outer product.
          do j=1,nm
            truej=res%cache(c)%order(j)
            if(mt(i,truej)==mt_bdl) mis(j)=min(mis(j),limit(j))
            xlr(truej)=mis(j)
          end do
          res%xlr_imputed(i,:)=xlr
          sumv=sumv+xlr; sum2=sum2+outer_product(xlr,xlr)
        end if
      end if
      deallocate(obs,pg,pm,mis,draw,limit)
    end do

    res%mean_clr=sumv/real(n,dp)
    res%covariance_source=sum2/real(n,dp)-outer_product(res%mean_clr,res%mean_clr)
    res%mean_clr=res%mean_clr-sum(res%mean_clr)/real(d,dp)
    allocate(h(d,d)); h=-1.0_dp/real(d,dp)
    do j=1,d; h(j,j)=h(j,j)+1.0_dp; end do
    res%covariance_clr=matmul(h,matmul(res%covariance_source,h))
    res%covariance_clr=0.5_dp*(res%covariance_clr+transpose(res%covariance_clr))
    res%ok=.true.
  end function acomp_clr_expectation

  pure function outer_product(a,b) result(c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp) :: c(size(a),size(b))
    integer :: j
    do j=1,size(b); c(:,j)=a*b(j); end do
  end function outer_product

  function xlr_as_matrix(base,order,missing,nm) result(a)
    real(dp), intent(in) :: base(:),missing(:)
    integer, intent(in) :: order(:),nm
    real(dp), allocatable :: a(:,:)
    integer :: j
    allocate(a(size(base),1)); a(:,1)=base
    do j=1,nm; a(order(j),1)=missing(j); end do
  end function xlr_as_matrix

  subroutine safe_cholesky(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:)
    real(dp) :: ridge
    integer :: tries,j
    if(size(a,1)==0) then; allocate(l(0,0)); info=0; return; end if
    ridge=0.0_dp
    do tries=1,8
      aa=a
      if(ridge>0.0_dp) then
        do j=1,size(aa,1); aa(j,j)=aa(j,j)+ridge; end do
      end if
      call chol_try(aa,l,info)
      if(info==0) return
      if(ridge==0.0_dp) then; ridge=1.0e-12_dp; else; ridge=10.0_dp*ridge; end if
    end do
  end subroutine safe_cholesky

  subroutine chol_try(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: i,j,k,n
    real(dp) :: v
    n=size(a,1); allocate(l(n,n)); l=0.0_dp; info=0
    do i=1,n
      v=a(i,i)
      do k=1,i-1; v=v-l(i,k)*l(i,k); end do
      if(v<=0.0_dp) then; info=i; return; end if
      l(i,i)=sqrt(v)
      do j=i+1,n
        v=a(j,i)
        do k=1,i-1; v=v-l(j,k)*l(i,k); end do
        l(j,i)=v/l(i,i)
      end do
    end do
  end subroutine chol_try

end module compositions_imputation_cache
