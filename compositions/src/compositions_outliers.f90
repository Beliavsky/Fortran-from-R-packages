! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_outliers
  !! Numerical translation of the Mahalanobis calibration and single-component
  !! outlier workflows in FunctionOutliers.R.  R's in-session pStore cache is
  !! intentionally replaced by explicit simulation arrays/seed control.
  use compositions_kinds, only: dp
  use compositions_geometry, only: ilr_rows
  use compositions_linalg, only: covariance_matrix, invert_matrix
  use robustbase_detmcd, only: detmcd_result, cov_detmcd
  use bayesm_rng, only: rng_seed, randn
  implicit none
  private

  type, public :: outlier_classification_result
    real(dp), allocatable :: distance(:)
    logical, allocatable :: is_outlier(:)
    logical, allocatable :: is_extreme(:)
    logical, allocatable :: explained_by_component(:,:)
    integer, allocatable :: best_component(:)
    real(dp) :: critical_corrected = 0.0_dp
    real(dp) :: critical_empirical = 0.0_dp
  end type outlier_classification_result

  public :: acomp_mahalanobis, r_empirical_mahalanobis, r_max_mahalanobis
  public :: q_empirical_mahalanobis, q_max_mahalanobis
  public :: p_empirical_mahalanobis, p_max_mahalanobis
  public :: is_mahalanobis_outlier, outlier_classifier_best

contains

  function acomp_mahalanobis(x,robust,good_mask) result(distance)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: robust,good_mask(:)
    real(dp), allocatable :: distance(:)
    real(dp), allocatable :: z(:,:),fitz(:,:),center(:),cov(:,:),inv(:,:),u(:)
    logical :: rb
    integer :: i,nfit
    type(detmcd_result) :: mcd

    z=ilr_rows(x); rb=.false.; if(present(robust)) rb=robust
    if(present(good_mask)) then
      if(size(good_mask)/=size(z,1)) error stop 'acomp_mahalanobis: good_mask mismatch'
      nfit=count(good_mask); allocate(fitz(nfit,size(z,2))); nfit=0
      do i=1,size(z,1)
        if(good_mask(i)) then; nfit=nfit+1; fitz(nfit,:)=z(i,:); end if
      end do
    else
      fitz=z
    end if
    if(size(fitz,1)<max(2,size(fitz,2)+1)) error stop 'acomp_mahalanobis: insufficient fit rows'
    if(rb) then
      call cov_detmcd(fitz,mcd); center=mcd%estimate%center; cov=mcd%estimate%covariance
    else
      call covariance_matrix(fitz,cov,center)
    end if
    call invert_matrix(cov,inv); allocate(distance(size(z,1)),u(size(z,2)))
    do i=1,size(z,1)
      u=z(i,:)-center
      distance(i)=sqrt(max(0.0_dp,dot_product(u,matmul(inv,u))))
    end do
  end function acomp_mahalanobis

  function r_empirical_mahalanobis(replicates,nobs,dim,robust,seed) result(s)
    integer, intent(in) :: replicates,nobs,dim
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), allocatable :: s(:,:)
    real(dp), allocatable :: x(:,:),center(:),cov(:,:),inv(:,:),u(:)
    logical :: rb
    integer :: r,i,j
    type(detmcd_result) :: mcd
    if(replicates<1.or.nobs<2.or.dim<1) error stop 'r_empirical_mahalanobis: invalid dimensions'
    rb=.true.; if(present(robust)) rb=robust
    if(present(seed)) call rng_seed(seed)
    allocate(s(replicates,nobs),x(nobs,dim),u(dim))
    do r=1,replicates
      do j=1,dim; do i=1,nobs; x(i,j)=randn(); end do; end do
      if(rb) then
        call cov_detmcd(x,mcd); center=mcd%estimate%center; cov=mcd%estimate%covariance
      else
        call covariance_matrix(x,cov,center)
      end if
      call invert_matrix(cov,inv)
      do i=1,nobs
        u=x(i,:)-center
        s(r,i)=sqrt(max(0.0_dp,dot_product(u,matmul(inv,u))))
      end do
    end do
  end function r_empirical_mahalanobis

  function r_max_mahalanobis(replicates,nobs,dim,robust,seed) result(s)
    integer, intent(in) :: replicates,nobs,dim
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), allocatable :: s(:),all(:,:)
    integer :: r
    all=r_empirical_mahalanobis(replicates,nobs,dim,robust,seed)
    allocate(s(replicates))
    do r=1,replicates; s(r)=maxval(all(r,:)); end do
  end function r_max_mahalanobis

  real(dp) function q_empirical_mahalanobis(p,nobs,dim,replicates,robust,seed,pow) result(q)
    real(dp), intent(in) :: p
    integer, intent(in) :: nobs,dim,replicates
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: pow
    real(dp), allocatable :: s(:,:),v(:)
    real(dp) :: pw
    pw=1.0_dp; if(present(pow)) pw=pow
    s=r_empirical_mahalanobis(replicates,nobs,dim,robust,seed)
    v=reshape(s,[size(s)])
    call sort_real(v); q=v(prob_index(p,size(v)))**pw
  end function q_empirical_mahalanobis

  real(dp) function q_max_mahalanobis(p,nobs,dim,replicates,robust,seed,pow) result(q)
    real(dp), intent(in) :: p
    integer, intent(in) :: nobs,dim,replicates
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: pow
    real(dp), allocatable :: s(:)
    real(dp) :: pw
    pw=1.0_dp; if(present(pow)) pw=pow
    s=r_max_mahalanobis(replicates,nobs,dim,robust,seed)
    call sort_real(s); q=s(prob_index(p,size(s)))**pw
  end function q_max_mahalanobis

  real(dp) function p_empirical_mahalanobis(q,nobs,dim,replicates,robust,seed,pow) result(p)
    real(dp), intent(in) :: q
    integer, intent(in) :: nobs,dim,replicates
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: pow
    real(dp), allocatable :: s(:,:)
    real(dp) :: pw
    pw=1.0_dp; if(present(pow)) pw=pow
    s=r_empirical_mahalanobis(replicates,nobs,dim,robust,seed)
    p=(1.0_dp+real(count(s**pw<=q),dp))/real(size(s)+2,dp)
  end function p_empirical_mahalanobis

  real(dp) function p_max_mahalanobis(q,nobs,dim,replicates,robust,seed,pow) result(p)
    real(dp), intent(in) :: q
    integer, intent(in) :: nobs,dim,replicates
    logical, intent(in), optional :: robust
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: pow
    real(dp), allocatable :: s(:)
    real(dp) :: pw
    pw=1.0_dp; if(present(pow)) pw=pow
    s=r_max_mahalanobis(replicates,nobs,dim,robust,seed)
    p=(1.0_dp+real(count(s**pw<=q),dp))/real(size(s)+2,dp)
  end function p_max_mahalanobis

  function is_mahalanobis_outlier(x,alpha,replicates,corrected,robust,seed,critical) result(flag)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: alpha,critical
    integer, intent(in), optional :: replicates,seed
    logical, intent(in), optional :: corrected,robust
    logical, allocatable :: flag(:)
    real(dp), allocatable :: d(:)
    real(dp) :: a,crit
    integer :: nr,sd
    logical :: corr,rb
    a=0.05_dp; if(present(alpha)) a=alpha
    nr=1000; if(present(replicates)) nr=replicates
    corr=.true.; if(present(corrected)) corr=corrected
    rb=.true.; if(present(robust)) rb=robust
    sd=910241; if(present(seed)) sd=seed
    if(present(critical)) then
      crit=critical
    else if(corr) then
      crit=q_max_mahalanobis(1.0_dp-a,size(x,1),size(x,2)-1,nr,rb,sd)
    else
      crit=q_empirical_mahalanobis(1.0_dp-a,size(x,1),size(x,2)-1,nr,rb,sd)
    end if
    d=acomp_mahalanobis(x,rb); allocate(flag(size(d))); flag=(d>crit)
  end function is_mahalanobis_outlier

  function outlier_classifier_best(x,alpha,replicates,robust,seed) result(res)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: replicates,seed
    logical, intent(in), optional :: robust
    type(outlier_classification_result) :: res
    real(dp), allocatable :: xr(:,:),dr(:)
    logical, allocatable :: rf(:)
    real(dp) :: a,best
    integer :: nr,sd,n,d,j,i,k
    logical :: rb
    a=0.05_dp; if(present(alpha)) a=alpha
    nr=1000; if(present(replicates)) nr=replicates
    rb=.true.; if(present(robust)) rb=robust
    sd=135791; if(present(seed)) sd=seed
    n=size(x,1); d=size(x,2)
    res%distance=acomp_mahalanobis(x,rb)
    res%critical_corrected=q_max_mahalanobis(1.0_dp-a,n,d-1,nr,rb,sd)
    res%critical_empirical=q_empirical_mahalanobis(1.0_dp-a,n,d-1,nr,rb,sd+1)
    allocate(res%is_outlier(n),res%is_extreme(n),res%explained_by_component(n,d),res%best_component(n))
    res%is_outlier=(res%distance>res%critical_corrected)
    res%is_extreme=(res%distance>res%critical_empirical)
    res%explained_by_component=.false.; res%best_component=0
    if(d<=2) return
    allocate(xr(n,d-1))
    do j=1,d
      k=0
      do i=1,d
        if(i==j) cycle
        k=k+1; xr(:,k)=x(:,i)
      end do
      rf=is_mahalanobis_outlier(xr,a,nr,.false.,rb,sd+17*j)
      dr=acomp_mahalanobis(xr,rb)
      do i=1,n
        if(res%is_outlier(i).and..not.rf(i)) res%explained_by_component(i,j)=.true.
        if(res%is_outlier(i).and.res%explained_by_component(i,j)) then
          if(res%best_component(i)==0) then
            res%best_component(i)=j
          else
            best=reduced_distance_for_component(x,i,res%best_component(i),rb)
            if(dr(i)<best) res%best_component(i)=j
          end if
        end if
      end do
    end do
    where(res%is_outlier .and. res%best_component==0) res%best_component=-1
  end function outlier_classifier_best

  real(dp) function reduced_distance_for_component(x,row,component,robust) result(v)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: row,component
    logical, intent(in) :: robust
    real(dp), allocatable :: xr(:,:),d(:)
    integer :: j,k
    allocate(xr(size(x,1),size(x,2)-1)); k=0
    do j=1,size(x,2)
      if(j==component) cycle
      k=k+1; xr(:,k)=x(:,j)
    end do
    d=acomp_mahalanobis(xr,robust); v=d(row)
  end function reduced_distance_for_component

  pure integer function prob_index(p,n) result(idx)
    real(dp), intent(in) :: p
    integer, intent(in) :: n
    idx=1+nint(max(0.0_dp,min(1.0_dp,p))*real(n-1,dp))
    idx=max(1,min(n,idx))
  end function prob_index

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

end module compositions_outliers
